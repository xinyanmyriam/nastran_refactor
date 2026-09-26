"""Build held-out reference oracles from the UNMODIFIED NASTRAN-95 kernels.

Purpose (JSS plan section 8.5, P1 pilot):
  For each held-out input ``holdout_inputs/B{id}_h{k}.json`` we generate a
  parameterised fixed-form Fortran driver that (a) loads that *single* input
  record, (b) calls the **frozen, unmodified** mis/<kernel>.f routine, and
  (c) prints a machine-readable reference value.  The kernel source is never
  edited; only the generated driver changes per (bench, variant).

  Two invariants are enforced:
    * h0 regression: the generated h0 driver's output is compared against the
      existing frozen reference (e.g. b9_inverd/reference_values.json).  A
      mismatch aborts - it would mean the input mapping misrepresents the kernel.
    * kernel untouched: the kernel file bytes are hashed before and reused, and
      the oracle JSON records that source-hash + driver source itself.

Only benchmark B10 (INVERD) is implemented in this pilot; the other four
(B1 KROD, B3 KTRMEM, B11 QRITER, B13 SROD1) are stubbed with real templates to
follow.  The structure is a per-benchmark generator registry so adding a kernel
means adding one ``_gen_*`` + ``_parse_*`` pair.

Usage:
    E:\\...\\refactor\\.venv\\Scripts\\python.exe build_holdout_oracles.py    # all implemented benches
    ... build_holdout_oracles.py B10                                        # only B10
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
IN_DIR = ROOT / "holdout_inputs"
OUT_DIR = ROOT / "holdout_oracles"
BUILD_DIR = ROOT / "holdout_build"
MIS = ROOT / "nastran" / "NASTRAN-95" / "mis"
GFORTRAN = r"C:\Users\lxy\w64devkit\bin\gfortran.exe"
FCFLAGS = ["-ffixed-form", "-std=legacy", "-w", "-fallow-argument-mismatch"]

STANDARD_TOL = 1e-6  # double-precision oracle comparison tolerance (h0 regression)

# Set of benchmark ids already wired up.  Others print a clear "not yet" line.
BENCH_ORDER = ["B10", "B2", "B1", "B3", "B4", "B5", "B6", "B7", "B9", "B11", "B12", "B13"]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def _freal(v: float) -> str:
    """Best-effort double-precision literal from a Python float."""
    if v is None:
        raise ValueError("cannot emit None as a Fortran literal")
    return f"{v:.17g}D0"


def _ffloat(v: float) -> str:
    """Single-precision (REAL) literal from a Python float, no 'D' suffix."""
    if v is None:
        raise ValueError("cannot emit None as a Fortran literal")
    return f"{v:.10g}"


def _idxs(n: int) -> str:
    """Fortran list of '1..n' for DO-loop WRITE lists."""
    return ",".join(str(i) for i in range(1, n + 1))


def _dispatch_gen(bench: str, rec: dict[str, Any]) -> dict[str, str]:
    """Return {filename: source} for a benchmark's generated driver + stubs."""
    return GENERATORS[bench]["gen"](rec)


def _gen_b10(rec: dict[str, Any]) -> dict[str, str]:
    K = rec["K"]
    n = len(K)
    if rec.get("b") is not None:
        b = list(rec["b"])
    else:
        x = rec["x"]
        b = [sum(K[i][j] * x[j] for j in range(n)) for i in range(n)]
    lines: list[str] = []
    lines.append("      PROGRAM HOLDSOLVE")
    lines.append("      IMPLICIT DOUBLE PRECISION (A-H,O-Z)")
    lines.append(f"      DIMENSION A({n},{n}), B({n},1), INDEX({n},3)")
    lines.append("      INTEGER ISING")
    lines.append("      DOUBLE PRECISION DETERM")
    lines.append("      COMMON /MACHIN/ MACH")
    lines.append("      MACH=1")
    for i in range(n):
        for j in range(n):
            lines.append(f"      A({i+1},{j+1})={_freal(K[i][j])}")
    for i in range(n):
        lines.append(f"      B({i+1},1)={_freal(b[i])}")
    lines.append("      M=1")
    lines.append("      ISING=0")
    lines.append(f"      CALL INVERD({n},A,{n},B,1,DETERM,ISING,INDEX)")
    lines.append("      IF (ISING .NE. 1) WRITE(6,*) 'ISING',ISING")
    lines.append("      WRITE(6,'(A1)') 'S'")
    lines.append(f"      WRITE(6,'({n}(ES25.17,1X))') (B(I,1),I=1,{n})")
    lines.append("      WRITE(6,*) 'DET', DETERM")
    lines.append("      STOP")
    lines.append("      END")
    lines.append("")
    return {"driver.f": "\n".join(lines)}


def _parse_b10(stdout: str) -> dict[str, Any]:
    sol = None
    det = None
    flag_s = False
    for line in stdout.splitlines():
        s = line.strip()
        if s == "S":
            flag_s = True
            continue
        if flag_s:
            toks = s.replace("D", "E").split()
            sol = [float(t) for t in toks]
            flag_s = False
            continue
        if s.startswith("ISING"):
            return {"missing_rhs": False, "ising": s, "solution": None}
        if s.startswith("DET"):
            det = float(s.split()[1].replace("D", "E"))
    return {"solution": sol, "determinant": det}


# comparison of h0 output vs frozen reference  -> return list of float abs diffs
def _compare_b10(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    sol = out.get("solution")
    refsol = ref.get("solution")
    if sol is None or refsol is None:
        return ["missing solution vector"]
    diffs = [abs(a - b) for a, b in zip(sol, refsol)]
    return [f"max|diff|={max(diffs):.3e}"] if max(diffs) > STANDARD_TOL else []


# =====================================================================
# B13 (SROD1) - scalar stress recovery (rod along x: A at (0,0,0), B at (L,0,0))
# Input: L,A,E,nu,J,C, node_a[6], node_b[6], apply_torsion
# SROD1 runs on ECPT(17) and emits SAT/SBT (axial) and SAR/SBR (torsion).
# We apply the displacement/rotation held-out values by hand to get
# axial_stress, axial_force, torsional_stress.
# =====================================================================
def _gen_b13(rec: dict[str, Any]) -> dict[str, str]:
    L = rec["L"]
    A = rec["A"]
    E = rec["E"]
    nu = rec["nu"]
    J = rec["J"]
    C = rec["C"]
    na = rec["node_a"]  # [ux,uy,uz, rx,ry,rz]
    nb = rec["node_b"]
    ua = na[:3]
    tha = na[3:6]
    ub = nb[:3]
    thb = nb[3:6]
    cjc = J / C if C else 0.0
    lines = []
    lines.append("      PROGRAM HOLDSROD1")
    lines.append("      IMPLICIT NONE")
    lines.append("C     /SDR2X5/ as srod1.f declares it")
    lines.append("      REAL ECPT(17), DUMMY1(83)")
    lines.append("      INTEGER IELID, ISILNO(2)")
    lines.append("      REAL SAT(3), SBT(3), SAR(3), SBR(3)")
    lines.append("      REAL ST, SDELTA, AREA, FJOVRC")
    lines.append("      REAL TSUBC0, SIGMAT, SIGMAC, SIGMAS")
    lines.append("      REAL SIGVEC(77), FORVEC(25)")
    lines.append("      COMMON /SDR2X5/ ECPT, DUMMY1, IELID, ISILNO,")
    lines.append("     1  SAT, SBT, SAR, SBR, ST, SDELTA, AREA, FJOVRC,")
    lines.append("     2  TSUBC0, SIGMAT, SIGMAC, SIGMAS, SIGVEC, FORVEC")
    lines.append("      REAL XN(6), TI(9)")
    lines.append("      REAL XL, EOVERL")
    lines.append("      INTEGER IBASE")
    lines.append("      COMMON /SDR2X6/ XN, TI, XL, EOVERL, IBASE")
    lines.append("      INTEGER IECPT(17)")
    lines.append("      EQUIVALENCE (IECPT(1), ECPT(1))")
    lines.append("      REAL AXS, AXF, TORS")
    lines.append("      INTEGER I")
    lines.append("      DO 10 I = 1, 17")
    lines.append("        ECPT(I) = 0.0")
    lines.append("   10 CONTINUE")
    lines.append("      IECPT(1)=100")
    lines.append("      IECPT(2)=1")
    lines.append("      IECPT(3)=2")
    lines.append("      IECPT(4)=1")
    lines.append(f"      ECPT(5)= {_ffloat(A)}")
    lines.append(f"      ECPT(6)= {_ffloat(J)}")
    lines.append(f"      ECPT(7)= {_ffloat(cjc)}")
    lines.append("      ECPT(8)= 0.0")
    lines.append("      IECPT(9)=0")
    lines.append("      ECPT(10)= 0.0")
    lines.append("      ECPT(11)= 0.0")
    lines.append("      ECPT(12)= 0.0")
    lines.append("      IECPT(13)=0")
    lines.append(f"      ECPT(14)= {_ffloat(L)}")
    lines.append("      ECPT(15)= 0.0")
    lines.append("      ECPT(16)= 0.0")
    lines.append("      ECPT(17)= 0.0")
    lines.append("      CALL SROD1")
    lines.append("      AXS = SAT(1)*" + _ffloat(ua[0]) + "+SAT(2)*" + _ffloat(ua[1]) +
                 "+SAT(3)*" + _ffloat(ua[2]) +
                 "+SBT(1)*" + _ffloat(ub[0]) + "+SBT(2)*" + _ffloat(ub[1]) +
                 "+SBT(3)*" + _ffloat(ub[2]))
    lines.append("      AXF = AXS * AREA")
    lines.append("      TORS = SAR(1)*" + _ffloat(tha[0]) + "+SAR(2)*" + _ffloat(tha[1]) +
                 "+SAR(3)*" + _ffloat(tha[2]) +
                 "+SBR(1)*" + _ffloat(thb[0]) + "+SBR(2)*" + _ffloat(thb[1]) +
                 "+SBR(3)*" + _ffloat(thb[2]))
    if not rec.get("apply_torsion", False):
        lines.append("      TORS = 0.0")
    lines.append("      WRITE(6,'(A1)') 'S'")
    lines.append("      WRITE(6,'(3(ES25.17,1X))') AXS, AXF, TORS")
    lines.append("      STOP")
    lines.append("      END")
    # ---- MAT stub: E/nu -> /MATOUT/ + TRANSS + GMMATS (CSID=0 so unused but must link)
    stubs = []
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATIDC, MATFLG")
    stubs.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      REAL E, G, NU, RHO, ALPHA, TSUB0, GSUBE, SIGT, SIGC, SIGS")
    stubs.append("      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUB0, GSUBE,")
    stubs.append("     1  SIGT, SIGC, SIGS")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      NU = {_ffloat(nu)}")
    g = E / (2.0 * (1.0 + nu))
    stubs.append(f"      G = {_ffloat(g)}")
    stubs.append("      RHO = 0.0")
    stubs.append("      ALPHA = 0.0")
    stubs.append("      TSUB0 = 0.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGT = 0.0")
    stubs.append("      SIGC = 0.0")
    stubs.append("      SIGS = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSS (ICSID, TI)")
    stubs.append("      INTEGER ICSID")
    stubs.append("      REAL TI(9)")
    stubs.append("      TI(1)=1.0; TI(2)=0.0; TI(3)=0.0")
    stubs.append("      TI(4)=0.0; TI(5)=1.0; TI(6)=0.0")
    stubs.append("      TI(7)=0.0; TI(8)=0.0; TI(9)=1.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE GMMATS (A, NROWA, NCOLA, ITRANSA,")
    stubs.append("     1  B, NROWB, NCOLB, ITRANSB, C)")
    stubs.append("      INTEGER NROWA, NCOLA, ITRANSA, NROWB, NCOLB, ITRANSB")
    stubs.append("      REAL A(NROWA,NCOLA), B(NROWB,NCOLB), C(*)")
    stubs.append("      INTEGER I, J, K, IDX")
    stubs.append("      REAL SUM")
    stubs.append("      IDX = 0")
    stubs.append("      DO 100 J = 1, NCOLB")
    stubs.append("        DO 90 I = 1, NROWA")
    stubs.append("          SUM = 0.0")
    stubs.append("          DO 80 K = 1, NCOLA")
    stubs.append("            IF (ITRANSA .EQ. 0) THEN")
    stubs.append("              IF (ITRANSB .EQ. 0) SUM = SUM + A(I,K)*B(K,J)")
    stubs.append("              IF (ITRANSB .NE. 0) SUM = SUM + A(I,K)*B(J,K)")
    stubs.append("            ELSE")
    stubs.append("              IF (ITRANSB .EQ. 0) SUM = SUM + A(K,I)*B(K,J)")
    stubs.append("              IF (ITRANSB .NE. 0) SUM = SUM + A(K,I)*B(J,K)")
    stubs.append("            END IF")
    stubs.append("   80      CONTINUE")
    stubs.append("          IDX = IDX + 1")
    stubs.append("          C(IDX) = SUM")
    stubs.append("   90    CONTINUE")
    stubs.append("  100 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _parse_b13(stdout: str) -> dict[str, Any]:
    flag = False
    vals = []
    for line in stdout.splitlines():
        s = line.strip()
        if s == "S":
            flag = True
            continue
        if flag:
            vals = [float(t.replace("D", "E")) for t in s.split()]
            break
    if len(vals) < 3:
        return {"missing": True}
    return {
        "axial_stress": vals[0],
        "axial_force": vals[1],
        "torsional_stress": vals[2],
    }


def _compare_b13(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    issues = []
    for key in ("axial_stress", "axial_force", "torsional_stress"):
        a = out.get(key)
        b = ref.get(key)
        if a is None or b is None:
            issues.append(f"{key}: missing")
            continue
        if a == 0.0 and b == 0.0:
            continue
        rel = abs(a - b) / max(abs(b), 1e-300)
        if rel > 1e-4:
            issues.append(f"{key}: {a:.6e} vs ref {b:.6e} (rel {rel:.2e})")
    return issues


# =====================================================================
# B12 (TRD1C) - central-difference transient time response (timeseries)
# Driver mirrors source_fortran/test_drivers/b12_trd1c_native/
#   test_trd1c_native.f: sets up K/M/C, u0/v0, effective matrices, the
#   mini-GINO COMMON (gino_sim.f), calls TRD1C once (ILOOP=1, one group),
#   then prints displacement histories from the UDV (file 201) vector store.
# Output contract: n_samples rows x n_dof displacement values, 'S'-marked.
# The naive/legacy semantics (form1.f ISTART<0 path) are preserved exactly
# as the native driver does (BLANK(6)=ISTART=-1).
# =====================================================================
def _gen_B12(rec: dict[str, Any]) -> dict[str, str]:
    K = rec["K"]
    M = rec["M"]
    C = rec.get("C", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    u0 = rec["u0"]
    v0 = rec["v0"]
    dt = rec["dt"]
    total_steps = rec["total_steps"]
    n_samples = rec["n_samples"]
    stride = rec.get("sample_stride", 10)
    n = len(K)  # n_dof (3 for all B12 variants)
    li: list[str] = []
    a = li.append
    a("      PROGRAM HOLDB12")
    a("      IMPLICIT NONE")
    a("      REAL MAT_STORE, VEC_STORE")
    a("      INTEGER MCB_STORE, FILE_POS, FILE_NCOL, NROW_G")
    a("      COMMON /MINIGINO/")
    a("     1  MAT_STORE(3,3,10),")
    a("     2  VEC_STORE(3,30),")
    a("     3  MCB_STORE(7,30),")
    a("     4  FILE_POS(30),")
    a("     5  FILE_NCOL(30),")
    a("     6  NROW_G")
    a("      REAL BLANK_R")
    a("      COMMON /BLANK/ BLANK_R")
    a("      INTEGER BLANK_I(8)")
    a("      EQUIVALENCE (BLANK_R, BLANK_I)")
    a("      INTEGER SYSTEM_I(100)")
    a("      COMMON /SYSTEM/ SYSTEM_I")
    a("      REAL Z(10000)")
    a("      INTEGER IZ(10000)")
    a("      COMMON /ZZZZZZ/ Z")
    a("      EQUIVALENCE (Z, IZ)")
    a("      INTEGER PACKX_IT1,PACKX_IT2,PACKX_II,PACKX_JJ,PACKX_INCR")
    a("      COMMON /PACKX/ PACKX_IT1,PACKX_IT2,PACKX_II,PACKX_JJ,")
    a("     1               PACKX_INCR")
    a("      INTEGER UNPAKX_IT3,UNPAKX_III,UNPAKX_JJJ,UNPAKX_INCR1")
    a("      COMMON /UNPAKX/ UNPAKX_IT3,UNPAKX_III,UNPAKX_JJJ,")
    a("     1                UNPAKX_INCR1")
    a("      INTEGER TRDXX(32)")
    a("      REAL TRDXX_R(32)")
    a("      COMMON /TRDXX/ TRDXX")
    a("      EQUIVALENCE (TRDXX, TRDXX_R)")
    a("      INTEGER NAMES(19)")
    a("      COMMON /NAMES/ NAMES")
    a("      INTEGER INFBSX(14)")
    a("      COMMON /INFBSX/ INFBSX")
    a("      INTEGER MACHIN_MACH")
    a("      COMMON /MACHIN/ MACHIN_MACH")
    a("C")
    a("      INTEGER I,J,N,NSTEPS,NOUT,NGROUP,ISING,INDEX(3,3)")
    a("      REAL DT,KEFF(3,3),A1(3,3),A2(3,3),KEFF_INV(3,3)")
    a("      DOUBLE PRECISION DKEFF(3,3),DB(3,1),DDETERM")
    a("      INTEGER IC_FILE,PD_FILE,UDV_FILE,NCOL5")
    a("      INTEGER SCR1_FILE,SCR2_FILE,SCR3_FILE,SCR4_FILE,SCR5_FILE")
    a("C")
    a("      N=3")
    a(f"      DT={_ffloat(dt)}")
    a(f"      NSTEPS={int(total_steps)}")
    a(f"      NOUT={int(stride)}")
    a("      NGROUP=1")
    a("      IC_FILE=101")
    a("      PD_FILE=102")
    a("      UDV_FILE=201")
    a("      SCR1_FILE=301")
    a("      SCR2_FILE=302")
    a("      SCR3_FILE=303")
    a("      SCR4_FILE=304")
    a("      SCR5_FILE=305")
    a("C")
    a("C     store K->slot7, M->slot8, C->slot10 (mini-GINO dense slots)")
    for i in range(n):
        for j in range(n):
            a(f"      MAT_STORE({i+1},{j+1},7)={_ffloat(K[i][j])}")
    for i in range(n):
        for j in range(n):
            a(f"      MAT_STORE({i+1},{j+1},8)={_ffloat(M[i][j])}")
    for i in range(n):
        for j in range(n):
            a(f"      MAT_STORE({i+1},{j+1},10)={_ffloat(C[i][j])}")
    a("C")
    a("C     effective matrices")
    a("C     KEFF = M/dt^2 ; A1 = 2M/dt^2 - K ; A2 = -M/dt^2")
    a("      DO 40 I=1,N")
    a("        DO 41 J=1,N")
    a("          KEFF(I,J)=MAT_STORE(I,J,8)/(DT*DT)")
    a("          A1(I,J)=2.0*MAT_STORE(I,J,8)/(DT*DT)")
    a("     1            -MAT_STORE(I,J,7)")
    a("          A2(I,J)=-MAT_STORE(I,J,8)/(DT*DT)")
    a(" 41     CONTINUE")
    a(" 40   CONTINUE")
    a("C")
    a("C     invert KEFF via INVERD")
    a("      DO 42 I=1,N")
    a("        DO 43 J=1,N")
    a("          DKEFF(I,J)=DBLE(KEFF(I,J))")
    a(" 43     CONTINUE")
    a(" 42   CONTINUE")
    a("      DB(1,1)=0.0D0")
    a("      ISING=-1")
    a("      CALL INVERD(3,DKEFF,3,DB,0,DDETERM,ISING,INDEX)")
    a("      IF (ISING.EQ.2) THEN")
    a("        WRITE(6,*) 'B12 SINGULAR'")
    a("        STOP")
    a("      ENDIF")
    a("      DO 44 I=1,N")
    a("        DO 45 J=1,N")
    a("          MAT_STORE(I,J,9)=REAL(DKEFF(I,J))")
    a(" 45     CONTINUE")
    a(" 44   CONTINUE")
    a("C")
    a("C     slot1=A1 (SCR1), slot4=A2 (SCR4)")
    for i in range(n):
        for j in range(n):
            a(f"      MAT_STORE({i+1},{j+1},1)=A1({i+1},{j+1})")
    for i in range(n):
        for j in range(n):
            a(f"      MAT_STORE({i+1},{j+1},4)=A2({i+1},{j+1})")
    a("C")
    a("C     initial conditions in IC file (vec slot 1: u0, v0)")
    for i in range(n):
        a(f"      VEC_STORE({i+1},1)={_ffloat(u0[i])}")
    for i in range(n):
        a(f"      VEC_STORE({i+1},2)={_ffloat(v0[i])}")
    a("      FILE_NCOL(1)=2")
    a("      FILE_POS(1)=1")
    a("C")
    a("C     load vectors (PD, all zero) vec slot 3")
    a("      DO 46 J=1,20")
    a("        DO 47 I=1,N")
    a("          VEC_STORE(I,2+J)=0.0")
    a(" 47     CONTINUE")
    a(" 46   CONTINUE")
    a("      FILE_NCOL(3)=20")
    a("      FILE_POS(3)=1")
    a("C")
    a("C     UDV output file empty")
    a("      FILE_NCOL(5)=0")
    a("      FILE_POS(5)=1")
    a("C")
    a("C     MCB for SCR2 (slot16) and SCR3 (slot17)")
    a("      MCB_STORE(1,16)=SCR2_FILE")
    a("      MCB_STORE(2,16)=N")
    a("      MCB_STORE(3,16)=N")
    a("      MCB_STORE(4,16)=6")
    a("      MCB_STORE(5,16)=1")
    a("      MCB_STORE(6,16)=0")
    a("      MCB_STORE(7,16)=0")
    a("      MCB_STORE(1,17)=SCR3_FILE")
    a("      MCB_STORE(2,17)=N")
    a("      MCB_STORE(3,17)=N")
    a("      MCB_STORE(4,17)=7")
    a("      MCB_STORE(5,17)=1")
    a("      MCB_STORE(6,17)=0")
    a("      MCB_STORE(7,17)=0")
    a("C")
    a("C     common blocks")
    a("      DO 48 I=1,8")
    a("        BLANK_I(I)=0")
    a(" 48   CONTINUE")
    a("      BLANK_I(5)=0")
    a("      BLANK_I(6)=-1")
    a("      DO 49 I=1,100")
    a("        SYSTEM_I(I)=0")
    a(" 49   CONTINUE")
    a("      SYSTEM_I(1)=100")
    a("      SYSTEM_I(2)=6")
    a("      SYSTEM_I(82)=0")
    a("      NAMES(1)=0")
    a("      NAMES(2)=1")
    a("      NAMES(3)=2")
    a("      NAMES(4)=3")
    a("      NAMES(5)=1")
    a("      NAMES(6)=0")
    a("      NAMES(7)=7")
    a("      NAMES(8)=1")
    a("      NAMES(9)=2")
    a("      NAMES(10)=3")
    a("      NAMES(11)=4")
    a("      NAMES(12)=1")
    a("      NAMES(13)=2")
    a("      NAMES(14)=6")
    a("      NAMES(15)=5")
    a("      NAMES(16)=7")
    a("      NAMES(17)=8")
    a("      NAMES(18)=9")
    a("      NAMES(19)=3")
    a("C     TRDXX file ids")
    a("      TRDXX(1)=401")
    a("      TRDXX(2)=N")
    a("      TRDXX(3)=N")
    a("      TRDXX(4)=1")
    a("      TRDXX(5)=1")
    a("      TRDXX(6)=0")
    a("      TRDXX(7)=0")
    a("      TRDXX(8)=402")
    a("      TRDXX(9)=N")
    a("      TRDXX(10)=N")
    a("      TRDXX(11)=1")
    a("      TRDXX(12)=1")
    a("      TRDXX(13)=0")
    a("      TRDXX(14)=0")
    a("      TRDXX(15)=403")
    a("      TRDXX(16)=N")
    a("      TRDXX(17)=N")
    a("      TRDXX(18)=1")
    a("      TRDXX(19)=1")
    a("      TRDXX(20)=0")
    a("      TRDXX(21)=0")
    a("      TRDXX(22)=SCR1_FILE")
    a("      TRDXX(23)=SCR2_FILE")
    a("      TRDXX(24)=SCR3_FILE")
    a("      TRDXX(25)=SCR4_FILE")
    a("      TRDXX(26)=SCR5_FILE")
    a("      TRDXX(27)=306")
    a("      TRDXX(28)=0")
    a("      TRDXX(29)=1")
    a("      TRDXX_R(30)=0.0")
    a("      TRDXX(31)=0")
    a("      TRDXX(32)=0")
    a("      MACHIN_MACH=1")
    a("C")
    a("C     time-step group at end of Z core")
    a("      IZ(9998)=NSTEPS")
    a("      Z(9999)=DT")
    a("      IZ(10000)=NOUT")
    a("C")
    a("C     call TRD1C")
    a("      CALL TRD1C(IC_FILE,PD_FILE,NGROUP,0,UDV_FILE,1,")
    a("     1           SCR5_FILE,0,0,0,0,0,0)")
    a("C")
    a("C     print displacement history from UDV (each sample = 3 recs)")
    a("      NCOL5=FILE_NCOL(5)")
    a("      WRITE(6,'(A1)') 'S'")
    a("      J=0")
    a("      DO 50 I=1,NCOL5,3")
    a("        WRITE(6,'(3(ES25.17,1X))')")
    a("     1    VEC_STORE(1,4+I),VEC_STORE(2,4+I),VEC_STORE(3,4+I)")
    a("        J=J+1")
    a(" 50   CONTINUE")
    a("      STOP")
    a("      END")
    driver = "\n".join(li)
    # stubs: full mini-GINO layer from native driver + KORSZ
    gino_path = ROOT / "source_fortran/test_drivers/b12_trd1c_native/gino_sim.f"
    gino = gino_path.read_text(encoding="utf-8")
    stubs = gino + """
      INTEGER FUNCTION KORSZ (Z)
      REAL Z(1)
      KORSZ = 10000
      RETURN
      END
"""
    return {"driver.f": driver, "stubs.f": stubs}


def _parse_B12(stdout: str) -> dict[str, Any]:
    rows: list[list[float]] = []
    flag = False
    for line in stdout.splitlines():
        s = line.strip()
        if s == "S":
            flag = True
            continue
        if flag:
            if not s:
                continue
            toks = s.replace("D", "E").split()
            rows.append([float(t) for t in toks])
    return {"timeseries": rows}


def _compare_B12(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    t = out.get("timeseries")
    rt = ref.get("timeseries") if isinstance(ref, dict) else ref
    if not t:
        return ["missing timeseries"]
    if rt is None:
        return ["reference has no timeseries"]
    # rt may be the raw list (JSON loader) -> normalize
    if isinstance(rt, dict):
        rt = rt.get("timeseries")
    issues: list[str] = []
    if len(t) != len(rt):
        issues.append(f"row-count {len(t)} != ref {len(rt)}")
        return issues
    for i, (row, rrow) in enumerate(zip(t, rt)):
        if len(row) != len(rrow):
            issues.append(f"row {i}: cols {len(row)} != ref {len(rrow)}")
            continue
        for j, (a, b) in enumerate(zip(row, rrow)):
            if a == 0.0 and b == 0.0:
                continue
            # reference was formatted with F12.6 (native driver) so tiny
            # terms carry up to 0.5e-6 absolute rounding; use an abs floor
            # in addition to the 1e-3 relative tolerance.
            if abs(a - b) <= 1e-6:
                continue
            rel = abs(a - b) / max(abs(b), 1e-300)
            if rel > 1e-3:
                issues.append(f"[{i},{j}] {a:.6e} vs {b:.6e} rel {rel:.2e}")
    return issues


# =====================================================================
# B11 (QRITER) - tridiagonal eigenvalues (pure algorithm; no stubs)
# =====================================================================
def _gen_b11(rec: dict[str, Any]) -> dict[str, str]:
    n = int(rec["n"])
    diag = rec["diag"]
    off = rec["offdiag"]
    lines = []
    lines.append("      PROGRAM HOLDQRITER")
    lines.append("      IMPLICIT NONE")
    lines.append(f"      INTEGER NDIM")
    lines.append(f"      PARAMETER (NDIM = {n})")
    lines.append("      DOUBLE PRECISION VAL(NDIM), O(NDIM)")
    lines.append("      INTEGER LOC(NDIM), I")
    lines.append("      INTEGER IDUM0(100), IDUM3, IDUM4, IDUM9, IDUM11, IDUM12, IDUM13")
    lines.append("      INTEGER NN, LAMA, NV, NE, NFOUND, NEVER, NMAX")
    lines.append("      REAL LFREQ, HFREQ")
    lines.append("      COMMON /GIVN  / IDUM0, NN, LFREQ, IDUM3, IDUM4, HFREQ, LAMA,")
    lines.append("     1  NV, NE, IDUM9, NFOUND, IDUM11, IDUM12, IDUM13, NEVER, NMAX")
    lines.append("      INTEGER SYSBUF, NOUT")
    lines.append("      COMMON /SYSTEM/ SYSBUF, NOUT")
    lines.append("      INTEGER IOPTN")
    lines.append("      COMMON /REIGKR/ IOPTN")
    lines.append("      DOUBLE PRECISION DLMDAS")
    lines.append("      COMMON /MGIVXX/ DLMDAS")
    lines.append("      CHARACTER UFM*23, UWM*25, UIM*29")
    lines.append("      COMMON /XMSSG / UFM, UWM, UIM")
    for i in range(100):
        lines.append("      IDUM0(%d) = 0" % (i + 1))
    lines.append("      NN = NDIM")
    lines.append("      LFREQ = 0.0")
    lines.append("      HFREQ = 99.0")
    lines.append("      LAMA = 1")
    lines.append("      NV = 0")
    lines.append("      NE = 0")
    lines.append("      NFOUND = 0")
    lines.append("      NEVER = 0")
    lines.append("      IDUM3 = 0")
    lines.append("      IDUM4 = 0")
    lines.append("      IDUM9 = 0")
    lines.append("      IDUM11 = 0")
    lines.append("      IDUM12 = 0")
    lines.append("      IDUM13 = 0")
    lines.append("      NMAX = 100*NDIM")
    lines.append("      SYSBUF = 100")
    lines.append("      NOUT = 6")
    lines.append("      IOPTN = 0")
    lines.append("      DLMDAS = 0.0D0")
    lines.append("      UFM = 'U'")
    lines.append("      UWM = 'W'")
    lines.append("      UIM = 'I'")
    for i in range(n):
        lines.append(f"      VAL({i+1}) = {_freal(diag[i])}")
    for i in range(n - 1):
        lines.append(f"      O({i+1}) = {_freal(off[i]*off[i])}")
    lines.append(f"      O({n}) = 0.0D0")
    lines.append("      CALL QRITER(VAL,O,LOC,0)")
    lines.append("      WRITE(6,'(A1)') 'S'")
    lines.append(f"      WRITE(6,'({n}(ES25.17,1X))') (VAL(I),I=1,{n})")
    lines.append("      STOP")
    lines.append("      END")
    # QRITER's runtime deps (GOPEN/CLOSE/WRITE/WRTTRL/MESAGE/KORSZ) come from
    # the existing infrastructure stub bundle used by the b10_qriter driver.
    stub_path = ROOT / "source_fortran" / "test_drivers" / "b10_qriter" / "stubs_qriter.f"
    stubs = stub_path.read_text(encoding="utf-8")
    return {"driver.f": "\n".join(lines), "stubs.f": stubs}


def _parse_b11(stdout: str) -> dict[str, Any]:
    flag = False
    vals = []
    for line in stdout.splitlines():
        s = line.strip()
        if s == "S":
            flag = True
            continue
        if flag:
            vals = [float(t.replace("D", "E")) for t in s.split()]
            break
    return {"eigenvalues": vals}


def _compare_b11(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    o = out.get("eigenvalues") or []
    r = ref.get("eigenvalues") or []
    if len(o) != len(r):
        return [f"count mismatch: {len(o)} vs {len(r)}"]
    diffs = [abs(a - b) for a, b in zip(sorted(o), sorted(r))]
    return [f"max|diff|={max(diffs):.3e}"] if (diffs and max(diffs) > 1e-6) else []


# =====================================================================
# B1 (KROD) - 12x12 rod stiffness matrix.  Kernel krod.f + real gmmatd.f.
# SMA1B captures the 12x12 (two pivots) into /KCAP/. Material from a
# parameterised MAT stub.
# =====================================================================
def _gen_b1(rec: dict[str, Any]) -> dict[str, str]:
    A = rec["A"]
    J = rec["J"]
    E = rec["E"]
    G = rec["G"]
    na = rec["nodeA"]
    nb = rec["nodeB"]
    lines = []
    lines.append("      PROGRAM HOLDKROD")
    lines.append("      IMPLICIT NONE")
    lines.append("      INTEGER ICOM")
    lines.append("      COMMON /BLANK/ ICOM")
    lines.append("      INTEGER ISYS")
    lines.append("      COMMON /SYSTEM/ ISYS")
    lines.append("      INTEGER IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT")
    lines.append("      INTEGER IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG")
    lines.append("      INTEGER IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW")
    lines.append("      INTEGER CLSNRW, CLSRW, NEOR, EOR, MCBKGG(7), MCB4GG(7)")
    lines.append("      COMMON /SMA1IO/ IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT,")
    lines.append("     1  IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG,")
    lines.append("     2  IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW,")
    lines.append("     3  CLSNRW, CLSRW, NEOR, EOR, MCBKGG, MCB4GG")
    lines.append("      INTEGER ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT")
    lines.append("      INTEGER I6X6K, N6X6K, I6X64, N6X64")
    lines.append("      COMMON /SMA1BK/ ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT,")
    lines.append("     1  I6X6K, N6X6K, I6X64, N6X64")
    lines.append("      INTEGER IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC")
    lines.append("      INTEGER NROWSC, TNROWS, JMAX, NLINKS, LINK(10), IDETCK")
    lines.append("      INTEGER DODET, NOGO")
    lines.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC,")
    lines.append("     1  NROWSC, TNROWS, JMAX, NLINKS, LINK, IDETCK, DODET, NOGO")
    lines.append("      LOGICAL HEAT")
    lines.append("      COMMON /SMA1HT/ HEAT")
    lines.append("      REAL SECPT(100)")
    lines.append("      COMMON /SMA1ET/ SECPT")
    lines.append("      INTEGER MATIDC, MATFLG")
    lines.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    lines.append("      COMMON /MATIN/ MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH")
    lines.append("      REAL E, G2, NU, RHO, ALPHA, TSUBO, GSUBE")
    lines.append("      REAL SIGT, SIGC, SIGS")
    lines.append("      COMMON /MATOUT/ E, G2, NU, RHO, ALPHA, TSUBO, GSUBE,")
    lines.append("     1  SIGT, SIGC, SIGS")
    lines.append("      REAL FK")
    lines.append("      COMMON /HMTOUT/ FK")
    lines.append("      DOUBLE PRECISION KSTORE(12,12)")
    lines.append("      INTEGER NPVTG")
    lines.append("      COMMON /KCAP/ KSTORE, NPVTG")
    lines.append("      INTEGER I, J, IPVT")
    for i in range(1, 13):
        for j in range(1, 13):
            lines.append(f"      KSTORE({i},{j}) = 0.0D0")
    lines.append("      HEAT = .FALSE.")
    lines.append("      NOGO = 0")
    lines.append("      IOPT4 = 0")
    lines.append("      K4GGSW = 0")
    lines.append("      IFKGG = 11")
    lines.append("      IF4GG = 13")
    lines.append("      ISYS = 1024")
    lines.append("      ICOM = 0")
    lines.append("      MATFLG = 1")
    lines.append("      STRESS = 0.0")
    lines.append("      SINTH = 0.0")
    lines.append("      COSTH = 1.0")
    lines.append("      ELTEMP = 0.0")
    lines.append("      FK = 0.0")
    lines.append("      DODET = 0")
    lines.append("      IDETCK = 0")
    lines.append("      DO 100 IPVT = 1, 2")
    lines.append("        CALL SETECP")
    lines.append("        NPVT = 100 + IPVT")
    lines.append("        NPVTG = NPVT")
    lines.append("        CALL KROD")
    lines.append("  100 CONTINUE")
    lines.append("      WRITE(6,'(A1)') 'S'")
    for i in range(1, 13):
        lines.append(f"      WRITE(6,'({_idxs_inline(12)})') (KSTORE({i},J),J=1,12)")
    lines.append("      STOP")
    lines.append("      END")
    lines.append("      SUBROUTINE SETECP")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL SECPT(100)")
    lines.append("      INTEGER IECPT(100)")
    lines.append("      COMMON /SMA1ET/ SECPT")
    lines.append("      EQUIVALENCE (SECPT(1), IECPT(1))")
    lines.append("      INTEGER I")
    lines.append("      DO 10 I = 1, 100")
    lines.append("        SECPT(I) = 0.0")
    lines.append("   10 CONTINUE")
    lines.append("      IECPT(1)=1001")
    lines.append("      IECPT(2)=101")
    lines.append("      IECPT(3)=102")
    lines.append("      IECPT(4)=1")
    lines.append(f"      SECPT(5)= {_ffloat(A)}")
    lines.append(f"      SECPT(6)= {_ffloat(J)}")
    lines.append("      SECPT(7)= 0.0")
    lines.append("      SECPT(8)= 0.0")
    lines.append("      IECPT(9)=0")
    lines.append(f"      SECPT(10)= {_ffloat(na[0])}")
    lines.append(f"      SECPT(11)= {_ffloat(na[1])}")
    lines.append(f"      SECPT(12)= {_ffloat(na[2])}")
    lines.append("      IECPT(13)=0")
    lines.append(f"      SECPT(14)= {_ffloat(nb[0])}")
    lines.append(f"      SECPT(15)= {_ffloat(nb[1])}")
    lines.append(f"      SECPT(16)= {_ffloat(nb[2])}")
    lines.append("      SECPT(17)= 0.0")
    lines.append("      RETURN")
    lines.append("      END")
    # ---- stubs: MAT (param), SMA1B, TRANSD, MESAGE, HMAT
    stubs = []
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATIDC, MATFLG")
    stubs.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      COMMON /MATIN/ MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      REAL E, G, NU, RHO, ALPHA, TSUBO, GSUBE")
    stubs.append("      REAL SIGT, SIGC, SIGS")
    stubs.append("      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUBO, GSUBE,")
    stubs.append("     1  SIGT, SIGC, SIGS")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      G = {_ffloat(G)}")
    stubs.append("      NU = 0.3")
    stubs.append("      RHO = 0.0")
    stubs.append("      ALPHA = 0.0")
    stubs.append("      TSUBO = 0.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGT = 0.0")
    stubs.append("      SIGC = 0.0")
    stubs.append("      SIGS = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KE, IGRID, ICODE, IFFILE, FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KE(36), FACTOR")
    stubs.append("      INTEGER IGRID, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(12,12)")
    stubs.append("      INTEGER NPVTG")
    stubs.append("      COMMON /KCAP/ KSTORE, NPVTG")
    stubs.append("      INTEGER I, J, IR, IC, IROW, ICOL")
    stubs.append("      IROW = NPVTG - 100")
    stubs.append("      ICOL = IGRID - 100")
    stubs.append("      IF (IROW .LT. 1 .OR. IROW .GT. 2) RETURN")
    stubs.append("      IF (ICOL .LT. 1 .OR. ICOL .GT. 2) RETURN")
    stubs.append("      DO 20 I = 1, 6")
    stubs.append("        DO 10 J = 1, 6")
    stubs.append("          IR = 6*(IROW-1) + I")
    stubs.append("          IC = 6*(ICOL-1) + J")
    stubs.append("          KSTORE(IR,IC) = KSTORE(IR,IC) + KE(6*(I-1) + J)")
    stubs.append("   10   CONTINUE")
    stubs.append("   20 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (RCSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      REAL RCSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      INTEGER I")
    stubs.append("      DO 10 I = 1, 9")
    stubs.append("        T(I) = 0.0D0")
    stubs.append("   10 CONTINUE")
    stubs.append("      T(1)=1.0D0")
    stubs.append("      T(5)=1.0D0")
    stubs.append("      T(9)=1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG")
    stubs.append("      REAL IPARM")
    stubs.append("      WRITE(6,900) ICODE, IMSG")
    stubs.append("  900 FORMAT(' MESAGE called: code=', I6, ' msg=', I6)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE HMAT (NECPT)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER NECPT(4)")
    stubs.append("      REAL FK")
    stubs.append("      COMMON /HMTOUT/ FK")
    stubs.append("      FK = 50.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _idxs_inline(n: int) -> str:
    f = "(" + "1X,E16.9" + ")" if False else "12(1X,E16.9)"
    return f


def _parse_matrix(stdout: str, rows: int, cols: int) -> list[list[float]]:
    flag = False
    mat = []
    for line in stdout.splitlines():
        s = line.strip()
        if s == "S":
            flag = True
            continue
        if flag:
            toks = s.replace("D", "E").split()
            row = [float(t) for t in toks]
            mat.append(row)
            if len(mat) >= rows:
                break
    return mat[:rows]


def _compare_matrix(out: dict[str, Any], ref: dict[str, Any], size: int) -> list[str]:
    a = out.get("matrix")
    bm = ref.get("stiffness_matrix")
    if not a or not bm:
        return ["missing matrix"]
    worst = None
    for i in range(size):
        for j in range(size):
            av, bv = a[i][j], bm[i][j]
            scale = max(abs(bv), 1e-300)
            rel = abs(av - bv) / scale
            if worst is None or rel > worst:
                worst = rel
    return [f"max rel diff={worst:.3e}"] if (worst is not None and worst > 1e-4) else []


def _parse_b1(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 12, 12)}


def _compare_b1(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 12)


# =====================================================================
# B2 (KBAR) - 12x12 beam stiffness matrix.  Kernel kbar.f + real gmmatd.f.
# SMA1B captures the 12x12 (two pivots x 6 DOF) into /KCAP/. Material from
# a parameterised MAT stub.  KBAR's /SMA1ET/ layout is the 42-word CBAR
# ECPT; /MATIN/ is only 4 words (no SINTH/COSTH, unlike KROD).
#
# ECPT (kbar.f lines 6-49):
#   1 IELID | 2,3 ISILNO | 4-6 SMALLV | 7 ICSSV | 8,9 IPINFL
#   10-12 ZA | 13-15 ZB | 16 IMATID | 17 A | 18 I1 | 19 I2 | 20 FJ
#   21 NSM | 22 FE | 23-30 C1..G2 | 31 K1 | 32 K2 | 33 I12
#   34 MCSIDA | 35-37 GPA | 38 MCSIDB | 39-41 GPB | 42 TEMPEL
#   INTEGER words: 1, 2, 3, 7, 8, 9, 16, 34, 38.
# =====================================================================
def _gen_B2(rec: dict[str, Any]) -> dict[str, str]:
    A = rec["A"]
    Iy = rec["Iy"]
    Iz = rec["Iz"]
    J = rec["J"]
    E = rec["E"]
    G = rec["G"]
    na = rec["nodeA"]
    nb = rec["nodeB"]
    refv = rec["ref_vec"]
    lines = []
    lines.append("      PROGRAM HOLDKBAR")
    lines.append("      IMPLICIT NONE")
    lines.append("      INTEGER ISYS")
    lines.append("      COMMON /SYSTEM/ ISYS")
    lines.append("      INTEGER IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT")
    lines.append("      INTEGER IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG")
    lines.append("      INTEGER IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW")
    lines.append("      INTEGER CLSNRW, CLSRW, NEOR, EOR, MCBKGG(7), MCB4GG(7)")
    lines.append("      COMMON /SMA1IO/ IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT,")
    lines.append("     1  IGECPT, IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG,")
    lines.append("     2  IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW,")
    lines.append("     3  CLSNRW, CLSRW, NEOR, EOR, MCBKGG, MCB4GG")
    lines.append("      INTEGER ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT")
    lines.append("      INTEGER I6X6K, N6X6K, I6X64, N6X64")
    lines.append("      COMMON /SMA1BK/ ICSTM, NCSTM, IGPCT, NGPCT, IPOINT,")
    lines.append("     1  NPOINT, I6X6K, N6X6K, I6X64, N6X64")
    lines.append("      INTEGER IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC")
    lines.append("      INTEGER NROWSC, TNROWS, JMAX, NLINKS, LINK(10), IDETCK")
    lines.append("      INTEGER DODET, NOGO")
    lines.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, LEFT, FROWIC,")
    lines.append("     1  LROWIC, NROWSC, TNROWS, JMAX, NLINKS, LINK, IDETCK,")
    lines.append("     2  DODET, NOGO")
    lines.append("      LOGICAL HEAT")
    lines.append("      COMMON /SMA1HT/ HEAT")
    lines.append("      REAL SECPT(100)")
    lines.append("      COMMON /SMA1ET/ SECPT")
    lines.append("      INTEGER MATIDC, MATFLG")
    lines.append("      REAL ELTEMP, STRESS")
    lines.append("      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS")
    lines.append("      REAL E, G, NU, RHO, ALPHA, TSUBO, GSUBE")
    lines.append("      REAL SIGT, SIGC, SIGS")
    lines.append("      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUBO, GSUBE,")
    lines.append("     1  SIGT, SIGC, SIGS")
    lines.append("      REAL FK")
    lines.append("      COMMON /HMTOUT/ FK")
    lines.append("      DOUBLE PRECISION KSTORE(12,12)")
    lines.append("      INTEGER NPVTG")
    lines.append("      COMMON /KCAP / KSTORE, NPVTG")
    lines.append("      INTEGER I, J, IPVT")
    for i in range(1, 13):
        for j in range(1, 13):
            lines.append(f"      KSTORE({i},{j}) = 0.0D0")
    lines.append("      HEAT = .FALSE.")
    lines.append("      NOGO = 0")
    lines.append("      IOPT4 = 0")
    lines.append("      K4GGSW = 0")
    lines.append("      IFKGG = 11")
    lines.append("      IF4GG = 13")
    lines.append("      ISYS = 1024")
    lines.append("      MATFLG = 1")
    lines.append("      STRESS = 0.0")
    lines.append("      ELTEMP = 0.0")
    lines.append("      FK = 0.0")
    lines.append("      DODET = 0")
    lines.append("      IDETCK = 0")
    lines.append("      DO 100 IPVT = 1, 2")
    lines.append("        CALL SETECP")
    lines.append("        NPVT = 100 + IPVT")
    lines.append("        NPVTG = NPVT")
    lines.append("        CALL KBAR")
    lines.append("  100 CONTINUE")
    lines.append("      WRITE(6,'(A1)') 'S'")
    for i in range(1, 13):
        lines.append(f"      WRITE(6,920) (KSTORE({i},J), J = 1, 12)")
    lines.append("  920 FORMAT (12(1X,E16.9))")
    lines.append("      STOP")
    lines.append("      END")
    lines.append("      SUBROUTINE SETECP")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL SECPT(100)")
    lines.append("      INTEGER IECPT(100)")
    lines.append("      COMMON /SMA1ET/ SECPT")
    lines.append("      EQUIVALENCE (SECPT(1), IECPT(1))")
    lines.append("      INTEGER I")
    lines.append("      DO 10 I = 1, 100")
    lines.append("        SECPT(I) = 0.0")
    lines.append("   10 CONTINUE")
    lines.append("      IECPT(1) = 1001")
    lines.append("      IECPT(2) = 101")
    lines.append("      IECPT(3) = 102")
    lines.append(f"      SECPT(4) = {_ffloat(refv[0])}")
    lines.append(f"      SECPT(5) = {_ffloat(refv[1])}")
    lines.append(f"      SECPT(6) = {_ffloat(refv[2])}")
    lines.append("      IECPT(7) = 0")
    lines.append("      IECPT(8) = 0")
    lines.append("      IECPT(9) = 0")
    lines.append("      IECPT(16) = 1")
    lines.append(f"      SECPT(17) = {_ffloat(A)}")
    lines.append(f"      SECPT(18) = {_ffloat(Iy)}")
    lines.append(f"      SECPT(19) = {_ffloat(Iz)}")
    lines.append(f"      SECPT(20) = {_ffloat(J)}")
    lines.append("      SECPT(21) = 0.0")
    lines.append("      SECPT(22) = 0.0")
    lines.append("      SECPT(31) = 0.0")
    lines.append("      SECPT(32) = 0.0")
    lines.append("      SECPT(33) = 0.0")
    lines.append("      IECPT(34) = 0")
    lines.append(f"      SECPT(35) = {_ffloat(na[0])}")
    lines.append(f"      SECPT(36) = {_ffloat(na[1])}")
    lines.append(f"      SECPT(37) = {_ffloat(na[2])}")
    lines.append("      IECPT(38) = 0")
    lines.append(f"      SECPT(39) = {_ffloat(nb[0])}")
    lines.append(f"      SECPT(40) = {_ffloat(nb[1])}")
    lines.append(f"      SECPT(41) = {_ffloat(nb[2])}")
    lines.append("      SECPT(42) = 0.0")
    lines.append("      RETURN")
    lines.append("      END")
    # ---- stubs: MAT (param E/G, 4-word /MATIN/, 10-word /MATOUT/),
    #              SMA1B (12x12), TRANSD, MESAGE, HMAT
    stubs = []
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATIDC, MATFLG")
    stubs.append("      REAL ELTEMP, STRESS")
    stubs.append("      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS")
    stubs.append("      REAL E, G, NU, RHO, ALPHA, TSUBO, GSUBE")
    stubs.append("      REAL SIGT, SIGC, SIGS")
    stubs.append("      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUBO, GSUBE,")
    stubs.append("     1  SIGT, SIGC, SIGS")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      G = {_ffloat(G)}")
    stubs.append("      NU = 0.3")
    stubs.append("      RHO = 0.0")
    stubs.append("      ALPHA = 0.0")
    stubs.append("      TSUBO = 0.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGT = 0.0")
    stubs.append("      SIGC = 0.0")
    stubs.append("      SIGS = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KEP, IGRID, ICODE, IFFILE, FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KEP(36), FACTOR")
    stubs.append("      INTEGER IGRID, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(12,12)")
    stubs.append("      INTEGER NPVTG")
    stubs.append("      COMMON /KCAP / KSTORE, NPVTG")
    stubs.append("      INTEGER I, J, IR, IC, IROW, ICOL")
    stubs.append("      IROW = NPVTG - 100")
    stubs.append("      ICOL = IGRID - 100")
    stubs.append("      IF (IROW .LT. 1 .OR. IROW .GT. 2) RETURN")
    stubs.append("      IF (ICOL .LT. 1 .OR. ICOL .GT. 2) RETURN")
    stubs.append("      DO 20 I = 1, 6")
    stubs.append("        DO 10 J = 1, 6")
    stubs.append("          IR = 6*(IROW-1) + I")
    stubs.append("          IC = 6*(ICOL-1) + J")
    stubs.append("          KSTORE(IR,IC) = KSTORE(IR,IC) + KEP(6*(I-1) + J)")
    stubs.append("   10   CONTINUE")
    stubs.append("   20 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (RCSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      REAL RCSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      INTEGER I")
    stubs.append("      DO 10 I = 1, 9")
    stubs.append("        T(I) = 0.0D0")
    stubs.append("   10 CONTINUE")
    stubs.append("      T(1)=1.0D0")
    stubs.append("      T(5)=1.0D0")
    stubs.append("      T(9)=1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG")
    stubs.append("      REAL IPARM")
    stubs.append("      WRITE(6,900) ICODE, IMSG")
    stubs.append("  900 FORMAT(' MESAGE called: code=', I6, ' msg=', I6)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE HMAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      REAL FK")
    stubs.append("      COMMON /HMTOUT/ FK")
    stubs.append("      FK = 50.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _parse_B2(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 12, 12)}


def _compare_B2(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 12)


# =====================================================================
# B3 (KTRMEM) - in-plane 6x6 membrane stiffness.  Kernel ktrmem.f + gmmatd.f.
# Output printed as the in-plane 6x6 block.
# =====================================================================
def _gen_b3(rec: dict[str, Any]) -> dict[str, str]:
    nodes = rec["nodes"]  # [[x,y],...]
    E = rec["E"]
    nu = rec["nu"]
    t = rec["t"]
    x = [n[0] for n in nodes]
    y = [n[1] for n in nodes]
    lines = []
    lines.append("      PROGRAM HOLDKTRMEM")
    lines.append("      IMPLICIT NONE")
    lines.append("      INTEGER MECPT(1), NGRID(3), MATID1")
    lines.append("      REAL ANGLE, T, FMU, DUMMY1")
    lines.append("      REAL X1, Y1, Z1, DUMMY2, X2, Y2, Z2, DUMMY3")
    lines.append("      REAL X3, Y3, Z3, DUMB(80)")
    lines.append("      COMMON /SMA1ET/ MECPT, NGRID, ANGLE, MATID1, T, FMU, DUMMY1,")
    lines.append("     1  X1, Y1, Z1, DUMMY2, X2, Y2, Z2, DUMMY3,")
    lines.append("     2  X3, Y3, Z3, DUMB")
    lines.append("      LOGICAL HEAT")
    lines.append("      COMMON /SMA1HT/ HEAT")
    lines.append("      INTEGER MATID, INFLAG")
    lines.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    lines.append("      COMMON /MATIN/ MATID, INFLAG, ELTEMP, STRESS, SINTH, COSTH")
    lines.append("      DOUBLE PRECISION KSTORE(9,9)")
    lines.append("      INTEGER NPVTG")
    lines.append("      COMMON /KCAP/ KSTORE, NPVTG")
    lines.append("      INTEGER IOPT4, K4GGSW, NPVT, ISKP, NOGOO")
    lines.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO")
    lines.append("      INTEGER DUM1, IFKGG, DUM2, IF4GG, DUM3")
    lines.append("      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3")
    lines.append("      INTEGER SYSBUF, IOUT, NOGO")
    lines.append("      COMMON /SYSTEM/ SYSBUF, IOUT, NOGO")
    lines.append("      REAL CONSTS(5)")
    lines.append("      COMMON /CONDAS/ CONSTS")
    lines.append("      REAL MATBUF(4)")
    lines.append("      COMMON /HMTOUT/ MATBUF")
    lines.append("      INTEGER I, J, IPVT, NTYPE")
    lines.append("      DOUBLE PRECISION KSMALL(6,6)")
    for i in range(1, 10):
        for j in range(1, 10):
            lines.append(f"      KSTORE({i},{j}) = 0.0D0")
    lines.append("      HEAT = .FALSE.")
    lines.append("      NOGO = 0")
    lines.append("      NOGOO = 0")
    lines.append("      IOPT4 = 0")
    lines.append("      K4GGSW = 0")
    lines.append("      ISKP = 0")
    lines.append("      IFKGG = 1")
    lines.append("      IF4GG = 2")
    lines.append("      SYSBUF = 1024")
    lines.append("      IOUT = 6")
    lines.append("      INFLAG = 2")
    lines.append("      STRESS = 0.0")
    lines.append("      SINTH = 0.0")
    lines.append("      COSTH = 1.0")
    lines.append("      ELTEMP = 0.0")
    lines.append("      CONSTS(1)=0.0")
    lines.append("      CONSTS(2)=0.0")
    lines.append("      CONSTS(3)=0.0")
    lines.append("      CONSTS(4)= 3.1415926535897932 / 180.0")
    lines.append("      CONSTS(5)=0.0")
    lines.append("      DO 6 I = 1, 4")
    lines.append("        MATBUF(I) = 0.0")
    lines.append("    6 CONTINUE")
    lines.append("      MECPT(1)=1001")
    lines.append("      NGRID(1)=101")
    lines.append("      NGRID(2)=102")
    lines.append("      NGRID(3)=103")
    lines.append("      ANGLE = 0.0")
    lines.append("      MATID1 = 1")
    lines.append(f"      T = {_ffloat(t)}")
    lines.append("      FMU = 0.0")
    lines.append(f"      X1 = {_ffloat(x[0])}")
    lines.append(f"      Y1 = {_ffloat(y[0])}")
    lines.append("      Z1 = 0.0")
    lines.append(f"      X2 = {_ffloat(x[1])}")
    lines.append(f"      Y2 = {_ffloat(y[1])}")
    lines.append("      Z2 = 0.0")
    lines.append(f"      X3 = {_ffloat(x[2])}")
    lines.append(f"      Y3 = {_ffloat(y[2])}")
    lines.append("      Z3 = 0.0")
    lines.append("      NTYPE = 0")
    lines.append("      DO 100 IPVT = 1, 3")
    lines.append("        NPVT = 100 + IPVT")
    lines.append("        NPVTG = NPVT")
    lines.append("        CALL KTRMEM(NTYPE)")
    lines.append("  100 CONTINUE")
    for i in range(1, 4):
        for j in range(1, 4):
            lines.append(f"      KSMALL(2*{i}-1,2*{j}-1)= KSTORE(3*{i}-2,3*{j}-2)")
            lines.append(f"      KSMALL(2*{i}-1,2*{j}  )= KSTORE(3*{i}-2,3*{j}-1)")
            lines.append(f"      KSMALL(2*{i}  ,2*{j}-1)= KSTORE(3*{i}-1,3*{j}-2)")
            lines.append(f"      KSMALL(2*{i}  ,2*{j}  )= KSTORE(3*{i}-1,3*{j}-1)")
    lines.append("      WRITE(6,'(A1)') 'S'")
    for i in range(1, 7):
        lines.append(f"      WRITE(6,'(6(1X,E16.9))') (KSMALL({i},J),J=1,6)")
    lines.append("      STOP")
    lines.append("      END")
    # ---- stubs: MAT (param G matrix), SMA1B (9x9), TRANSD, MESAGE, HMAT
    stubs = []
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATID, INFLAG")
    stubs.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      COMMON /MATIN/ MATID, INFLAG, ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      REAL G11, G12, G13, G22, G23, G33, RHO")
    stubs.append("      REAL ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE")
    stubs.append("      REAL SIGTEN, SIGCOM, SIGSHE")
    stubs.append("      REAL G2X211, G2X212, G2X222")
    stubs.append("      COMMON /MATOUT/ G11, G12, G13, G22, G23, G33, RHO,")
    stubs.append("     1  ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE,")
    stubs.append("     2  SIGTEN, SIGCOM, SIGSHE, G2X211, G2X212, G2X222")
    stubs.append("      REAL E, NU, FAC")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      NU = {_ffloat(nu)}")
    stubs.append("      FAC = E / (1.0 - NU*NU)")
    stubs.append("      G11 = FAC")
    stubs.append("      G12 = NU * FAC")
    stubs.append("      G13 = 0.0")
    stubs.append("      G22 = FAC")
    stubs.append("      G23 = 0.0")
    stubs.append("      G33 = FAC * (1.0 - NU) / 2.0")
    stubs.append("      RHO = 0.0")
    stubs.append("      ALPHA1 = 0.0")
    stubs.append("      ALPHA2 = 0.0")
    stubs.append("      ALP12 = 0.0")
    stubs.append("      TSUB0 = 0.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGTEN = 0.0")
    stubs.append("      SIGCOM = 0.0")
    stubs.append("      SIGSHE = 0.0")
    stubs.append("      G2X211 = 0.0")
    stubs.append("      G2X212 = 0.0")
    stubs.append("      G2X222 = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KIJ(36), FACTOR")
    stubs.append("      INTEGER IGRID, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(9,9)")
    stubs.append("      INTEGER NPVTG")
    stubs.append("      COMMON /KCAP/ KSTORE, NPVTG")
    stubs.append("      INTEGER I, J, IR, IC, IROW, ICOL")
    stubs.append("      IROW = NPVTG - 100")
    stubs.append("      ICOL = IGRID - 100")
    stubs.append("      IF (IROW .LT. 1 .OR. IROW .GT. 3) RETURN")
    stubs.append("      IF (ICOL .LT. 1 .OR. ICOL .GT. 3) RETURN")
    stubs.append("      DO 20 I = 1, 3")
    stubs.append("        DO 10 J = 1, 3")
    stubs.append("          IR = 3*(IROW-1) + I")
    stubs.append("          IC = 3*(ICOL-1) + J")
    stubs.append("          KSTORE(IR,IC) = KSTORE(IR,IC) + KIJ(6*(I-1) + J)")
    stubs.append("   10   CONTINUE")
    stubs.append("   20 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (ICSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      INTEGER I")
    stubs.append("      DO 10 I = 1, 9")
    stubs.append("        T(I) = 0.0D0")
    stubs.append("   10 CONTINUE")
    stubs.append("      T(1)=1.0D0")
    stubs.append("      T(5)=1.0D0")
    stubs.append("      T(9)=1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG")
    stubs.append("      REAL IPARM")
    stubs.append("      WRITE(6,900) ICODE, IMSG")
    stubs.append("  900 FORMAT(' MESAGE called: code=', I6, ' msg=', I6)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE HMAT (NECPT)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER NECPT(4)")
    stubs.append("      REAL MATBUF(7)")
    stubs.append("      COMMON /HMTOUT/ MATBUF")
    stubs.append("      MATBUF(1)=50.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _parse_b3(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 6, 6)}


def _compare_b3(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 6)


# =====================================================================
# B4 (KQDMEM) - 8x8 in-plane quadrilateral membrane stiffness.
# Kernel kqdmem.f (calls KTRMEM + SMA1B) + real ktrmem.f + real gmmatd.f.
# KQDMEM does NOT split the quad into two CSTs; it uses four overlapping
# corner triangles at t/2, calling KTRMEM once per pivot for the three
# sub-triangles that contain that pivot (JNOT skips one).  Each call emits
# four 3x3 blocks via SMA1B into a 12x12 (4 nodes x 3 DOF); we then drop
# the out-of-plane (z) DOF rows/cols 3,6,9,12 to get the in-plane 8x8.
#
# ECPT layout required by KQDMEM (kqdmem.f comments):
#   1 EL.ID | 2-5 GRID A,B,C,D | 6 THETA | 7 MAT.ID | 8 T | 9 NSM
#   10 CSID1 | 11-13 X1,Y1,Z1 | 14 CSID2 | 15-17 X2,Y2,Z2
#   18 CSID3 | 19-21 X3,Y3,Z3 | 22 CSID4 | 23-25 X4,Y4,Z4 | 26 ELTEMP
#   INTEGER words: 1,2,3,4,5,7,10,14,18,22.
# KQDMEM rewrites ECPT(5..8) and 21, plus NECPT(2..4) and coordinates,
# for each sub-triangle, so SETECP must be re-run before every pivot.
# The ECPT/comms here mirror b4_ktriqd/test_kqdmem.f exactly.
# =====================================================================
def _gen_B4(rec: dict[str, Any]) -> dict[str, str]:
    nodes = rec["nodes"]  # [[x,y],...] 4 nodes
    E = rec["E"]
    nu = rec["nu"]
    t = rec["t"]
    x = [n[0] for n in nodes]
    y = [n[1] for n in nodes]
    lines = []
    lines.append("      PROGRAM HOLDKQDMEM")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL CONSTS(5)")
    lines.append("      COMMON /CONDAS/ CONSTS")
    lines.append("      LOGICAL HEAT")
    lines.append("      COMMON /SMA1HT/ HEAT")
    lines.append("      REAL SECPT(100)")
    lines.append("      COMMON /SMA1ET/ SECPT")
    lines.append("      INTEGER DUM1(10), IFKGG, DUM2(1), IF4GG, DUM3(23)")
    lines.append("      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3")
    lines.append("      INTEGER IOPT4, K4GGSW, NPVT, DUMCL(7), LINK(10)")
    lines.append("      INTEGER IDETCK, DODET, NOGO")
    lines.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, DUMCL, LINK,")
    lines.append("     1  IDETCK, DODET, NOGO")
    lines.append("      DOUBLE PRECISION DPBLK(600)")
    lines.append("      COMMON /SMA1DP/ DPBLK")
    lines.append("      INTEGER MATID, INFLAG")
    lines.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    lines.append("      COMMON /MATIN/ MATID, INFLAG, ELTEMP, STRESS, SINTH, COSTH")
    lines.append("      REAL G11, G12, G13, G22, G23, G33, RHO")
    lines.append("      REAL ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE")
    lines.append("      REAL SIGTEN, SIGCOM, SIGSHE")
    lines.append("      COMMON /MATOUT/ G11, G12, G13, G22, G23, G33, RHO,")
    lines.append("     1  ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE,")
    lines.append("     2  SIGTEN, SIGCOM, SIGSHE")
    lines.append("      REAL MATBUF(7)")
    lines.append("      COMMON /HMTOUT/ MATBUF")
    lines.append("      INTEGER SYSBUF, IOUT, NOGOS")
    lines.append("      COMMON /SYSTEM/ SYSBUF, IOUT, NOGOS")
    lines.append("      REAL BSKIP, VOLUME, SURFAC")
    lines.append("      COMMON /BLANK/ BSKIP, VOLUME, SURFAC")
    lines.append("      LOGICAL HYDRO")
    lines.append("      COMMON /HYDROE/ HYDRO")
    lines.append("      DOUBLE PRECISION KSTORE(12,12)")
    lines.append("      INTEGER NPVTG")
    lines.append("      COMMON /KCAP/ KSTORE, NPVTG")
    lines.append("      INTEGER I, J, IPVT")
    lines.append("      DOUBLE PRECISION KSMALL(8,8)")
    lines.append("      DO 5 I = 1, 12")
    lines.append("        DO 4 J = 1, 12")
    lines.append("          KSTORE(I,J) = 0.0D0")
    lines.append("    4   CONTINUE")
    lines.append("    5 CONTINUE")
    lines.append("      DO 6 I = 1, 600")
    lines.append("        DPBLK(I) = 0.0D0")
    lines.append("    6 CONTINUE")
    lines.append("      HEAT = .FALSE.")
    lines.append("      HYDRO = .FALSE.")
    lines.append("      NOGO = 0")
    lines.append("      NOGOS = 0")
    lines.append("      IOPT4 = 0")
    lines.append("      K4GGSW = 0")
    lines.append("      IFKGG = 1")
    lines.append("      IF4GG = 2")
    lines.append("      SYSBUF = 1024")
    lines.append("      IOUT = 6")
    lines.append("      INFLAG = 2")
    lines.append("      STRESS = 0.0")
    lines.append("      SINTH = 0.0")
    lines.append("      COSTH = 1.0")
    lines.append("      ELTEMP = 0.0")
    lines.append("      IDETCK = 0")
    lines.append("      DODET = 0")
    lines.append("      BSKIP = 0.0")
    lines.append("      VOLUME = 0.0")
    lines.append("      SURFAC = 0.0")
    lines.append("      CONSTS(1) = 0.0")
    lines.append("      CONSTS(2) = 0.0")
    lines.append("      CONSTS(3) = 0.0")
    lines.append("      CONSTS(4) = 3.1415926535897932 / 180.0")
    lines.append("      CONSTS(5) = 0.0")
    lines.append("      DO 7 I = 1, 7")
    lines.append("        MATBUF(I) = 0.0")
    lines.append("    7 CONTINUE")
    lines.append("      DO 100 IPVT = 1, 4")
    lines.append("        CALL SETECP")
    lines.append("        NPVT = 100 + IPVT")
    lines.append("        NPVTG = NPVT")
    lines.append("        CALL KQDMEM")
    lines.append("  100 CONTINUE")
    lines.append("      DO 320 I = 1, 4")
    lines.append("        DO 310 J = 1, 4")
    lines.append("          KSMALL(2*I-1,2*J-1) = KSTORE(3*I-2,3*J-2)")
    lines.append("          KSMALL(2*I-1,2*J)   = KSTORE(3*I-2,3*J-1)")
    lines.append("          KSMALL(2*I,  2*J-1) = KSTORE(3*I-1,3*J-2)")
    lines.append("          KSMALL(2*I,  2*J)   = KSTORE(3*I-1,3*J-1)")
    lines.append("  310   CONTINUE")
    lines.append("  320 CONTINUE")
    lines.append("      WRITE(6,'(A1)') 'S'")
    lines.append("      DO 400 I = 1, 8")
    lines.append("        WRITE(6,930) (KSMALL(I,J), J = 1, 8)")
    lines.append("  400 CONTINUE")
    lines.append("  930 FORMAT (8(1X,E16.9))")
    lines.append("      STOP")
    lines.append("      END")
    # SETECP fills the kqdmem /SMA1ET/ layout, re-run per pivot
    lines.append("      SUBROUTINE SETECP")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL SECPT(100)")
    lines.append("      INTEGER IECPT(100)")
    lines.append("      COMMON /SMA1ET/ SECPT")
    lines.append("      EQUIVALENCE (SECPT(1), IECPT(1))")
    lines.append("      INTEGER I")
    lines.append("      DO 10 I = 1, 100")
    lines.append("        SECPT(I) = 0.0")
    lines.append("   10 CONTINUE")
    lines.append("      IECPT(1) = 1001")
    lines.append("      IECPT(2) = 101")
    lines.append("      IECPT(3) = 102")
    lines.append("      IECPT(4) = 103")
    lines.append("      IECPT(5) = 104")
    lines.append("      SECPT(6) = 0.0")
    lines.append("      IECPT(7) = 1")
    lines.append(f"      SECPT(8) = {_ffloat(t)}")
    lines.append("      SECPT(9) = 0.0")
    lines.append("      IECPT(10) = 0")
    lines.append(f"      SECPT(11) = {_ffloat(x[0])}")
    lines.append(f"      SECPT(12) = {_ffloat(y[0])}")
    lines.append("      SECPT(13) = 0.0")
    lines.append("      IECPT(14) = 0")
    lines.append(f"      SECPT(15) = {_ffloat(x[1])}")
    lines.append(f"      SECPT(16) = {_ffloat(y[1])}")
    lines.append("      SECPT(17) = 0.0")
    lines.append("      IECPT(18) = 0")
    lines.append(f"      SECPT(19) = {_ffloat(x[2])}")
    lines.append(f"      SECPT(20) = {_ffloat(y[2])}")
    lines.append("      SECPT(21) = 0.0")
    lines.append("      IECPT(22) = 0")
    lines.append(f"      SECPT(23) = {_ffloat(x[3])}")
    lines.append(f"      SECPT(24) = {_ffloat(y[3])}")
    lines.append("      SECPT(25) = 0.0")
    lines.append("      SECPT(26) = 0.0")
    lines.append("      RETURN")
    lines.append("      END")
    # ---- stubs: MAT (plane-stress G, 18-word /MATOUT/), SMA1B (12x12,
    #              3x3 blocks), TRANSD, MESAGE, HMAT
    stubs = []
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATID, INFLAG")
    stubs.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      COMMON /MATIN/ MATID, INFLAG, ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      REAL G11, G12, G13, G22, G23, G33, RHO")
    stubs.append("      REAL ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE")
    stubs.append("      REAL SIGTEN, SIGCOM, SIGSHE")
    stubs.append("      REAL G2X211, G2X212, G2X222")
    stubs.append("      COMMON /MATOUT/ G11, G12, G13, G22, G23, G33, RHO,")
    stubs.append("     1  ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE,")
    stubs.append("     2  SIGTEN, SIGCOM, SIGSHE, G2X211, G2X212, G2X222")
    stubs.append("      REAL E, NU, FAC")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      NU = {_ffloat(nu)}")
    stubs.append("      FAC = E / (1.0 - NU*NU)")
    stubs.append("      G11 = FAC")
    stubs.append("      G12 = NU * FAC")
    stubs.append("      G13 = 0.0")
    stubs.append("      G22 = FAC")
    stubs.append("      G23 = 0.0")
    stubs.append("      G33 = FAC * (1.0 - NU) / 2.0")
    stubs.append("      RHO = 0.0")
    stubs.append("      ALPHA1 = 0.0")
    stubs.append("      ALPHA2 = 0.0")
    stubs.append("      ALP12 = 0.0")
    stubs.append("      TSUB0 = 0.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGTEN = 0.0")
    stubs.append("      SIGCOM = 0.0")
    stubs.append("      SIGSHE = 0.0")
    stubs.append("      G2X211 = 0.0")
    stubs.append("      G2X212 = 0.0")
    stubs.append("      G2X222 = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KIJ(36), FACTOR")
    stubs.append("      INTEGER IGRID, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(12,12)")
    stubs.append("      INTEGER NPVTG")
    stubs.append("      COMMON /KCAP/ KSTORE, NPVTG")
    stubs.append("      INTEGER I, J, IR, IC, IROW, ICOL")
    stubs.append("      IROW = NPVTG - 100")
    stubs.append("      ICOL = IGRID - 100")
    stubs.append("      IF (IROW .LT. 1 .OR. IROW .GT. 4) RETURN")
    stubs.append("      IF (ICOL .LT. 1 .OR. ICOL .GT. 4) RETURN")
    stubs.append("      DO 20 I = 1, 3")
    stubs.append("        DO 10 J = 1, 3")
    stubs.append("          IR = 3*(IROW-1) + I")
    stubs.append("          IC = 3*(ICOL-1) + J")
    stubs.append("          KSTORE(IR,IC) = KSTORE(IR,IC) + KIJ(6*(I-1) + J)")
    stubs.append("   10   CONTINUE")
    stubs.append("   20 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (ICSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      INTEGER I")
    stubs.append("      DO 10 I = 1, 9")
    stubs.append("        T(I) = 0.0D0")
    stubs.append("   10 CONTINUE")
    stubs.append("      T(1) = 1.0D0")
    stubs.append("      T(5) = 1.0D0")
    stubs.append("      T(9) = 1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG")
    stubs.append("      REAL IPARM")
    stubs.append("      WRITE (6,900) ICODE, IMSG")
    stubs.append("  900 FORMAT (' MESAGE called: code=', I6, ' msg=', I6)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE HMAT (NECPT)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER NECPT(4)")
    stubs.append("      REAL MATBUF(7)")
    stubs.append("      COMMON /HMTOUT/ MATBUF")
    stubs.append("      MATBUF(1) = 50.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _parse_B4(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 8, 8)}


def _compare_B4(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 8)


# =====================================================================
# B7 (KTETRA) - 4-node 3D tetrahedron solid, 12x12 stiffness matrix.
# Kernel ktetra.f (calls INVERD for the 4x4 shape-diff H inverse, GMMATD
# for the C/G products) + real inverd.f + real gmmatd.f.
#
# /SMA1ET/ ECPT (ktetra.f lines 10-34):
#   1 IELID | 2 MAT.ID | 3-6 SIL(1..4) | 7-10 CSID/X1,Y1,Z1
#   11-14 CSID/X2,Y2,Z2 | 15-18 CSID/X3,Y3,Z3 | 19-22 CSID/X4,Y4,Z4 | 23 ELTEMP
#   INTEGER words: 1,2,3,4,5,6,7,11,15,19.
# KTETRA is called once per grid point as NPVT (4 pivots); each call emits
# a 3x3 KIJ sub-block per column grid via SMA1B.  We capture the 3x3 blocks
# into /CAPTURE/ (as test_drivers/stubs.f) and GET_KMATRIX assembles 12x12.
# Material comes from a parameterised MAT stub (3-word /MATIN/, 10-word
# /MATOUT/).  CSID all 0 => TRANSD never called (kept as identity stub).
# =====================================================================
def _gen_B7(rec: dict[str, Any]) -> dict[str, str]:
    nodes = rec["nodes3d"]  # [[x,y,z],...] 4 nodes
    E = rec["E"]
    nu = rec["nu"]
    GG = E / (2.0 * (1.0 + nu))
    lines = []
    lines.append("      PROGRAM HOLDTETRA")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL ECPT(100)")
    lines.append("      INTEGER NECPT(100)")
    lines.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    lines.append("      COMMON /SMA1ET/ ECPT")
    lines.append("      DOUBLE PRECISION C(72), G(36), H(16), TEMP(12), T(9)")
    lines.append("      DOUBLE PRECISION CT(18), GCT(18), KIJ(36)")
    lines.append("      DOUBLE PRECISION HDETER, TEMP1, TVOL")
    lines.append("      INTEGER NGPT, DIREC, KOUNT")
    lines.append("      COMMON /SMA1DP/ C, G, H, TEMP, T, CT, GCT, KIJ,")
    lines.append("     1                HDETER, TEMP1, NGPT, DIREC, KOUNT, TVOL")
    lines.append("      LOGICAL HEAT")
    lines.append("      COMMON /SMA1HT/ HEAT")
    lines.append("      LOGICAL HYDRO")
    lines.append("      COMMON /HYDROE/ HYDRO")
    lines.append("      INTEGER MATID, INFLAG")
    lines.append("      REAL ELTEMP")
    lines.append("      COMMON /MATIN / MATID, INFLAG, ELTEMP")
    lines.append("      REAL E, GG, NU, RHO, ALPHA, TSUB0, GSUBE")
    lines.append("      REAL SIGT, SIGC, SIGS")
    lines.append("      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUB0, GSUBE,")
    lines.append("     1                SIGT, SIGC, SIGS")
    lines.append("      REAL MATBUF(7)")
    lines.append("      COMMON /HMTOUT/ MATBUF")
    lines.append("      INTEGER IOPT4, K4GGSW, NPVT, ISKP(17)")
    lines.append("      LOGICAL NOGOO")
    lines.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO")
    lines.append("      REAL DUM1(10), IFKGG, DUM2(1), IF4GG, DUM3(23)")
    lines.append("      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3")
    lines.append("      INTEGER SYSBUF, OUT")
    lines.append("      LOGICAL NOGO")
    lines.append("      COMMON /SYSTEM/ SYSBUF, OUT, NOGO")
    lines.append("      INTEGER MACH")
    lines.append("      COMMON /MACHIN/ MACH")
    lines.append("      REAL SKIP(16), VOLUME, SURFAC")
    lines.append("      COMMON /BLANK / SKIP, VOLUME, SURFAC")
    lines.append("      DOUBLE PRECISION KFULL(12,12)")
    lines.append("      INTEGER IPVT, I, J")
    lines.append("      MACH = 1")
    lines.append("      OUT = 6")
    lines.append("      SYSBUF = 1")
    lines.append("      NOGO = .FALSE.")
    lines.append("      HEAT = .FALSE.")
    lines.append("      HYDRO = .FALSE.")
    lines.append("      NOGOO = .FALSE.")
    lines.append("      IOPT4 = 0")
    lines.append("      K4GGSW = 0")
    lines.append("      VOLUME = 0.0")
    lines.append("      SURFAC = 0.0")
    lines.append("      DO 5 I = 1, 10")
    lines.append("        DUM1(I) = 0.0")
    lines.append("    5 CONTINUE")
    lines.append("      IFKGG = 0.0")
    lines.append("      DUM2(1) = 0.0")
    lines.append("      IF4GG = 0.0")
    lines.append("      DO 6 I = 1, 23")
    lines.append("        DUM3(I) = 0.0")
    lines.append("    6 CONTINUE")
    lines.append("      DO 10 I = 1, 100")
    lines.append("        ECPT(I) = 0.0")
    lines.append("   10 CONTINUE")
    lines.append("      NECPT(1) = 1")
    lines.append("      NECPT(2) = 1")
    lines.append("      NECPT(3) = 101")
    lines.append("      NECPT(4) = 102")
    lines.append("      NECPT(5) = 103")
    lines.append("      NECPT(6) = 104")
    for i in range(4):
        base = 7 + 4 * i
        lines.append(f"      NECPT({base}) = 0")
        lines.append(f"      ECPT({base+1}) = {_ffloat(nodes[i][0])}")
        lines.append(f"      ECPT({base+2}) = {_ffloat(nodes[i][1])}")
        lines.append(f"      ECPT({base+3}) = {_ffloat(nodes[i][2])}")
    lines.append("      ECPT(23) = 0.0")
    lines.append("      DO 20 I = 1, 12")
    lines.append("        DO 15 J = 1, 12")
    lines.append("          KFULL(I,J) = 0.0D0")
    lines.append("   15   CONTINUE")
    lines.append("   20 CONTINUE")
    lines.append("      DO 100 IPVT = 1, 4")
    lines.append("        NPVT = NECPT(IPVT + 2)")
    lines.append("        NECPT(1) = IPVT")
    lines.append("        CALL KTETRA(0, 0)")
    lines.append("  100 CONTINUE")
    lines.append("      CALL GET_KMATRIX(KFULL)")
    lines.append("      WRITE(6,'(A1)') 'S'")
    lines.append("      DO 200 I = 1, 12")
    lines.append("        WRITE(6,300) (KFULL(I,J), J = 1, 12)")
    lines.append("  200 CONTINUE")
    lines.append("  300 FORMAT (12(1X,E16.9))")
    lines.append("      STOP")
    lines.append("      END")
    # ---- stubs: CAPTURE block, MAT (param E/NU/GG), SMA1B (3x3 blocks),
    #              GET_KMATRIX, TRANSD, MESAGE, WRITE, HMAT
    stubs = []
    stubs.append("      BLOCK DATA CAPTUREDATA")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KSTORE(3,3,4,4)")
    stubs.append("      INTEGER NBLOCKS")
    stubs.append("      COMMON /CAPTURE/ KSTORE, NBLOCKS")
    stubs.append("      DATA NBLOCKS /0/")
    stubs.append("      DATA KSTORE /144*0.0D0/")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATID, INFLAG")
    stubs.append("      REAL ELTEMP")
    stubs.append("      COMMON /MATIN / MATID, INFLAG, ELTEMP")
    stubs.append("      REAL E, GG, NU, RHO, ALPHA, TSUB0, GSUBE")
    stubs.append("      REAL SIGT, SIGC, SIGS")
    stubs.append("      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUB0, GSUBE,")
    stubs.append("     1                SIGT, SIGC, SIGS")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      NU = {_ffloat(nu)}")
    stubs.append(f"      GG = {_ffloat(GG)}")
    stubs.append("      RHO = 0.0")
    stubs.append("      ALPHA = 0.0")
    stubs.append("      TSUB0 = 0.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGT = 0.0")
    stubs.append("      SIGC = 0.0")
    stubs.append("      SIGS = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KIJ(6,6), FACTOR")
    stubs.append("      INTEGER IGRID, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(3,3,4,4)")
    stubs.append("      INTEGER NBLOCKS")
    stubs.append("      COMMON /CAPTURE/ KSTORE, NBLOCKS")
    stubs.append("      REAL ECPT(100)")
    stubs.append("      INTEGER NECPT(100)")
    stubs.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    stubs.append("      COMMON /SMA1ET/ ECPT")
    stubs.append("      INTEGER IOPT4, K4GGSW, NPVT, ISKP(17)")
    stubs.append("      LOGICAL NOGOO")
    stubs.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO")
    stubs.append("      INTEGER I, J, IPVT, IGRP")
    stubs.append("      IF (FACTOR .NE. 0.0D0) RETURN")
    stubs.append("      DO 10 IPVT = 1, 4")
    stubs.append("        IF (NPVT .EQ. NECPT(IPVT+2)) GO TO 20")
    stubs.append("   10 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("   20 DO 30 IGRP = 1, 4")
    stubs.append("        IF (IGRID .EQ. NECPT(IGRP+2)) GO TO 40")
    stubs.append("   30 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("   40 DO 50 I = 1, 3")
    stubs.append("        DO 50 J = 1, 3")
    stubs.append("          KSTORE(I, J, IPVT, IGRP) = KIJ(I, J)")
    stubs.append("   50 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE GET_KMATRIX (KFULL)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KFULL(12,12)")
    stubs.append("      DOUBLE PRECISION KSTORE(3,3,4,4)")
    stubs.append("      INTEGER NBLOCKS")
    stubs.append("      COMMON /CAPTURE/ KSTORE, NBLOCKS")
    stubs.append("      INTEGER I, J, IPVT, IGRP, IROW, ICOL")
    stubs.append("      DO 100 IPVT = 1, 4")
    stubs.append("        DO 90 IGRP = 1, 4")
    stubs.append("          DO 80 I = 1, 3")
    stubs.append("            DO 80 J = 1, 3")
    stubs.append("              IROW = (IPVT-1)*3 + I")
    stubs.append("              ICOL = (IGRP-1)*3 + J")
    stubs.append("              KFULL(IROW, ICOL) = KSTORE(I, J, IPVT, IGRP)")
    stubs.append("   80         CONTINUE")
    stubs.append("   90       CONTINUE")
    stubs.append("  100 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (ICSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      T(1)=1.0D0; T(2)=0.0D0; T(3)=0.0D0")
    stubs.append("      T(4)=0.0D0; T(5)=1.0D0; T(6)=0.0D0")
    stubs.append("      T(7)=0.0D0; T(8)=0.0D0; T(9)=1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG, IPARM(2)")
    stubs.append("      INTEGER OUT, SYSBUF")
    stubs.append("      LOGICAL NOGO")
    stubs.append("      COMMON /SYSTEM/ SYSBUF, OUT, NOGO")
    stubs.append("      WRITE(OUT,900) ICODE, IMSG")
    stubs.append("  900 FORMAT(' MESAGE code=',I6,' msg=',I6)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE WRITE (IFILE, BUF, NWORDS, IFLG)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IFILE, NWORDS, IFLG")
    stubs.append("      REAL BUF(1)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE HMAT (NECPT)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER NECPT(4)")
    stubs.append("      REAL MATBUF(7)")
    stubs.append("      COMMON /HMTOUT/ MATBUF")
    stubs.append("      MATBUF(1) = 50.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _parse_B7(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 12, 12)}


def _compare_B7(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 12)


def _load_b7_reference(text: str) -> dict[str, Any]:
    """Load the frozen ktetra reference from the plain-text output file.

    Format (test_drivers/reference_output.txt): a line 'Full 12x12
    Stiffness Matrix' is followed by 'Row/Col:...', a '====' separator,
    then 12 rows of 12 E14.7 values.
    """
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if "Full 12x12" in line:
            start = i + 1
            break
    if start is None:
        raise ValueError("B7 reference: no 'Full 12x12' header")
    idx = start
    while idx < len(lines) and "=====" not in lines[idx]:
        idx += 1
    idx += 1  # skip the '=====' separator
    mat = []
    while idx < len(lines) and len(mat) < 12:
        toks = lines[idx].strip().split()
        if len(toks) >= 12:
            mat.append([float(t) for t in toks[:12]])
        idx += 1
    if len(mat) != 12:
        raise ValueError(f"B7 reference: parsed {len(mat)} rows, expected 12")
    return {"stiffness_matrix": mat}


# =====================================================================
# B5 / B6 (KTRPLT / KQDPLT) - triangular & quadrilateral plate bending.
# Kernel ktrplt.f / kqdplt.f + real ktrbsc.f + real inverd.f + real
# gmmatd.f.  KTRPLT/KQDPLT call KTRBSC (which emits 6x6 SMA1B blocks
# whose bending DOF live at rows/cols 3,4,5 = w,thx,thy).  The driver
# captures each 6x6 block per (pivot,grid), then extracts the 9x9 / 12x12
# bending matrix by taking columns 3,4,5 of every 6-block.
#
# ECPT layout (ktrplt.f): as documented in mis/ktrplt.f:
#   1 IELID | 2,3,4 NGRID | 5 THETA | 6 MATID1 | 7 I | 8 MATID2
#   9 T2 | 10 NSM | 11 Z1 | 12 Z2 | CSID/XYZ triples | 25 ELTEMP
#   node A: NECPT(13)+ECPT(14,15,16); B: NECPT(17)+ECPT(18,19,20);
#   C: NECPT(21)+ECPT(22,23,24); ELTEMP at ECPT(25).
# kqdplt.f uses the SAME layout shifted right by one word (its ECPT
# comment 'ELEMENT' column): ESLID at 2..5, THETA=6, ... nodes at
#   A: NECPT(14)+ECPT(15..17); B: NECPT(18)+19..21; C: NECPT(22)+23..25;
#   D: NECPT(26)+27..29; ELTEMP at ECPT(30).
# KQDPLT performs an internal ECPT shift on each call, so the driver must
# restore ECPT from a saved copy before every pivot.
# =====================================================================
def _plate_stubs(E: float, nu: float) -> str:
    stubs = []
    stubs.append("      BLOCK DATA CAPTUREDATAP")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KSTORE(6,6,4,4)")
    stubs.append("      INTEGER CURPVT, NPVTS")
    stubs.append("      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS")
    stubs.append("      DATA CURPVT /0/")
    stubs.append("      DATA NPVTS /0/")
    stubs.append("      DATA KSTORE /576*0.0D0/")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATID, INFLAG")
    stubs.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      COMMON /MATIN / MATID, INFLAG, ELTEMP, STRESS,")
    stubs.append("     1  SINTH, COSTH")
    stubs.append("      REAL G11, G12, G13, G22, G23, G33, RHO")
    stubs.append("      REAL ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE")
    stubs.append("      REAL SIGTEN, SIGCOM, SIGSHE")
    stubs.append("      REAL G2X211, G2X212, G2X222")
    stubs.append("      COMMON /MATOUT/ G11, G12, G13, G22, G23, G33, RHO,")
    stubs.append("     1  ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE,")
    stubs.append("     2  SIGTEN, SIGCOM, SIGSHE, G2X211, G2X212, G2X222")
    stubs.append("      REAL EE, NU, FAC")
    stubs.append(f"      EE = {_ffloat(E)}")
    stubs.append(f"      NU = {_ffloat(nu)}")
    stubs.append("      IF (INFLAG .EQ. 3) THEN")
    stubs.append("        G2X211 = 0.0")
    stubs.append("        G2X212 = 0.0")
    stubs.append("        G2X222 = 0.0")
    stubs.append("      ELSE")
    stubs.append("        FAC = EE / (1.0 - NU*NU)")
    stubs.append("        G11 = FAC")
    stubs.append("        G12 = NU * FAC")
    stubs.append("        G13 = 0.0")
    stubs.append("        G22 = FAC")
    stubs.append("        G23 = 0.0")
    stubs.append("        G33 = FAC * (1.0 - NU) / 2.0")
    stubs.append("      END IF")
    stubs.append("      RHO = 7800.0")
    stubs.append("      ALPHA1 = 12.0E-6")
    stubs.append("      ALPHA2 = 12.0E-6")
    stubs.append("      ALP12 = 0.0")
    stubs.append("      TSUB0 = 20.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGTEN = 0.0")
    stubs.append("      SIGCOM = 0.0")
    stubs.append("      SIGSHE = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE,")
    stubs.append("     1  FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KIJ(6,6), FACTOR")
    stubs.append("      INTEGER IGRID, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(6,6,4,4)")
    stubs.append("      INTEGER CURPVT, NPVTS")
    stubs.append("      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS")
    stubs.append("      REAL ECPT(100)")
    stubs.append("      INTEGER NECPT(100)")
    stubs.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    stubs.append("      COMMON /SMA1ET/ ECPT")
    stubs.append("      INTEGER IOPT4, K4GGSW, NPVT")
    stubs.append("      REAL DUMCL(7)")
    stubs.append("      INTEGER LINK(10), IDETCK")
    stubs.append("      REAL DODET")
    stubs.append("      INTEGER NOGO")
    stubs.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, DUMCL, LINK,")
    stubs.append("     1  IDETCK, DODET, NOGO")
    stubs.append("      INTEGER I, J, IPVT, IGRP")
    stubs.append("      IF (FACTOR .NE. 0.0D0) RETURN")
    stubs.append("      DO 10 IPVT = 1, NPVTS")
    stubs.append("        IF (NPVT .EQ. NECPT(IPVT+1)) GO TO 20")
    stubs.append("   10 CONTINUE")
    stubs.append("      WRITE(6,*) 'SMA1B WARNING pivot', NPVT")
    stubs.append("      RETURN")
    stubs.append("   20 DO 30 IGRP = 1, NPVTS")
    stubs.append("        IF (IGRID .EQ. NECPT(IGRP+1)) GO TO 40")
    stubs.append("   30 CONTINUE")
    stubs.append("      WRITE(6,*) 'SMA1B WARNING grid', IGRID")
    stubs.append("      RETURN")
    stubs.append("   40 DO 50 I = 1, 6")
    stubs.append("        DO 50 J = 1, 6")
    stubs.append("          KSTORE(I, J, IPVT, IGRP) = KIJ(I, J)")
    stubs.append("   50 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (ICSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      T(1)=1.0D0; T(2)=0.0D0; T(3)=0.0D0")
    stubs.append("      T(4)=0.0D0; T(5)=1.0D0; T(6)=0.0D0")
    stubs.append("      T(7)=0.0D0; T(8)=0.0D0; T(9)=1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG, IPARM(2)")
    stubs.append("      WRITE(6,*) 'MESAGE code=', ICODE, ' msg=', IMSG")
    stubs.append("      RETURN")
    stubs.append("      END")
    return "\n".join(stubs)


def _gen_B5(rec: dict[str, Any]) -> dict[str, str]:
    nodes = rec["nodes3d"]  # 3 nodes, [x,y,z]
    E = rec["E"]
    nu = rec["nu"]
    t = rec["t"]
    Ipl = t * t * t / 12.0
    z1 = -t / 2.0
    z2 = t / 2.0
    L = []
    L.append("      PROGRAM HOLDKTRPLT")
    L.append("      IMPLICIT NONE")
    L.append("      REAL ECPT(100)")
    L.append("      INTEGER NECPT(100)")
    L.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    L.append("      COMMON /SMA1ET/ ECPT")
    L.append("      REAL CONSTS(5)")
    L.append("      COMMON /CONDAS/ CONSTS")
    L.append("      REAL DUM1(10)")
    L.append("      INTEGER IFKGG")
    L.append("      REAL DUM2(1)")
    L.append("      INTEGER IF4GG")
    L.append("      REAL DUM3(23)")
    L.append("      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3")
    L.append("      INTEGER IOPT4, K4GGSW, NPVT")
    L.append("      REAL DUMCL(7)")
    L.append("      INTEGER LINK(10), IDETCK")
    L.append("      REAL DODET")
    L.append("      INTEGER NOGO")
    L.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, DUMCL, LINK,")
    L.append("     1  IDETCK, DODET, NOGO")
    L.append("      DOUBLE PRECISION DPWORK(400)")
    L.append("      COMMON /SMA1DP/ DPWORK")
    L.append("      DOUBLE PRECISION KSTORE(6,6,4,4)")
    L.append("      INTEGER CURPVT, NPVTS")
    L.append("      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS")
    L.append("      DOUBLE PRECISION K9(9,9)")
    L.append("      INTEGER IPVT, IGRP, I, J, II, JJ")
    L.append("      CONSTS(4) = 3.14159265358979D0 / 180.0D0")
    L.append("      IOPT4 = 0")
    L.append("      K4GGSW = 0")
    L.append("      NOGO = 0")
    L.append("      IDETCK = 0")
    L.append("      DODET = 0.0")
    L.append("      IFKGG = 0")
    L.append("      IF4GG = 0")
    L.append("      NPVTS = 3")
    L.append("      DO 10 I = 1, 100")
    L.append("        ECPT(I) = 0.0")
    L.append("   10 CONTINUE")
    L.append("      NECPT(1) = 1001")
    L.append("      NECPT(2) = 101")
    L.append("      NECPT(3) = 102")
    L.append("      NECPT(4) = 103")
    L.append("      ECPT(5) = 0.0")
    L.append("      NECPT(6) = 1")
    L.append(f"      ECPT(7) = {_ffloat(Ipl)}")
    L.append("      NECPT(8) = 0")
    L.append("      ECPT(9) = 0.0")
    L.append("      ECPT(10) = 0.0")
    L.append(f"      ECPT(11) = {_ffloat(z1)}")
    L.append(f"      ECPT(12) = {_ffloat(z2)}")
    for i in range(3):
        base = 13 + 4 * i
        L.append(f"      NECPT({base}) = 0")
        L.append(f"      ECPT({base+1}) = {_ffloat(nodes[i][0])}")
        L.append(f"      ECPT({base+2}) = {_ffloat(nodes[i][1])}")
        L.append(f"      ECPT({base+3}) = {_ffloat(nodes[i][2])}")
    L.append("      ECPT(25) = 20.0")
    L.append("      DO 100 IPVT = 1, 3")
    L.append("        NPVT = NECPT(IPVT + 1)")
    L.append("        CALL KTRPLT")
    L.append("  100 CONTINUE")
    L.append("      DO 200 IPVT = 1, 3")
    L.append("        DO 200 IGRP = 1, 3")
    L.append("          DO 200 II = 1, 3")
    L.append("            DO 200 JJ = 1, 3")
    L.append("              K9(3*(IPVT-1)+II, 3*(IGRP-1)+JJ) =")
    L.append("     1          KSTORE(II+2, JJ+2, IPVT, IGRP)")
    L.append("  200 CONTINUE")
    L.append("      WRITE(6,'(A1)') 'S'")
    L.append("      DO 300 I = 1, 9")
    L.append("        WRITE(6,310) (K9(I,J), J = 1, 9)")
    L.append("  300 CONTINUE")
    L.append("  310 FORMAT (9(1X,E16.9))")
    L.append("      STOP")
    L.append("      END")
    return {"driver.f": "\n".join(L), "stubs.f": _plate_stubs(E, nu)}


def _parse_B5(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 9, 9)}


def _compare_B5(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 9)


def _gen_B6(rec: dict[str, Any]) -> dict[str, str]:
    nodes = rec["nodes3d"]  # 4 nodes
    E = rec["E"]
    nu = rec["nu"]
    t = rec["t"]
    Ipl = t * t * t / 12.0
    z1 = -t / 2.0
    z2 = t / 2.0
    L = []
    L.append("      PROGRAM HOLDKQDPLT")
    L.append("      IMPLICIT NONE")
    L.append("      REAL ECPT(100)")
    L.append("      INTEGER NECPT(100)")
    L.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    L.append("      COMMON /SMA1ET/ ECPT")
    L.append("      REAL CONSTS(5)")
    L.append("      COMMON /CONDAS/ CONSTS")
    L.append("      REAL DUM1(10)")
    L.append("      INTEGER IFKGG")
    L.append("      REAL DUM2(1)")
    L.append("      INTEGER IF4GG")
    L.append("      REAL DUM3(23)")
    L.append("      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3")
    L.append("      INTEGER IOPT4, K4GGSW, NPVT")
    L.append("      REAL DUMCL(7)")
    L.append("      INTEGER LINK(10), IDETCK")
    L.append("      REAL DODET")
    L.append("      INTEGER NOGO")
    L.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, DUMCL, LINK,")
    L.append("     1  IDETCK, DODET, NOGO")
    L.append("      DOUBLE PRECISION DPWORK(1000)")
    L.append("      COMMON /SMA1DP/ DPWORK")
    L.append("      DOUBLE PRECISION KSTORE(6,6,4,4)")
    L.append("      INTEGER CURPVT, NPVTS")
    L.append("      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS")
    L.append("      DOUBLE PRECISION K12(12,12)")
    L.append("      REAL ECPT_SAVE(100)")
    L.append("      INTEGER IPVT, IGRP, I, J, II, JJ")
    L.append("      CONSTS(4) = 3.14159265358979D0 / 180.0D0")
    L.append("      IOPT4 = 0")
    L.append("      K4GGSW = 0")
    L.append("      NOGO = 0")
    L.append("      IDETCK = 0")
    L.append("      DODET = 0.0")
    L.append("      IFKGG = 0")
    L.append("      IF4GG = 0")
    L.append("      NPVTS = 4")
    L.append("      DO 10 I = 1, 100")
    L.append("        ECPT(I) = 0.0")
    L.append("   10 CONTINUE")
    L.append("      NECPT(1) = 1001")
    L.append("      NECPT(2) = 101")
    L.append("      NECPT(3) = 102")
    L.append("      NECPT(4) = 103")
    L.append("      NECPT(5) = 104")
    L.append("      ECPT(6) = 0.0")
    L.append("      NECPT(7) = 1")
    L.append(f"      ECPT(8) = {_ffloat(Ipl)}")
    L.append("      NECPT(9) = 0")
    L.append("      ECPT(10) = 0.0")
    L.append("      ECPT(11) = 0.0")
    L.append(f"      ECPT(12) = {_ffloat(z1)}")
    L.append(f"      ECPT(13) = {_ffloat(z2)}")
    for i in range(4):
        base = 14 + 4 * i
        L.append(f"      NECPT({base}) = 0")
        L.append(f"      ECPT({base+1}) = {_ffloat(nodes[i][0])}")
        L.append(f"      ECPT({base+2}) = {_ffloat(nodes[i][1])}")
        L.append(f"      ECPT({base+3}) = {_ffloat(nodes[i][2])}")
    L.append("      ECPT(30) = 20.0")
    L.append("      DO 15 I = 1, 100")
    L.append("        ECPT_SAVE(I) = ECPT(I)")
    L.append("   15 CONTINUE")
    L.append("      DO 100 IPVT = 1, 4")
    L.append("        DO 20 I = 1, 100")
    L.append("          ECPT(I) = ECPT_SAVE(I)")
    L.append("   20   CONTINUE")
    L.append("        NPVT = NECPT(IPVT + 1)")
    L.append("        CALL KQDPLT")
    L.append("  100 CONTINUE")
    L.append("      DO 200 IPVT = 1, 4")
    L.append("        DO 200 IGRP = 1, 4")
    L.append("          DO 200 II = 1, 3")
    L.append("            DO 200 JJ = 1, 3")
    L.append("              K12(3*(IPVT-1)+II, 3*(IGRP-1)+JJ) =")
    L.append("     1          KSTORE(II+2, JJ+2, IPVT, IGRP)")
    L.append("  200 CONTINUE")
    L.append("      WRITE(6,'(A1)') 'S'")
    L.append("      DO 300 I = 1, 12")
    L.append("        WRITE(6,310) (K12(I,J), J = 1, 12)")
    L.append("  300 CONTINUE")
    L.append("  310 FORMAT (12(1X,E16.9))")
    L.append("      STOP")
    L.append("      END")
    return {"driver.f": "\n".join(L), "stubs.f": _plate_stubs(E, nu)}


def _parse_B6(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 12, 12)}


def _compare_B6(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 12)


def _load_plate_reference(text: str, header: str, n: int) -> dict[str, Any]:
    """Load NxN plate reference from the plain-text driver output."""
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if header in line:
            start = i + 1
            break
    if start is None:
        raise ValueError(f"plate reference: no '{header}' header")
    idx = start
    while idx < len(lines) and "=====" not in lines[idx]:
        idx += 1
    idx += 1  # skip separator
    mat = []
    while idx < len(lines) and len(mat) < n:
        toks = lines[idx].strip().split()
        if len(toks) >= n:
            mat.append([float(t) for t in toks[:n]])
        idx += 1
    if len(mat) != n:
        raise ValueError(f"plate reference: parsed {len(mat)} rows, expected {n}")
    return {"stiffness_matrix": mat}


def _load_b5_reference(text: str) -> dict[str, Any]:
    return _load_plate_reference(text, "9x9 Plate Bending", 9)


def _load_b6_reference(text: str) -> dict[str, Any]:
    return _load_plate_reference(text, "12x12 Plate Bending", 12)


# =====================================================================
# B8 (KSOLID ITYPE=1) - 6-node 3D wedge solid, 18x18 stiffness matrix.
# Kernel ksolid.f (decomposes the wedge into 3 tetrahedra and calls ktetra.f)
# + real ktetra.f + real inverd.f + real gmmatd.f.
#
# Reference is the FULL 18x18 wedge matrix (reference_matrix_wedge.json),
# stored as JSON under the 'stiffness_matrix' key -> default json.loads in the
# h0 gate works, so no custom load_ref is needed.
#
# /SMA1ET/ wedge ECPT (ksolid.f comments):
#   1 ELID | 2 MATID | 3-8 GRID(1..6) | per node i (1-based):
#     CSID at NECPT(5+4i), XYZ at ECPT(6+4i .. 8+4i)  (i.e. node 1: 9,10,11,12)
#   33 ELTEMP.  INTEGER words: 1,2,3..8,9,13,17,21,25,29.
# KSOLID IS CALLED ONCE PER PIVOT GRID POINT (6 pivots -> NPVT grid).  Each
# call rebuilds the per-tetra ECPT (thus the driver must restore ECPT from a
# saved copy before every pivot) and emits 3x3 KIJ sub-blocks via SMA1B.  We
# capture the 3x3 blocks into /CAPTURW/ (KSTORE(3,3,6,6)) and an
# GET_KMATRIX_WEDGE routine assembles the 18x18 (mirrors b8_ksolid stubs).
# Material comes from a parameterised MAT stub (3-word /MATIN/, 10-word
# /MATOUT/).  CSID all 0 => TRANSD leaves identity.  SAXB/SADOTB/KPLTST are
# required by ksolid.f (geometry checks).
# =====================================================================
def _gen_B8(rec: dict[str, Any]) -> dict[str, str]:
    nodes = rec["nodes3d"]  # 6 nodes, [[x,y,z],...]
    E = rec["E"]
    nu = rec["nu"]
    lines = []
    lines.append("      PROGRAM HOLDKSOLID")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL ECPT(100)")
    lines.append("      INTEGER NECPT(100)")
    lines.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    lines.append("      COMMON /SMA1ET/ ECPT")
    lines.append("      DOUBLE PRECISION C(72), G(36), H(16), TEMP(12), T(9)")
    lines.append("      DOUBLE PRECISION CT(18), GCT(18), KIJ(36), HDETER,")
    lines.append("     1  TEMP1")
    lines.append("      INTEGER NGPT, DIREC, KOUNT")
    lines.append("      DOUBLE PRECISION TVOL")
    lines.append("      COMMON /SMA1DP/ C, G, H, TEMP, T, CT, GCT, KIJ,")
    lines.append("     1  HDETER, TEMP1, NGPT, DIREC, KOUNT, TVOL")
    lines.append("      LOGICAL HEAT")
    lines.append("      COMMON /SMA1HT/ HEAT")
    lines.append("      LOGICAL HYDRO")
    lines.append("      COMMON /HYDROE/ HYDRO")
    lines.append("      INTEGER MATID, INFLAG")
    lines.append("      REAL ELTEMP")
    lines.append("      COMMON /MATIN / MATID, INFLAG, ELTEMP")
    lines.append("      REAL E, GG, NU, RHO, ALPHA, TSUB0, GSUBE")
    lines.append("      REAL SIGT, SIGC, SIGS")
    lines.append("      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUB0,")
    lines.append("     1  GSUBE, SIGT, SIGC, SIGS")
    lines.append("      REAL MATBUF(7)")
    lines.append("      COMMON /HMTOUT/ MATBUF")
    lines.append("      INTEGER IOPT4, K4GGSW, NPVT, ISKP(19)")
    lines.append("      LOGICAL NOGOO")
    lines.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO")
    lines.append("      REAL DUM1(10), IFKGG, DUM2(1), IF4GG, DUM3(23)")
    lines.append("      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3")
    lines.append("      INTEGER SYSBUF, OUT")
    lines.append("      LOGICAL NOGO")
    lines.append("      COMMON /SYSTEM/ SYSBUF, OUT, NOGO")
    lines.append("      INTEGER MACH")
    lines.append("      COMMON /MACHIN/ MACH")
    lines.append("      REAL BSKIP(16), VOLUME, SURFAC")
    lines.append("      COMMON /BLANK / BSKIP, VOLUME, SURFAC")
    lines.append("      CHARACTER UFM*23")
    lines.append("      COMMON /XMSSG / UFM")
    lines.append("      DOUBLE PRECISION KFULL(18,18)")
    lines.append("      INTEGER I, J, IPVT")
    lines.append("      REAL ECPT_SAVE(100)")
    lines.append("      UFM = ' *** USER FATAL MESSAGE'")
    lines.append("      MACH = 1")
    lines.append("      OUT = 6")
    lines.append("      SYSBUF = 1")
    lines.append("      NOGO = .FALSE.")
    lines.append("      HEAT = .FALSE.")
    lines.append("      HYDRO = .FALSE.")
    lines.append("      NOGOO = .FALSE.")
    lines.append("      IOPT4 = 0")
    lines.append("      K4GGSW = 0")
    lines.append("      VOLUME = 0.0")
    lines.append("      SURFAC = 0.0")
    lines.append("      DO 5 I = 1, 10")
    lines.append("        DUM1(I) = 0.0")
    lines.append("    5 CONTINUE")
    lines.append("      IFKGG = 0.0")
    lines.append("      DUM2(1) = 0.0")
    lines.append("      IF4GG = 0.0")
    lines.append("      DO 6 I = 1, 23")
    lines.append("        DUM3(I) = 0.0")
    lines.append("    6 CONTINUE")
    lines.append("      DO 10 I = 1, 100")
    lines.append("        ECPT(I) = 0.0")
    lines.append("   10 CONTINUE")
    lines.append("      NECPT(1) = 1")
    lines.append("      NECPT(2) = 1")
    for i in range(6):
        lines.append(f"      NECPT({3+i}) = {101+i}")
    for i, n in enumerate(nodes):
        base = 9 + 4 * i
        lines.append(f"      NECPT({base}) = 0")
        lines.append(f"      ECPT({base+1}) = {_ffloat(n[0])}")
        lines.append(f"      ECPT({base+2}) = {_ffloat(n[1])}")
        lines.append(f"      ECPT({base+3}) = {_ffloat(n[2])}")
    lines.append("      ECPT(33) = 20.0")
    lines.append("      DO 20 I = 1, 18")
    lines.append("        DO 15 J = 1, 18")
    lines.append("          KFULL(I,J) = 0.0D0")
    lines.append("   15   CONTINUE")
    lines.append("   20 CONTINUE")
    lines.append("      DO 22 I = 1, 100")
    lines.append("        ECPT_SAVE(I) = ECPT(I)")
    lines.append("   22 CONTINUE")
    lines.append("      DO 100 IPVT = 1, 6")
    lines.append("        DO 24 I = 1, 100")
    lines.append("          ECPT(I) = ECPT_SAVE(I)")
    lines.append("   24   CONTINUE")
    lines.append("        NPVT = 100 + IPVT")
    lines.append("        CALL KSOLID(1)")
    lines.append("  100 CONTINUE")
    lines.append("      CALL GET_KMATRIX_WEDGE(KFULL)")
    lines.append("      WRITE(6,'(A1)') 'S'")
    lines.append("      DO 200 I = 1, 18")
    lines.append("        WRITE(6,300) (KFULL(I,J), J = 1, 18)")
    lines.append("  200 CONTINUE")
    lines.append("  300 FORMAT (18(1X,E16.9))")
    lines.append("      STOP")
    lines.append("      END")
    # ---- stubs: CAPTURE block (18x18 via /CAPTURW/), MAT (param E/NU),
    #              SMA1B (3x3 node blocks), GET_KMATRIX_WEDGE, SAXB, SADOTB,
    #              KPLTST, TRANSD, MESAGE, WRITE, HMAT
    stubs = []
    stubs.append("      BLOCK DATA CAPTUREDATAW")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KSTORE(3,3,6,6)")
    stubs.append("      INTEGER NCAPTURE")
    stubs.append("      COMMON /CAPTURW/ KSTORE, NCAPTURE")
    stubs.append("      DATA NCAPTURE /0/")
    stubs.append("      DATA KSTORE /324*0.0D0/")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATID, INFLAG")
    stubs.append("      REAL ELTEMP")
    stubs.append("      COMMON /MATIN / MATID, INFLAG, ELTEMP")
    stubs.append("      REAL E, GG, NU, RHO, ALPHA, TSUB0, GSUBE, SIGT,")
    stubs.append("     1  SIGC, SIGS")
    stubs.append("      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUB0,")
    stubs.append("     1  GSUBE, SIGT, SIGC, SIGS")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      NU = {_ffloat(nu)}")
    stubs.append("      GG = E / (2.0*(1.0+NU))")
    stubs.append("      RHO = 7800.0")
    stubs.append("      ALPHA = 12.0E-6")
    stubs.append("      TSUB0 = 20.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGT = 0.0")
    stubs.append("      SIGC = 0.0")
    stubs.append("      SIGS = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE,")
    stubs.append("     1  FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KIJ(6,6), FACTOR")
    stubs.append("      INTEGER IGRID, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(3,3,6,6)")
    stubs.append("      INTEGER NCAPTURE")
    stubs.append("      COMMON /CAPTURW/ KSTORE, NCAPTURE")
    stubs.append("      REAL ECPT(100)")
    stubs.append("      INTEGER NECPT(100)")
    stubs.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    stubs.append("      COMMON /SMA1ET/ ECPT")
    stubs.append("      INTEGER IOPT4, K4GGSW, NPVT, ISKP(19)")
    stubs.append("      LOGICAL NOGOO")
    stubs.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO")
    stubs.append("      INTEGER I, J, IPVT, IGRP")
    stubs.append("      INTEGER GRIDS(6)")
    stubs.append("      DATA GRIDS /101,102,103,104,105,106/")
    stubs.append("      IF (FACTOR .NE. 0.0D0) RETURN")
    stubs.append("      DO 10 IPVT = 1, 6")
    stubs.append("        IF (NPVT .EQ. GRIDS(IPVT)) GO TO 20")
    stubs.append("   10 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("   20 DO 30 IGRP = 1, 6")
    stubs.append("        IF (IGRID .EQ. GRIDS(IGRP)) GO TO 40")
    stubs.append("   30 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("   40 DO 50 I = 1, 3")
    stubs.append("        DO 50 J = 1, 3")
    stubs.append("          KSTORE(I, J, IPVT, IGRP) = KSTORE(I, J, IPVT,")
    stubs.append("     1      IGRP) + KIJ(I, J)")
    stubs.append("   50 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE GET_KMATRIX_WEDGE (KFULL)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KFULL(18,18)")
    stubs.append("      DOUBLE PRECISION KSTORE(3,3,6,6)")
    stubs.append("      INTEGER NCAPTURE")
    stubs.append("      COMMON /CAPTURW/ KSTORE, NCAPTURE")
    stubs.append("      INTEGER I, J, IPVT, IGRP, IROW, ICOL")
    stubs.append("      DO 100 IPVT = 1, 6")
    stubs.append("        DO 90 IGRP = 1, 6")
    stubs.append("          DO 80 I = 1, 3")
    stubs.append("            DO 80 J = 1, 3")
    stubs.append("              IROW = (IPVT-1)*3 + I")
    stubs.append("              ICOL = (IGRP-1)*3 + J")
    stubs.append("              KFULL(IROW, ICOL) = KSTORE(I, J, IPVT,")
    stubs.append("     1          IGRP)")
    stubs.append("   80       CONTINUE")
    stubs.append("   90     CONTINUE")
    stubs.append("  100 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SAXB (A, B, C)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      REAL A(3), B(3), C(3)")
    stubs.append("      C(1) = A(2)*B(3) - A(3)*B(2)")
    stubs.append("      C(2) = A(3)*B(1) - A(1)*B(3)")
    stubs.append("      C(3) = A(1)*B(2) - A(2)*B(1)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      FUNCTION SADOTB (A, B)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      REAL SADOTB")
    stubs.append("      REAL A(3), B(3)")
    stubs.append("      SADOTB = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE KPLTST (P1, P2, P3, P4)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      REAL P1(3), P2(3), P3(3), P4(3)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (ICSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      T(1)=1.0D0; T(2)=0.0D0; T(3)=0.0D0")
    stubs.append("      T(4)=0.0D0; T(5)=1.0D0; T(6)=0.0D0")
    stubs.append("      T(7)=0.0D0; T(8)=0.0D0; T(9)=1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG, IPARM(2)")
    stubs.append("      INTEGER OUT, SYSBUF")
    stubs.append("      LOGICAL NOGO")
    stubs.append("      COMMON /SYSTEM/ SYSBUF, OUT, NOGO")
    stubs.append("      WRITE(OUT,*) 'MESAGE code=', ICODE, ' msg=', IMSG")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE WRITE (IFILE, BUF, NWORDS, IFLG)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IFILE, NWORDS, IFLG")
    stubs.append("      REAL BUF(1)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE HMAT (NECPT)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER NECPT(4)")
    stubs.append("      REAL MATBUF(7)")
    stubs.append("      COMMON /HMTOUT/ MATBUF")
    stubs.append("      MATBUF(1) = 50.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _parse_B8(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 18, 18)}


def _compare_B8(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 18)


# =====================================================================
# B9 (KELBOW) - 12x12 curved-bar / elbow element stiffness matrix.
# Kernel kelbow.f (calls GMMATD, INVERD, TRANSD, MAT, SMA1B) + real
# gmmatd.f + real inverd.f.  KELBOW does NOT call GMMATS.
#
# /SMA1ET/ layout (kelbow.f comments, field-named EQUIVALENCE):
#   1 IELID | 2,3 ISILNO | 4-6 SMALLV | 7 ICSSV | 8 IMATID | 9 A
#   10 I1 | 11 I2 | 12 FJ | 13 NSM | 14 FE | 15-22 R1..T4
#   23 K1 | 24 K2 | 25 C | 26 KX | 27 KY | 28 KZ | 29 R | 30 BETAR
#   31 MCSIDA | 32-34 GPA | 35 MCSIDB | 36-38 GPB | 39 TEMPEL
#   INTEGER words: 1, 2, 3, 7, 8, 31, 35.
#
# KELBOW is called twice (once per pivot SIL).  Per pivot it emits two 6x6
# sub-blocks via SMA1B: first the self-coupling K(NPVT,NPVT), then the
# coupling K(NPVT,J).  Unlike KBAR, both calls carry INDEX = pivot SIL
# (not the other grid), so we capture by call-count per pivot into
# KSTORE(6,6,2,2) and assemble via GET_KMATRIX_ELBOW (mirrors the frozen
# b9_kelbow stubs).  A/I1/I2/FJ are null in the inputs -> derived from the
# circular-pipe geometry outer_radius/wall_thickness exactly as the frozen
# test_kelbow.f driver does.
# =====================================================================
def _gen_B9(rec: dict[str, Any]) -> dict[str, str]:
    import math
    ga = rec["gridA"]
    gb = rec["gridB"]
    refv = rec["ref_vec"]
    R = rec["R"]
    ang = rec["angle_deg"]
    E = rec["E"]
    nu = rec["nu"]
    ro = rec["outer_radius"]
    tw = rec["wall_thickness"]
    ri = ro - tw
    A = math.pi * (ro**2 - ri**2)
    I1 = math.pi / 4.0 * (ro**4 - ri**4)
    I2 = I1
    FJ = math.pi / 2.0 * (ro**4 - ri**4)
    lines = []
    lines.append("      PROGRAM HOLDB9")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL ECPT(100)")
    lines.append("      INTEGER NECPT(100)")
    lines.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    lines.append("      COMMON /SMA1ET/ ECPT")
    lines.append("      INTEGER IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC")
    lines.append("      INTEGER NROWSC, TNROWS, JMAX, NLINKS, LINK(10)")
    lines.append("      INTEGER IDETCK, DODET, NOGO_CL")
    lines.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, LEFT, FROWIC,")
    lines.append("     1  LROWIC, NROWSC, TNROWS, JMAX, NLINKS, LINK, IDETCK,")
    lines.append("     2  DODET, NOGO_CL")
    lines.append("      LOGICAL HEAT")
    lines.append("      COMMON /SMA1HT/ HEAT")
    lines.append("      INTEGER MATIDC, MATFLG")
    lines.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    lines.append("      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS,")
    lines.append("     1  SINTH, COSTH")
    lines.append("      REAL E, GG, NU, RHO, ALPHA, TSUBO, GSUBE")
    lines.append("      REAL SIGT, SIGC, SIGS")
    lines.append("      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUBO, GSUBE,")
    lines.append("     1  SIGT, SIGC, SIGS")
    lines.append("      REAL FK")
    lines.append("      COMMON /HMTOUT/ FK")
    lines.append("      INTEGER IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT")
    lines.append("      INTEGER IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG")
    lines.append("      INTEGER IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW")
    lines.append("      INTEGER CLSNRW, CLSRW, NEOR, EOR, MCBKGG(7), MCB4GG(7)")
    lines.append("      COMMON /SMA1IO/ IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT,")
    lines.append("     1  IGECPT, IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG,")
    lines.append("     2  IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW,")
    lines.append("     3  CLSNRW, CLSRW, NEOR, EOR, MCBKGG, MCB4GG")
    lines.append("      INTEGER ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT")
    lines.append("      INTEGER I6X6K, N6X6K, I6X64, N6X64")
    lines.append("      COMMON /SMA1BK/ ICSTM, NCSTM, IGPCT, NGPCT, IPOINT,")
    lines.append("     1  NPOINT, I6X6K, N6X6K, I6X64, N6X64")
    lines.append("      DOUBLE PRECISION KE(144), KEP(144), DELA(6), DELB(6)")
    lines.append("      COMMON /SMA1DP/ KE, KEP, DELA, DELB")
    lines.append("      INTEGER SYSBUF, NOUT")
    lines.append("      COMMON /SYSTEM/ SYSBUF, NOUT")
    lines.append("      INTEGER MACH")
    lines.append("      COMMON /MACHIN/ MACH")
    lines.append("      DOUBLE PRECISION Z(1000)")
    lines.append("      COMMON /ZZZZZZ/ Z")
    lines.append("      DOUBLE PRECISION KFULL(12,12)")
    lines.append("      INTEGER I, J")
    lines.append("      MACH = 1")
    lines.append("      NOUT = 6")
    lines.append("      SYSBUF = 1")
    lines.append("      HEAT = .FALSE.")
    lines.append("      IOPT4 = 0")
    lines.append("      K4GGSW = 0")
    lines.append("      NOGO_CL = 0")
    lines.append("      IFCSTM=0; IFMPT=0; IFDIT=0; IDUM1=0; IFECPT=0")
    lines.append("      IGECPT=0; IFGPCT=0; IGGPCT=0; IFGEI=0; IGGEI=0")
    lines.append("      IFKGG=11; IGKGG=0; IF4GG=13; IG4GG=0")
    lines.append("      IFGPST=0; IGGPST=0; INRW=0; OUTRW=0")
    lines.append("      CLSNRW=0; CLSRW=0; NEOR=0; EOR=0")
    lines.append("      DO 5 I = 1, 7")
    lines.append("        MCBKGG(I) = 0")
    lines.append("        MCB4GG(I) = 0")
    lines.append("    5 CONTINUE")
    lines.append("      ICSTM=0; NCSTM=0; IGPCT=0; NGPCT=0; IPOINT=0")
    lines.append("      NPOINT=0; I6X6K=0; N6X6K=0; I6X64=0; N6X64=0")
    lines.append("      DO 10 I = 1, 1000")
    lines.append("        Z(I) = 0.0D0")
    lines.append("   10 CONTINUE")
    lines.append("      DO 20 I = 1, 12")
    lines.append("        DO 15 J = 1, 12")
    lines.append("          KFULL(I,J) = 0.0D0")
    lines.append("   15   CONTINUE")
    lines.append("   20 CONTINUE")
    lines.append("      CALL SETECP")
    lines.append("      NPVT = NECPT(2)")
    lines.append("      CALL KELBOW")
    lines.append("      NPVT = NECPT(3)")
    lines.append("      CALL KELBOW")
    lines.append("      CALL GET_KMATRIX_ELBOW(KFULL)")
    lines.append("      WRITE(6,'(A1)') 'S'")
    lines.append("      DO 200 I = 1, 12")
    lines.append("        WRITE(6,920) (KFULL(I,J), J = 1, 12)")
    lines.append("  200 CONTINUE")
    lines.append("  920 FORMAT (12(1X,E16.9))")
    lines.append("      STOP")
    lines.append("      END")
    # ---- SETECP: parameterised kelbow ECPT (all -basic- coord, CSID=0) ----
    lines.append("      SUBROUTINE SETECP")
    lines.append("      IMPLICIT NONE")
    lines.append("      REAL ECPT(100)")
    lines.append("      INTEGER NECPT(100)")
    lines.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    lines.append("      COMMON /SMA1ET/ ECPT")
    lines.append("      INTEGER I")
    lines.append("      DO 10 I = 1, 100")
    lines.append("        ECPT(I) = 0.0")
    lines.append("   10 CONTINUE")
    lines.append("      NECPT(1) = 101")
    lines.append("      NECPT(2) = 1")
    lines.append("      NECPT(3) = 2")
    lines.append(f"      ECPT(4) = {_ffloat(refv[0])}")
    lines.append(f"      ECPT(5) = {_ffloat(refv[1])}")
    lines.append(f"      ECPT(6) = {_ffloat(refv[2])}")
    lines.append("      NECPT(7) = 0")
    lines.append("      NECPT(8) = 1")
    lines.append(f"      ECPT(9) = {_ffloat(A)}")
    lines.append(f"      ECPT(10) = {_ffloat(I1)}")
    lines.append(f"      ECPT(11) = {_ffloat(I2)}")
    lines.append(f"      ECPT(12) = {_ffloat(FJ)}")
    lines.append("      ECPT(23) = 0.9")
    lines.append("      ECPT(24) = 0.9")
    lines.append("      ECPT(25) = 1.0")
    lines.append("      ECPT(26) = 1.0")
    lines.append("      ECPT(27) = 1.0")
    lines.append("      ECPT(28) = 1.0")
    lines.append(f"      ECPT(29) = {_ffloat(R)}")
    lines.append(f"      ECPT(30) = {_ffloat(ang)}")
    lines.append("      NECPT(31) = 0")
    lines.append(f"      ECPT(32) = {_ffloat(ga[0])}")
    lines.append(f"      ECPT(33) = {_ffloat(ga[1])}")
    lines.append(f"      ECPT(34) = {_ffloat(ga[2])}")
    lines.append("      NECPT(35) = 0")
    lines.append(f"      ECPT(36) = {_ffloat(gb[0])}")
    lines.append(f"      ECPT(37) = {_ffloat(gb[1])}")
    lines.append(f"      ECPT(38) = {_ffloat(gb[2])}")
    lines.append("      ECPT(39) = 20.0")
    lines.append("      RETURN")
    lines.append("      END")
    # ---- stubs: CAPTURE block, MAT (param E/G/NU), SMA1B (call-count per
    #              pivot into KSTORE(6,6,2,2)), GET_KMATRIX_ELBOW, TRANSD,
    #              MESAGE, HMAT.  KELBOW never calls GMMATS.
    stubs = []
    stubs.append("      BLOCK DATA CAPTUREDATA")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KSTORE(6,6,2,2)")
    stubs.append("      INTEGER NCALLS, CPIVOT")
    stubs.append("      COMMON /CAPTURE/ KSTORE, NCALLS, CPIVOT")
    stubs.append("      DATA NCALLS /0/")
    stubs.append("      DATA CPIVOT /0/")
    stubs.append("      DATA KSTORE /144*0.0D0/")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      INTEGER MATIDC, MATFLG")
    stubs.append("      REAL ELTEMP, STRESS, SINTH, COSTH")
    stubs.append("      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS,")
    stubs.append("     1  SINTH, COSTH")
    stubs.append("      REAL E, GG, NU, RHO, ALPHA, TSUBO, GSUBE, SIGT,")
    stubs.append("     1  SIGC, SIGS")
    stubs.append("      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUBO, GSUBE,")
    stubs.append("     1  SIGT, SIGC, SIGS")
    stubs.append(f"      E = {_ffloat(E)}")
    stubs.append(f"      NU = {_ffloat(nu)}")
    stubs.append(f"      GG = {_ffloat(E / (2.0*(1.0+nu)))}")
    stubs.append("      RHO = 7800.0")
    stubs.append("      ALPHA = 12.0E-6")
    stubs.append("      TSUBO = 20.0")
    stubs.append("      GSUBE = 0.0")
    stubs.append("      SIGT = 0.0")
    stubs.append("      SIGC = 0.0")
    stubs.append("      SIGS = 0.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE SMA1B (KIJ, INDEX, ICODE, IFFILE,")
    stubs.append("     1  FACTOR)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KIJ(6,6), FACTOR")
    stubs.append("      INTEGER INDEX, ICODE, IFFILE")
    stubs.append("      DOUBLE PRECISION KSTORE(6,6,2,2)")
    stubs.append("      INTEGER NCALLS, CPIVOT")
    stubs.append("      COMMON /CAPTURE/ KSTORE, NCALLS, CPIVOT")
    stubs.append("      REAL ECPT(100)")
    stubs.append("      INTEGER NECPT(100)")
    stubs.append("      EQUIVALENCE (NECPT(1), ECPT(1))")
    stubs.append("      COMMON /SMA1ET/ ECPT")
    stubs.append("      INTEGER IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC")
    stubs.append("      INTEGER NROWSC, TNROWS, JMAX, NLINKS, LINK(10)")
    stubs.append("      INTEGER IDETCK, DODET, NOGO")
    stubs.append("      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, LEFT, FROWIC,")
    stubs.append("     1  LROWIC, NROWSC, TNROWS, JMAX, NLINKS, LINK, IDETCK,")
    stubs.append("     2  DODET, NOGO")
    stubs.append("      INTEGER I, J, IPVT, IBLK")
    stubs.append("      IF (FACTOR .NE. 0.0D0) RETURN")
    stubs.append("      IPVT = 0")
    stubs.append("      IF (NPVT .EQ. NECPT(2)) IPVT = 1")
    stubs.append("      IF (NPVT .EQ. NECPT(3)) IPVT = 2")
    stubs.append("      IF (IPVT .EQ. 0) RETURN")
    stubs.append("      IF (IPVT .NE. CPIVOT) THEN")
    stubs.append("        CPIVOT = IPVT")
    stubs.append("        NCALLS = 0")
    stubs.append("      END IF")
    stubs.append("      NCALLS = NCALLS + 1")
    stubs.append("      IBLK = NCALLS")
    stubs.append("      IF (IBLK .GT. 2) IBLK = 2")
    stubs.append("      DO 50 I = 1, 6")
    stubs.append("        DO 50 J = 1, 6")
    stubs.append("          KSTORE(I, J, IPVT, IBLK) = KIJ(I, J)")
    stubs.append("   50 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE GET_KMATRIX_ELBOW (KFULL)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      DOUBLE PRECISION KFULL(12,12)")
    stubs.append("      DOUBLE PRECISION KSTORE(6,6,2,2)")
    stubs.append("      INTEGER NCALLS, CPIVOT")
    stubs.append("      COMMON /CAPTURE/ KSTORE, NCALLS, CPIVOT")
    stubs.append("      INTEGER I, J")
    stubs.append("      DO 100 I = 1, 6")
    stubs.append("        DO 100 J = 1, 6")
    stubs.append("          KFULL(I, J) = KSTORE(I, J, 1, 1)")
    stubs.append("          KFULL(I, J+6) = KSTORE(I, J, 1, 2)")
    stubs.append("          KFULL(I+6, J) = KSTORE(I, J, 2, 2)")
    stubs.append("          KFULL(I+6, J+6) = KSTORE(I, J, 2, 1)")
    stubs.append("  100 CONTINUE")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE TRANSD (ICSID, T)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICSID")
    stubs.append("      DOUBLE PRECISION T(9)")
    stubs.append("      T(1)=1.0D0; T(2)=0.0D0; T(3)=0.0D0")
    stubs.append("      T(4)=0.0D0; T(5)=1.0D0; T(6)=0.0D0")
    stubs.append("      T(7)=0.0D0; T(8)=0.0D0; T(9)=1.0D0")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER ICODE, IMSG, IPARM")
    stubs.append("      WRITE(6,900) ICODE, IMSG")
    stubs.append("  900 FORMAT(' MESAGE called: code=', I6, ' msg=', I6)")
    stubs.append("      RETURN")
    stubs.append("      END")
    stubs.append("      SUBROUTINE HMAT (IELID)")
    stubs.append("      IMPLICIT NONE")
    stubs.append("      INTEGER IELID")
    stubs.append("      REAL FK")
    stubs.append("      COMMON /HMTOUT/ FK")
    stubs.append("      FK = 50.0")
    stubs.append("      RETURN")
    stubs.append("      END")
    return {"driver.f": "\n".join(lines), "stubs.f": "\n".join(stubs)}


def _parse_B9(stdout: str) -> dict[str, Any]:
    return {"matrix": _parse_matrix(stdout, 12, 12)}


def _compare_B9(out: dict[str, Any], ref: dict[str, Any]) -> list[str]:
    return _compare_matrix(out, ref, 12)


GENERATORS: dict[str, dict[str, Any]] = {
    "B10": {
        "kernel_files": ("inverd.f",),
        "gen": _gen_b10,
        "parse": _parse_b10,
        "compare_ref": _compare_b10,
        "ref_path": Path("source_fortran/test_drivers/b9_inverd/reference_values.json"),
    },
    "B12": {
        "kernel_files": ("trd1c.f", "step.f", "form1.f", "form2.f", "inverd.f"),
        "gen": _gen_B12,
        "parse": _parse_B12,
        "compare_ref": _compare_B12,
        "ref_path": Path("source_fortran/test_drivers/b12_trd1c_native/reference_timeseries.json"),
    },
    "B13": {
        "kernel_files": ("srod1.f",),
        "gen": _gen_b13,
        "parse": _parse_b13,
        "compare_ref": _compare_b13,
        "ref_path": Path("source_fortran/test_drivers/b12_srod1/reference_values.json"),
    },
    "B11": {
        "kernel_files": ("qriter.f",),
        "gen": _gen_b11,
        "parse": _parse_b11,
        "compare_ref": _compare_b11,
        "ref_path": Path("source_fortran/test_drivers/b10_qriter/reference_eigenvalues.json"),
    },
    "B2": {
        "kernel_files": ("kbar.f", "gmmatd.f"),
        "gen": _gen_B2,
        "parse": _parse_B2,
        "compare_ref": _compare_B2,
        "ref_path": Path("source_fortran/test_drivers/b2_kbar/reference_matrix.json"),
    },
    "B1": {
        "kernel_files": ("krod.f", "gmmatd.f"),
        "gen": _gen_b1,
        "parse": _parse_b1,
        "compare_ref": _compare_b1,
        "ref_path": Path("source_fortran/test_drivers/b1_krod/reference_matrix.json"),
    },
    "B3": {
        "kernel_files": ("ktrmem.f", "gmmatd.f"),
        "gen": _gen_b3,
        "parse": _parse_b3,
        "compare_ref": _compare_b3,
        "ref_path": Path("source_fortran/test_drivers/b3_ktrmem/reference_matrix.json"),
    },
    "B4": {
        "kernel_files": ("kqdmem.f", "ktrmem.f", "gmmatd.f"),
        "gen": _gen_B4,
        "parse": _parse_B4,
        "compare_ref": _compare_B4,
        "ref_path": Path("source_fortran/test_drivers/b4_ktriqd/reference_matrix.json"),
    },
    "B7": {
        "kernel_files": ("ktetra.f", "inverd.f", "gmmatd.f"),
        "gen": _gen_B7,
        "parse": _parse_B7,
        "compare_ref": _compare_B7,
        "load_ref": _load_b7_reference,
        "ref_path": Path("source_fortran/test_drivers/reference_output.txt"),
    },
    "B5": {
        "kernel_files": ("ktrplt.f", "ktrbsc.f", "inverd.f", "gmmatd.f"),
        "gen": _gen_B5,
        "parse": _parse_B5,
        "compare_ref": _compare_B5,
        "load_ref": _load_b5_reference,
        "ref_path": Path("source_fortran/test_drivers/b5_ktrplt/reference_output.txt"),
    },
    "B6": {
        "kernel_files": ("kqdplt.f", "ktrbsc.f", "inverd.f", "gmmatd.f"),
        "gen": _gen_B6,
        "parse": _parse_B6,
        "compare_ref": _compare_B6,
        "load_ref": _load_b6_reference,
        "ref_path": Path("source_fortran/test_drivers/b6_kqdplt/reference_output.txt"),
    },
    "B8": {
        "kernel_files": ("ksolid.f", "ktetra.f", "inverd.f", "gmmatd.f"),
        "gen": _gen_B8,
        "parse": _parse_B8,
        "compare_ref": _compare_B8,
        "ref_path": Path("source_fortran/test_drivers/b8_ksolid/reference_matrix_wedge.json"),
    },
    "B9": {
        "kernel_files": ("kelbow.f", "gmmatd.f", "inverd.f"),
        "gen": _gen_B9,
        "parse": _parse_B9,
        "compare_ref": _compare_B9,
        "ref_path": Path("source_fortran/test_drivers/b9_kelbow/reference_matrix.json"),
    },
}


def _build(work: Path, bench: str) -> Path:
    spec = GENERATORS[bench]
    exe = work / "holdoracle.exe"
    files = [str(work / "driver.f")]
    if (work / "stubs.f").is_file():
        files.append(str(work / "stubs.f"))
    files += [str(work / f) for f in spec["kernel_files"]]
    cmd = [GFORTRAN, "-o", str(exe), *files, *FCFLAGS]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"compile failed [{bench}]:\n{proc.stdout}\n{proc.stderr}"
        )
    return exe


def process(rec: dict[str, Any], bench: str) -> dict[str, Any]:
    k = rec["_variant"]
    spec = GENERATORS[bench]
    work = BUILD_DIR / f"{bench}_{k}"
    work.mkdir(parents=True, exist_ok=True)
    files = spec["gen"](rec)  # dict filename->source
    srcs = {}
    for name, content in files.items():
        (work / name).write_bytes(content.encode("utf-8"))
        srcs[name] = content
    kernel_hashes = {}
    for f in spec["kernel_files"]:
        kernel_hashes[f] = sha256_file(MIS / f)
        shutil.copy(MIS / f, work / f)
    exe = _build(work, bench)
    proc = subprocess.run([str(exe)], capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        raise RuntimeError(f"run failed [{bench} {k}]:\n{proc.stderr}")
    parsed = spec["parse"](proc.stdout)
    result: dict[str, Any] = {
        "_schema_version": "jss-holdout-oracle-1",
        "_benchmark": bench,
        "_variant": k,
        "_generated_utc": datetime.now(timezone.utc).isoformat(),
        "_method": "parameterised Fortran driver + unmodified NASTRAN kernel",
        "_toolchain": {"gfortran": GFORTRAN, "flags": FCFLAGS},
        "_kernel_sha256": kernel_hashes,
        "_driver_source": {
            "driver.f": srcs.get("driver.f", ""),
            "stubs.f": srcs.get("stubs.f", ""),
        },
        "_result": parsed,
    }
    return result


def freeze_manifest(records: list[dict[str, Any]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    entries = {}
    for r in records:
        path = OUT_DIR / f"{r['_benchmark']}_{r['_variant']}.json"
        text = json.dumps(r, ensure_ascii=False, indent=2) + "\n"
        path.write_bytes(text.encode("utf-8"))
        entries[str(path.relative_to(ROOT)).replace("\\", "/")] = sha256_bytes(text.encode("utf-8"))
    manifest = {
        "schema_version": "jss-holdout-oracles-manifest-1",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "n_files": len(records),
        "files": entries,
    }
    mtext = json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    (ROOT / "holdout_oracle_manifest.json").write_bytes(mtext.encode("utf-8"))
    print(f"oracle manifest: {len(entries)} files -> holdout_oracle_manifest.json")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("benches", nargs="*", help="bench ids to build (default: all implemented)")
    args = ap.parse_args()
    targets = args.benches or list(GENERATORS.keys())
    records = []
    for bid in targets:
        if bid not in GENERATORS:
            print(f"[skip] {bid}: generator not implemented yet (P1 pending)")
            continue
        for k in ("h0", "h1", "h2", "h3"):
            rec_path = IN_DIR / f"{bid}_{k}.json"
            if not rec_path.is_file():
                print(f"[skip] {bid}_{k}: no input file")
                continue
            rec = json.loads(rec_path.read_text(encoding="utf-8"))
            # ---- h0 regression gate ----
            if k == "h0":
                spec = GENERATORS[bid]
                loader = spec.get("load_ref", json.loads)
                reftxt = (ROOT / spec["ref_path"]).read_text(encoding="utf-8")
                ref = loader(reftxt) if callable(loader) else loader
                r0 = process(rec, bid)
                probs = spec["compare_ref"](r0["_result"], ref)
                print(f"[h0-regression] {bid}_h0 issues={probs}")
                if probs:
                    print(f"  !! ORACLE DOES NOT MATCH FROZEN REFERENCE: {probs}")
                    shutil.rmtree(BUILD_DIR / f"{bid}_h0", ignore_errors=True)
                    raise SystemExit(1)
                records.append(r0)
                continue
            r = process(rec, bid)
            records.append(r)
            print(f"[ok] {bid}_{k}")
    freeze_manifest(records)
    # clean build dirs (provenance kept in oracle json _driver_source)
    for p in sorted(BUILD_DIR.glob("*")):
        shutil.rmtree(p, ignore_errors=True)


if __name__ == "__main__":
    main()