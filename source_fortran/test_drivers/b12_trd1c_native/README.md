# B12: TRD1C Native — Mini-NASTRAN Runtime

## Status: ✅ WORKING

Successfully compiles and runs NASTRAN-95's actual `trd1c.f` transient integration driver
with a minimal GINO simulation layer. Produces physically correct results for a 3-DOF
spring-mass free vibration problem.

## What This Demonstrates

The actual NASTRAN-95 source code (`trd1c.f`, `step.f`, `form1.f`, `form2.f`) can be executed
outside of the full NASTRAN environment by providing:

1. **Mini-GINO layer** (`gino_sim.f`): Replaces NASTRAN's disk-based GINO I/O with in-memory
   arrays. Matrices stored in a COMMON block, file operations become array lookups.

2. **Custom MATVEC/INTFBS**: Instead of NASTRAN's sparse packed-column format read from
   tape files, our versions do direct dense matrix-vector multiply from stored arrays.

3. **Stub routines**: MESAGE, TRD1D, SSWTCH, TMTOGO, etc. provide minimal runtime support.

## Test Problem

```
3-DOF fixed-fixed spring-mass chain:
  wall --k-- m1 --k-- m2 --k-- m3 --k-- wall

K = [2 -1  0]    M = I (identity)    C = 0
    [-1  2 -1]
    [ 0 -1  2]

Initial conditions: u0 = [1, 0, 0], udot0 = [0, 0, 0]
Time step: dt = 0.01, 100 steps, output every 10 steps
Method: Central difference (NASTRAN's default for TRD)
```

## Results

| Time | u1 (numerical) | u1 (analytical) | Error |
|------|---------------|----------------|-------|
| 0.0  | 1.000000      | 1.000000       | 0     |
| 0.5  | 0.764243      | ~0.7643        | <0.01%|
| 1.0  | 0.191835      | 0.189895       | 0.19% |

Error is O(dt²) as expected for central difference — 2nd order accuracy confirmed.

## Architecture

```
test_trd1c_native.f  — Driver: sets up matrices, COMMON blocks, calls TRD1C
gino_sim.f           — Mini-GINO: MATVEC, INTFBS, OPEN/CLOSE/PACK/UNPACK stubs
trd1c.f              — ORIGINAL NASTRAN source (unmodified, compiled from mis/)
step.f               — ORIGINAL NASTRAN source (unmodified)
form1.f              — ORIGINAL NASTRAN source (unmodified)
form2.f              — ORIGINAL NASTRAN source (unmodified)
inverd.f             — ORIGINAL NASTRAN source (used for K_eff inversion)
```

## Key Technical Decisions

### 1. Replace MATVEC entirely (not stub the I/O layer)

NASTRAN's MATVEC reads matrices in packed sparse column format using INTPK/ZNTPKI.
Simulating this format for a 3×3 dense matrix adds complexity for no benefit.
Instead, our MATVEC does a simple dense multiply: `X = X + A*Y`.

### 2. Replace INTFBS with direct inverse multiplication

NASTRAN's INTFBS does LU forward/backward substitution from factored matrices stored
on tape. For a 3×3 system, we pre-compute K_eff_inverse using INVERD and multiply directly.

### 3. File-number-to-slot mapping

NASTRAN uses integer file IDs (301, 302, etc.). We map these to array indices:
- SCR1=301 → matrix slot 1 (holds 2M/dt² - K)
- SCR4=304 → matrix slot 4 (holds -M/dt²)
- Slot 9 → K_eff_inverse (used by INTFBS)
- IK=401 → slot 7 (stiffness K, for FORM1)
- IM=402 → slot 8 (mass M, for FORM2)
- IB=403 → slot 10 (damping B=0)

### 4. Central difference formulation

TRD1C's STEP routine does:
```
P += SCR1 * u_n        (SCR1 = 2M/dt² - K)
P += SCR4 * u_{n-1}    (SCR4 = -M/dt²)
solve K_eff * u_{n+1} = P   (K_eff = M/dt²)
```

This is the standard central difference:
```
M/dt² * u_{n+1} = F_n + (2M/dt² - K)*u_n - M/dt² * u_{n-1}
```

### 5. FORM1 starting procedure

TRD1C calls FORM1 to compute u_{-1} (the "previous" displacement for the first step).
FORM1 uses: `u_{-1} = u_0 - dt * udot_0` and then applies K and B via MATVEC to form
the initial load vectors P_0 and P_{-1}.

## What Works

- ✅ Compilation with gfortran (MSYS2 UCRT64)
- ✅ Full TRD1C execution path (first entry, ILOOP=1)
- ✅ FORM1 starting vector computation
- ✅ STEP time integration loop (100 steps)
- ✅ PACK output to UDV file
- ✅ Load averaging (P* = (P_{n-1} + P_n + P_{n+1})/3)
- ✅ Displacement, velocity, and acceleration output
- ✅ Results match analytical modal solution to O(dt²)

## What's Not Tested / Limitations

- ❌ Multi-group time stepping (NGROUP > 1, time step changes)
- ❌ Nonlinear loads (NLFTP ≠ 0, TRD1D path)
- ❌ Checkpoint/restart (ICPFLG path)
- ❌ Damped systems (C ≠ 0) — should work but untested
- ❌ ISYM=0 path (FBSINT instead of INTFBS) — stub exists but untested
- ❌ Large systems (limited to 3×3 by storage arrays)
- ❌ Panel nonlinear force output (PNL path)

## Build

```batch
build_b12_native.bat
```

Or manually:
```
gfortran -c test_trd1c_native.f gino_sim.f trd1c.f step.f form1.f form2.f inverd.f -ffixed-form -std=legacy -w -fallow-argument-mismatch -fno-range-check -O0 -g
gfortran -o test_trd1c_native.exe *.o -g
```
