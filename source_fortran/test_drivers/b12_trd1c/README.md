# B12: TRD1C Standalone Newmark-β Time Integration

## Why TRD1C.F Cannot Be Run Directly

NASTRAN-95's `trd1c.f` (658 lines) is the main transient analysis time-stepping driver.
It depends heavily on:

1. **GINO Tape I/O System** — All matrices (K, M, C, effective stiffness, decomposed L/L^T)
   are stored on "scratch files" (SCR1–SCR6) managed by NASTRAN's GINO file system.
   Operations like `OPEN`, `CLOSE`, `FREAD`, `WRITE`, `SKPREC`, `RDTRL`, `WRTTRL` 
   are GINO-specific and have no standard Fortran equivalent.

2. **MATVEC** — Matrix-vector product reading matrices from tape files.

3. **INTFBS/FBSINT** — Forward/backward substitution using pre-decomposed L/L^T 
   factors stored on tape.

4. **PACK/UNPACK** — Sparse vector packing/unpacking routines.

5. **Memory Manager (KORSZ, /ZZZZZZ/)** — NASTRAN's dynamic memory allocation through 
   a single large blank common array with pointer arithmetic.

6. **Checkpoint/Restart Logic** — Complex state save/restore for job continuation.

7. **FORM1/FORM2** — Pre-computation of effective matrices stored on scratch files.

Simulating all of these would require reimplementing a significant fraction of the 
NASTRAN runtime environment — effectively building a "mini-NASTRAN" just to run one 
subroutine.

## What This Standalone Implementation Does

Instead of attempting to simulate the full GINO infrastructure, this test extracts the
**mathematical algorithm** from `trd1c.f` and `step.f` and implements it directly with
explicit arrays:

### Algorithm (Newmark-β, Average Acceleration)

The Newmark-β method with β=1/4, γ=1/2 (average acceleration, unconditionally stable):

```
K_eff = K + (4/(dt²)) · M + (2/(dt)) · C

At each time step:
  RHS = F(t+dt) + M·(a0·u + a2·v + a3·a) + C·(a1·u + a4·v + a5·a)
  Solve: K_eff · u_{n+1} = RHS
  Update: a_{n+1} = a0·(u_{n+1} - u_n) - a2·v_n - a3·a_n
          v_{n+1} = v_n + a6·a_n + a7·a_{n+1}

where: a0 = 1/(β·dt²), a2 = 1/(β·dt), a3 = 1/(2β) - 1
       a6 = dt·(1-γ), a7 = γ·dt
```

### Correspondence to NASTRAN's TRD1C/STEP

| NASTRAN Operation | Standalone Equivalent |
|---|---|
| `FORM1` pre-computes effective matrices on SCR1, SCR4 | We form `K_eff` directly |
| `STEP → MATVEC(SCR1, U1)` | Matrix-vector: `K_eff_partial · u_n` |
| `STEP → MATVEC(SCR4, U0)` | Matrix-vector for `u_{n-1}` contribution |
| `STEP → INTFBS(P, U2)` | INVERD solve of `K_eff · u = RHS` |
| GINO scratch file I/O | Direct array operations |
| NASTRAN memory management | Fortran static arrays |

### NASTRAN Component Used

The linear system solve uses **INVERD.F** — the actual NASTRAN-95 Gauss-Jordan solver
with partial pivoting. This is the same routine that would ultimately be called (through
INTFBS) during the decomposition phase.

## Test Problem

- **System**: 3-DOF spring-mass chain, both ends fixed
- **M** = I (unit masses)
- **K** = [2 -1 0; -1 2 -1; 0 -1 2] (spring stiffness = 1)
- **C** = 0 (no damping)
- **Initial**: u₀ = [1, 0, 0], v₀ = [0, 0, 0]
- **Time step**: dt = 0.01, 100 steps (T_final = 1.0 s)

### Analytical Solution

Eigenvalues: λ_k = 2 - 2·cos(kπ/4), k = 1,2,3
- ω₁ = 0.7654 rad/s, ω₂ = 1.4142 rad/s, ω₃ = 1.8478 rad/s

Modal superposition: u(t) = Σ φ_k · q_k(0) · cos(ω_k · t)

## Results

- **Max L2 error**: 3.0×10⁻⁵ (PASS: within 10⁻⁴ tolerance)
- **Energy conservation**: drift < 10⁻¹² (machine precision)
- The average acceleration method exactly preserves energy for undamped systems

## Build & Run

```batch
build_b12_trd1c.bat
```

Or manually:
```
F:\msys64\ucrt64\bin\gfortran.exe -o test_newmark.exe test_newmark.f inverd.f -ffixed-form -std=legacy -w -fallow-argument-mismatch
test_newmark.exe
```

## Files

| File | Description |
|---|---|
| `test_newmark.f` | Standalone Newmark-β implementation |
| `inverd.f` | NASTRAN-95 linear equation solver (copied from mis/) |
| `build_b12_trd1c.bat` | Build and run script |
| `reference_output.txt` | Captured output for verification |
| `README.md` | This file |
