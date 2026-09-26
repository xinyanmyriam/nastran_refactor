# NASTRAN-95 Benchmark Source Files

Legacy Fortran subroutines extracted from NASA NASTRAN-95 (`mis/` directory) for use as benchmarks in the legacy code refactoring paper. These files represent real-world structural analysis routines spanning element stiffness generation, numerical solvers, and stress recovery.

**Source:** `NASTRAN-95/mis/` — Machine Independent Subroutines  
**Language:** FORTRAN 77 (fixed-form, 6-column indent)

---

## Benchmark Index

| ID  | File        | Element/Function                  | Lines | Matrix/DOF       |
|-----|-------------|-----------------------------------|-------|------------------|
| B1  | krod.f      | Rod element stiffness             | ~50   | 2×2 (axial)      |
| B2  | kbar.f      | Beam element stiffness            | ~250  | 12×12 (6 DOF/node)|
| B3  | ktrmem.f    | Triangular membrane stiffness     | ~200  | 6×6 (2 DOF/node) |
| B4  | ktriqd.f    | Quadrilateral membrane (driver)   | ~150  | 8×8 (2 DOF/node) |
| B4+ | kqdmem.f    | Quadrilateral membrane (kernel)   | 288   | 8×8 (2 DOF/node) |
| B5  | ktrplt.f    | Triangular plate bending          | 566   | 9×9 (3 DOF/node) |
| B6  | kqdplt.f    | Quadrilateral plate bending       | 342   | 12×12 (3 DOF/node)|
| B7  | ktetra.f    | Tetrahedron solid element         | 384   | 12×12 (3 DOF/node)|
| B8  | ksolid.f    | Hexahedron solid (wedge/hexa)     | 236   | 24×24 (3 DOF/node)|
| B9  | solve.f     | Linear equation solver (AX=B)     | 220   | N×N system       |
| B10 | invpwr.f    | Inverse power eigenvalue method   | 235   | N×N generalized  |
| B11 | trd1c.f     | Transient integration (Newmark)   | 658   | N time steps     |
| B12 | srod1.f     | Rod stress recovery               | 118   | per element      |

---

## Detailed File Descriptions

### B5 — ktrplt.f (Triangular Plate Bending Element, 566 lines)

```fortran
SUBROUTINE KTRPLT
```

**Purpose:** Generates 3 sets of 6×6 stiffness matrices with respect to one pivot point for a triangular plate bending element (TRIA1/TRIA2). Based on Kirchhoff thin plate theory.

**Key parameters (via ECPT COMMON):**
- Grid points A, B, C (3 nodes)
- Material properties (bending rigidity D matrix)
- Plate thickness, moment of inertia

**Calls:** KTRBSC (basic bending triangle), TRANSD (3×3 transformations), INVERD (matrix inversion), SMA1B (assembly insertion), GMMATD (general matrix multiply)

**Reference:** FMMS-55, November 1, 1967

---

### B6 — kqdplt.f (Quadrilateral Plate Bending Element, 342 lines)

```fortran
SUBROUTINE KQDPLT
```

**Purpose:** Generates 4 sets of 6×6 stiffness matrices with respect to one pivot point for a quadrilateral plate bending element (QUAD1/QUAD2). Decomposes the quadrilateral into triangles and assembles.

**Key parameters (via ECPT COMMON):**
- Grid points A, B, C, D (4 nodes)
- Material properties (G matrix), plate thickness
- THETA (material angle)

**Calls:** KTRBSC, TRANSD, SMA1B, GMMATD, MESAGE

**References:** FMMS-44 (July 18, 1967), FMMS-48 (August 1, 1967)

---

### B7 — ktetra.f (Tetrahedron Solid Element, 384 lines)

```fortran
SUBROUTINE KTETRA (IOPT, JTYPE)
```

**Purpose:** Element stiffness matrix generator for the tetrahedron solid element. Used as the building block for wedge (3 tetrahedra) and hexahedron (5 or 10 tetrahedra) elements.

**Key parameters:**
- `IOPT` — operation option
- `JTYPE` — element type: 0=tetra, 1=wedge, 2=hexa1, 3=hexa2
- 4 grid points with coordinates (ECPT 3–22)
- MAT1 material (Young's modulus, Poisson's ratio)

**ECPT layout:** Element ID, Material ID, 4 grid SIL points, coordinate system IDs, coordinates (X,Y,Z for each), element temperature

**DOF:** 3 translations per node → 12×12 stiffness matrix

---

### B8 — ksolid.f (Hexahedron/Wedge Solid Element Driver, 236 lines)

```fortran
SUBROUTINE KSOLID (ITYPE)
```

**Purpose:** Driver routine for 3D solid element stiffness generation. Decomposes wedge (6-node) and hexahedron (8-node) solids into tetrahedra and calls KTETRA for each sub-element.

**Key parameters:**
- `ITYPE` — 1=wedge (3 tetras), 2=hexa (5 tetras), 3=hexa (10 tetras)
- ECPT varies by element type (23–43 words)

**DOF:** Up to 8 nodes × 3 DOF = 24×24 stiffness matrix

---

### B9 — solve.f (Linear Equation Solver, 220 lines)

```fortran
SUBROUTINE SOLVE
```

**Purpose:** DMAP driver to solve the matrix equation AX=B. Performs symmetric or unsymmetric decomposition followed by forward/backward substitution (FBS).

**DMAP call syntax:** `SOLVE A,B/X/SYM/SIGN/PREC/TYPE $`

**Key parameters:**
- `SYM` — 1=symmetric decomp, 0=auto-detect, -1=unsymmetric
- `SIGN` — 1 for AX=B, -1 for AX=-B
- `PREC` — precision for FBS pass
- `TYPE` — desired output matrix type

**Method:** LU decomposition with forward/backward substitution. Handles both real and complex matrices.

---

### B10 — invpwr.f (Inverse Power Eigenvalue Method, 235 lines)

```fortran
SUBROUTINE INVPWR
```

**Purpose:** Given a real symmetric stiffness matrix K and mass matrix M, solves for all eigenvalues and eigenvectors within a specified range using the inverse power (inverse iteration) method.

**Key parameters:**
- `FILEK(7)` — MCB for stiffness matrix K
- `FILEM(7)` — MCB for mass matrix M
- `FILELM(7)` — MCB for output eigenvalues
- `FILEVC(7)` — MCB for output eigenvectors
- `LAMMIN/LAMMAX` — eigenvalue search range
- `EPS` — convergence criterion
- `NOEST` — estimated number of eigenvalues in range

**Method:** Inverse iteration with Sturm sequence check. Uses 7 scratch files for intermediate LU factors and orthogonalization vectors.

---

### B11 — trd1c.f (Transient Time Integration — Newmark Method, 658 lines)

```fortran
SUBROUTINE TRD1C(IC, PD, NGROUP, NLFTP, UDV, ILOOP, SCR1, DIT, NLFT, NOUE, MODAL, PNL, ISKIP)
```

**Purpose:** Steps the time integration procedure for transient response analysis (direct or modal). Implements the Newmark-beta method for second-order ODEs (M·ü + C·u̇ + K·u = P(t)).

**Key parameters:**
- `IC` — initial condition flag
- `PD` — load file
- `NGROUP` — number of groups
- `NLFTP` — nonlinear force type
- `UDV` — displacement/velocity/acceleration file
- `ILOOP` — current time step loop counter
- `MODAL` — modal vs. direct flag
- `DELTAT` — time step size (from COMMON /TRDD1/)
- `NROW` — system size (from IK(3))

**COMMON blocks:** /TRDXX/, /TRDD1/, /PACKX/, /UNPAKX/, /SYSTEM/

---

### B12 — srod1.f (Rod Stress Recovery, 118 lines)

```fortran
SUBROUTINE SROD1
```

**Purpose:** Phase I of stress data recovery for the CROD element. Computes axial stress, torsional stress, and safety margins from element displacements.

**Key parameters (via COMMON /SDR2X5/):**
- `ECPT(17)` — element connectivity/property table
- `IELID` — element ID
- `ISILNO(2)` — SIL numbers for the 2 grid points
- `AREA` — cross-sectional area
- `FJOVRC` — J/C for torsion
- `SIGMAT/SIGMAC/SIGMAS` — allowable stresses (tension/compression/shear)

**Output:** Axial force, torque, axial stress, torsional stress, margins of safety (tension, compression)

---

### B4+ — kqdmem.f (Quadrilateral Membrane Kernel, 288 lines)

```fortran
SUBROUTINE KQDMEM
```

**Purpose:** Quadrilateral membrane element stiffness matrix generation. Decomposes the quad into 4 overlapping triangles and calls KTRMEM for each, then assembles using area-weighted averaging.

**Key parameters (via ECPT COMMON):**
- Grid points A, B, C, D (4 nodes, in-plane DOF)
- Material ID, thickness T, non-structural mass
- Coordinate systems and grid coordinates

**Calls:** KTRMEM (triangular membrane kernel), SMA1B (assembly), MESAGE

---

## Notes

- All files use FORTRAN 77 fixed-form format (columns 1–6 for labels/continuation, 7–72 for code)
- Heavy use of COMMON blocks for inter-routine communication (no explicit argument passing for element data)
- ECPT = Element Connection and Property Table — the primary data structure for element information
- Double precision arithmetic used throughout for numerical stability
- Original NASA documentation references (FMMS series) cited in source comments
- B1–B4 source files are located in the parent `source_fortran/` directory (krod.f, kbar.f, ktrmem.f, ktriqd.f)
