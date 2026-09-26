# NASTRAN-95 Source Fortran Files

This folder contains Fortran 77 source files extracted from NASA NASTRAN-95
for refactoring experiments.

## How to obtain source files

```bash
git clone https://github.com/nasa/NASTRAN-95.git
```

The relevant files are in the `mis/` directory.

## Target files for Phase 1 (CROD element)

The CROD (rod) element in NASTRAN-95 involves these key subroutines:
- Element stiffness matrix computation
- Stress recovery
- Element force computation

## NASTRAN-95 Repository Structure

```
NASTRAN-95/
├── mis/          # Core computational subroutines (Fortran 77)
├── mds/          # Module data sets
├── inp/          # Input test files
├── demoout/      # Demo output files (reference results)
├── um/           # User manual text files
├── bin/          # Binary/executable related
├── rf/           # Rigid format DMAP sequences
└── README.md
```

## Notes on NASTRAN-95 Fortran Code

- Fixed-format Fortran 77 (columns 1-72)
- Column 1: 'C' for comments
- Column 6: continuation character
- Columns 7-72: statements
- Variable names ≤ 6 characters
- Implicit typing: I-N are integers, others are real
- Heavy use of COMMON blocks for data sharing
- EQUIVALENCE statements for memory reuse
