@echo off
REM Build script for SOLVE test driver (B9 - Linear Equation Solver)
REM Uses gfortran from MSYS2 UCRT64
REM
REM NOTE: This is a standalone reference implementation.
REM The NASTRAN solve.f uses GINO file I/O and is not suitable
REM for standalone testing. This implements the same Cholesky
REM decomposition algorithm independently.

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe

echo Compiling SOLVE test driver (standalone Cholesky)...
"%GFORTRAN%" -o test_solve.exe test_solve.f -ffixed-form -std=legacy -w

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Compilation successful!
    echo Running test...
    echo.
    test_solve.exe > reference_output.txt 2>&1
    type reference_output.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
