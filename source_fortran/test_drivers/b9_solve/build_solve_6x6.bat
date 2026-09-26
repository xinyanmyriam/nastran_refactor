@echo off
REM Build script for SOLVE 6x6 test driver (B9 - Linear Equation Solver)
REM Uses gfortran from MSYS2 UCRT64
REM
REM This is the 6x6 version with known solution x=[1,2,3,4,5,6]
REM implementing Cholesky decomposition as a standalone reference.

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe

echo ============================================================
echo  B9: SOLVE - 6x6 Cholesky Decomposition Reference
echo ============================================================
echo.

echo Compiling SOLVE 6x6 test driver...
"%GFORTRAN%" -o test_solve_6x6.exe test_solve_6x6.f -ffixed-form -std=legacy -w

if %ERRORLEVEL% EQU 0 (
    echo Compilation successful!
    echo.
    echo Running test...
    echo ============================================================
    test_solve_6x6.exe > reference_output_6x6.txt 2>&1
    type reference_output_6x6.txt
    echo.
    echo Output saved to reference_output_6x6.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
