@echo off
REM Build script for INVPWR test driver (B10 - Inverse Power Eigenvalue)
REM Uses gfortran from MSYS2 UCRT64
REM
REM NOTE: This is a standalone reference implementation.
REM The NASTRAN invpwr.f uses GINO file I/O and is not suitable
REM for standalone testing. This implements the same inverse power
REM iteration algorithm independently.

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe

echo Compiling INVPWR test driver (standalone inverse iteration)...
"%GFORTRAN%" -o test_invpwr.exe test_invpwr.f -ffixed-form -std=legacy -w

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Compilation successful!
    echo Running test...
    echo.
    test_invpwr.exe > reference_output.txt 2>&1
    type reference_output.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
