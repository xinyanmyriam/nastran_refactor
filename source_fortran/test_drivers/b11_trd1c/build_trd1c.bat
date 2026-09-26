@echo off
REM Build script for TRD1C test driver (B11 - Newmark-Beta Integration)
REM Uses gfortran from MSYS2 UCRT64
REM
REM NOTE: This is a standalone reference implementation.
REM The NASTRAN trd1c.f is deeply embedded in the transient response
REM module. This implements the same Newmark-beta algorithm independently
REM for both 1-DOF and 2-DOF systems.

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe

echo Compiling TRD1C test driver (standalone Newmark-beta)...
"%GFORTRAN%" -o test_trd1c.exe test_trd1c.f -ffixed-form -std=legacy -w

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Compilation successful!
    echo Running test...
    echo.
    test_trd1c.exe > reference_output.txt 2>&1
    type reference_output.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
