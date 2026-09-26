@echo off
REM ============================================================
REM Build script for B12: Mini-NASTRAN runtime for TRD1C
REM 
REM This builds the actual NASTRAN trd1c.f, step.f, form1.f 
REM with a mini-GINO simulation layer for 3-DOF transient analysis
REM
REM VERIFIED WORKING: 2025 - compiles and runs correctly
REM ============================================================

set FC=F:\msys64\ucrt64\bin\gfortran.exe
set SRC=f:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis
set FFLAGS=-ffixed-form -std=legacy -w -fallow-argument-mismatch -fno-range-check -O0 -g -fbacktrace

echo.
echo === Building B12: TRD1C Native Mini-NASTRAN ===
echo.

REM Step 1: Compile all object files
echo [1/2] Compiling source files...
%FC% -c test_trd1c_native.f %FFLAGS%
if %ERRORLEVEL% NEQ 0 goto :fail
%FC% -c gino_sim.f %FFLAGS%
if %ERRORLEVEL% NEQ 0 goto :fail
%FC% -c %SRC%\trd1c.f %FFLAGS%
if %ERRORLEVEL% NEQ 0 goto :fail
%FC% -c %SRC%\step.f %FFLAGS%
if %ERRORLEVEL% NEQ 0 goto :fail
%FC% -c %SRC%\form1.f %FFLAGS%
if %ERRORLEVEL% NEQ 0 goto :fail
%FC% -c %SRC%\form2.f %FFLAGS%
if %ERRORLEVEL% NEQ 0 goto :fail
%FC% -c %SRC%\inverd.f %FFLAGS%
if %ERRORLEVEL% NEQ 0 goto :fail

REM Step 2: Link
echo [2/2] Linking...
%FC% -o test_trd1c_native.exe test_trd1c_native.o gino_sim.o trd1c.o step.o form1.o form2.o inverd.o -g
if %ERRORLEVEL% NEQ 0 goto :fail

echo.
echo === Build successful! ===
echo.
echo Running test...
echo ================================================
echo.

test_trd1c_native.exe

echo.
echo ================================================
echo === Done ===
goto :end

:fail
echo.
echo *** BUILD FAILED ***
echo.

:end
pause
