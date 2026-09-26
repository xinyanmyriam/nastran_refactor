@echo off
REM ================================================================
REM Build script for B12_TRD1C: Standalone Newmark-Beta test
REM Uses NASTRAN-95's INVERD for linear solve
REM ================================================================

set GFC=F:\msys64\ucrt64\bin\gfortran.exe
set OPTS=-ffixed-form -std=legacy -w -fallow-argument-mismatch

echo Building test_newmark.exe ...
%GFC% -o test_newmark.exe test_newmark.f inverd.f %OPTS%

if %ERRORLEVEL% NEQ 0 (
    echo BUILD FAILED
    exit /b 1
)

echo Build successful.
echo Running test_newmark.exe ...
test_newmark.exe > reference_output.txt
type reference_output.txt

echo.
echo Output saved to reference_output.txt
