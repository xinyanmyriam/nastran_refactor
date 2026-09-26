@echo off
REM Build script for B9 INVERD test driver
REM NASTRAN-95 matrix inversion and linear equation solver
REM

echo ============================================
echo  Building B9 INVERD Test Driver
echo ============================================
echo.

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe

REM Check if gfortran exists
if not exist "%GFORTRAN%" (
    echo ERROR: gfortran not found at %GFORTRAN%
    exit /b 1
)

REM Compile
echo Compiling test_inverd.f + inverd.f ...
%GFORTRAN% -o test_inverd.exe test_inverd.f inverd.f -ffixed-form -std=legacy -w

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ERROR: Compilation failed!
    exit /b 1
)

echo.
echo Build successful: test_inverd.exe
echo.

REM Run
echo Running test_inverd.exe ...
echo ============================================
test_inverd.exe
echo.
echo ============================================
echo Done.
