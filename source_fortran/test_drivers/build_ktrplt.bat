@echo off
REM Build script for KTRPLT (B5 triangular plate bending) test driver
REM Uses gfortran from MSYS2 UCRT64

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe
set SRCDIR=F:\ma_sim_paper\nastran_refactor\source_fortran\benchmarks
set NASDIR=F:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis

echo Copying source files...
copy "%SRCDIR%\ktrplt.f" . >nul 2>&1
copy "%NASDIR%\ktrbsc.f" . >nul 2>&1

echo.
echo Compiling KTRPLT test driver...
echo   Sources: test_ktrplt.f stubs_plate.f ktrplt.f ktrbsc.f inverd.f gmmatd.f
echo.

"%GFORTRAN%" -o test_ktrplt.exe ^
    test_ktrplt.f ^
    stubs_plate.f ^
    ktrplt.f ^
    ktrbsc.f ^
    inverd.f ^
    gmmatd.f ^
    -ffixed-form -std=legacy -w -fallow-argument-mismatch

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Compilation successful!
    echo Running test...
    echo.
    test_ktrplt.exe
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
    echo.
    echo Common issues:
    echo   - COMMON block size mismatch: check SMA1DP workspace size
    echo   - Argument type mismatch: use -fallow-argument-mismatch
    echo   - Missing subroutines: check stubs_plate.f has all needed stubs
)
