@echo off
REM Build script for B12 SROD1 (Rod Stress Recovery) test driver
REM Uses gfortran from MSYS2 UCRT64

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe
set NASDIR=F:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis

echo ============================================================
echo  B12: SROD1 - Rod Stress Recovery Test Driver
echo ============================================================
echo.

echo Copying source files from NASTRAN...
copy "%NASDIR%\srod1.f" . >nul 2>&1

echo Compiling...
echo   Sources: test_srod1.f stubs_srod1.f srod1.f
echo.

"%GFORTRAN%" -o test_srod1.exe ^
    test_srod1.f ^
    stubs_srod1.f ^
    srod1.f ^
    -ffixed-form -std=legacy -w -fallow-argument-mismatch

if %ERRORLEVEL% EQU 0 (
    echo Compilation successful!
    echo.
    echo Running test...
    echo ============================================================
    test_srod1.exe > reference_output.txt 2>&1
    type reference_output.txt
    echo.
    echo Output saved to reference_output.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
