@echo off
REM Build script for B4 KQDMEM (CQUAD4 quadrilateral membrane) test driver
REM
REM NOTE: the MSYS2 bin directory MUST be on PATH. gfortran's driver resolves
REM f951.exe, as.exe and collect2.exe through PATH; without it the driver runs,
REM prints nothing, exits 0, and produces NO output file. Diagnosed 2026-08-19.
set MSYS2BIN=F:\msys64\ucrt64\bin
set PATH=%MSYS2BIN%;%PATH%
set GFORTRAN=%MSYS2BIN%\gfortran.exe
set NASDIR=F:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis

echo ============================================================
echo  B4: KQDMEM + KTRMEM - CQUAD4 quadrilateral membrane
echo ============================================================
echo.
echo Copying unmodified sources from NASTRAN...
copy "%NASDIR%\kqdmem.f" . >nul 2>&1
copy "%NASDIR%\ktrmem.f" . >nul 2>&1
copy "%NASDIR%\gmmatd.f" . >nul 2>&1
echo Compiling...
echo   Sources: test_kqdmem.f stubs_kqdmem.f kqdmem.f ktrmem.f gmmatd.f
echo.
"%GFORTRAN%" -o test_kqdmem.exe ^
    test_kqdmem.f ^
    stubs_kqdmem.f ^
    kqdmem.f ^
    ktrmem.f ^
    gmmatd.f ^
    -ffixed-form -std=legacy -w -fallow-argument-mismatch
if %ERRORLEVEL% EQU 0 (
    echo Compilation successful!
    echo.
    echo Running test...
    echo ============================================================
    test_kqdmem.exe > reference_output.txt 2>&1
    type reference_output.txt
    echo.
    echo Output saved to reference_output.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
