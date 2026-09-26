@echo off
REM Build script for B8 KSOLID WEDGE (6-node triangular prism) test driver
REM Uses gfortran from MSYS2 UCRT64

REM NOTE: the MSYS2 bin directory MUST be on PATH. gfortran's driver resolves
REM f951.exe, as.exe and collect2.exe through PATH; without it the driver runs,
REM prints nothing, exits 0, and produces NO output file. Diagnosed 2026-08-19.
set MSYS2BIN=F:\msys64\ucrt64\bin
set PATH=%MSYS2BIN%;%PATH%
set GFORTRAN=%MSYS2BIN%\gfortran.exe
set NASDIR=F:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis

echo ============================================================
echo  B8: KSOLID - Wedge Element (ITYPE=1, 3 tetrahedra)
echo ============================================================
echo.

echo Copying source files from NASTRAN...
copy "%NASDIR%\ksolid.f" . >nul 2>&1
copy "%NASDIR%\ktetra.f" . >nul 2>&1
copy "%NASDIR%\inverd.f" . >nul 2>&1
copy "%NASDIR%\gmmatd.f" . >nul 2>&1

echo Compiling...
echo   Sources: test_ksolid_wedge.f stubs_ksolid_wedge.f ksolid.f ktetra.f inverd.f gmmatd.f
echo.

"%GFORTRAN%" -o test_ksolid_wedge.exe ^
    test_ksolid_wedge.f ^
    stubs_ksolid_wedge.f ^
    ksolid.f ^
    ktetra.f ^
    inverd.f ^
    gmmatd.f ^
    -ffixed-form -std=legacy -w -fallow-argument-mismatch

if %ERRORLEVEL% EQU 0 (
    echo Compilation successful!
    echo.
    echo Running test...
    echo ============================================================
    test_ksolid_wedge.exe > reference_output_wedge.txt 2>&1
    type reference_output_wedge.txt
    echo.
    echo Output saved to reference_output_wedge.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
