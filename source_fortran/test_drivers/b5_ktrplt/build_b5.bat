@echo off
REM Build script for B5 KTRPLT (Triangular Plate Bending) test driver
REM Uses gfortran from MSYS2 UCRT64
REM Test geometry: right-angled triangle (0,0,0)-(1,0,0)-(0,1,0)

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe
set NASDIR=F:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis

echo ============================================================
echo  B5: KTRPLT - Triangular Plate Bending Element Test Driver
echo  Geometry: Right-angled triangle in XY plane
echo ============================================================
echo.

echo Copying source files from NASTRAN...
copy "%NASDIR%\ktrplt.f" . >nul 2>&1
copy "%NASDIR%\ktrbsc.f" . >nul 2>&1
copy "%NASDIR%\inverd.f" . >nul 2>&1
copy "%NASDIR%\gmmatd.f" . >nul 2>&1

echo Compiling...
echo   Sources: test_ktrplt.f stubs_ktrplt.f ktrplt.f ktrbsc.f inverd.f gmmatd.f
echo.

"%GFORTRAN%" -o test_ktrplt.exe ^
    test_ktrplt.f ^
    stubs_ktrplt.f ^
    ktrplt.f ^
    ktrbsc.f ^
    inverd.f ^
    gmmatd.f ^
    -ffixed-form -std=legacy -w -fallow-argument-mismatch

if %ERRORLEVEL% EQU 0 (
    echo Compilation successful!
    echo.
    echo Running test...
    echo ============================================================
    test_ktrplt.exe > reference_output.txt 2>&1
    type reference_output.txt
    echo.
    echo Output saved to reference_output.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
