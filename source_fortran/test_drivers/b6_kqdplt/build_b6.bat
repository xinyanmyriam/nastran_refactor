@echo off
REM Build script for B6 KQDPLT (Quadrilateral Plate Bending) test driver
REM Uses gfortran from MSYS2 UCRT64

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe
set NASDIR=F:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis

echo ============================================================
echo  B6: KQDPLT - Quadrilateral Plate Bending Element Test Driver
echo ============================================================
echo.

echo Copying source files from NASTRAN...
copy "%NASDIR%\kqdplt.f" . >nul 2>&1
copy "%NASDIR%\ktrbsc.f" . >nul 2>&1
copy "%NASDIR%\inverd.f" . >nul 2>&1
copy "%NASDIR%\gmmatd.f" . >nul 2>&1

echo Compiling...
echo   Sources: test_kqdplt.f stubs_kqdplt.f kqdplt.f ktrbsc.f inverd.f gmmatd.f
echo.

"%GFORTRAN%" -o test_kqdplt.exe ^
    test_kqdplt.f ^
    stubs_kqdplt.f ^
    kqdplt.f ^
    ktrbsc.f ^
    inverd.f ^
    gmmatd.f ^
    -ffixed-form -std=legacy -w -fallow-argument-mismatch

if %ERRORLEVEL% EQU 0 (
    echo Compilation successful!
    echo.
    echo Running test...
    echo ============================================================
    test_kqdplt.exe > reference_output.txt 2>&1
    type reference_output.txt
    echo.
    echo Output saved to reference_output.txt
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
