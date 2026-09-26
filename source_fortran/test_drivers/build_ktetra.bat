@echo off
REM Build script for KTETRA test driver
REM Uses gfortran from MSYS2 UCRT64

set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe
set SRCDIR=F:\ma_sim_paper\nastran_refactor\nastran\NASTRAN-95\mis

echo Copying source files...
copy "%SRCDIR%\ktetra.f" . >nul 2>&1
copy "%SRCDIR%\inverd.f" . >nul 2>&1
copy "%SRCDIR%\gmmatd.f" . >nul 2>&1

echo Compiling KTETRA test driver...
"%GFORTRAN%" -o test_ktetra.exe test_ktetra.f stubs.f ktetra.f inverd.f gmmatd.f -ffixed-form -std=legacy -w

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Compilation successful!
    echo Running test...
    echo.
    test_ktetra.exe
) else (
    echo.
    echo Compilation FAILED with error code %ERRORLEVEL%
)
