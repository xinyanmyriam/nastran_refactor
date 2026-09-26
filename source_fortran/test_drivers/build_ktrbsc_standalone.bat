@echo off
set GFORTRAN=F:\msys64\ucrt64\bin\gfortran.exe

echo Compiling KTRBSC standalone test...
"%GFORTRAN%" -o test_ktrbsc_standalone.exe test_ktrbsc_standalone.f stubs_plate.f ktrbsc.f inverd.f gmmatd.f -ffixed-form -std=legacy -w -fallow-argument-mismatch

if %ERRORLEVEL% EQU 0 (
    echo Compilation OK. Running...
    test_ktrbsc_standalone.exe
) else (
    echo FAILED
)
