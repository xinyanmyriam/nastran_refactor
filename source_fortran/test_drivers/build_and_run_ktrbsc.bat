@echo off
cd /d %~dp0
F:\msys64\ucrt64\bin\gfortran.exe -o test_ktrbsc_standalone.exe test_ktrbsc_standalone.f stubs_plate.f ktrbsc.f inverd.f gmmatd.f -ffixed-form -std=legacy -w -fallow-argument-mismatch
if %ERRORLEVEL% NEQ 0 (
    echo COMPILE FAILED
    exit /b 1
)
echo COMPILE OK
test_ktrbsc_standalone.exe
echo EXIT CODE: %ERRORLEVEL%
