@echo off
setlocal
cd /d "%~dp0"
if "%~1"=="" (
  echo Drag an AutoChem/MKS CSV or Excel export onto this file.
  echo.
  pause
  exit /b 1
)
python ms_deconvolution.py "%~1"
if errorlevel 1 (
  echo.
  echo Analysis failed. Read the error above.
  pause
  exit /b 1
)
echo.
echo Analysis finished. Results are next to the input file.
pause
