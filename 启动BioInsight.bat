@echo off
setlocal

set "APP_DIR=%~dp0"
set "R_SCRIPT="

if exist "%APP_DIR%R\bin\Rscript.exe" (
  set "R_SCRIPT=%APP_DIR%R\bin\Rscript.exe"
)

if not defined R_SCRIPT if exist "%APP_DIR%R\bin\x64\Rscript.exe" (
  set "R_SCRIPT=%APP_DIR%R\bin\x64\Rscript.exe"
)

if not defined R_SCRIPT if exist "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" (
  set "R_SCRIPT=C:\Program Files\R\R-4.5.3\bin\Rscript.exe"
)

if not defined R_SCRIPT (
  for /f "delims=" %%R in ('where Rscript 2^>nul') do (
    if not defined R_SCRIPT set "R_SCRIPT=%%R"
  )
)

if not defined R_SCRIPT (
  echo Cannot find Rscript.exe. Please install R or add Rscript to PATH.
  pause
  exit /b 1
)

echo Starting BioInsight one-click bioinformatics platform...
"%R_SCRIPT%" -e "shiny::runApp('%APP_DIR:\=/%', launch.browser = TRUE, host = '127.0.0.1', port = 3838)"
pause


