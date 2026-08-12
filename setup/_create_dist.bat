@echo off

REM This script builds a working Python environment into ..\dist of the current file location.
REM Modeled after FinalCif's scripts\_create_dist.bat.

REM Accept PYTHON_VERSION as first argument, see make_win_release.bat
set PYTHON_VERSION=%~1

set PYTHON_URL=https://www.python.org/ftp/python/%PYTHON_VERSION%/python-%PYTHON_VERSION%-embed-amd64.zip

for %%A in ("%~dp0.") do set "SCRIPT_DIR=%%~fA"
set BUILD_DIR=%SCRIPT_DIR%\..\dist
set PACKAGE_DIR=%BUILD_DIR%\python_dist

REM Check if uv is available
where uv >NUL 2>&1
if %errorlevel% neq 0 (
    echo uv not found, installing...
    powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
    set "PATH=%USERPROFILE%\.local\bin;%PATH%;%PACKAGE_DIR%\Scripts"
)

setlocal enabledelayedexpansion
for %%a in (!PYTHON_VERSION!) do (
    set "NEW_PYTHON_VERSION=%%~na"
)
set "SHORT_PYTHON_VERSION=!NEW_PYTHON_VERSION:.=!"

mkdir %BUILD_DIR% 2>NUL
cd %BUILD_DIR%

curl %PYTHON_URL% -o python-%PYTHON_VERSION%.zip

del /S /Q /F %PACKAGE_DIR%\*.* >NUL 2>NUL
rmdir /s /q %PACKAGE_DIR% > NUL 2>NUL
mkdir %PACKAGE_DIR%

REM The Windows embeddable package already ships without tkinter/tcl, Doc, Tools,
REM the test suite or IDLE -- there is nothing extra to strip out.
tar -xf python-%PYTHON_VERSION%.zip -C %PACKAGE_DIR%
del python-%PYTHON_VERSION%.zip

echo python%SHORT_PYTHON_VERSION%.zip > %PACKAGE_DIR%\python%SHORT_PYTHON_VERSION%._pth
echo . >> %PACKAGE_DIR%\python%SHORT_PYTHON_VERSION%._pth
echo import site >> %PACKAGE_DIR%\python%SHORT_PYTHON_VERSION%._pth
endlocal

del vc_redist.x64.exe 2>NUL

curl -L https://aka.ms/vs/17/release/vc_redist.x64.exe -o vc_redist.x64.exe

cd %SCRIPT_DIR%\..

REM Install the project's runtime dependencies (if any are ever added to
REM requirements.txt) straight into the embedded interpreter. Currently a no-op
REM since DSR ships all its dependencies (mpmath, networkx) vendored in-tree, but
REM this keeps the release script future-proof without hardcoding "no deps".
uv pip install --python "%PACKAGE_DIR%\python.exe" -r requirements.txt --link-mode=copy
if %errorlevel% neq 0 (
    echo uv pip install failed. Stopping now.
    exit /b %errorlevel%
)

echo - finished environment install!

CALL "%PACKAGE_DIR%\python.exe" "%SCRIPT_DIR%\make_win_release.py"
