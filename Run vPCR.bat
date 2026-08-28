@echo off
rem Launch the virtualPCR Studio GUI. Double-click this file.
rem
rem The console stays open behind the window on purpose: it is where Python
rem reports anything the GUI cannot. Close it and the app closes with it.
rem
rem The interpreter is found by searching PATH inside the FOR below, not by
rem running `py --version`: py.exe is a launcher, so every such check spawns two
rem processes and flashes a window. Dependency checks live in vpcr/__main__.py.
rem (Do not write the FOR substitution token in a comment: cmd expands it here
rem too, and the file stops parsing.)
setlocal
cd /d "%~dp0"

set "PY="
for %%I in (python.exe) do if not defined PY set "PY=%%~$PATH:I"

if not defined PY (
    echo Python was not found.
    echo.
    echo Install Python 3.10 or newer from https://www.python.org/downloads/
    echo and tick "Add python.exe to PATH" during setup.
    goto :fail
)

echo Starting virtualPCR Studio...
echo Keep this window open. Errors, if any, appear here.
echo.
"%PY%" -m vpcr
if errorlevel 1 goto :fail
exit /b 0

:fail
echo.
echo Press any key to close.
pause >nul
exit /b 1
