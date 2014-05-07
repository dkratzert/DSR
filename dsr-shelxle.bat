@ECHO OFF

rem DSR startup script for ShelXle in windows
rem echo %DSRDIR%
rem set DSRDIR="C:\Program Files (x86)\DSR\"

SET cmd=%1
SET args=%*

rem deletes '"' from arguments for compatibility with ShelXle:
SET args=%args:"=%


IF DEFINED args (GOTO main) ELSE (GOTO help)


:main
    "%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py %args%
GOTO end

:help
    "%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -h
GOTO end


:end
pause
