@ECHO OFF

rem DSR startup script for ShelXle in windows
rem echo %DSRDIR%
rem set DSRDIR="C:\Program Files (x86)\DSR\"

SET cmd=%1
SET args=%*

rem cls

IF "%args%"=="" (GOTO help)


GOTO main

:main
    "%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py %args%
GOTO end

:help
    "%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -h
GOTO end


:end
pause
