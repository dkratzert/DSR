@ECHO OFF

rem DSR startup script for windows
rem echo %DSRDIR%
rem set DSRDIR="C:\Program Files (x86)\DSR\"

TITLE "DSR - Disordered Structure Refinement"

SET cmd=%1
SET args=%*

rem cls

IF NOT DEFINED DSRDIR (GOTO setdsrdir) ELSE (GOTO arguments)


:setdsrdir
    echo -----------------------------------------------------------------
    echo You should define the DSRDIR environment variable for your system
    echo pointing to the DSR install directory.
    echo -----------------------------------------------------------------
    SET DSRDIR="."
    goto arguments

:arguments
    IF DEFINED args (GOTO main) ELSE (GOTO help)    
    goto main

:main
    "%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py %args%
    GOTO end

:help
    "%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -h
    GOTO end


:end
rem This is for ShelXle. Otherwise the DSR output would
rem vanish immediately in detached mode.
IF /i %0 NEQ dsr PAUSE