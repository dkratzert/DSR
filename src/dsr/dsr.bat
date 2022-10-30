@ECHO OFF

rem DSR startup script for windows
rem echo %DSR_DIR%
rem set DSR_DIR="C:\Program Files (x86)\DSR\"

TITLE "DSR - Disordered Structure Refinement"

SET args=%*

rem cls

IF NOT DEFINED DSR_DIR (GOTO setdsrdir) ELSE (GOTO arguments)


:setdsrdir
    SET DSR_DIR="."
    goto arguments

:arguments
    IF DEFINED args (GOTO main) ELSE (GOTO help)    
    goto main

:main
    "%DSR_DIR%"\python3\python.exe "%DSR_DIR%"\dsr.py %args%
    GOTO end

:help
    "%DSR_DIR%"\python3\python.exe "%DSR_DIR%"\dsr.py -h
    GOTO end

:end