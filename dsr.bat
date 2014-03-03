@ECHO OFF

rem DSR startup script for windows
rem echo %DSRDIR%
rem set DSRDIR="C:\Program Files (x86)\DSR\"
SET cmd=%1
SET args=%2
cls
IF "%cmd%"=="" GOTO help

:argloop
IF NOT "%3"=="" SET args=%args% %3
SHIFT
IF NOT "%3"=="" GOTO argloop
IF "%cmd%"=="-r" GOTO resfile
IF "%cmd%"=="-e" GOTO export
IF "%cmd%"=="-ea" GOTO export-all
IF "%cmd%"=="-h" GOTO help
IF "%cmd%"=="-l" GOTO list
IF "%cmd%"=="-n" GOTO noref
IF "%cmd%"=="-i" GOTO import
IF "%cmd%"=="--help" GOTO help
GOTO end

:resfile
rem IF "%args%" == "" GOTO help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -r %args%
GOTO end

:export
rem IF "%args%" == "" GOTO help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -e "%args%"
GOTO end

:export-all
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -ea
GOTO end

:help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -h
GOTO end

:list
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -l
GOTO end

:noref
IF "%args%" == "" GOTO help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -n %args%
GOTO end

:import
rem IF "%args%" == "" GOTO help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -i %args%
GOTO end

:end
