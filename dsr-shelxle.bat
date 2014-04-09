@ECHO OFF

rem DSR startup script for ShelXle in windows
rem echo %DSRDIR%

SET cmd=%1
SET args=%2

IF "%cmd%"=="" GOTO help

:argloop
IF NOT "%3"=="" SET args=%args% %3
SHIFT
IF NOT "%3"=="" GOTO argloop
IF "%cmd%"=="-r" GOTO resfile
IF "%cmd%"=="-re" GOTO resfile_ext
IF "%cmd%"=="-e" GOTO export
IF "%cmd%"=="-o" GOTO clip
IF "%cmd%"=="-h" GOTO help
IF "%cmd%"=="-l" GOTO list
IF "%cmd%"=="-n" GOTO noref
IF "%cmd%"=="--help" GOTO help
GOTO end

:resfile
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -r %args%
GOTO end

:resfile_ext
rem IF "%args%" == "" GOTO help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -re %args%
GOTO end

:export
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -e %args%
GOTO end

:clip
rem IF "%args%" == "" GOTO help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -o "%args%"
GOTO end

:help
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -h
GOTO end

:list
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -l
GOTO end

:noref
"%DSRDIR%"\python27\python.exe "%DSRDIR%"\dsr.py -n %args%
GOTO end

:end
pause
