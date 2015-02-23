@ECHO OFF

rem DSR startup script for windows
echo %DSRDIR%
set DSRDIR="d:\Programme\DSR"
rem cd C:\Users\ntrapp\Documents\GitHub\DSR
set GIT="D:\GitHub\DSR"
rem set GIT="."

xcopy /Y %GIT%\afix.py %DSRDIR%
xcopy /Y %GIT%\atomhandling.py %DSRDIR%
xcopy /Y %GIT%\atoms.py %DSRDIR%
xcopy /Y %GIT%\constants.py %DSRDIR%
xcopy /Y %GIT%\dbfile.py %DSRDIR%
xcopy /Y %GIT%\dsr.py %DSRDIR%
xcopy /Y %GIT%\dsrparse.py %DSRDIR%
xcopy /Y %GIT%\elements.py %DSRDIR%
xcopy /Y %GIT%\export.py %DSRDIR%
xcopy /Y %GIT%\misc.py %DSRDIR%
xcopy /Y %GIT%\options.py %DSRDIR%
xcopy /Y %GIT%\resfile.py %DSRDIR%
xcopy /Y %GIT%\resi.py %DSRDIR%
xcopy /Y %GIT%\pyperclip.py %DSRDIR%
xcopy /Y %GIT%\refine.py %DSRDIR%
xcopy /Y %GIT%\restraints.py %DSRDIR%
xcopy /Y %GIT%\dsr.bat %DSRDIR%
xcopy /Y %GIT%\terminalsize.py %DSRDIR%
xcopy /Y %GIT%\dsr_db.txt %DSRDIR%
xcopy /Y %GIT%\manuals\DSR-manual.pdf %DSRDIR%
xcopy /Y %GIT%\example\p21c.res %DSRDIR%\example

rem xcopy /Y %GIT%\update-dsr.bat %DSRDIR%

del %DSRDIR%\*.pyc

rem copy %GIT%\dsr.bat %DSRDIR%
rem sleep 2s
rem pause
