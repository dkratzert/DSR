@ECHO ON

rem DSR startup script for windows
echo %DSR_DIR%

rem set DSRDIR="d:\Programme\DSR"
rem cd C:\Users\ntrapp\Documents\GitHub\DSR
set GIT=D:\GitHub\DSR
rem set GIT="."

cd %GIT%
rem python %GIT%\rss-feed.py

xcopy /Y "%GIT%\afix.py" "%DSR_DIR%"
xcopy /Y "%GIT%\atomhandling.py" "%DSR_DIR%"
xcopy /Y "%GIT%\atoms.py" "%DSR_DIR%"
xcopy /Y "%GIT%\constants.py" "%DSR_DIR%"
xcopy /Y "%GIT%\dbfile.py" "%DSR_DIR%"
xcopy /Y "%GIT%\dsr.py" "%DSR_DIR%"
xcopy /Y "%GIT%\dsrparse.py" "%DSR_DIR%"
xcopy /Y "%GIT%\elements.py" "%DSR_DIR%"
xcopy /Y "%GIT%\export.py" "%DSR_DIR%"
xcopy /Y "%GIT%\misc.py" "%DSR_DIR%"
xcopy /Y "%GIT%\options.py" "%DSR_DIR%"
xcopy /Y "%GIT%\resfile.py" "%DSR_DIR%"
xcopy /Y "%GIT%\resi.py" "%DSR_DIR%"
xcopy /Y "%GIT%\pyperclip.py" "%DSR_DIR%"
xcopy /Y "%GIT%\refine.py" "%DSR_DIR%"
xcopy /Y "%GIT%\restraints.py" "%DSR_DIR%"
xcopy /Y "%GIT%\cf3fit.py" "%DSR_DIR%"
xcopy /Y "%GIT%\selfupdate.py" "%DSR_DIR%"
xcopy /Y "%GIT%\dsr.bat" "%DSR_DIR%"
xcopy /Y "%GIT%\terminalsize.py" "%DSR_DIR%"
xcopy /Y "%GIT%\dsr_db.txt" "%DSR_DIR%"
xcopy /Y "%GIT%\manuals\DSR-manual.pdf" "%DSR_DIR%\manuals"
xcopy /Y "%GIT%\example\p21c.res" "%DSR_DIR%\example"
xcopy /Y "%GIT%\example\p21c.hkl" "%DSR_DIR%\example"
xcopy "%GIT%\fit" "%DSR_DIR%\fit" /e /v /r /y /f

rem xcopy /Y %GIT%\update-dsr.bat %DSRDIR%

del "%DSR_DIR%\*.pyc"

rem copy %GIT%\dsr.bat %DSRDIR%
sleep 1s
rem pause
