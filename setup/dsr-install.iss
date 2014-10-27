; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "DSR - Disordered Structure Refinement"
#define MyAppVersion "1.5.13"
#define MyAppPublisher "Daniel Kratzert"

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{D8BF1E9B-3EFD-408E-91CB-A69D39445C7F}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={pf}\DSR
OutputBaseFilename=DSR-setup-{#MyAppVersion}
Compression=lzma/fast
SolidCompression=yes
SetupLogging=True
CloseApplications=False
RestartApplications=False
ShowLanguageDialog=no
ChangesAssociations=True
RestartIfNeededByRun=False
ChangesEnvironment=True
DisableFinishedPage=True
DisableReadyPage=True
DisableReadyMemo=True
DisableWelcomePage=True
AlwaysShowDirOnReadyPage=True
InternalCompressLevel=fast
EnableDirDoesntExistWarning=True
DirExistsWarning=no
UninstallLogMode=new
VersionInfoVersion={#MyAppVersion}
MinVersion=0,5.01
DefaultGroupName=DSR
DisableProgramGroupPage=yes
AppendDefaultGroupName=True
AppContact=dkratzert@gmx.de
AppCopyright=Daniel Kratzert
AppSupportPhone=+49 761 201 6156
VersionInfoProductName=DSR - Disordered Structure Refinement
AlwaysShowComponentsList=False
ShowComponentSizes=False

[Run]
;Filename: "msiexec.exe"; Parameters: "TARGETDIR=""{app}\python27"" ADDLOCAL=Extensions,SharedCRT /i ""{app}\python-2.7.5.msi"""; Description: "Install Python"

; why does the quiet installation of python with parameter /qn not install the visual C++ redistributable?
;Filename: "msiexec.exe"; Parameters: "/qn TARGETDIR=""{app}\python27"" /qn ADDLOCAL=Extensions,SharedCRT /i ""{app}\python-2.7.5.msi"""; Description: "Install Python"

[UninstallRun]
;Filename: "msiexec.exe"; Parameters: "/qn /x{{DBDD570E-0952-475f-9453-AB88F3DD5659}"; WorkingDir: "{app}"

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

; adds a new page to the setup where you can choose if the path should be added
[Tasks]
Name: "modifypath"; Description: "Add application directory to your environmental path"; GroupDescription: "DSR needs to be in the PATH variable to find its components."

[Files]
Source: "..\dsr_db.txt"; DestDir: "{app}"; Flags: ignoreversion; Permissions: users-modify
Source: "..\dsr.bat"; DestDir: "{app}"; Flags: ignoreversion; Permissions: users-modify
Source: "..\dsr-shelxle.bat"; DestDir: "{app}"; Flags: ignoreversion; Permissions: users-modify
;Source: "C:\Users\ntrapp\Downloads\python-2.7.5.msi"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\dsr.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\manuals\DSR-manual.pdf"; DestDir: "{app}\manual"; Flags: ignoreversion
Source: "..\atoms.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\afix.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\atomhandling.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\constants.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\dbfile.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\dsrparse.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\export.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\terminalsize.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\misc.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\options.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\refine.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\restraints.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\resfile.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\resi.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\pyperclip.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\setup\OlexDSR.py"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\setup\custom.xld"; DestDir: "{app}"; Flags: ignoreversion
;Source: "C:\Users\ntrapp\Documents\Python27\*"; DestDir: "{app}\Python27"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "C:\Python27-dsr\*"; DestDir: "{app}\Python27"; Flags: ignoreversion createallsubdirs recursesubdirs
;Source: "..\setup\vcredist_x86.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\example\*"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: ".\Output\changelog.txt"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\dsr_user_db.txt"; DestDir: "{app}"; Flags: ignoreversion confirmoverwrite uninsneveruninstall onlyifdoesntexist; Permissions: users-modify
Source: "..\networkx\*"; DestDir: "{app}\networkx"; Flags: ignoreversion createallsubdirs recursesubdirs

[Icons]
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{group}\DSR manual"; Filename: "{app}\manual\DSR-manual.pdf"; WorkingDir: "{app}"
Name: "{group}\DSR database"; Filename: "{app}\dsr_db.txt"; WorkingDir: "{app}"
Name: "{group}\DSR user-database"; Filename: "{app}\dsr_user_db.txt"

[Registry]
Root: "HKLM"; Subkey: "Software\DSR"; ValueType: string; ValueName: "dsr_directory"; ValueData: "{app}"; Flags: uninsdeletekey
Root: "HKLM"; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "DSRDIR"; ValueData: "{app}"; Flags: deletevalue
Root: "HKLM"; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "DSR_DB_DIR"; ValueData: "{app}"; Flags: deletevalue

[UninstallDelete]
Type: files; Name: "{app}\*.pyc"
Type: filesandordirs; Name: "{app}\python27"

[Dirs]
Name: "{app}\example"; Permissions: authusers-full
Name: "{app}\manual"

[Code]
const 
    ModPathName = 'modifypath'; 
    ModPathType = 'system'; 

function ModPathDir(): TArrayOfString; 
begin 
    setArrayLength(Result, 1) 
    Result[0] := ExpandConstant('{app}'); 
end; 
#include "modpath.iss"
