; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "DSR"
;#define MyAppVersion "237" ; Now using commandline argument
#define MyAppPublisher "Daniel Kratzert"

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{D8BF1E9B-3EFD-408E-91CB-A69D39445C7F}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={commonpf}\DSR
OutputBaseFilename=DSR-setup-{#MyAppVersion}
Compression=lzma2/fast
SolidCompression=yes
SetupLogging=True
CloseApplications=False
RestartApplications=False
ShowLanguageDialog=no
ChangesAssociations=True
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
MinVersion=10.0.10240
DefaultGroupName=DSR
DisableProgramGroupPage=yes
AppendDefaultGroupName=True
AppContact=dkratzert@gmx.de
AppCopyright=Daniel Kratzert
AppSupportPhone=+49 761 203 6156
VersionInfoProductName=DSR - Disordered Structure Refinement
AlwaysShowComponentsList=False
ShowComponentSizes=False
;SignTool=signtool

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
Name: "modifypath"; Description: "Add application directory to your system environment path variable"; GroupDescription: "DSR needs to be in the PATH variable to find its components."

[Files]
Source: "..\..\python311\*"; DestDir: "{app}\python311"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\example\p21c.hkl"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\example\p21c.res"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs; Permissions: everyone-full
Source: "..\example\p21c_step0.res"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\example\p21c_step1.res"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\example\p21c_step3.res"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\example\p21c_step2.ins"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\example\p21c_final.res"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\example\p21n_cf3.hkl"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: "..\example\p21n_cf3.res"; DestDir: "{app}\example"; Flags: ignoreversion createallsubdirs recursesubdirs
Source: ".\Output\changelog.txt"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\manuals\DSR-manual.pdf"; DestDir: "{app}\manuals\"; Flags: ignoreversion
Source: "..\src\dsr_shelx\*"; DestDir: "{app}"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"

[Icons]
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{group}\DSR manual"; Filename: "{app}\manuals\DSR-manual.pdf"; WorkingDir: "{app}"
Name: "{group}\DSR database"; Filename: "{app}\dsr_db.txt"; WorkingDir: "{app}"
Name: "{group}\DSR user-database"; Filename: "{%USERPROFILE}\dsr_user_db.txt"

[Registry]
Root: "HKLM"; Subkey: "Software\DSR"; ValueType: string; ValueName: "dsr_directory"; ValueData: "{app}"; Flags: uninsdeletekey
Root: "HKLM"; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "DSR_DIR"; ValueData: "{app}"; Flags: deletevalue
Root: "HKLM"; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: string; ValueName: "DSR_DB_DIR"; ValueData: "{app}"; Flags: deletevalue

[InstallDelete]
Type: filesandordirs; Name: "{app}\networkx"
Type: filesandordirs; Name: "{app}\mpmath"
Type: filesandordirs; Name: "{app}\python27"

[UninstallDelete]
Type: files; Name: "{app}\*.pyc"
Type: files; Name: "{app}\dsr"
Type: filesandordirs; Name: "{app}\python311"
Type: filesandordirs; Name: "{app}\networkx"
Type: filesandordirs; Name: "{app}\mpmath"
Type: filesandordirs; Name: "{app}\example"
Type: filesandordirs; Name: "{app}\fit"
Type: filesandordirs; Name: "{app}\manuals"
Type: filesandordirs; Name: "{app}\setup"
Type: filesandordirs; Name: "{app}"

[Dirs]
Name: "{app}\example"; Permissions: everyone-full
Name: "{app}\manuals"; Permissions: everyone-full
Name: "{app}\."; Permissions: everyone-full

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
