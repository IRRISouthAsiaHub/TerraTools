; installer.iss — Inno Setup script for TerraTools.
;
; Produces a single Setup.exe that installs:
;   1. The PyInstaller-built app (dist/TerraTools/*)
;   2. A complete, pre-packaged portable R runtime with all required
;      packages already installed (r_runtime_build/*)
;
; End users run Setup.exe and get a working app with zero other installs
; — no Python, no R, no conda, nothing else needed.
;
; Prerequisites before compiling this script:
;   - Download Inno Setup: https://jrsoftware.org/isinfo.php (free)
;   - Have already run `pyinstaller terratools.spec` (produces dist/TerraTools/)
;   - Have already prepared a portable R build in r_runtime_build/
;     (see BUILD_GUIDE.md for exactly how to prepare this folder)
;
; Compile with the Inno Setup Compiler (ISCC.exe) or the Inno Setup IDE.

#define MyAppName "TerraTools"
#define MyAppVersion "1.0.0"
#define MyAppPublisher "IRRI-SAH"
#define MyAppExeName "TerraTools.exe"

[Setup]
AppId={{95CA17BB-50CA-4B87-984A-F67CA81664D1}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
OutputDir=installer_output
OutputBaseFilename=TerraTools_Setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"

[Files]
; The full PyInstaller output folder (app exe + all bundled Python deps)
Source: "dist\TerraTools\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

; The prepared portable R runtime, placed at {app}\r_runtime so
; r_locator.py's frozen-mode lookup (<exe_dir>/r_runtime/bin/Rscript.exe)
; finds it automatically.
Source: "r_runtime_build\*"; DestDir: "{app}\r_runtime"; Flags: ignoreversion recursesubdirs createallsubdirs

[Dirs]
; Program Files is write-protected for standard (non-admin) users by
; default. Several tools write their results/ folders relative to their
; own install location — without this, every tool's output would
; silently fail to save for anyone not running as admin. This grants
; standard users modify (read/write/create/delete) access to the whole
; install folder, which is the simplest fix that doesn't require
; touching any tool's Python code.
Name: "{app}"; Permissions: users-modify

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch {#MyAppName}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
Type: filesandordirs; Name: "{localappdata}\TerraTools"
