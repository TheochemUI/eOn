# Activate MSVC include/lib for meson without putting MSVC bin first on PATH
# (MSVC bin first breaks torch DLL load: WinError 127 on shm.dll).

$kill = @(
  'VSINSTALLDIR','VS_VERSION','VS_MAJOR','VS_YEAR','VCVARSBAT',
  'NEWER_VS_WITH_OLDER_VC','WindowsSdkDir','WindowsSDKVer','MSSdk',
  'DISTUTILS_USE_SDK','VCINSTALLDIR','VCToolsInstallDir','VCToolsVersion',
  'CMAKE_GENERATOR','CMAKE_GENERATOR_TOOLSET','CMAKE_GENERATOR_PLATFORM',
  'CMAKE_GEN','CMAKE_PLAT','USE_NEW_CMAKE_GEN_SYNTAX','VSCMD_VER',
  'VSCMD_ARG_TGT_ARCH','VSCMD_ARG_HOST_ARCH','UniversalCRTSdkDir',
  'UCRTVersion','WindowsLibPath','WindowsSdkLibPath','WindowsSdkVerBinPath',
  'WindowsSDKVersion','ExtensionSdkDir','VCToolsRedistDir','FrameworkDir',
  'FrameworkDir64','FrameworkVersion','FrameworkVersion64','Framework40Version'
)
foreach ($v in $kill) { Remove-Item "Env:$v" -ErrorAction SilentlyContinue }

$vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $vswhere)) { throw "vswhere not found" }
$installPath = & $vswhere -latest -products * `
  -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
  -property installationPath
if (-not $installPath) { throw "vswhere: no VC tools" }
$vcvars = Join-Path $installPath "VC\Auxiliary\Build\vcvars64.bat"
if (-not (Test-Path $vcvars)) { throw "missing $vcvars" }
Write-Host "Using $vcvars"

$afterLines = cmd.exe /c "`"$vcvars`" >nul 2>&1 && set"
$after = @{}
foreach ($line in $afterLines) {
  if ($line -match '^([^=]+)=(.*)$') { $after[$Matches[1]] = $Matches[2] }
}

# Include/lib and VC locations only — NOT PATH (torch needs conda DLL order)
$exportNames = @(
  'INCLUDE','LIB','LIBPATH','VCINSTALLDIR','VCToolsInstallDir','VCToolsVersion',
  'WindowsSdkDir','WindowsSDKVersion','WindowsSdkVerBinPath','UniversalCRTSdkDir',
  'UCRTVersion','WindowsLibPath','WindowsSdkLibPath','VSINSTALLDIR','VSCMD_VER',
  'VSCMD_ARG_TGT_ARCH','VSCMD_ARG_HOST_ARCH','FrameworkDir','FrameworkDir64',
  'FrameworkVersion','FrameworkVersion64','DevEnvDir','ExtensionSdkDir',
  'VCToolsRedistDir','WindowsSdkBinPath','NETFXSDKDir','WindowsSDK_ExecutablePath_x64',
  'WindowsSDK_ExecutablePath_x86','WindowsSDKLibVersion','VCIDEInstallDir',
  'VSCMD_ARG_app_plat','IFCPATH','ExternalCompilerOptions'
)
foreach ($name in $exportNames) {
  if ($after.ContainsKey($name)) {
    Add-Content -Path $env:GITHUB_ENV -Value "$name=$($after[$name])"
  }
}

# Absolute tool paths for later steps (no PATH mutation in GITHUB_ENV)
$cl = $null
$link = $null
if ($after.ContainsKey('PATH')) {
  foreach ($dir in $after['PATH'].Split(';')) {
    if (-not $cl -and (Test-Path (Join-Path $dir 'cl.exe'))) { $cl = Join-Path $dir 'cl.exe' }
    if (-not $link -and (Test-Path (Join-Path $dir 'link.exe'))) { $link = Join-Path $dir 'link.exe' }
  }
}
if (-not $cl) { throw "cl.exe not found after vcvars" }
if (-not $link) { throw "link.exe not found after vcvars" }
Add-Content -Path $env:GITHUB_ENV -Value "EON_MSVC_CL=$cl"
Add-Content -Path $env:GITHUB_ENV -Value "EON_MSVC_LINK=$link"
Add-Content -Path $env:GITHUB_ENV -Value "EON_MSVC_BIN=$(Split-Path $cl -Parent)"
Write-Host "EON_MSVC_CL=$cl"
Write-Host "EON_MSVC_LINK=$link"
