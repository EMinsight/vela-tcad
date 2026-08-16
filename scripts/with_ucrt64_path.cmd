@echo off
setlocal EnableExtensions

if "%~1"=="" (
    >&2 echo with_ucrt64_path.cmd: missing compiler command
    exit /b 2
)

rem GCC loads cc1/cc1plus and their DLL dependencies from the toolchain bin
rem directory.  An absolute path to gcc.exe alone is not sufficient when the
rem parent process has not inherited the MSYS2 UCRT64 PATH.
for %%I in ("%~1") do set "VELA_UCRT64_BIN=%%~dpI"
for %%I in ("%VELA_UCRT64_BIN%..\..\usr\bin") do set "VELA_MSYS_USR_BIN=%%~fI"
set "PATH=%VELA_UCRT64_BIN%;%VELA_MSYS_USR_BIN%;%PATH%"

%*
exit /b %ERRORLEVEL%
