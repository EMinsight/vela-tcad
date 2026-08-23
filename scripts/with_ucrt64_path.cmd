@echo off
setlocal EnableExtensions

if "%~1"=="" (
    >&2 echo with_ucrt64_path.cmd: missing compiler command
    exit /b 2
)

rem UCRT64-built compilers and applications need their runtime DLLs on PATH.
rem Use the repository's documented MSYS2 installation instead of inferring
rem the toolchain from the target executable (which may live in build-release).
set "VELA_UCRT64_BIN=D:\msys64\ucrt64\bin"
set "VELA_MSYS_USR_BIN=D:\msys64\usr\bin"
if not exist "%VELA_UCRT64_BIN%\libgcc_s_seh-1.dll" (
    >&2 echo with_ucrt64_path.cmd: UCRT64 runtime not found at %VELA_UCRT64_BIN%
    exit /b 3
)
set "PATH=%VELA_UCRT64_BIN%;%VELA_MSYS_USR_BIN%;%PATH%"

%*
exit /b %ERRORLEVEL%
