@echo off
setlocal
pushd "%~dp0" || exit /b 1
docker compose -p ai_primetta_v33_windows -f compose.yaml stop
set "RESULT=%ERRORLEVEL%"
if not "%RESULT%"=="0" pause
popd
exit /b %RESULT%
