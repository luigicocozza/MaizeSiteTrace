@echo off
setlocal
pushd "%~dp0" || exit /b 1
docker info --format "{{.OSType}}" 2>nul | findstr /x linux >nul
if errorlevel 1 goto errore
docker compose -p ai_primetta_v33_windows -f compose.yaml config -q
if errorlevel 1 goto errore
for /f "usebackq delims=" %%I in ("images.txt") do (
    docker image inspect "%%I" >nul 2>nul
    if errorlevel 1 (
        docker image load -i "docker_images.tar.gz"
        if errorlevel 1 goto errore
    )
)
docker compose -p ai_primetta_v33_windows -f compose.yaml up -d --no-build --pull never
if errorlevel 1 goto errore
powershell.exe -NoProfile -Command "$ErrorActionPreference='Stop'; $urls=@('http://127.0.0.1:8503','http://127.0.0.1:8502'); foreach ($url in $urls) { $ok=$false; for ($i=0; $i -lt 45; $i++) { try { $r=Invoke-WebRequest ($url+'/_stcore/health') -UseBasicParsing -TimeoutSec 2; if ($r.StatusCode -eq 200 -and $r.Content.Trim() -eq 'ok') { $ok=$true; break } } catch {}; Start-Sleep -Seconds 2 }; if (-not $ok) { Write-Host ('GUI non pronta: '+$url); exit 1 } }; foreach ($url in $urls) { Start-Process $url }"
if errorlevel 1 goto errore
popd
exit /b 0
:errore
echo [ERRORE] Avvia Docker Desktop in modalita Linux e controlla i log.
docker compose -p ai_primetta_v33_windows -f compose.yaml ps
docker compose -p ai_primetta_v33_windows -f compose.yaml logs --tail 30
pause
popd
exit /b 1
