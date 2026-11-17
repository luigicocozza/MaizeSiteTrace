#!/usr/bin/env bash
set -e

cd "$(dirname "$0")"

stop_pid_file() {
  local file="$1"
  local name="$2"

  if [ -f "$file" ]; then
    local pid
    pid=$(cat "$file")

    if kill -0 "$pid" 2>/dev/null; then
      echo "Fermo $name PID=$pid"
      kill "$pid" || true
      sleep 1
    else
      echo "$name PID=$pid non attivo"
    fi

    rm -f "$file"
  else
    echo "Nessun PID file per $name"
  fi
}

stop_pid_file "logs/api.pid" "API"
stop_pid_file "logs/gui.pid" "GUI"

echo "Stop completato."
