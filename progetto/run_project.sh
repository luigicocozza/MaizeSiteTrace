#!/usr/bin/env bash
set -e

cd "$(dirname "$0")"

if [ -d ".venv" ]; then
  source .venv/bin/activate
fi

mkdir -p logs

is_url_up() {
  python - "$1" <<'PY'
import sys
import urllib.request

url = sys.argv[1]
try:
    with urllib.request.urlopen(url, timeout=2) as r:
        sys.exit(0 if r.status < 500 else 1)
except Exception:
    sys.exit(1)
PY
}

echo "============================================================"
echo "Avvio progetto maize_drought_ncrna_gnn"
echo "============================================================"

if is_url_up "http://localhost:8000/health"; then
  echo "API già attiva su http://localhost:8000"
else
  echo "Avvio API FastAPI su http://localhost:8000 ..."
  nohup uvicorn api.main:app --host 0.0.0.0 --port 8000 > logs/api.log 2>&1 &
  echo $! > logs/api.pid
fi

if is_url_up "http://localhost:8501"; then
  echo "GUI già attiva su http://localhost:8501"
else
  echo "Avvio GUI Streamlit su http://localhost:8501 ..."
  nohup streamlit run app/streamlit_app.py --server.address 0.0.0.0 --server.port 8501 > logs/gui.log 2>&1 &
  echo $! > logs/gui.pid
fi

echo ""
echo "Attesa servizi..."

for i in $(seq 1 30); do
  api_ok=0
  gui_ok=0

  if is_url_up "http://localhost:8000/health"; then
    api_ok=1
  fi

  if is_url_up "http://localhost:8501"; then
    gui_ok=1
  fi

  if [ "$api_ok" = "1" ] && [ "$gui_ok" = "1" ]; then
    echo "Servizi pronti."
    echo ""
    echo "API: http://localhost:8000/docs"
    echo "GUI: http://localhost:8501"
    echo ""
    echo "Log API: logs/api.log"
    echo "Log GUI: logs/gui.log"
    exit 0
  fi

  sleep 1
done

echo "ATTENZIONE: uno o più servizi non risultano pronti."
echo "Controlla:"
echo "  tail -n 80 logs/api.log"
echo "  tail -n 80 logs/gui.log"
exit 1
