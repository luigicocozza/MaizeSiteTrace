#!/usr/bin/env bash
set -e

cd /app

mkdir -p logs

echo "Avvio FastAPI su porta 8000..."
uvicorn api.main:app --host 0.0.0.0 --port 8000 > logs/docker_api.log 2>&1 &

echo "Avvio Streamlit su porta 8501..."
streamlit run app/streamlit_app.py \
  --server.address 0.0.0.0 \
  --server.port 8501 \
  --server.headless true \
  > logs/docker_gui.log 2>&1 &

echo "Servizi avviati."
echo "API: http://localhost:8000/docs"
echo "GUI: http://localhost:8501"

tail -f logs/docker_api.log logs/docker_gui.log
