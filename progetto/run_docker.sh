#!/usr/bin/env bash
set -e

cd "$(dirname "$0")"

mkdir -p logs

docker rm -f maize_ncrna_gnn 2>/dev/null || true

docker build -t maize_ncrna_gnn:latest .

docker run -d \
  --name maize_ncrna_gnn \
  -p 8000:8000 \
  -p 8501:8501 \
  -v "$PWD/outputs:/app/outputs" \
  -v "$PWD/data/processed:/app/data/processed" \
  -v "$PWD/logs:/app/logs" \
  maize_ncrna_gnn:latest

echo ""
echo "Docker avviato."
echo "GUI: http://localhost:8501"
echo "API: http://localhost:8000/docs"
