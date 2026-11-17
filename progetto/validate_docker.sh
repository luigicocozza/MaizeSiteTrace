#!/usr/bin/env bash
set -e

echo "Controllo container..."
docker ps --filter "name=maize_ncrna_gnn"

echo ""
echo "Attesa avvio servizi..."
sleep 8

echo ""
echo "Test API health..."
curl -s http://localhost:8000/health | python -m json.tool

echo ""
echo "Test API summary..."
curl -s http://localhost:8000/summary | python -m json.tool

echo ""
echo "Test RAG API..."
curl -s -X POST http://localhost:8000/rag/ask \
  -H "Content-Type: application/json" \
  -d '{"question":"Quanti ncRNA sono stati analizzati?","provider":"extractive","top_k":6}' \
  | python -m json.tool

echo ""
echo "Test GUI HTTP..."
python - <<'PY'
import urllib.request

url = "http://localhost:8501"
with urllib.request.urlopen(url, timeout=10) as r:
    print("GUI status:", r.status)
    print("GUI OK:", r.status < 500)
PY

echo ""
echo "Validazione Docker completata."
