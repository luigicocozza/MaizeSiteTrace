#!/usr/bin/env bash
set -e

docker rm -f maize_ncrna_gnn 2>/dev/null || true

echo "Docker fermato."
