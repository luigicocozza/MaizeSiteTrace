#!/usr/bin/env python3
from pathlib import Path
import torch

PYG_PATH = Path("data/processed/pyg/maize_ncrna_graph.pt")

try:
    obj = torch.load(PYG_PATH, weights_only=False)
except TypeError:
    obj = torch.load(PYG_PATH)

print("ISPEZIONE FILE PT")
print("=" * 80)
print("Path:", PYG_PATH)
print("Tipo oggetto:", type(obj))

if isinstance(obj, dict):
    print("")
    print("Chiavi dizionario:")
    for k, v in obj.items():
        print("-" * 80)
        print("Chiave:", k)
        print("Tipo:", type(v))

        if hasattr(v, "shape"):
            print("Shape:", v.shape)

        if hasattr(v, "num_nodes"):
            print("num_nodes:", v.num_nodes)

        if hasattr(v, "edge_index"):
            print("edge_index shape:", v.edge_index.shape)

        if hasattr(v, "x") and v.x is not None:
            print("x shape:", v.x.shape)

        if isinstance(v, dict):
            print("Sotto-chiavi:", list(v.keys())[:30])

        if isinstance(v, list):
            print("Lista lunghezza:", len(v))
            print("Primi elementi:", v[:5])

else:
    print("")
    print("Attributi principali:")
    for attr in ["num_nodes", "edge_index", "x", "y", "train_mask"]:
        if hasattr(obj, attr):
            val = getattr(obj, attr)
            print(attr, type(val))
            if hasattr(val, "shape"):
                print("  shape:", val.shape)

print("")
print("FINE ISPEZIONE")
