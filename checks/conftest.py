"""Lets `pytest checks` import the shared helpers the same way run_all.py does."""
import sys
from pathlib import Path

HERE = str(Path(__file__).resolve().parent)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
