#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Top-level wrapper to run 432 coils simulation
"""
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = ROOT_DIR / "variants" / "gabbia_sferica_doppio_gruppo_alta_densita_432coils" / "scripts" / "run_doppio_gruppo_432coils_simulation.py"

if __name__ == "__main__":
    import runpy
    sys.path.insert(0, str(SCRIPT_PATH.parent))
    runpy.run_path(str(SCRIPT_PATH), run_name="__main__")
