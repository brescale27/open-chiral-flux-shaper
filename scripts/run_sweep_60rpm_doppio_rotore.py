#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convenience Runner per la root del progetto:
Esegue lo Sweep di Frequenza a 60 RPM per Assetto 2 (Gabbia Sferica Doppio Rotore Ortogonale a 90°).
"""

import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = ROOT_DIR / "variants" / "gabbia_sferica_doppio_rotore_90deg" / "scripts" / "run_sweep_60rpm_doppio_rotore.py"

if __name__ == "__main__":
    import importlib.util
    spec = importlib.util.spec_from_file_location("run_sweep_60rpm_doppio_rotore", str(SCRIPT_PATH))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["run_sweep_60rpm_doppio_rotore"] = mod
    spec.loader.exec_module(mod)
    mod.run_sweep()
