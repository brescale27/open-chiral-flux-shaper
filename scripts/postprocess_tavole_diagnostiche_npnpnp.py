#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convenience Runner per la Generazione delle Tavole Diagnostiche a 300 DPI:
Gabbia Sferica Metamateriale a Doppio Rotore Ortogonale a 90° (NPNPNP a 120°).
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
TARGET_SCRIPT = SCRIPT_DIR / "variants" / "gabbia_sferica_doppio_rotore_90deg" / "scripts" / "postprocess_tavole_diagnostiche_npnpnp.py"

if __name__ == '__main__':
    # Esegui lo script principale
    import subprocess
    cmd = [sys.executable, str(TARGET_SCRIPT)]
    sys.exit(subprocess.call(cmd))
