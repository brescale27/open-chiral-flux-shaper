#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner per la Simulazione Elmer FEM:
Variante con Rotore Centrato a Z=0 e 6 Bobine a Polarità Alternate Specchiate (N-S-N-S-N-S)
con pilotaggio a Semionde Pulsate (40 timestep, T=20 ms, dt=0.5 ms).
"""

import os
import sys
import shutil
import subprocess
import time
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
VARIANT_DIR = ROOT_DIR / "variants" / "rotore_centrato_poli_alternati_semionda"
CONFIG_SIF = VARIANT_DIR / "config" / "case_poli_alternati_semionda.sif"
RESULTS_DIR = VARIANT_DIR / "results"
EXPECTED_VTUS = 40

def find_elmersolver():
    cmd = shutil.which("ElmerSolver")
    if cmd:
        return cmd
    candidates = [
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerSolver.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerSolver.exe"
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return "ElmerSolver"

def main():
    print("=" * 80)
    print("Simulazione Elmer FEM: Variante Polarita Alternate Specchiate a Semionde Pulsate")
    print(f"Directory Root: {ROOT_DIR}")
    print(f"SIF Config:     {CONFIG_SIF}")
    print(f"Output:         {RESULTS_DIR}")
    print(f"Timesteps:      {EXPECTED_VTUS} (dt = 0.5 ms, T = 20 ms)")
    print("=" * 80)

    elmer_bin = find_elmersolver()
    print(f"Solutore individuato: {elmer_bin}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    existing_vtus = list(RESULTS_DIR.glob("*.vtu"))
    if len(existing_vtus) == EXPECTED_VTUS and "--force" not in sys.argv:
        print(f"[SKIP] Trovati già {EXPECTED_VTUS} file VTU completi in {RESULTS_DIR}. Usa --force per rieseguire.")
        return 0

    for f in existing_vtus:
        try:
            f.unlink()
        except Exception:
            pass

    startinfo_path = ROOT_DIR / "ELMERSOLVER_STARTINFO"
    rel_sif = CONFIG_SIF.relative_to(ROOT_DIR)
    with open(startinfo_path, "w", encoding="utf-8") as f:
        f.write(f"{rel_sif.as_posix()}\n1\n")

    print(f"\n[AVVIO] Esecuzione ElmerSolver {rel_sif.as_posix()}...")
    t0 = time.time()
    
    # Esecuzione con cattura o streaming
    process = subprocess.Popen(
        [elmer_bin, str(rel_sif.as_posix())],
        cwd=str(ROOT_DIR),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    )

    last_step = 0
    for line in iter(process.stdout.readline, ''):
        if "Time:" in line or "MAIN: Current time" in line or "Timestep:" in line:
            print(line.strip())
            sys.stdout.flush()
        elif "WhitneyAVSolver" in line and "Solve done" in line:
            last_step += 1
            elapsed = time.time() - t0
            print(f"  -> Timestep {last_step}/{EXPECTED_VTUS} completato ({elapsed:.1f}s trascorsi)")
            sys.stdout.flush()

    process.stdout.close()
    return_code = process.wait()
    total_time = time.time() - t0

    if return_code != 0:
        print(f"\n[ERRORE] Simulazione fallita con codice {return_code} dopo {total_time:.1f}s")
        return return_code

    vtus = list(RESULTS_DIR.glob("*.vtu"))
    print("\n" + "=" * 80)
    print(f"[SUCCESSO] Simulazione terminata in {total_time:.1f}s ({total_time/60:.2f} min).")
    print(f"File VTU generati: {len(vtus)}/{EXPECTED_VTUS}")
    print("=" * 80)

    if len(vtus) < EXPECTED_VTUS:
        print(f"[ATTENZIONE] Trovati {len(vtus)} VTU invece dei {EXPECTED_VTUS} attesi.")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
