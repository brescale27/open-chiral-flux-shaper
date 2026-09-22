#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline di esecuzione numerica Elmer FEM per lo Sweep di Sfasamento Spazio-Temporale e Modulazione AM (CERN-OHL-S v2):
1. Regime A (Sincrono / In-Phase): config/case_sweep_regime_A.sif -> results_regime_A/ (10 timesteps, dt=1.0ms)
2. Regime B (Progressivo Co-rotante 60°): config/case_sweep_regime_B.sif -> results_regime_B/ (10 timesteps, dt=1.0ms)
3. Regime C (Coppie in Quadratura 90°/180°): config/case_sweep_regime_C.sif -> results_regime_C/ (10 timesteps, dt=1.0ms)
4. Regime D (Progressivo Contro-rotante 60°): config/case_sweep_regime_D.sif -> results_regime_D/ (10 timesteps, dt=1.0ms)
5. Caso AM (Modulazione 10 Hz su portante 100 Hz): config/case_am_modulation.sif -> results_am_modulation/ (10 timesteps, dt=5.0ms)
"""

import os
import sys
import shutil
import subprocess
import time
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent

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

def run_simulation(sif_rel_path, out_dir_name, expected_steps):
    elmer_bin = find_elmersolver()
    print(f"\n{'='*75}")
    print(f"Avvio simulazione: {sif_rel_path} -> {out_dir_name}/ ({expected_steps} timesteps)")
    print(f"{'='*75}")

    out_dir = ROOT_DIR / out_dir_name
    out_dir.mkdir(parents=True, exist_ok=True)
    existing_vtus = list(out_dir.glob("*.vtu"))
    if len(existing_vtus) == expected_steps and "--force" not in sys.argv:
        print(f"[SKIP] {out_dir_name} contiene già i {expected_steps} file VTU completi. (Usa --force per rieseguire)")
        return True
    for f in existing_vtus:
        try:
            f.unlink()
        except OSError:
            pass

    startinfo_path = ROOT_DIR / "ELMERSOLVER_STARTINFO"
    with open(startinfo_path, "w", encoding="utf-8") as f:
        f.write(f"{sif_rel_path}\n1\n")

    t0 = time.time()
    res = subprocess.run([elmer_bin, sif_rel_path], cwd=ROOT_DIR, capture_output=True, text=True)
    dt = time.time() - t0

    if res.returncode != 0:
        print(f"[ERRORE] Simulazione {sif_rel_path} fallita con codice {res.returncode}")
        print("STDERR:\n", res.stderr[-1000:] if res.stderr else "Nessun output stderr")
        print("STDOUT:\n", res.stdout[-1000:] if res.stdout else "Nessun output stdout")
        return False

    vtus = list(out_dir.glob("*.vtu"))
    print(f"[SUCCESSO] {sif_rel_path} completata in {dt:.1f}s | Generati {len(vtus)} file VTU in {out_dir_name}/")
    if len(vtus) < expected_steps:
        print(f"[ATTENZIONE] Attesi {expected_steps} file VTU, trovati {len(vtus)}")
        return False
    return True

def main():
    sims = [
        (os.path.join("config", "case_sweep_regime_A.sif"), "results_regime_A", 10),
        (os.path.join("config", "case_sweep_regime_B.sif"), "results_regime_B", 10),
        (os.path.join("config", "case_sweep_regime_C.sif"), "results_regime_C", 10),
        (os.path.join("config", "case_sweep_regime_D.sif"), "results_regime_D", 10),
        (os.path.join("config", "case_am_modulation.sif"), "results_am_modulation", 10)
    ]

    for sif, out, steps in sims:
        ok = run_simulation(sif, out, steps)
        if not ok:
            sys.exit(1)

    print("\n" + "="*75)
    print("TUTTI I 4 REGIMI DI SFASAMENTO E IL CASO AM COMPLETATI CON SUCCESSO!")
    print("="*75)

if __name__ == "__main__":
    main()
