#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline di esecuzione numerica Elmer FEM da zero (CERN-OHL-S v2):
1. Caso Nominale 30°: config/case_mesh_stirata.sif -> results_nominal_30deg/
2. Caso Sensitività 25°: config/case_mesh_25deg.sif -> results_25deg/
3. Caso Sensitività 35°: config/case_mesh_35deg.sif -> results_35deg/
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

def run_simulation(sif_rel_path, out_dir_name):
    elmer_bin = find_elmersolver()
    print(f"\n{'='*75}")
    print(f"Avvio simulazione: {sif_rel_path} -> {out_dir_name}/")
    print(f"{'='*75}")
    
    out_dir = ROOT_DIR / out_dir_name
    out_dir.mkdir(parents=True, exist_ok=True)
    existing_vtus = list(out_dir.glob("*.vtu"))
    if len(existing_vtus) == 10 and "--force" not in sys.argv:
        print(f"[SKIP] {out_dir_name} contiene già i 10 file VTU completi. (Usa --force per rieseguire)")
        return True
    for f in existing_vtus:
        f.unlink()
    
    # Aggiorna ELMERSOLVER_STARTINFO nella root
    startinfo_path = ROOT_DIR / "ELMERSOLVER_STARTINFO"
    with open(startinfo_path, "w", encoding="utf-8") as f:
        f.write(f"{sif_rel_path}\n1\n")
        
    t0 = time.time()
    res = subprocess.run([elmer_bin, sif_rel_path], cwd=ROOT_DIR, capture_output=True, text=True)
    dt = time.time() - t0
    
    if res.returncode != 0:
        print(f"[ERRORE] Simulazione {sif_rel_path} fallita con codice {res.returncode}")
        print("STDERR:\n", res.stderr[-1000:])
        print("STDOUT:\n", res.stdout[-1000:])
        return False
        
    vtus = list(out_dir.glob("*.vtu"))
    print(f"[SUCCESSO] {sif_rel_path} completata in {dt:.1f}s | Generati {len(vtus)} file VTU in {out_dir_name}/")
    if len(vtus) < 10:
        print(f"[ATTENZIONE] Attesi 10 file VTU, trovati {len(vtus)}")
        return False
    return True

def main():
    sims = [
        (os.path.join("config", "case_mesh_stirata.sif"), "results_nominal_30deg"),
        (os.path.join("config", "case_mesh_25deg.sif"), "results_25deg"),
        (os.path.join("config", "case_mesh_35deg.sif"), "results_35deg")
    ]
    
    for sif, out in sims:
        ok = run_simulation(sif, out)
        if not ok:
            sys.exit(1)
            
    print("\n" + "="*75)
    print("TUTTE LE 3 SIMULAZIONI FEM SONO STATE COMPLETATE CON SUCCESSO (CODICE 0)!")
    print("="*75)

if __name__ == "__main__":
    main()
