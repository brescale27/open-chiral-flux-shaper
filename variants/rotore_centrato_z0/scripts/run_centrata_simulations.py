#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline di esecuzione numerica Elmer FEM per la Variante con Rotore Centrato a Z=0:
1. Simulazione Sincrona (Baseline Centrata, 10 timestep, T=10ms, dt=1.0ms):
   config/case_centrata_sincrono.sif -> results_sincrono/
2. Simulazione Regime B (Co-rotante 60°, 20 timestep, T=10ms, dt=0.5ms):
   config/case_centrata_regime_B.sif -> results_regime_B/
"""

import os
import sys
import shutil
import subprocess
import time
from pathlib import Path

# ROOT_DIR punta a simulazione/
ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
VARIANTS_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0"

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

def run_simulation(sif_rel_path, out_dir_name, expected_vtus=10):
    elmer_bin = find_elmersolver()
    print(f"\n{'='*75}")
    print(f"Avvio simulazione: {sif_rel_path} -> {out_dir_name}/")
    print(f"{'='*75}")
    
    out_dir = ROOT_DIR / out_dir_name
    out_dir.mkdir(parents=True, exist_ok=True)
    existing_vtus = list(out_dir.glob("*.vtu"))
    if len(existing_vtus) == expected_vtus and "--force" not in sys.argv:
        print(f"[SKIP] {out_dir_name} contiene già i {expected_vtus} file VTU completi. (Usa --force per rieseguire)")
        return True
        
    for f in existing_vtus:
        try:
            f.unlink()
        except Exception:
            pass
            
    # Aggiorna ELMERSOLVER_STARTINFO nella root di simulazione
    startinfo_path = ROOT_DIR / "ELMERSOLVER_STARTINFO"
    with open(startinfo_path, "w", encoding="utf-8") as f:
        f.write(f"{sif_rel_path}\n1\n")
        
    t0 = time.time()
    res = subprocess.run([elmer_bin, sif_rel_path], cwd=ROOT_DIR, capture_output=True, text=True)
    dt = time.time() - t0
    
    if res.returncode != 0:
        print(f"[ERRORE] Simulazione {sif_rel_path} fallita con codice {res.returncode}")
        print("STDERR:\n", res.stderr[-1500:] if res.stderr else "None")
        print("STDOUT:\n", res.stdout[-1500:] if res.stdout else "None")
        return False
        
    vtus = list(out_dir.glob("*.vtu"))
    print(f"[COMPLETATO] {sif_rel_path} in {dt:.1f}s - File VTU generati: {len(vtus)}/{expected_vtus}")
    if len(vtus) < expected_vtus:
        print(f"[ATTENZIONE] Attesi {expected_vtus} VTU, trovati {len(vtus)}")
        return False
        
    return True

def main():
    print("===========================================================================")
    print("Suite di Simulazione FEM: Variante Rotore Centrato a Z=0 (CERN-OHL-S v2)")
    print("===========================================================================")
    
    runs = [
        ("variants/rotore_centrato_z0/config/case_centrata_sincrono.sif", "variants/rotore_centrato_z0/results_sincrono", 10),
        ("variants/rotore_centrato_z0/config/case_centrata_regime_B.sif", "variants/rotore_centrato_z0/results_regime_B", 20),
    ]
    
    success = True
    for sif, out_dir, count in runs:
        ok = run_simulation(sif, out_dir, count)
        if not ok:
            print(f"\n[INTERRUZIONE] Fallimento su {sif}")
            success = False
            break
            
    if success:
        print("\n" + "="*75)
        print("TUTTE LE SIMULAZIONI DELLA VARIANTE CENTRATA SONO STATE COMPLETATE!")
        print("===========================================================================")
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
