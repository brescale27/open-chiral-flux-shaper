#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico Comparativo: Variante con Mantello Chiuso a Barattolo.
Include coperchio superiore (Z = +H/2) e coperchio inferiore (Z = -H/2) in rete stirata anisotropa.

Esegue:
1. Run Nominale (+30 gradi, 100 Hz, 1200 RPM, 20 timestep)
2. Run Speculare (-30 gradi, 100 Hz, 1200 RPM, 20 timestep)
3. Decoupling rigoroso del Mesh Bias (F_bias) e ricavo del vero Lift Chirale Netto (F_z,chiral)
4. Scomposizione delle perdite Joule P_J e delle forze F_z nei 4 settori:
   - Parete laterale (r in [47, 50] mm, |z| <= 50 mm)
   - Coperchio superiore (z in [50, 53] mm, r <= 50 mm)
   - Coperchio inferiore (z in [-53, -50] mm, r <= 50 mm)
   - Rotore e nucleo interno (r < 47 mm, |z| <= 50 mm)
5. Confronto 1:1 con il benchmark a tubo aperto.

Autore: Alessandro Brescacin
Licenza: CERN-OHL-S-2.0
"""

import os
import sys
import json
import time
import shutil
import subprocess
from pathlib import Path
import numpy as np
import meshio

SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_DIR = SCRIPT_DIR.parent
ROOT_DIR = VARIANT_DIR.parent.parent
CONFIG_DIR = VARIANT_DIR / "config"
DATA_DIR = VARIANT_DIR / "data"
FIGURES_DIR = VARIANT_DIR / "figures"
WORK_DIR = VARIANT_DIR / "work_dirs"
MESH_DIR = VARIANT_DIR / "mesh"

def find_elmersolver():
    cmd = os.environ.get("ELMER_SOLVER") or shutil.which("ElmerSolver")
    if cmd:
        return cmd
    candidates = [
        r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe",
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerSolver.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerSolver.exe",
        "/usr/local/bin/ElmerSolver",
        "/usr/bin/ElmerSolver",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return "ElmerSolver"

ELMER_SOLVER = find_elmersolver()
TIMESTEPS = 20
DT = 0.0005
F_HZ = 100
RPM = 1200

def run_simulation(case_name, sif_file, work_sub):
    work_sub.mkdir(parents=True, exist_ok=True)
    res_dir = work_sub / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepara SIF per esecuzione
    sif_src = CONFIG_DIR / sif_file
    sif_dest = work_sub / "case.sif"
    
    # Adegua percorso Mesh DB e Results Directory al run_dir
    mesh_db_rel = os.path.relpath(str(MESH_DIR), str(work_sub)).replace('\\', '/')
    res_dir_rel = "results"
    
    lines = sif_src.read_text(encoding="utf-8").splitlines()
    new_lines = []
    for line in lines:
        if 'Mesh DB' in line:
            new_lines.append(f'  Mesh DB "{mesh_db_rel}" "macchina_mantello_chiuso"')
        elif 'Results Directory' in line:
            new_lines.append(f'  Results Directory "{res_dir_rel}"')
        else:
            new_lines.append(line)
            
    sif_dest.write_text("\n".join(new_lines), encoding="utf-8")
    (work_sub / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")
    
    # Controlla se i VTU esistono gia
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] {case_name}: trovati {len(vtus)} VTU gia calcolati in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]
        
    print(f"  [SOLVER] Avvio ElmerSolver per {case_name} ({TIMESTEPS} timestep, 100 Hz, 1200 RPM)...")
    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    proc = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(work_sub), env=env, capture_output=True, text=True)
    duration = time.time() - t0
    
    if proc.returncode != 0:
        print(f"  [ERRORE] ElmerSolver fallito per {case_name} (codice {proc.returncode}):")
        print(proc.stderr[-1000:])
        sys.exit(1)
        
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    print(f"  [SOLVER] {case_name} completata con successo in {duration:.1f}s ({len(vtus)} VTU).")
    return vtus[:TIMESTEPS]

def evaluate_vtu_series(vtus, mesh_data):
    cells = mesh_data["cells"]
    elem_vols = mesh_data["elem_vols"]
    mask_lat = mesh_data["mask_lat"]
    mask_top = mesh_data["mask_top"]
    mask_bot = mesh_data["mask_bot"]
    mask_inner = mesh_data["mask_inner"]
    mask_can = mask_lat | mask_top | mask_bot
    mask_total = mask_can | mask_inner
    
    fz_series = {
        "lateral": [], "top_cap": [], "bottom_cap": [],
        "can_total": [], "inner_core": [], "total": []
    }
    pj_series = {
        "lateral": [], "top_cap": [], "bottom_cap": [],
        "can_total": [], "inner_core": [], "total": []
    }
    
    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        jxb_z = m.point_data['jxb'][:, 2]
        pj = m.point_data['joule heating'].ravel()
        
        jxb_z_elem = np.mean(jxb_z[cells], axis=1)
        pj_elem = np.mean(pj[cells], axis=1)
        
        # Integrazione nei 4 settori
        # 1. Parete laterale
        fz_series["lateral"].append(float(np.sum(elem_vols[mask_lat] * jxb_z_elem[mask_lat])))
        pj_series["lateral"].append(float(np.sum(elem_vols[mask_lat] * pj_elem[mask_lat])))
        
        # 2. Coperchio superiore
        fz_series["top_cap"].append(float(np.sum(elem_vols[mask_top] * jxb_z_elem[mask_top])))
        pj_series["top_cap"].append(float(np.sum(elem_vols[mask_top] * pj_elem[mask_top])))
        
        # 3. Coperchio inferiore
        fz_series["bottom_cap"].append(float(np.sum(elem_vols[mask_bot] * jxb_z_elem[mask_bot])))
        pj_series["bottom_cap"].append(float(np.sum(elem_vols[mask_bot] * pj_elem[mask_bot])))
        
        # Totale barattolo (alluminio)
        fz_series["can_total"].append(float(np.sum(elem_vols[mask_can] * jxb_z_elem[mask_can])))
        pj_series["can_total"].append(float(np.sum(elem_vols[mask_can] * pj_elem[mask_can])))
        
        # 4. Nucleo e rotore interno
        fz_series["inner_core"].append(float(np.sum(elem_vols[mask_inner] * jxb_z_elem[mask_inner])))
        pj_series["inner_core"].append(float(np.sum(elem_vols[mask_inner] * pj_elem[mask_inner])))
        
        # Totale intera macchina attiva
        fz_series["total"].append(float(np.sum(elem_vols[mask_total] * jxb_z_elem[mask_total])))
        pj_series["total"].append(float(np.sum(elem_vols[mask_total] * pj_elem[mask_total])))
        
    res = {}
    for k in fz_series:
        arr_fz = np.array(fz_series[k])
        arr_pj = np.array(pj_series[k])
        res[k] = {
            "mean_Fz_uN": float(np.mean(arr_fz) * 1e6),
            "peak_Fz_uN": float(np.max(arr_fz) * 1e6),
            "min_Fz_uN": float(np.min(arr_fz) * 1e6),
            "Fz_series_uN": (arr_fz * 1e6).tolist(),
            "mean_Pj_mW": float(np.mean(arr_pj) * 1000.0),
            "peak_Pj_mW": float(np.max(arr_pj) * 1000.0),
            "Pj_series_mW": (arr_pj * 1000.0).tolist()
        }
    return res

def main():
    print("=" * 85)
    print("STUDIO ELETTRODINAMICO: VARIANTE CON MANTELLO CHIUSO A BARATTOLO")
    print("Sanity Check Mesh Bias (+-30 gradi) e Scomposizione Forze e Perdite nei 4 Settori")
    print("=" * 85)
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    
    # Caricamento e preparazione maschere mesh
    mesh_msh = MESH_DIR / "macchina_mantello_chiuso.msh"
    if not mesh_msh.is_file():
        print("[ERRORE] Mesh macchina_mantello_chiuso.msh non trovata. Esegui prima build_mesh.py")
        sys.exit(1)
        
    m = meshio.read(str(mesh_msh))
    cells = m.cells_dict['tetra']
    pts = m.points
    v0, v1, v2, v3 = pts[cells[:, 0]], pts[cells[:, 1]], pts[cells[:, 2]], pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0
    elem_com = (v0 + v1 + v2 + v3) / 4.0
    elem_r = np.hypot(elem_com[:, 0], elem_com[:, 1])
    elem_z = elem_com[:, 2]
    
    mesh_data = {
        "cells": cells,
        "elem_vols": elem_vols,
        "mask_lat": (elem_r >= 0.0465) & (elem_r <= 0.0505) & (np.abs(elem_z) <= 0.050),
        "mask_top": (elem_z > 0.050) & (elem_z <= 0.0535) & (elem_r <= 0.0505),
        "mask_bot": (elem_z < -0.050) & (elem_z >= -0.0535) & (elem_r <= 0.0505),
        "mask_inner": (elem_r < 0.047) & (np.abs(elem_z) <= 0.050)
    }
    
    # 1. Esecuzione Run Nominale (+30 gradi)
    vtus_plus = run_simulation("Nominale (+30°)", "case_plus30deg.sif", WORK_DIR / "run_plus30deg")
    
    # 2. Esecuzione Run Speculare (-30 gradi)
    vtus_minus = run_simulation("Speculare (-30°)", "case_minus30deg.sif", WORK_DIR / "run_minus30deg")
    
    # 3. Post-processing di entrambe le serie
    print("\n[POST-PROCESSING] Integrazione settori per Nominale (+30°)...")
    res_plus = evaluate_vtu_series(vtus_plus, mesh_data)
    
    print("[POST-PROCESSING] Integrazione settori per Speculare (-30°)...")
    res_minus = evaluate_vtu_series(vtus_minus, mesh_data)
    
    # 4. Decoupling del Mesh Bias settore per settore
    decoupled = {}
    sectors = ["lateral", "top_cap", "bottom_cap", "can_total", "inner_core", "total"]
    
    for s in sectors:
        f_p = res_plus[s]["mean_Fz_uN"]
        f_m = res_minus[s]["mean_Fz_uN"]
        pj_p = res_plus[s]["mean_Pj_mW"]
        pj_m = res_minus[s]["mean_Pj_mW"]
        
        f_bias = (f_p + f_m) / 2.0
        f_chiral = (f_p - f_m) / 2.0
        pj_avg = (pj_p + pj_m) / 2.0
        
        decoupled[s] = {
            "F_plus30_uN": f_p,
            "F_minus30_uN": f_m,
            "F_bias_uN": f_bias,
            "F_chiral_uN": f_chiral,
            "Pj_plus30_mW": pj_p,
            "Pj_minus30_mW": pj_m,
            "Pj_mean_mW": pj_avg,
            "efficiency_uN_per_W": (f_chiral / (pj_p * 1e-3 + 1e-12)) if pj_p > 0 else 0.0
        }
        
    # 5. Dati benchmark aperto da validazione_mst_chiral_bias.json
    open_benchmark_file = ROOT_DIR / "variants" / "rotore_centrato_z0_resonance_sweep" / "data" / "validazione_mst_chiral_bias.json"
    open_data = None
    if open_benchmark_file.is_file():
        try:
            with open(open_benchmark_file, "r", encoding="utf-8") as f:
                open_data = json.load(f)
        except Exception:
            pass
            
    # Stampa a video dei risultati
    print("\n" + "=" * 95)
    print("SINTESI QUANTITATIVA: SCOMPOSIZIONE NEI 4 SETTORI E DECOUPLING BIAS (MANTELLO CHIUSO)")
    print("=" * 95)
    header = f"{'Settore Geometria':<22} | {'F(+30°) [µN]':>12} | {'F(-30°) [µN]':>12} | {'F_bias [µN]':>12} | {'F_chiral [µN]':>13} | {'P_J (+30°) [mW]':>15}"
    print(header)
    print("-" * 95)
    
    sector_labels = {
        "lateral": "1. Parete Laterale",
        "top_cap": "2. Coperchio Superiore",
        "bottom_cap": "3. Coperchio Inferiore",
        "can_total": "-> Totale Barattolo",
        "inner_core": "4. Rotore/Nucleo Int.",
        "total": "==> TOTALE COMPLESSIVO"
    }
    
    for s in sectors:
        d = decoupled[s]
        lbl = sector_labels[s]
        print(f"{lbl:<22} | {d['F_plus30_uN']:+12.3f} | {d['F_minus30_uN']:+12.3f} | {d['F_bias_uN']:+12.3f} | {d['F_chiral_uN']:+13.3f} | {d['Pj_plus30_mW']:15.3f}")
    print("=" * 95)
    
    # Confronto con Mantello Aperto
    print("\n" + "=" * 85)
    print("CONFRONTO DIRETTO: MANTELLO APERTO (TUBO) vs MANTELLO CHIUSO (BARATTOLO)")
    print("Regime Nominale: 100 Hz, 1200 RPM, J0 = 1.0e5 A/m^2")
    print("=" * 85)
    
    # Dati tubo aperto
    f_open_raw = +4.669
    f_open_bias = +4.917
    f_open_chiral = -0.247
    pj_open_mW = 1.524
    
    f_closed_raw = decoupled["total"]["F_plus30_uN"]
    f_closed_bias = decoupled["total"]["F_bias_uN"]
    f_closed_chiral = decoupled["total"]["F_chiral_uN"]
    pj_closed_mW = decoupled["can_total"]["Pj_plus30_mW"]
    pj_closed_tot_mW = decoupled["total"]["Pj_plus30_mW"]
    
    delta_pj_pct = ((pj_closed_tot_mW - pj_open_mW) / pj_open_mW) * 100.0
    
    print(f"Metrica / Parametro                 | Tubo Aperto (Baseline) | Barattolo Chiuso       | Variazione Delta")
    print(f"------------------------------------+------------------------+------------------------+-----------------")
    print(f"Spinta Grezza <F_z(+30°)>           | {f_open_raw:+10.3f} µN        | {f_closed_raw:+10.3f} µN        | {f_closed_raw - f_open_raw:+10.3f} µN")
    print(f"Bias Numerico Mesh F_bias           | {f_open_bias:+10.3f} µN        | {f_closed_bias:+10.3f} µN        | {f_closed_bias - f_open_bias:+10.3f} µN")
    print(f"VERO LIFT CHIRALE F_z,chiral        | {f_open_chiral:+10.3f} µN        | {f_closed_chiral:+10.3f} µN        | {f_closed_chiral - f_open_chiral:+10.3f} µN")
    print(f"Perdite Joule Mantello (mW)         | {pj_open_mW:10.3f} mW        | {pj_closed_mW:10.3f} mW        | {pj_closed_mW - pj_open_mW:+10.3f} mW")
    print(f"Perdite Joule Totali (mW)           | {pj_open_mW:10.3f} mW        | {pj_closed_tot_mW:10.3f} mW        | {delta_pj_pct:+10.1f} %")
    print(f"Contributo Coperchio Sup. (F_z)     |       N/A              | {decoupled['top_cap']['F_chiral_uN']:+10.3f} µN        | coperchio sup.")
    print(f"Contributo Coperchio Inf. (F_z)     |       N/A              | {decoupled['bottom_cap']['F_chiral_uN']:+10.3f} µN        | coperchio inf.")
    print("=" * 85)
    
    # Salvataggio JSON completo
    consolidated = {
        "metadata": {
            "title": "Caratterizzazione Elettrodinamica Variante con Mantello Chiuso a Barattolo",
            "author": "Alessandro Brescacin",
            "date": "2026-09-23",
            "license": "CERN-OHL-S-2.0",
            "geometry": "Guscio cilindrico chiuso in alluminio (parete + coperchi sup/inf t=3mm)",
            "frequency_Hz": F_HZ,
            "rpm": RPM,
            "dt_s": DT,
            "timesteps": TIMESTEPS
        },
        "nominal_plus30deg": res_plus,
        "specular_minus30deg": res_minus,
        "decoupling_by_sector": decoupled,
        "comparison_closed_vs_open": {
            "open_tube_baseline": {
                "mean_Fz_raw_uN": f_open_raw,
                "F_bias_uN": f_open_bias,
                "F_chiral_pure_uN": f_open_chiral,
                "mean_Pj_mW": pj_open_mW
            },
            "closed_can_variant": {
                "mean_Fz_raw_uN": f_closed_raw,
                "F_bias_uN": f_closed_bias,
                "F_chiral_pure_uN": f_closed_chiral,
                "mean_Pj_mantle_can_mW": pj_closed_mW,
                "mean_Pj_total_mW": pj_closed_tot_mW,
                "delta_Pj_pct": delta_pj_pct
            }
        },
        "time_series_ms": (np.arange(1, TIMESTEPS + 1) * DT * 1e3).tolist()
    }
    
    json_path = DATA_DIR / "risultati_mantello_chiuso_bias_chiral.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(consolidated, f, indent=2)
    print(f"\nDataset completo salvato in: {json_path}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
