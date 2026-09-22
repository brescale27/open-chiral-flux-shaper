#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Harvesting e post-processing dei risultati parziali dello sweep 2D (Frequenza vs RPM)
interrotto dall'utente. Raccoglie tutti i dati completi e parziali simulati.
"""

import os
import sys
import json
from pathlib import Path
import numpy as np
import meshio

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
SWEEP_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0_resonance_sweep"
WORK_BASE = SWEEP_DIR / "work_dirs"
DATA_DIR = SWEEP_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

def main():
    print("=" * 85)
    print("RAGGRUPPAMENTO RISULTATI DELLO SWEEP 2D (FREQUENZA VS RPM)")
    print("=" * 85)
    
    ref_vtu = ROOT_DIR / "variants" / "rotore_centrato_z0" / "results_regime_B" / "macchina_out_t0001.vtu"
    m_ref = meshio.read(str(ref_vtu))
    pts = m_ref.points
    cells = m_ref.cells_dict['tetra']
    v0 = pts[cells[:, 0]]
    v1 = pts[cells[:, 1]]
    v2 = pts[cells[:, 2]]
    v3 = pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0
    elem_com = (v0 + v1 + v2 + v3) / 4.0
    elem_r = np.hypot(elem_com[:, 0], elem_com[:, 1])
    elem_z = np.abs(elem_com[:, 2])

    active_mask = (elem_r <= 0.051) & (elem_z <= 0.051)
    mantle_mask = (elem_r >= 0.046) & (elem_r <= 0.051) & (elem_z <= 0.051)

    results = []

    for run_dir in sorted(WORK_BASE.glob("run_f*_rpm*")):
        res_json = run_dir / "result.json"
        if res_json.is_file():
            with open(res_json, "r", encoding="utf-8") as f:
                data = json.load(f)
                data["status"] = "COMPLETO (20/20 step)"
                data["steps_computed"] = 20
                results.append(data)
        else:
            vtus = sorted(list((run_dir / "results").glob("macchina_out_t*.vtu")))
            if len(vtus) > 0:
                name = run_dir.name
                parts = name.split("_")
                f_hz = int(parts[1].replace("f", ""))
                rpm = int(parts[2].replace("rpm", ""))
                f_slip = abs(f_hz - (3.0 * rpm / 60.0))
                
                fz_series = []
                pj_series = []
                for vf in vtus:
                    m = meshio.read(str(vf))
                    jxb_z = m.point_data['jxb'][:, 2]
                    pj = m.point_data['joule heating'].ravel()
                    jxb_z_elem = np.mean(jxb_z[cells], axis=1)
                    pj_elem = np.mean(pj[cells], axis=1)
                    fz_tot = float(np.sum(elem_vols[active_mask] * jxb_z_elem[active_mask]))
                    pj_tot = float(np.sum(elem_vols[mantle_mask] * pj_elem[mantle_mask]))
                    fz_series.append(fz_tot)
                    pj_series.append(pj_tot)
                    
                fz_series = np.array(fz_series)
                pj_series = np.array(pj_series)
                mean_fz = float(np.mean(fz_series))
                peak_fz = float(np.max(fz_series))
                min_fz = float(np.min(fz_series))
                mean_pj = float(np.mean(pj_series))
                
                data = {
                    "frequency_Hz": f_hz,
                    "rpm": rpm,
                    "f_slip_Hz": f_slip,
                    "mean_Fz_mN": mean_fz * 1e3,
                    "mean_Fz_uN": mean_fz * 1e6,
                    "peak_Fz_mN": peak_fz * 1e3,
                    "peak_Fz_uN": peak_fz * 1e6,
                    "min_Fz_mN": min_fz * 1e3,
                    "min_Fz_uN": min_fz * 1e6,
                    "mean_Poule_W": mean_pj,
                    "efficiency_uN_per_W": (mean_fz * 1e6) / (mean_pj + 1e-12),
                    "steps_computed": len(vtus),
                    "status": f"AVANZATO ({len(vtus)}/20 step)"
                }
                results.append(data)

    results.sort(key=lambda r: (r["frequency_Hz"], r["rpm"]))

    # Salvataggio JSON
    out_path = DATA_DIR / "sweep_risonanza_parziale.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"Salvataggio effettuato con successo in: {out_path}\n")
    print(f"{'Frequenza [Hz]':<15} | {'RPM':<8} | {'f_slip [Hz]':<12} | {'<F_z> [uN]':<12} | {'P_Joule [W]':<12} | {'Stato':<20}")
    print("-" * 88)
    for r in results:
        print(f"{r['frequency_Hz']:<15} | {r['rpm']:<8} | {r['f_slip_Hz']:<12.1f} | {r['mean_Fz_uN']:<+12.2f} | {r['mean_Poule_W']:<12.4f} | {r['status']:<20}")
    print("-" * 88)

if __name__ == "__main__":
    main()
