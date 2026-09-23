#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sweep Parametrico 2D (Frequenza Elettrica vs RPM Meccanici) su Griglia 5x5:
Massimizzazione della spinta assiale ponderomotrice di Lorentz (Lift Elettromagnetico)
e mappatura della superficie di risonanza dello scorrimento (Slip Resonance).

Configurazione:
- Variante Rotore Centrato a Z=0 (Mantello Anisotropo a 30°, Nucleo Equatoriale)
- Regime B co-rotante a 60° (s_k = +1, unipolare omogenea)
- Frequenze: [50, 100, 200, 400, 800] Hz
- RPM: [0, 600, 1200, 2400, 4800] RPM
- Griglia: 25 combinazioni
- Esecuzione parallela multi-core headless ElmerSolver
"""

import os
import sys
import json
import time
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import meshio

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
SWEEP_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0_resonance_sweep"
WORK_BASE = SWEEP_DIR / "work_dirs"
DATA_DIR = SWEEP_DIR / "data"
CONFIG_DIR = SWEEP_DIR / "config"
MESH_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0" / "mesh" / "macchina_centrata"

FREQUENCIES = [50, 100, 200, 400, 800]
RPMS = [0, 600, 1200, 2400, 4800]
P_POLE_PAIRS = 3
TIMESTEPS_PER_CYCLE = 20

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

def generate_sif_content(f_hz, rpm, mesh_db_dir, out_dir):
    w_e = 2.0 * np.pi * f_hz
    w_m = 2.0 * np.pi * (rpm / 60.0)
    dt = 1.0 / (TIMESTEPS_PER_CYCLE * f_hz)
    
    mesh_db_posix = mesh_db_dir.as_posix()
    out_dir_posix = out_dir.as_posix()
    
    sif = f"""! Elmer FEM: Sweep Risonanza f={f_hz} Hz, n={rpm} RPM
! Regime B co-rotante a 60 gradi
$ function rotating_current_regime_B(tx) {{\\
  x = tx(0); y = tx(1); z = tx(2); t = tx(3);\\
  res = 0.0;\\
  if (abs(z) > 0.05) {{ res = 0.0; }} else {{\\
    rsq = x*x + y*y;\\
    if (rsq < 0.0007) {{ res = 0.0; }} else {{ if (rsq > 0.0019) {{ res = 0.0; }} else {{\\
      w = {w_m:.14f};\\
      rc = 0.035; rwsq = 4.9e-05; j0 = 1.0e+05;\\
      th0 = w * t;\\
      c0 = cos(th0); s0 = sin(th0);\\
      s60 = 0.8660254037844386;\\
      x1 = rc * c0; y1 = rc * s0;\\
      x2 = rc * (0.5*c0 - s60*s0); y2 = rc * (0.5*s0 + s60*c0);\\
      x3 = rc * (-0.5*c0 - s60*s0); y3 = rc * (-0.5*s0 + s60*c0);\\
      se = sin({w_e:.14f} * t); ce = cos({w_e:.14f} * t);\\
      s1 = se; s2 = 0.5*se - s60*ce; s3 = -0.5*se - s60*ce;\\
      d1 = (x-x1)*(x-x1) + (y-y1)*(y-y1);\\
      d2 = (x-x2)*(x-x2) + (y-y2)*(y-y2);\\
      d3 = (x-x3)*(x-x3) + (y-y3)*(y-y3);\\
      d4 = (x+x1)*(x+x1) + (y+y1)*(y+y1);\\
      d5 = (x+x2)*(x+x2) + (y+y2)*(y+y2);\\
      d6 = (x+x3)*(x+x3) + (y+y3)*(y+y3);\\
      if (d1 < rwsq) {{ res = j0 * 0.5 * (s1 + abs(s1)); }} else {{\\
      if (d2 < rwsq) {{ res = j0 * 0.5 * (s2 + abs(s2)); }} else {{\\
      if (d3 < rwsq) {{ res = j0 * 0.5 * (s3 + abs(s3)); }} else {{\\
      if (d4 < rwsq) {{ res = j0 * 0.5 * (abs(s1) - s1); }} else {{\\
      if (d5 < rwsq) {{ res = j0 * 0.5 * (abs(s2) - s2); }} else {{\\
      if (d6 < rwsq) {{ res = j0 * 0.5 * (abs(s3) - s3); }} else {{ res = 0.0; }}; }}; }}; }}; }}; }};\\
    }}; }};\\
  }};\\
  _rotating_current_regime_B = res;\\
}}

$ function sigma_cyl(tx) {{\\
  x = tx(0); y = tx(1);\\
  r = sqrt(x*x + y*y);\\
  if (r < 1.0e-5) {{ r = 1.0e-5; }} else {{ r = r; }};\\
  c = x / r; s = y / r;\\
  s_rr = 1.75e6; s_tt = 1.75e6; s_zz = 1.22e7; s_tz = 3.031089e6;\\
  sxz = -s_tz * s;\\
  syz =  s_tz * c;\\
  a = zeros(3,3);\\
  a(0,0) = s_rr; a(1,1) = s_tt; a(2,2) = s_zz;\\
  a(0,2) = sxz;  a(2,0) = sxz;\\
  a(1,2) = syz;  a(2,1) = syz;\\
  _sigma_cyl = a;\\
}}

Header
  CHECK KEYWORDS Warn
  Mesh DB "{mesh_db_posix}" "macchina_centrata"
  Include Path ""
  Results Directory "{out_dir_posix}"
End

Simulation
  Max Output Level = 3
  Coordinate System = Cartesian 3D
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Transient
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
  Timestep Intervals(1) = {TIMESTEPS_PER_CYCLE}
  Timestep Sizes(1) = {dt:.12e}
End

Constants
  Permittivity of Vacuum = 8.8542e-12
  Permeability of Vacuum = 1.2566370614e-6
End

Equation 1
  Name = "Coupled Electromagnetics"
  Active Solvers(3) = 1 2 3
End

Solver 1
  Equation = "MGDynamics"
  Variable = "AV"
  Procedure = "MagnetoDynamics" "WhitneyAVSolver"
  Fix Input Current Density = Logical False
  Use Piola Transform = Logical False
  Linear System Solver = "Direct"
  Linear System Direct Method = Umfpack
End

Solver 2
  Equation = "MGDynamicsCalc"
  Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"
  Potential Variable = String "AV"
  Calculate Elemental Fields = Logical True
  Calculate Nodal Fields = Logical True
  Calculate Magnetic Flux Density = Logical True
  Calculate Magnetic Vector Potential = Logical True
  Calculate Current Density = Logical True
  Calculate Electric Field = Logical True
  Calculate Joule Heating = Logical True
  Calculate JxB = Logical True
  Linear System Solver = "Direct"
  Linear System Direct Method = Umfpack
End

Solver 3
  Equation = "ResultOutput"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = "macchina_out"
  Output Format = "vtu"
  Vtu Format = Logical True
  Save Geometry Ids = Logical True
  Discontinuous Bodies = Logical True
End

Material 1
  Name = "Air"
  Relative Permeability = 1.0
  Electric Conductivity = 0.0
End

Material 2
  Name = "ExpandedMeshMantle"
  Relative Permeability = 1.0
  Electric Conductivity(3,3) = Variable Coordinate 1, Coordinate 2
    Real MATC "sigma_cyl(tx)"
End

Material 3
  Name = "FerromagneticCore"
  Relative Permeability = 1000.0
  Electric Conductivity = 0.0
End

Body Force 1
  Current Density 3 = Variable Coordinate 1, Coordinate 2, Coordinate 3, time
    Real MATC "rotating_current_regime_B(tx)"
End

Body 1
  Target Bodies(1) = 1
  Name = "AirExterior"
  Equation = 1
  Material = 1
End

Body 2
  Target Bodies(1) = 2
  Name = "Aluminum"
  Equation = 1
  Material = 2
End

Body 3
  Target Bodies(1) = 3
  Name = "RotorAir"
  Equation = 1
  Material = 1
  Body Force = 1
End

Body 4
  Target Bodies(1) = 4
  Name = "FerromagneticCore"
  Equation = 1
  Material = 3
  Body Force = 1
End

Boundary Condition 1
  Target Boundaries(1) = 1
  Name = "FarFieldInfinity"
  AV {{e}} = Real 0.0
  AV = Real 0.0
End
"""
    return sif

def evaluate_point(params):
    f_hz, rpm, elmer_bin, mesh_data, f_bias, cached_entry = params
    f_slip = abs(f_hz - (P_POLE_PAIRS * rpm / 60.0))
    dt_val = 1.0 / (TIMESTEPS_PER_CYCLE * f_hz)
    
    run_dir = WORK_BASE / f"run_f{f_hz}_rpm{rpm}"
    res_cached = run_dir / "result.json"
    
    # 1. Se gia presente nel file result.json del run_dir
    if res_cached.is_file():
        try:
            with open(res_cached, "r", encoding="utf-8") as f:
                c_res = json.load(f)
            raw_val = c_res.get("mean_Fz_raw_uN", c_res.get("mean_Fz_uN", 0.0))
            print(f"  [CACHE] f={f_hz:3d} Hz | n={rpm:4d} RPM -> <Fz_raw> = {raw_val:+6.2f} uN")
            return c_res
        except Exception:
            pass

    # 2. Se presente nei dati consolidati parziali
    if cached_entry is not None:
        mean_raw = cached_entry.get("mean_Fz_uN", 0.0)
        mean_chiral = mean_raw - f_bias
        pj = cached_entry.get("mean_Poule_W", cached_entry.get("mean_Joule_W", 0.0))
        eff = mean_chiral / (pj + 1e-12)
        res = {
            "frequency_Hz": f_hz,
            "rpm": rpm,
            "f_slip_Hz": f_slip,
            "omega_e_rad_s": 2.0 * np.pi * f_hz,
            "omega_m_rad_s": 2.0 * np.pi * (rpm / 60.0),
            "dt_s": cached_entry.get("dt_s", dt_val),
            "mean_Fz_raw_mN": mean_raw * 1e-3,
            "mean_Fz_raw_uN": mean_raw,
            "mean_Fz_chiral_mN": mean_chiral * 1e-3,
            "mean_Fz_chiral_uN": mean_chiral,
            "mean_Fz_mN": mean_chiral * 1e-3,
            "mean_Fz_uN": mean_chiral,
            "peak_Fz_mN": cached_entry.get("peak_Fz_mN", 0.0),
            "peak_Fz_uN": cached_entry.get("peak_Fz_uN", 0.0),
            "min_Fz_mN": cached_entry.get("min_Fz_mN", 0.0),
            "min_Fz_uN": cached_entry.get("min_Fz_uN", 0.0),
            "mean_Poule_W": pj,
            "efficiency_uN_per_W": eff,
            "time_series_ms": cached_entry.get("time_series_ms", []),
            "Fz_series_mN": cached_entry.get("Fz_series_mN", []),
            "solve_duration_s": cached_entry.get("solve_duration_s", 0.0),
            "status": cached_entry.get("status", "REUSED_PARTIAL"),
            "steps_computed": cached_entry.get("steps_computed", 20)
        }
        print(f"  [REUSE] f={f_hz:3d} Hz | n={rpm:4d} RPM | f_slip={f_slip:5.1f} Hz -> <Fz_raw> = {mean_raw:+6.2f} uN | <Fz_chiral> = {mean_chiral:+6.2f} uN | P_J = {pj*1e3:6.2f} mW")
        return res

    # 3. Altrimenti esegui la simulazione con ElmerSolver
    run_dir.mkdir(parents=True, exist_ok=True)
    res_dir = run_dir / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    
    mesh_db_parent = ROOT_DIR / "variants" / "rotore_centrato_z0" / "mesh"
    sif_text = generate_sif_content(f_hz, rpm, mesh_db_parent, res_dir)
    sif_path = run_dir / "case.sif"
    with open(sif_path, "w", encoding="utf-8") as f:
        f.write(sif_text)
        
    startinfo = run_dir / "ELMERSOLVER_STARTINFO"
    with open(startinfo, "w", encoding="utf-8") as f:
        f.write("case.sif\n1\n")

    # Se nominale 100 Hz, 1200 RPM, riusa i 20 VTU gia calcolati se presenti
    if f_hz == 100 and rpm == 1200:
        reg_b_vtus = sorted(list((ROOT_DIR / "variants" / "rotore_centrato_z0" / "results_regime_B").glob("*.vtu")))
        if len(reg_b_vtus) >= TIMESTEPS_PER_CYCLE:
            for vf in reg_b_vtus[:TIMESTEPS_PER_CYCLE]:
                dest = res_dir / vf.name
                if not dest.is_file():
                    shutil.copy2(vf, dest)
            print(f"  [REUSE] f=100 Hz, n=1200 RPM: riutilizzati i 20 VTU di benchmark da results_regime_B.")

    vtus_exist = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    solve_duration = 0.0
    if len(vtus_exist) < TIMESTEPS_PER_CYCLE:
        t0 = time.time()
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "2"
        p = subprocess.run([elmer_bin, "case.sif"], cwd=str(run_dir), env=env, capture_output=True, text=True)
        solve_duration = time.time() - t0
        
        if p.returncode != 0:
            print(f"[ERRORE] f={f_hz} Hz, n={rpm} RPM fallito con codice {p.returncode}")
            return None
        
    # Elaborazione VTU generati
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    if len(vtus) < TIMESTEPS_PER_CYCLE:
        print(f"[ERRORE] Generati solo {len(vtus)}/{TIMESTEPS_PER_CYCLE} VTU per f={f_hz}, n={rpm}")
        return None
        
    cells = mesh_data["cells"]
    elem_vols = mesh_data["elem_vols"]
    active_mask = mesh_data["active_mask"]
    mantle_mask = mesh_data["mantle_mask"]
    rotor_mask = mesh_data["rotor_mask"]
    
    fz_series = []
    pj_series = []
    
    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        jxb = m.point_data['jxb']
        jxb_z = jxb[:, 2]
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
    mean_raw_uN = mean_fz * 1e6
    mean_chiral_uN = mean_raw_uN - f_bias
    peak_fz = float(np.max(fz_series))
    min_fz = float(np.min(fz_series))
    mean_pj = float(np.mean(pj_series))
    
    efficiency = mean_chiral_uN / (mean_pj + 1e-12)  # uN / W
    
    # Pulisci i VTU intermedi per risparmiare spazio su disco
    is_nominal = (f_hz == 100 and rpm == 1200)
    if not is_nominal:
        for vf in vtus:
            try:
                vf.unlink()
            except Exception:
                pass
                
    time_series = np.linspace(dt_val, TIMESTEPS_PER_CYCLE * dt_val, TIMESTEPS_PER_CYCLE)
    
    result = {
        "frequency_Hz": f_hz,
        "rpm": rpm,
        "f_slip_Hz": f_slip,
        "omega_e_rad_s": 2.0 * np.pi * f_hz,
        "omega_m_rad_s": 2.0 * np.pi * (rpm / 60.0),
        "dt_s": dt_val,
        "mean_Fz_raw_mN": mean_raw_uN * 1e-3,
        "mean_Fz_raw_uN": mean_raw_uN,
        "mean_Fz_chiral_mN": mean_chiral_uN * 1e-3,
        "mean_Fz_chiral_uN": mean_chiral_uN,
        "mean_Fz_mN": mean_chiral_uN * 1e-3,
        "mean_Fz_uN": mean_chiral_uN,
        "peak_Fz_mN": peak_fz * 1e3,
        "peak_Fz_uN": peak_fz * 1e6,
        "min_Fz_mN": min_fz * 1e3,
        "min_Fz_uN": min_fz * 1e6,
        "mean_Poule_W": mean_pj,
        "efficiency_uN_per_W": efficiency,
        "time_series_ms": (time_series * 1e3).tolist(),
        "Fz_series_mN": (fz_series * 1e3).tolist(),
        "solve_duration_s": solve_duration,
        "status": f"COMPLETO ({TIMESTEPS_PER_CYCLE}/{TIMESTEPS_PER_CYCLE} step)",
        "steps_computed": TIMESTEPS_PER_CYCLE
    }
    
    try:
        with open(res_cached, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
    except Exception:
        pass
        
    print(f"  [DONE] f={f_hz:3d} Hz | n={rpm:4d} RPM | f_slip={f_slip:5.1f} Hz -> <Fz_raw> = {mean_raw_uN:+6.2f} uN | <Fz_chiral> = {mean_chiral_uN:+6.2f} uN | P_J = {mean_pj*1e3:6.2f} mW | {solve_duration:.1f}s")
    return result

def main():
    print("=" * 85)
    print("FASE 2: SWEEP PARAMETRICO 2D (FREQUENZA VS RPM) SU GRIGLIA 5x5")
    print("Variante Rotore Centrato a Z=0 - Decoupling Mesh Bias e Mappatura Risonanza")
    print("=" * 85)
    
    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    WORK_BASE.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    
    elmer_bin = find_elmersolver()
    print(f"Solutore individuato: {elmer_bin}")
    
    # 1. Carica Mesh Bias da validazione Fase 1
    bias_json = DATA_DIR / "validazione_mst_chiral_bias.json"
    f_bias = 0.0
    if bias_json.is_file():
        try:
            with open(bias_json, "r", encoding="utf-8") as f:
                b_data = json.load(f)
            f_bias = b_data.get("bias_decoupling", {}).get("total_active_domain", {}).get("F_bias_uN", 0.0)
            print(f"Mesh Bias caricato da validazione Fase 1: F_bias = {f_bias:+.3f} µN")
        except Exception as e:
            print(f"[WARN] Impossibile leggere F_bias da {bias_json}: {e}")
    else:
        print("[WARN] validazione_mst_chiral_bias.json non trovato, F_bias impostato a 0.0 µN")
        
    # 2. Carica punti parziali gia noti
    partial_json = DATA_DIR / "sweep_risonanza_parziale.json"
    partial_map = {}
    if partial_json.is_file():
        try:
            with open(partial_json, "r", encoding="utf-8") as f:
                p_data = json.load(f)
            for c in p_data.get("configurations", []):
                partial_map[(c["frequency_Hz"], c["rpm"])] = c
            print(f"Caricati {len(partial_map)} punti pre-calcolati da sweep_risonanza_parziale.json")
        except Exception as e:
            print(f"[WARN] Impossibile leggere dati parziali: {e}")
            
    # 3. Pre-caricamento mesh
    ref_vtu = ROOT_DIR / "variants" / "rotore_centrato_z0" / "results_regime_B" / "macchina_out_t0001.vtu"
    if not ref_vtu.is_file():
        ref_vtu = ROOT_DIR / "variants" / "rotore_centrato_poli_alternati_semionda" / "results" / "macchina_out_t0001.vtu"
        
    print(f"Caricamento topologia mesh da: {ref_vtu}")
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
    
    mesh_data = {
        "cells": cells,
        "elem_vols": elem_vols,
        "active_mask": (elem_r <= 0.051) & (elem_z <= 0.051),
        "mantle_mask": (elem_r >= 0.046) & (elem_r <= 0.051) & (elem_z <= 0.051),
        "rotor_mask": (elem_r < 0.046) & (elem_z <= 0.051)
    }
    
    # 4. Preparazione task della matrice 5x5
    tasks = []
    for f in FREQUENCIES:
        for n in RPMS:
            c_entry = partial_map.get((f, n), None)
            tasks.append((f, n, elmer_bin, mesh_data, f_bias, c_entry))
            
    n_cached = sum(1 for t in tasks if t[5] is not None)
    n_to_run = len(tasks) - n_cached
    print(f"\nPunti della griglia 5x5: {len(tasks)} totali ({n_cached} riutilizzati da cache, {n_to_run} da calcolare)")
    print(f"Avvio elaborazione parallela con 4 worker...")
    t_start = time.time()
    
    matrix_results = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(evaluate_point, t): (t[0], t[1]) for t in tasks}
        for fut in as_completed(futures):
            f_val, n_val = futures[fut]
            try:
                res = fut.result()
                if res is not None:
                    matrix_results.append(res)
            except Exception as e:
                print(f"[ECCEZIONE] Punto f={f_val}, n={n_val}: {e}")
                
    total_time = time.time() - t_start
    print("\n" + "=" * 85)
    print(f"Sweep 2D completato in {total_time:.1f}s ({total_time/60:.2f} minuti).")
    print(f"Combinazioni calcolate con successo: {len(matrix_results)}/25")
    print("=" * 85)
    
    matrix_results.sort(key=lambda r: (r["frequency_Hz"], r["rpm"]))
    
    # Individua punto di picco e benchmark nominale
    peak_point = max(matrix_results, key=lambda r: r["mean_Fz_chiral_uN"])
    nominal_point = next((r for r in matrix_results if r["frequency_Hz"] == 100 and r["rpm"] == 1200), None)
    
    boost_pct = 0.0
    if nominal_point and nominal_point["mean_Fz_chiral_uN"] != 0:
        boost_pct = ((peak_point["mean_Fz_chiral_uN"] - nominal_point["mean_Fz_chiral_uN"]) / abs(nominal_point["mean_Fz_chiral_uN"])) * 100.0
        
    consolidated_output = {
        "metadata": {
            "title": "Sweep Parametrico 2D Frequenza vs RPM - Matrice Completa Consolidata",
            "author": "Alessandro Brescacin",
            "date": "2026-09-23",
            "license": "CERN-OHL-S-2.0",
            "variant": "rotore_centrato_z0",
            "mantle": "rete stirata romboidale alluminio 30 gradi MATC cilindrico",
            "core": "nucleo ferromagnetico equatoriale Z=0 (t=1.5cm)",
            "mesh_bias_correction_uN": f_bias,
            "total_points": len(matrix_results)
        },
        "key_findings": {
            "peak_operating_point": {
                "frequency_Hz": peak_point["frequency_Hz"],
                "rpm": peak_point["rpm"],
                "f_slip_Hz": peak_point["f_slip_Hz"],
                "mean_Fz_raw_uN": peak_point["mean_Fz_raw_uN"],
                "mean_Fz_chiral_uN": peak_point["mean_Fz_chiral_uN"],
                "peak_Fz_uN": peak_point["peak_Fz_uN"],
                "mean_Joule_W": peak_point["mean_Poule_W"],
                "efficiency_uN_per_W": peak_point["efficiency_uN_per_W"],
                "boost_vs_nominal_percentage": boost_pct
            },
            "nominal_benchmark": {
                "frequency_Hz": 100,
                "rpm": 1200,
                "f_slip_Hz": 40.0,
                "mean_Fz_raw_uN": nominal_point["mean_Fz_raw_uN"] if nominal_point else None,
                "mean_Fz_chiral_uN": nominal_point["mean_Fz_chiral_uN"] if nominal_point else None,
                "mean_Joule_W": nominal_point["mean_Poule_W"] if nominal_point else None
            }
        },
        "full_scale_1e7_projection": {
            "scale_factor_current": 100.0,
            "scale_factor_Force": 10000.0,
            "projected_peak_mean_Fz_mN": peak_point["mean_Fz_chiral_uN"] * 1e-3 * 10000.0 / 1e3,
            "projected_nominal_mean_Fz_mN": (nominal_point["mean_Fz_chiral_uN"] * 1e-3 * 10000.0 / 1e3) if nominal_point else None
        },
        "configurations": matrix_results
    }
    
    matrice_path = DATA_DIR / "sweep_risonanza_matrice.json"
    with open(matrice_path, "w", encoding="utf-8") as f:
        json.dump(consolidated_output, f, indent=2)
    print(f"\nMatrice consolidata salvata in: {matrice_path}")
    
    # Stampa tabella riepilogativa ASCII
    print("\n" + "=" * 95)
    print(f"{'f [Hz]':>8} | {'RPM':>6} | {'f_slip [Hz]':>12} | {'<Fz_raw> [µN]':>14} | {'<Fz_chiral> [µN]':>16} | {'P_J [mW]':>10} | {'Efficienza [µN/W]':>18}")
    print("-" * 95)
    for r in matrix_results:
        print(f"{r['frequency_Hz']:8d} | {r['rpm']:6d} | {r['f_slip_Hz']:12.1f} | {r['mean_Fz_raw_uN']:+14.2f} | {r['mean_Fz_chiral_uN']:+16.2f} | {r['mean_Poule_W']*1000.0:10.2f} | {r['efficiency_uN_per_W']:18.1f}")
    print("=" * 95)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
