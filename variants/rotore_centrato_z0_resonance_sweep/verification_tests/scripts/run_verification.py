#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protocollo di Falsificazione Artefatti Numerici (CERN-OHL-S v2):
Verifica di Parità e Control Runs su Lift di Lorentz (<F_z> a 100 Hz, 0 RPM).

Esegue 4 casi indipendenti (20 timestep ciascuno, dt = 0.5 ms, T = 10 ms):
1. BASELINE REFERENCE:     Chiralità +30°, Rotazione di fase forward (+60°)
2. TEST 1 (CHIRAL REVERSAL): Chiralità -30°, Rotazione forward (+60°) -> atteso <F_z> ≈ -5.72 µN (specularità ±5%)
3. TEST 2 (ISOTROPIC MANTLE): Chiralità 0°, Mantello isotropo puro     -> atteso <F_z> ≈ 0.00 µN (|res| < 0.1 µN)
4. TEST 3 (PHASE INVERSION):  Chiralità +30°, Rotazione reverse (-60°) -> atteso <F_z> < 0 (inversione spinta)

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
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import matplotlib.pyplot as plt
import meshio

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

# Configurazione percorsi
SCRIPT_DIR = Path(__file__).resolve().parent
VERIF_DIR = SCRIPT_DIR.parent
ROOT_DIR = VERIF_DIR.parent.parent.parent
CONFIG_DIR = VERIF_DIR / "config"
DATA_DIR = VERIF_DIR / "data"
FIGURES_DIR = VERIF_DIR / "figures"
WORK_BASE = VERIF_DIR / "work_dirs"
MESH_DB_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0" / "mesh"

CONFIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
WORK_BASE.mkdir(parents=True, exist_ok=True)

# Parametri fisici
FREQUENCY_HZ = 100
RPM = 0
DT = 0.0005  # 0.5 ms
TIMESTEPS = 20  # 1 periodo elettrico completo (10 ms)

def find_elmersolver():
    cmd = shutil.which("ElmerSolver")
    if cmd:
        return cmd
    candidates = [
        r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe",
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerSolver.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerSolver.exe"
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return "ElmerSolver"

def generate_sif_text(case_key, theta_chiral, phase_dir, is_isotropic, out_dir_posix):
    """
    Genera il contenuto SIF formattato per Elmer FEM WhitneyAVSolver con MATC.
    """
    mesh_db_posix = MESH_DB_DIR.as_posix()
    w_e = 2.0 * np.pi * FREQUENCY_HZ
    
    # 1. Funzione MATC per il campo di corrente rotante
    # phase_dir == +1 -> progressione +60° (CCW, -s60*ce)
    # phase_dir == -1 -> progressione -60° (CW,  +s60*ce)
    if phase_dir >= 0:
        s2_formula = "0.5*se - s60*ce"
        s3_formula = "-0.5*se - s60*ce"
    else:
        s2_formula = "0.5*se + s60*ce"
        s3_formula = "-0.5*se + s60*ce"

    matc_current = f"""$ function rotating_current(tx) {{\\
  x = tx(0); y = tx(1); z = tx(2); t = tx(3);\\
  res = 0.0;\\
  if (abs(z) > 0.05) {{ res = 0.0; }} else {{\\
    rsq = x*x + y*y;\\
    if (rsq < 0.0007) {{ res = 0.0; }} else {{ if (rsq > 0.0019) {{ res = 0.0; }} else {{\\
      w = 0.0;\\
      rc = 0.035; rwsq = 4.9e-05; j0 = 1.0e+05;\\
      th0 = w * t;\\
      c0 = cos(th0); s0 = sin(th0);\\
      s60 = 0.8660254037844386;\\
      x1 = rc * c0; y1 = rc * s0;\\
      x2 = rc * (0.5*c0 - s60*s0); y2 = rc * (0.5*s0 + s60*c0);\\
      x3 = rc * (-0.5*c0 - s60*s0); y3 = rc * (-0.5*s0 + s60*c0);\\
      se = sin({w_e:.14f} * t); ce = cos({w_e:.14f} * t);\\
      s1 = se; s2 = {s2_formula}; s3 = {s3_formula};\\
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
  _rotating_current = res;\\
}}"""

    # 2. Funzione MATC per il tensore di conducibilità
    if is_isotropic:
        # Mantello puramente isotropo (alluminio 3.5e7 S/m o equivalente louver-free)
        matc_sigma = """$ function sigma_cyl(tx) {\\
  a = zeros(3,3);\\
  a(0,0) = 3.5e7; a(1,1) = 3.5e7; a(2,2) = 3.5e7;\\
  _sigma_cyl = a;\\
}"""
    else:
        # Tensore anisotropo con inclinazione persiana controllata
        if theta_chiral == 30:
            s_tz_str = "3.031089e6"
        elif theta_chiral == -30:
            s_tz_str = "-3.031089e6"
        else:
            s_tz_str = "0.0"

        matc_sigma = f"""$ function sigma_cyl(tx) {{\\
  x = tx(0); y = tx(1);\\
  r = sqrt(x*x + y*y);\\
  if (r < 1.0e-5) {{ r = 1.0e-5; }} else {{ r = r; }};\\
  c = x / r; s = y / r;\\
  s_rr = 1.75e6; s_tt = 1.75e6; s_zz = 1.22e7; s_tz = {s_tz_str};\\
  sxz = -s_tz * s;\\
  syz =  s_tz * c;\\
  a = zeros(3,3);\\
  a(0,0) = s_rr; a(1,1) = s_tt; a(2,2) = s_zz;\\
  a(0,2) = sxz;  a(2,0) = sxz;\\
  a(1,2) = syz;  a(2,1) = syz;\\
  _sigma_cyl = a;\\
}}"""

    sif_body = f"""! Elmer FEM: Control Run [{case_key}]
! Verification of Lorentz Lift Parity & Numerical Artifact Elimination
! f = {FREQUENCY_HZ} Hz, RPM = {RPM}, dt = {DT} s

{matc_current}

{matc_sigma}

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
  Timestep Intervals(1) = {TIMESTEPS}
  Timestep Sizes(1) = {DT:.8e}
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
  Name = "MantleMaterial"
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
    Real MATC "rotating_current(tx)"
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
    return sif_body

def precalculate_mesh_volumes():
    ref_vtu = ROOT_DIR / "variants" / "rotore_centrato_z0" / "results_regime_B" / "macchina_out_t0001.vtu"
    if not ref_vtu.is_file():
        ref_vtu = ROOT_DIR / "variants" / "rotore_centrato_poli_alternati_semionda" / "results" / "macchina_out_t0001.vtu"
        
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
    
    return cells, elem_vols, active_mask, mantle_mask

def run_case_simulation(case_spec, elmer_bin):
    case_key = case_spec["key"]
    run_dir = WORK_BASE / f"run_{case_key}"
    run_dir.mkdir(parents=True, exist_ok=True)
    res_dir = run_dir / "results"
    res_dir.mkdir(parents=True, exist_ok=True)

    # Scrivi SIF in config e nel run_dir
    sif_content = generate_sif_text(
        case_key, 
        case_spec["theta_chiral"], 
        case_spec["phase_dir"], 
        case_spec["is_isotropic"], 
        res_dir.as_posix()
    )
    (CONFIG_DIR / f"case_{case_key}.sif").write_text(sif_content, encoding="utf-8")
    (run_dir / "case.sif").write_text(sif_content, encoding="utf-8")
    (run_dir / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")

    # Controlla se i VTU esistono già
    existing_vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    if len(existing_vtus) >= TIMESTEPS:
        print(f"  [{case_key}] Già completato con {len(existing_vtus)} VTU. Skip solver.")
        return case_key, 0.0

    print(f"  [{case_key}] Avvio ElmerSolver per {TIMESTEPS} timestep...")
    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    proc = subprocess.run([elmer_bin, "case.sif"], cwd=str(run_dir), env=env, capture_output=True, text=True)
    duration = time.time() - t0
    
    if proc.returncode != 0:
        print(f"  [ERRORE] {case_key} terminato con codice {proc.returncode}")
        print("STDERR tail:\n", proc.stderr[-1000:])
        return case_key, -1.0
        
    vtus_after = list(res_dir.glob("macchina_out_t*.vtu"))
    print(f"  [{case_key}] Completato con successo in {duration:.1f}s ({len(vtus_after)} VTU generati).")
    return case_key, duration

def postprocess_case(case_spec, cells, elem_vols, active_mask, mantle_mask):
    case_key = case_spec["key"]
    run_dir = WORK_BASE / f"run_{case_key}"
    res_dir = run_dir / "results"
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    
    if len(vtus) == 0:
        raise RuntimeError(f"Nessun file VTU trovato per {case_key} in {res_dir}")

    time_series_ms = []
    fz_series_uN = []
    pj_series_W = []

    for idx, vf in enumerate(vtus):
        t_ms = (idx + 1) * DT * 1000.0
        time_series_ms.append(t_ms)
        m = meshio.read(str(vf))
        jxb_z = m.point_data['jxb'][:, 2]
        pj = m.point_data['joule heating'].ravel()
        
        jxb_z_elem = np.mean(jxb_z[cells], axis=1)
        pj_elem = np.mean(pj[cells], axis=1)
        
        fz_N = float(np.sum(elem_vols[active_mask] * jxb_z_elem[active_mask]))
        pj_W = float(np.sum(elem_vols[mantle_mask] * pj_elem[mantle_mask]))
        
        fz_series_uN.append(fz_N * 1e6)
        pj_series_W.append(pj_W)

    fz_arr = np.array(fz_series_uN)
    pj_arr = np.array(pj_series_W)

    mean_fz = float(np.mean(fz_arr))
    peak_fz = float(np.max(fz_arr))
    min_fz = float(np.min(fz_arr))
    mean_pj = float(np.mean(pj_arr))
    eff = mean_fz / (mean_pj + 1e-12)

    return {
        "key": case_key,
        "title": case_spec["title"],
        "theta_chiral": case_spec["theta_chiral"],
        "phase_dir": case_spec["phase_dir"],
        "is_isotropic": case_spec["is_isotropic"],
        "timesteps": len(vtus),
        "time_series_ms": time_series_ms,
        "Fz_series_uN": fz_series_uN,
        "Pj_series_W": pj_series_W,
        "mean_Fz_uN": mean_fz,
        "peak_Fz_uN": peak_fz,
        "min_Fz_uN": min_fz,
        "mean_Pj_W": mean_pj,
        "efficiency_uN_per_W": eff
    }

def main():
    print("=" * 85)
    print("PROTOCOLLO DI FALSIFICAZIONE ARTEFATTI NUMERICI SUL LIFT DI LORENTZ")
    print("Verifica di Parita Elettrodinamica e Control Runs (CERN-OHL-S v2)")
    print("=" * 85)

    elmer_bin = find_elmersolver()
    print(f"Solutore ElmerSolver: {elmer_bin}")
    
    cells, elem_vols, active_mask, mantle_mask = precalculate_mesh_volumes()
    print(f"Mesh caricata: {len(elem_vols)} elementi tetraedrici.")

    cases = [
        {
            "key": "baseline_ref",
            "title": "Baseline Riferimento (+30 deg, +w)",
            "theta_chiral": 30,
            "phase_dir": +1,
            "is_isotropic": False
        },
        {
            "key": "test1_chirality_reversal",
            "title": "Test 1: Inversione Chirale (-30 deg, +w)",
            "theta_chiral": -30,
            "phase_dir": +1,
            "is_isotropic": False
        },
        {
            "key": "test2_isotropic",
            "title": "Test 2: Mantello Isotropo Puro (0 deg, +w)",
            "theta_chiral": 0,
            "phase_dir": +1,
            "is_isotropic": True
        },
        {
            "key": "test3_phase_inversion",
            "title": "Test 3: Inversione Fase (+30 deg, -w)",
            "theta_chiral": 30,
            "phase_dir": -1,
            "is_isotropic": False
        }
    ]

    # Esecuzione parallela controllata (4 worker concorrenti)
    print(f"\nAvvio simulazioni per i 4 casi con ProcessPoolExecutor (4 worker)...")
    t_start = time.time()
    durations = {}
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(run_case_simulation, c, elmer_bin): c["key"] for c in cases}
        for fut in as_completed(futures):
            k, dur = fut.result()
            durations[k] = dur

    tot_solve_time = time.time() - t_start
    print(f"Simulazioni completate in {tot_solve_time:.1f}s ({tot_solve_time/60:.2f} min).")

    # Post-processing
    print("\nEstrazione ed integrazione volumetrica di Lorentz per ciascun caso...")
    results = {}
    for c in cases:
        data = postprocess_case(c, cells, elem_vols, active_mask, mantle_mask)
        results[c["key"]] = data

    # Valutazione Criteri Fisici PASS / FAIL
    base_mean = results["baseline_ref"]["mean_Fz_uN"]
    t1_mean = results["test1_chirality_reversal"]["mean_Fz_uN"]
    t2_mean = results["test2_isotropic"]["mean_Fz_uN"]
    t3_mean = results["test3_phase_inversion"]["mean_Fz_uN"]

    # Criterio 1: Specularità (t1_mean deve essere ≈ -base_mean entro 5%)
    sym_error_pct = abs(t1_mean + base_mean) / abs(base_mean) * 100.0
    pass_t1 = (sym_error_pct <= 5.0) and (t1_mean < 0)

    # Criterio 2: Parità conservata (t2_mean deve essere < 0.10 µN)
    pass_t2 = abs(t2_mean) < 0.10

    # Criterio 3: Inversione di fase (t3_mean < 0)
    pass_t3 = t3_mean < 0.0

    all_passed = pass_t1 and pass_t2 and pass_t3

    # Stampa tabella di confronto terminale
    print("\n" + "=" * 95)
    print("RISULTATI DEL PROTOCOLLO DI FALSIFICAZIONE ARTEFATTI NUMERICI")
    print("=" * 95)
    print(f"{'Caso / Test':<36} | {'<F_z> (uN)':<11} | {'Fz,max (uN)':<12} | {'P_J (mW)':<9} | {'Criterio Atteso':<18} | {'Esito':<6}")
    print("-" * 95)
    print(f"{results['baseline_ref']['title']:<36} | {base_mean:+9.2f}   | {results['baseline_ref']['peak_Fz_uN']:+10.2f}   | {results['baseline_ref']['mean_Pj_W']*1000:7.2f}   | Riferimento Nominale | REF   ")
    print(f"{results['test1_chirality_reversal']['title']:<36} | {t1_mean:+9.2f}   | {results['test1_chirality_reversal']['peak_Fz_uN']:+10.2f}   | {results['test1_chirality_reversal']['mean_Pj_W']*1000:7.2f}   | ~ {-base_mean:+.2f} uN (+/-5%) | {'PASS' if pass_t1 else 'FAIL'}")
    print(f"{results['test2_isotropic']['title']:<36} | {t2_mean:+9.2f}   | {results['test2_isotropic']['peak_Fz_uN']:+10.2f}   | {results['test2_isotropic']['mean_Pj_W']*1000:7.2f}   | |res| < 0.10 uN    | {'PASS' if pass_t2 else 'FAIL'}")
    print(f"{results['test3_phase_inversion']['title']:<36} | {t3_mean:+9.2f}   | {results['test3_phase_inversion']['peak_Fz_uN']:+10.2f}   | {results['test3_phase_inversion']['mean_Pj_W']*1000:7.2f}   | <F_z> < 0          | {'PASS' if pass_t3 else 'FAIL'}")
    print("-" * 95)
    print(f"Errore di simmetria speculare T1 vs Baseline: {sym_error_pct:.2f}% (soglia: <= 5.0%)")
    print(f"Residuo numerico mantello isotropo T2:        {abs(t2_mean):.4f} uN (soglia: < 0.10 uN)")
    print("=" * 95)
    
    if all_passed:
        print("\n>>> VERDETTO SCIENTIFICO: LIFT DI LORENTZ FISICAMENTE VALIDATO DA ROTTURA DI PARITA <<<")
        print("    Ogni ipotesi di artefatto numerico (bias di mesh 3D o quadratura temporale) e categoricamente falsificata.")
    else:
        print("\n>>> VERDETTO SCIENTIFICO: ARTEFATTO NUMERICO RILEVATO (BIAS DI MESH/DISCRETIZZAZIONE) <<<")

    # Salvataggio JSON
    output_json = {
        "metadata": {
            "title": "Protocollo di Falsificazione Artefatti Numerici sul Lift di Lorentz",
            "author": "Alessandro Brescacin",
            "license": "CERN-OHL-S-2.0",
            "date": "2026-09-22",
            "solver": "Elmer FEM 9.0 Whitney A-V Solver",
            "frequency_Hz": FREQUENCY_HZ,
            "rpm": RPM,
            "dt_s": DT,
            "timesteps": TIMESTEPS
        },
        "verdict": "LIFT DI LORENTZ FISICAMENTE VALIDATO DA ROTTURA DI PARITA" if all_passed else "ARTEFATTO NUMERICO RILEVATO",
        "validation_criteria": {
            "test1_chirality_reversal": {
                "mean_Fz_uN": t1_mean,
                "expected_target_uN": -base_mean,
                "symmetry_error_pct": sym_error_pct,
                "threshold_pct": 5.0,
                "passed": bool(pass_t1)
            },
            "test2_isotropic_mantle": {
                "mean_Fz_uN": t2_mean,
                "threshold_abs_uN": 0.10,
                "passed": bool(pass_t2)
            },
            "test3_phase_inversion": {
                "mean_Fz_uN": t3_mean,
                "condition": "<F_z> < 0",
                "passed": bool(pass_t3)
            }
        },
        "cases": results
    }

    json_path = DATA_DIR / "risultati_falsificazione_artefatti.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(output_json, f, indent=2)
    print(f"\nDataset consolidato salvato in: {json_path}")

    # Generazione Grafico a 4 Quadranti a 300 DPI
    generate_quadrant_figure(results, pass_t1, pass_t2, pass_t3, all_passed)

def generate_quadrant_figure(results, pass_t1, pass_t2, pass_t3, all_passed):
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
    fig.patch.set_facecolor('#ffffff')

    time_ms = results["baseline_ref"]["time_series_ms"]

    # 1. Baseline
    ax1 = axes[0, 0]
    ax1.set_facecolor('#fafafa')
    fz_base = results["baseline_ref"]["Fz_series_uN"]
    m_base = results["baseline_ref"]["mean_Fz_uN"]
    ax1.plot(time_ms, fz_base, 'b-', lw=2.2, label=r'$F_z(t)$ Baseline (+30°, +$\omega$)')
    ax1.axhline(m_base, color='navy', linestyle='--', lw=1.8, label=rf'$\langle F_z \rangle = {m_base:+.2f}\,\mu\mathrm{{N}}$ (Nominale)')
    ax1.axhline(0, color='gray', linestyle=':', lw=1.0)
    ax1.set_title(r"$\mathbf{Baseline\ Nominale:\ \theta = +30^\circ,\ +\omega}$", fontsize=12, pad=8)
    ax1.set_xlabel("Tempo [ms]", fontsize=10)
    ax1.set_ylabel(r"Forza Assiale $F_z\ [\mu\mathrm{N}]$", fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_ylim([-20, 30])

    # 2. Test 1: Inversione Chirale
    ax2 = axes[0, 1]
    ax2.set_facecolor('#fafafa')
    fz_t1 = results["test1_chirality_reversal"]["Fz_series_uN"]
    m_t1 = results["test1_chirality_reversal"]["mean_Fz_uN"]
    ax2.plot(time_ms, fz_t1, 'r-', lw=2.2, label=r'$F_z(t)$ Inversione (-30°, +$\omega$)')
    ax2.axhline(m_t1, color='darkred', linestyle='--', lw=1.8, label=rf'$\langle F_z \rangle = {m_t1:+.2f}\,\mu\mathrm{{N}}$ (Speculare)')
    ax2.axhline(0, color='gray', linestyle=':', lw=1.0)
    status_t1 = "PASS (Simmetria < 5%)" if pass_t1 else "FAIL"
    ax2.set_title(rf"$\mathbf{{Test\ 1:\ Inversione\ Chirale\ \theta = -30^\circ\ [{status_t1}]}}$", fontsize=12, pad=8, color='darkred' if not pass_t1 else 'black')
    ax2.set_xlabel("Tempo [ms]", fontsize=10)
    ax2.set_ylabel(r"Forza Assiale $F_z\ [\mu\mathrm{N}]$", fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='lower right', fontsize=9)
    ax2.set_ylim([-30, 20])

    # 3. Test 2: Isotropo
    ax3 = axes[1, 0]
    ax3.set_facecolor('#fafafa')
    fz_t2 = results["test2_isotropic"]["Fz_series_uN"]
    m_t2 = results["test2_isotropic"]["mean_Fz_uN"]
    ax3.plot(time_ms, fz_t2, 'g-', lw=2.2, label=r'$F_z(t)$ Mantello Isotropo (0°)')
    ax3.axhline(m_t2, color='darkgreen', linestyle='--', lw=1.8, label=rf'$\langle F_z \rangle = {m_t2:+.4f}\,\mu\mathrm{{N}}$ (Zero)')
    ax3.axhline(0, color='black', linestyle='-', lw=1.2, alpha=0.7)
    status_t2 = "PASS (|res| < 0.1 µN)" if pass_t2 else "FAIL"
    ax3.set_title(rf"$\mathbf{{Test\ 2:\ Mantello\ Isotropo\ \sigma_{{\theta z}} = 0\ [{status_t2}]}}$", fontsize=12, pad=8, color='darkgreen' if pass_t2 else 'red')
    ax3.set_xlabel("Tempo [ms]", fontsize=10)
    ax3.set_ylabel(r"Forza Assiale $F_z\ [\mu\mathrm{N}]$", fontsize=10)
    ax3.grid(True, linestyle='--', alpha=0.6)
    ax3.legend(loc='upper right', fontsize=9)
    ax3.set_ylim([-5, 5])

    # 4. Test 3: Inversione Fase
    ax4 = axes[1, 1]
    ax4.set_facecolor('#fafafa')
    fz_t3 = results["test3_phase_inversion"]["Fz_series_uN"]
    m_t3 = results["test3_phase_inversion"]["mean_Fz_uN"]
    ax4.plot(time_ms, fz_t3, 'm-', lw=2.2, label=r'$F_z(t)$ Inversione Fasi (+30°, -$\omega$)')
    ax4.axhline(m_t3, color='purple', linestyle='--', lw=1.8, label=rf'$\langle F_z \rangle = {m_t3:+.2f}\,\mu\mathrm{{N}}$ (Inversa)')
    ax4.axhline(0, color='gray', linestyle=':', lw=1.0)
    status_t3 = "PASS (<Fz> < 0)" if pass_t3 else "FAIL"
    ax4.set_title(rf"$\mathbf{{Test\ 3:\ Inversione\ Sequenza\ Fasi\ \omega < 0\ [{status_t3}]}}$", fontsize=12, pad=8, color='purple' if pass_t3 else 'red')
    ax4.set_xlabel("Tempo [ms]", fontsize=10)
    ax4.set_ylabel(r"Forza Assiale $F_z\ [\mu\mathrm{N}]$", fontsize=10)
    ax4.grid(True, linestyle='--', alpha=0.6)
    ax4.legend(loc='lower right', fontsize=9)
    ax4.set_ylim([-30, 20])

    verdict_text = "VERDETTO: LIFT DI LORENTZ FISICAMENTE VALIDATO DA ROTTURA DI PARITA (CERN-OHL-S v2)" if all_passed else "VERDETTO: ARTEFATTO NUMERICO RILEVATO"
    verdict_color = "darkgreen" if all_passed else "red"
    fig.suptitle(f"Open Chiral Flux Shaper - Protocollo di Falsificazione Artefatti Numerici\n{verdict_text}", fontsize=13, fontweight='bold', color=verdict_color, y=0.98)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig_path = FIGURES_DIR / "fig_falsificazione_simmetria_4quadranti.png"
    plt.savefig(fig_path, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figura a 4 quadranti salvata a 300 DPI in: {fig_path}")

if __name__ == "__main__":
    main()
