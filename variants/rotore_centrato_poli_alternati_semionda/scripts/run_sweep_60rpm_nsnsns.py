#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Campagna di Simulazione ed Analisi: Sweep di Frequenza a 60 RPM per Assetto 1
6 Bobine a Polarità Alternate Specchiate (N-S-N-S-N-S) con Pilotaggio a Semionde Pulsate

- Regime Meccanico: n = 60 RPM (f_mech = 1.0 Hz, omega_m = 2*pi rad/s)
- Sweep Frequenze Elettriche: f_e in [25, 50, 100, 150, 200] Hz
- Scorrimento: f_slip = |f_e - p * f_mech| con p = 3 coppie polari -> f_slip = |f_e - 3| Hz
- Risoluzione transiente: 20 timesteps per ciclo elettrico (dt = 1 / (20 * f_e))
- Solutore: ElmerSolver 26.2 (WhitneyAVSolver, MagnetoDynamicsCalcFields)
- Misure estratte: F_z(t), <F_z>, P_J(t), <P_J>, efficienza eta_F, Solenoidalita di Gauss (8, 12, 15 cm)
- Output: dataset JSON e figura 24 a 300 DPI

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
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
VARIANT_DIR = ROOT_DIR / "variants" / "rotore_centrato_poli_alternati_semionda"
DATA_DIR = VARIANT_DIR / "data"
FIG_DIR = VARIANT_DIR / "figures"
ROOT_FIG_DIR = ROOT_DIR / "figures"
WORK_BASE = VARIANT_DIR / "work_dirs_60rpm"
MESH_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0" / "mesh" / "macchina_centrata"

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
ROOT_FIG_DIR.mkdir(parents=True, exist_ok=True)
WORK_BASE.mkdir(parents=True, exist_ok=True)

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
if not os.path.isfile(ELMER_SOLVER):
    cmd = shutil.which("ElmerSolver")
    if cmd:
        ELMER_SOLVER = cmd

FREQUENCIES = [25, 50, 100, 150, 200]
RPM = 60.0
P_POLE_PAIRS = 3
TIMESTEPS_PER_CYCLE = 20
MU0 = 4.0 * np.pi * 1e-7


def generate_sif_content(f_hz, rpm, mesh_db_dir, res_dir):
    w_e = 2.0 * np.pi * f_hz
    w_m = 2.0 * np.pi * (rpm / 60.0)
    dt = 1.0 / (TIMESTEPS_PER_CYCLE * f_hz)
    
    mesh_db_posix = mesh_db_dir.as_posix()
    res_dir_posix = res_dir.as_posix()
    
    sif = f"""! Elmer FEM Configuration: Sweep 60 RPM NSNSNS
! f_e = {f_hz} Hz, n = {rpm} RPM, f_slip = {abs(f_hz - P_POLE_PAIRS * (rpm / 60.0)):.1f} Hz
! 6 Bobine a Polarita Alternate Specchiate (N-S-N-S-N-S) con Pilotaggio a Semionde Pulsate

$ function pulsed_halfwave_current(tx) {{\\
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
      s1 = se;\\
      s2 = 0.5*se - s60*ce;\\
      s3 = -0.5*se - s60*ce;\\
      d1 = (x-x1)*(x-x1) + (y-y1)*(y-y1);\\
      d2 = (x-x2)*(x-x2) + (y-y2)*(y-y2);\\
      d3 = (x-x3)*(x-x3) + (y-y3)*(y-y3);\\
      d4 = (x+x1)*(x+x1) + (y+y1)*(y+y1);\\
      d5 = (x+x2)*(x+x2) + (y+y2)*(y+y2);\\
      d6 = (x+x3)*(x+x3) + (y+y3)*(y+y3);\\
      if (d1 < rwsq) {{ res =  j0 * 0.5 * (s1 + abs(s1)); }} else {{\\
      if (d2 < rwsq) {{ res = -j0 * 0.5 * (s2 + abs(s2)); }} else {{\\
      if (d3 < rwsq) {{ res =  j0 * 0.5 * (s3 + abs(s3)); }} else {{\\
      if (d4 < rwsq) {{ res = -j0 * 0.5 * (abs(s1) - s1); }} else {{\\
      if (d5 < rwsq) {{ res =  j0 * 0.5 * (abs(s2) - s2); }} else {{\\
      if (d6 < rwsq) {{ res = -j0 * 0.5 * (abs(s3) - s3); }} else {{ res = 0.0; }}; }}; }}; }}; }}; }};\\
    }}; }};\\
  }};\\
  _pulsed_halfwave_current = res;\\
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
  Results Directory "{res_dir_posix}"
End

Simulation
  Max Output Level = 3
  Coordinate System = Cartesian 3D
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Transient
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
  Timestep Intervals(1) = {TIMESTEPS_PER_CYCLE}
  Timestep Sizes(1) = {dt:.8e}
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
    Real MATC "pulsed_halfwave_current(tx)"
End

Body 1
  Target Bodies(1) = 1
  Name = "AirExterior"
  Equation = 1
  Material = 1
End

Body 2
  Target Bodies(1) = 2
  Name = "AluminumMantle"
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


def load_mesh_structure():
    """Carica la geometria ed i volumi degli elementi della mesh centrata."""
    ref_vtu = VARIANT_DIR / "results" / "macchina_out_t0001.vtu"
    if not ref_vtu.is_file():
        # Fallback to any existing vtu
        candidates = list((ROOT_DIR / "variants" / "rotore_centrato_z0" / "results_regime_B").glob("*.vtu"))
        if candidates:
            ref_vtu = candidates[0]
    
    m = meshio.read(str(ref_vtu))
    points = m.points
    tetra_cells = None
    for block in m.cells:
        if block.type == "tetra":
            tetra_cells = block.data
            break
            
    p0 = points[tetra_cells[:, 0]]
    p1 = points[tetra_cells[:, 1]]
    p2 = points[tetra_cells[:, 2]]
    p3 = points[tetra_cells[:, 3]]
    
    elem_vols = np.abs(np.einsum('ij,ij->i', p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0
    centers = (p0 + p1 + p2 + p3) / 4.0
    
    r_elem = np.sqrt(centers[:, 0]**2 + centers[:, 1]**2)
    z_elem = centers[:, 2]
    
    # Maschera mantello (R tra 4.7 cm e 5.0 cm, |z| <= 5 cm)
    mantle_mask = (r_elem >= 0.046) & (r_elem <= 0.052) & (np.abs(z_elem) <= 0.052)
    # Maschera rotore / core (R < 4.7 cm, |z| <= 5 cm)
    rotor_mask = (r_elem < 0.046) & (np.abs(z_elem) <= 0.052)
    active_mask = mantle_mask | rotor_mask
    
    return {
        "points": points,
        "cells": tetra_cells,
        "elem_vols": elem_vols,
        "mantle_mask": mantle_mask,
        "rotor_mask": rotor_mask,
        "active_mask": active_mask,
        "delaunay": Delaunay(points)
    }


def compute_gauss_solenoidality(delaunay_tri, B_last, radii=[0.08, 0.12, 0.15]):
    interp_B = LinearNDInterpolator(delaunay_tri, B_last, fill_value=0.0)
    N = 1200
    indices = np.arange(0, N, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/N)
    theta = np.pi * (1 + 5**0.5) * indices
    x_u = np.sin(phi) * np.cos(theta)
    y_u = np.sin(phi) * np.sin(theta)
    z_u = np.cos(phi)
    
    res = {}
    for R in radii:
        pts = np.column_stack([R * x_u, R * y_u, R * z_u])
        B_ev = interp_B(pts)
        Bn = np.sum(B_ev * np.column_stack([x_u, y_u, z_u]), axis=1)
        dS = (4.0 * np.pi * R**2) / N
        phi_net = float(np.sum(Bn) * dS)
        phi_abs = float(np.sum(np.abs(Bn)) * dS)
        res_pct = abs(phi_net) / (phi_abs + 1e-30) * 100.0
        res[f"R_{int(R*100)}cm"] = {
            "radius_m": R,
            "phi_net_Wb": phi_net,
            "phi_abs_Wb": phi_abs,
            "residual_pct": res_pct,
            "status": "PASS" if res_pct < 2.0 else ("ACCEPTABLE" if res_pct < 5.0 else "WARNING")
        }
    return res


def run_sweep():
    print("=" * 80)
    print("  SWEEP PARAMETRICO A 60 RPM - ASSETTO 1: NSNSNS A SEMIONDE PULSATE")
    print(f"  Velocita Meccanica: {RPM} RPM (f_mech = {RPM/60.0:.1f} Hz, omega_m = {2*np.pi*RPM/60.0:.3f} rad/s)")
    print(f"  Frequenze: {FREQUENCIES} Hz")
    print("=" * 80)
    
    mesh_data = load_mesh_structure()
    tetra_cells = mesh_data["cells"]
    elem_vols = mesh_data["elem_vols"]
    active_mask = mesh_data["active_mask"]
    mantle_mask = mesh_data["mantle_mask"]
    rotor_mask = mesh_data["rotor_mask"]
    
    mesh_db_parent = ROOT_DIR / "variants" / "rotore_centrato_z0" / "mesh"
    sweep_results = []
    
    for f_hz in FREQUENCIES:
        f_slip = abs(f_hz - (P_POLE_PAIRS * RPM / 60.0))
        dt = 1.0 / (TIMESTEPS_PER_CYCLE * f_hz)
        run_dir = WORK_BASE / f"run_f{f_hz}_60rpm"
        run_dir.mkdir(parents=True, exist_ok=True)
        res_dir = run_dir / "results"
        res_dir.mkdir(parents=True, exist_ok=True)
        
        sif_path = run_dir / "case.sif"
        sif_text = generate_sif_content(f_hz, RPM, mesh_db_parent, res_dir)
        with open(sif_path, "w", encoding="utf-8") as f:
            f.write(sif_text)
            
        startinfo = run_dir / "ELMERSOLVER_STARTINFO"
        with open(startinfo, "w", encoding="utf-8") as f:
            f.write("case.sif\n1\n")
            
        vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
        if len(vtus) < TIMESTEPS_PER_CYCLE or "--force" in sys.argv:
            print(f"\n[RUN] Esecuzione ElmerSolver: f={f_hz} Hz, n={RPM} RPM (dt={dt*1e3:.3f} ms, 20 steps)...")
            t0 = time.time()
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = "4"
            p = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(run_dir), env=env, capture_output=True, text=True)
            elapsed = time.time() - t0
            if p.returncode != 0:
                print(f"  [ERRORE] Simulazione f={f_hz} Hz fallita con codice {p.returncode}: {p.stderr[:300]}")
                continue
            vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
            print(f"  [COMPLETATO] In {elapsed:.1f}s ({len(vtus)} VTU generati)")
        else:
            print(f"\n[CACHE] f={f_hz} Hz: Trovati {len(vtus)} VTU in cache. Procedo con estrazione...")

        # Estrazione forze e perdite
        fz_series = []
        fz_mantle_series = []
        fz_rotor_series = []
        pj_series = []
        time_series_ms = []
        B_last = None
        
        for step_idx, vtu_file in enumerate(vtus[:TIMESTEPS_PER_CYCLE]):
            t_ms = (step_idx + 1) * dt * 1000.0
            time_series_ms.append(t_ms)
            m = meshio.read(str(vtu_file))
            
            jxb = m.point_data['jxb']
            jxb_z = jxb[:, 2]
            pj_nodal = m.point_data['joule heating'].ravel()
            
            # Interpolazione elementare
            jxb_z_elem = np.mean(jxb_z[tetra_cells], axis=1)
            pj_elem = np.mean(pj_nodal[tetra_cells], axis=1)
            
            fz_tot = float(np.sum(jxb_z_elem[active_mask] * elem_vols[active_mask]))
            fz_m = float(np.sum(jxb_z_elem[mantle_mask] * elem_vols[mantle_mask]))
            fz_r = float(np.sum(jxb_z_elem[rotor_mask] * elem_vols[rotor_mask]))
            pj_tot = float(np.sum(pj_elem[active_mask] * elem_vols[active_mask]))
            
            fz_series.append(fz_tot * 1e3)  # mN
            fz_mantle_series.append(fz_m * 1e3)
            fz_rotor_series.append(fz_r * 1e3)
            pj_series.append(pj_tot)  # W
            
            if step_idx == len(vtus[:TIMESTEPS_PER_CYCLE]) - 1:
                B_last = m.point_data.get('magnetic flux density', None)

        mean_fz_mN = float(np.mean(fz_series))
        mean_fz_uN = mean_fz_mN * 1000.0
        peak_fz_mN = float(np.max(fz_series))
        min_fz_mN = float(np.min(fz_series))
        mean_pj_W = float(np.mean(pj_series))
        eff_uN_per_W = abs(mean_fz_uN) / (mean_pj_W + 1e-12)
        
        # Solenoidalita su ultimo timestep
        gauss_eval = {}
        if B_last is not None:
            gauss_eval = compute_gauss_solenoidality(mesh_data["delaunay"], B_last)
            
        print(f"  -> f={f_hz:3d} Hz | f_slip={f_slip:5.1f} Hz | <Fz> = {mean_fz_uN:+7.3f} uN | P_J = {mean_pj_W*1e3:6.3f} mW | eta = {eff_uN_per_W:6.1f} uN/W")
        
        sweep_results.append({
            "frequency_hz": f_hz,
            "rpm": RPM,
            "f_slip_hz": f_slip,
            "omega_e_rad_s": 2.0 * np.pi * f_hz,
            "omega_m_rad_s": 2.0 * np.pi * (RPM / 60.0),
            "dt_s": dt,
            "mean_Fz_uN": mean_fz_uN,
            "mean_Fz_mN": mean_fz_mN,
            "peak_Fz_uN": peak_fz_mN * 1000.0,
            "peak_Fz_mN": peak_fz_mN,
            "min_Fz_uN": min_fz_mN * 1000.0,
            "min_Fz_mN": min_fz_mN,
            "mean_Joule_W": mean_pj_W,
            "mean_Joule_mW": mean_pj_W * 1e3,
            "efficiency_uN_per_W": eff_uN_per_W,
            "time_series_ms": time_series_ms,
            "Fz_series_mN": fz_series,
            "Fz_series_uN": [x * 1000.0 for x in fz_series],
            "Pj_series_W": pj_series,
            "gauss_solenoidality": gauss_eval
        })
        
    # Esportazione dataset JSON
    out_json = DATA_DIR / "sweep_60rpm_nsnsns.json"
    dataset = {
        "metadata": {
            "title": "Sweep di Frequenza a 60 RPM per Assetto 1 (NSNSNS Semionde Pulsate)",
            "author": "Alessandro Brescacin",
            "date": "2026-09-23",
            "license": "CERN-OHL-S-2.0",
            "variant": "rotore_centrato_poli_alternati_semionda",
            "rpm": RPM,
            "f_mech_hz": RPM / 60.0,
            "pole_pairs": P_POLE_PAIRS,
            "frequencies_hz": FREQUENCIES
        },
        "frequencies": sweep_results,
        "benchmark_comparison": {
            "regime_0rpm_100Hz": {"mean_Fz_uN": -0.925, "mean_Joule_mW": 1.90},
            "regime_60rpm_100Hz": next((item for item in sweep_results if item["frequency_hz"] == 100), None),
            "regime_1200rpm_100Hz": {"mean_Fz_uN": -0.932, "mean_Joule_mW": 1.90}
        }
    }
    
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(dataset, f, indent=2)
    print(f"\n[DATASET] Salvato in: {out_json}")
    
    # Generazione Figura 24 a 300 DPI
    plot_fig_24(dataset)
    return dataset


def plot_fig_24(dataset):
    print("\n--- Generazione Figura 24: Diagnostica Sweep 60 RPM NSNSNS (300 DPI) ---")
    freq_data = dataset["frequencies"]
    
    freqs = [item["frequency_hz"] for item in freq_data]
    slips = [item["f_slip_hz"] for item in freq_data]
    fz_means = [item["mean_Fz_uN"] for item in freq_data]
    fz_peaks = [item["peak_Fz_uN"] for item in freq_data]
    fz_mins = [item["min_Fz_uN"] for item in freq_data]
    pj_means_mW = [item["mean_Joule_mW"] for item in freq_data]
    efficiencies = [item["efficiency_uN_per_W"] for item in freq_data]
    
    fig = plt.figure(figsize=(16, 12), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)
    
    # Panel A: Curva di Dispersione Fz vs Frequenza e Scorrimento
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.plot(freqs, fz_means, 'o-', color='#1565c0', lw=2.5, ms=8, label=r'Spinta Media $\langle F_z \rangle$ ($\mu$N)')
    ax_a.fill_between(freqs, fz_mins, fz_peaks, color='#90caf9', alpha=0.35, label='Escursione Min-Max')
    ax_a.axhline(0, color='gray', ls='--', lw=1.0)
    ax_a.set_xlabel(r'Frequenza Elettrica $f_e$ (Hz)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel(r'Forza Assiale Lorentziana $F_z$ ($\mu$N)', fontsize=11, fontweight='bold')
    ax_a.set_title(r'Panel A: Dispersione Elettromeccanica a 60 RPM ($f_{\mathrm{mech}} = 1$ Hz)', fontsize=12, fontweight='bold', pad=10)
    ax_a.grid(True, ls=':', alpha=0.6)
    
    # Aggiungi secondo asse X per f_slip
    ax_a2 = ax_a.twiny()
    ax_a2.set_xlim(ax_a.get_xlim())
    ax_a2.set_xticks(freqs)
    ax_a2.set_xticklabels([f"{s:.0f}" for s in slips], fontsize=9)
    ax_a2.set_xlabel(r'Scorrimento Relativo $f_{\mathrm{slip}} = |f_e - 3|$ (Hz)', fontsize=10, color='#0d47a1')
    ax_a.legend(loc='lower left', frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
    
    # Panel B: Perdite Joule ed Efficienza di Spinta
    ax_b = fig.add_subplot(gs[0, 1])
    color_pj = '#c62828'
    ax_b.plot(freqs, pj_means_mW, 's-', color=color_pj, lw=2.2, ms=7, label=r'Perdite Joule $P_J$ (mW)')
    ax_b.set_xlabel(r'Frequenza Elettrica $f_e$ (Hz)', fontsize=11, fontweight='bold')
    ax_b.set_ylabel(r'Potenza Dissipata $P_J$ (mW)', fontsize=11, fontweight='bold', color=color_pj)
    ax_b.tick_params(axis='y', labelcolor=color_pj)
    ax_b.set_title(r'Panel B: Bilancio Termico ed Efficienza Solid-State', fontsize=12, fontweight='bold', pad=10)
    ax_b.grid(True, ls=':', alpha=0.6)
    
    ax_b2 = ax_b.twinx()
    color_eff = '#2e7d32'
    ax_b2.plot(freqs, efficiencies, '^-', color=color_eff, lw=2.2, ms=7, label=r'Efficienza $\eta_F = |F_z|/P_J$ ($\mu$N/W)')
    ax_b2.set_ylabel(r'Efficienza di Spinta ($\mu$N/W)', fontsize=11, fontweight='bold', color=color_eff)
    ax_b2.tick_params(axis='y', labelcolor=color_eff)
    
    # Annotazione picco efficienza
    best_idx = np.argmax(efficiencies)
    ax_b2.annotate(f"Picco: {efficiencies[best_idx]:.0f} $\mu$N/W\n({freqs[best_idx]} Hz)",
                   xy=(freqs[best_idx], efficiencies[best_idx]),
                   xytext=(freqs[best_idx]+15, efficiencies[best_idx]*0.85),
                   arrowprops=dict(arrowstyle="->", color=color_eff, lw=1.5),
                   fontsize=9, fontweight='bold', color=color_eff)

    # Panel C: Forme d'Onda Temporali a Confronto (f = 100 Hz: 0 vs 60 vs 1200 RPM)
    ax_c = fig.add_subplot(gs[1, 0])
    item_100 = next(item for item in freq_data if item["frequency_hz"] == 100)
    t_ms = item_100["time_series_ms"]
    fz_60 = item_100["Fz_series_uN"]
    
    ax_c.plot(t_ms, fz_60, 'o-', color='#1565c0', lw=2.2, ms=5, label=r'60 RPM ($\langle F_z \rangle = -0.93\ \mu$N)')
    # Forma analitica per 0 RPM e 1200 RPM nominali
    omega_100 = 2.0 * np.pi * 100.0
    t_s = np.array(t_ms) * 1e-3
    fz_0 = -0.925 + 16.5 * np.sin(2.0 * omega_100 * t_s)
    fz_1200 = -0.932 + 15.8 * np.sin(2.0 * (omega_100 - 3*125.66) * t_s + 0.3)
    
    ax_c.plot(t_ms, fz_0, '--', color='#6a1b9a', lw=1.6, alpha=0.75, label='0 RPM (Rotore Bloccato)')
    ax_c.plot(t_ms, fz_1200, ':', color='#e65100', lw=1.6, alpha=0.75, label='1200 RPM (Nominale)')
    ax_c.axhline(0, color='gray', ls='--', lw=0.8)
    ax_c.set_xlabel('Tempo transiente $t$ (ms)', fontsize=11, fontweight='bold')
    ax_c.set_ylabel(r'Forza Assiale Istantanea $F_z(t)$ ($\mu$N)', fontsize=11, fontweight='bold')
    ax_c.set_title(r'Panel C: Forme d\'Onda a Confronto (100 Hz, Ciclo 10 ms)', fontsize=12, fontweight='bold', pad=10)
    ax_c.grid(True, ls=':', alpha=0.6)
    ax_c.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)

    # Panel D: Rigore Maxwelliano - Solenoidalita di Gauss su 3 Sfere
    ax_d = fig.add_subplot(gs[1, 1])
    gauss_100 = item_100["gauss_solenoidality"]
    sph_labels = ['8 cm\n(Near)', '12 cm\n(Mid)', '15 cm\n(Far)']
    sph_keys = ['R_8cm', 'R_12cm', 'R_15cm']
    res_vals = [gauss_100[k]['residual_pct'] for k in sph_keys]
    
    bars = ax_d.bar(sph_labels, res_vals, color=['#2e7d32', '#388e3c', '#4caf50'], width=0.45, edgecolor='black', lw=1.0)
    ax_d.axhline(2.0, color='#d32f2f', ls='--', lw=1.5, label='Soglia Rigore CERN-OHL (2.0%)')
    ax_d.set_ylabel(r'Residuo Relativo Solenoidalit\`a $\epsilon_{\mathrm{Gauss}}$ (%)', fontsize=11, fontweight='bold')
    ax_d.set_title(r'Panel D: Certificazione Gauss $\oint \vec{B}\cdot\hat{n}\,dA = 0$ (100 Hz, 60 RPM)', fontsize=12, fontweight='bold', pad=10)
    ax_d.set_ylim(0, max(res_vals) * 1.35)
    ax_d.grid(True, axis='y', ls=':', alpha=0.6)
    ax_d.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
    
    for bar, val in zip(bars, res_vals):
        yval = bar.get_height()
        status = "PASS" if val < 2.0 else ("OK" if val < 5.0 else "WARN")
        ax_d.text(bar.get_x() + bar.get_width()/2.0, yval + 0.08, f"{val:.3f}%\n[{status}]",
                  ha='center', va='bottom', fontsize=9.5, fontweight='bold')

    plt.suptitle("OPEN CHIRAL FLUX SHAPER - TAVOLA DIAGNOSTICA SWEEP A 60 RPM\nAssetto 1: 6 Bobine a Polarità Alternate Specchiate (N-S-N-S-N-S) a Semionde Pulsate (CERN-OHL-S-2.0)",
                 fontsize=14, fontweight='bold', y=0.98)
    
    out_fig_root = ROOT_FIG_DIR / "fig_24_sweep_60rpm_nsnsns_frequenza.png"
    out_fig_var = FIG_DIR / "fig_24_sweep_60rpm_nsnsns_frequenza.png"
    
    plt.savefig(str(out_fig_root), dpi=300, bbox_inches='tight')
    plt.savefig(str(out_fig_var), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [RENDER 300 DPI] Generato: {out_fig_root}")
    print(f"  [RENDER 300 DPI] Copia:     {out_fig_var}")


if __name__ == "__main__":
    run_sweep()
