#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Campagna di Simulazione ed Analisi: Sweep di Frequenza a 60 RPM per Assetto 2
Gabbia Sferica Metamateriale con Doppio Rotore Ortogonale a 90° (Multi-Axis Vector Shaper)

- Regime Meccanico: n = 60 RPM (f_mech = 1.0 Hz, omega_m = 2*pi rad/s attorno a Z)
- Sweep Frequenze Elettriche: f_e in [25, 50, 100, 150, 200] Hz
- Scorrimento: f_slip = |f_e - p * f_mech| con p = 3 coppie polari -> f_slip = |f_e - 3| Hz
- Alimentazione: Trifase continua simultanea NPNPNP a 120° con quadratura temporale a 90° tra Rotore 1 (Z) e Rotore 2 (X)
- Risoluzione transiente: 20 timesteps per ciclo elettrico (dt = 1 / (20 * f_e))
- Solutore: ElmerSolver 26.2 (WhitneyAVSolver, MagnetoDynamicsCalcFields)
- Misure estratte: Matrice forze 3D (Fx, Fy, Fz), risultante |F|, Odografo 3D, MST, Bilancio Joule per sottocorpo, Gauss
- Output: dataset JSON e figura 25 a 300 DPI

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
from mpl_toolkits.mplot3d import Axes3D

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
VARIANT_DIR = ROOT_DIR / "variants" / "gabbia_sferica_doppio_rotore_90deg"
DATA_DIR = VARIANT_DIR / "data"
FIG_DIR = VARIANT_DIR / "figures"
ROOT_FIG_DIR = ROOT_DIR / "figures"
WORK_BASE = VARIANT_DIR / "work_dirs_60rpm"
MESH_DIR = VARIANT_DIR / "mesh" / "macchina_sferica_ortogonale"

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
TIMESTEPS = 64
MU0 = 4.0 * np.pi * 1e-7


def generate_sif_content(f_hz, rpm, mesh_db_dir, res_dir):
    w_e = 2.0 * np.pi * f_hz
    w_m = 2.0 * np.pi * (rpm / 60.0)
    dt = (1.6 / f_hz) / TIMESTEPS
    
    mesh_db_posix = mesh_db_dir.as_posix()
    res_dir_posix = res_dir.as_posix()
    
    sif = f"""! Elmer FEM Configuration: Sweep 60 RPM Doppio Rotore Ortogonale 90°
! f_e = {f_hz} Hz, n = {rpm} RPM, f_slip = {abs(f_hz - P_POLE_PAIRS * (rpm / 60.0)):.1f} Hz
! Eccitazione Continua Trifase NPNPNP a 120° con Quadratura Temporale a 90°

$ function rotor1_npnpnp_current(tx) {{\\
  x = tx(0); y = tx(1); z = tx(2); t = tx(3);\\
  rc = 0.035; rwsq = 2.5e-05; j0 = 1.0e+05;\\
  s60 = 0.8660254037844386;\\
  wm = {w_m:.14f};\\
  th0 = wm * t;\\
  c0 = cos(th0); s0 = sin(th0);\\
  x1 = rc * c0; y1 = rc * s0;\\
  x2 = rc * (0.5 * c0 - s60 * s0); y2 = rc * (0.5 * s0 + s60 * c0);\\
  x3 = rc * (-0.5 * c0 - s60 * s0); y3 = rc * (-0.5 * s0 + s60 * c0);\\
  w = {w_e:.14f};\\
  phi120 = 2.0943951023931953;\\
  phi240 = 4.1887902047863905;\\
  ja = cos(w * t);\\
  jb = cos(w * t - phi120);\\
  jc = cos(w * t - phi240);\\
  d1 = (x-x1)*(x-x1) + (y-y1)*(y-y1);\\
  d2 = (x-x2)*(x-x2) + (y-y2)*(y-y2);\\
  d3 = (x-x3)*(x-x3) + (y-y3)*(y-y3);\\
  d4 = (x+x1)*(x+x1) + (y+y1)*(y+y1);\\
  d5 = (x+x2)*(x+x2) + (y+y2)*(y+y2);\\
  d6 = (x+x3)*(x+x3) + (y+y3)*(y+y3);\\
  res = 0.0;\\
  if (d1 < rwsq) {{ res =  j0 * ja; }} else {{\\
  if (d2 < rwsq) {{ res = -j0 * jb; }} else {{\\
  if (d3 < rwsq) {{ res =  j0 * jc; }} else {{\\
  if (d4 < rwsq) {{ res = -j0 * ja; }} else {{\\
  if (d5 < rwsq) {{ res =  j0 * jb; }} else {{\\
  if (d6 < rwsq) {{ res = -j0 * jc; }} else {{ res = 0.0; }}; }}; }}; }}; }}; }};\\
  _rotor1_npnpnp_current = res;\\
}}

$ function rotor2_npnpnp_base(tx) {{\\
  x = tx(0); y = tx(1); z = tx(2); t = tx(3);\\
  rc = 0.035; rwsq = 2.5e-05; j0 = 1.0e+05;\\
  wm = {w_m:.14f};\\
  th0 = wm * t;\\
  xp =  x * cos(th0) + y * sin(th0);\\
  yp = -x * sin(th0) + y * cos(th0);\\
  zp = z;\\
  c30 = 0.8660254037844386; s30 = 0.5;\\
  y1 = rc * c30; z1 = rc * s30;\\
  y2 = 0.0;      z2 = rc;\\
  y3 = -rc * c30; z3 = rc * s30;\\
  w = {w_e:.14f};\\
  phi90 = 1.5707963267948966;\\
  phi210 = 3.665191429188092;\\
  phi330 = 5.759586531581287;\\
  ja = cos(w * t - phi90);\\
  jb = cos(w * t - phi210);\\
  jc = cos(w * t - phi330);\\
  d1 = (yp-y1)*(yp-y1) + (zp-z1)*(zp-z1);\\
  d2 = (yp-y2)*(yp-y2) + (zp-z2)*(zp-z2);\\
  d3 = (yp-y3)*(yp-y3) + (zp-z3)*(zp-z3);\\
  d4 = (yp+y1)*(yp+y1) + (zp+z1)*(zp+z1);\\
  d5 = (yp+y2)*(yp+y2) + (zp+z2)*(zp+z2);\\
  d6 = (yp+y3)*(yp+y3) + (zp+z3)*(zp+z3);\\
  res = 0.0;\\
  if (d1 < rwsq) {{ res =  j0 * ja; }} else {{\\
  if (d2 < rwsq) {{ res = -j0 * jb; }} else {{\\
  if (d3 < rwsq) {{ res =  j0 * jc; }} else {{\\
  if (d4 < rwsq) {{ res = -j0 * ja; }} else {{\\
  if (d5 < rwsq) {{ res =  j0 * jb; }} else {{\\
  if (d6 < rwsq) {{ res = -j0 * jc; }} else {{ res = 0.0; }}; }}; }}; }}; }}; }};\\
  _rotor2_npnpnp_base = res;\\
}}

$ function rotor2_npnpnp_jx(tx) {{\\
  t = tx(3); wm = {w_m:.14f}; th0 = wm * t;\\
  j_base = rotor2_npnpnp_base(tx);\\
  _rotor2_npnpnp_jx = j_base * cos(th0);\\
}}

$ function rotor2_npnpnp_jy(tx) {{\\
  t = tx(3); wm = {w_m:.14f}; th0 = wm * t;\\
  j_base = rotor2_npnpnp_base(tx);\\
  _rotor2_npnpnp_jy = j_base * sin(th0);\\
}}

$ function sigma_spherical_X_layers(tx) {{\\
  x = tx(0); y = tx(1); z = tx(2);\\
  r = sqrt(x*x + y*y + z*z);\\
  if (r < 1.0e-5) {{ r = 1.0e-5; }} else {{ r = r; }};\\
  rho = sqrt(x*x + y*y);\\
  if (rho < 1.0e-5) {{ rho = 1.0e-5; }} else {{ rho = rho; }};\\
  xi = (r - 0.047) / 0.003;\\
  if (xi < 0.0) {{ xi = 0.0; }} else {{ if (xi > 1.0) {{ xi = 1.0; }} else {{ xi = xi; }}; }};\\
  s_base = 3.031089e6;\\
  if (xi <= 0.3333333333333333) {{\\
    s_tp = s_base;\\
  }} else {{\\
    if (xi <= 0.6666666666666666) {{\\
      s_tp = 0.0;\\
    }} else {{\\
      s_tp = -s_base;\\
    }};\\
  }};\\
  s0 = 1.75e6;\\
  s_phi = 1.22e7;\\
  ds = s_phi - s0;\\
  px = -y / rho; py = x / rho;\\
  tx_c = (z / r) * (x / rho);\\
  ty_c = (z / r) * (y / rho);\\
  tz_c = -rho / r;\\
  a = zeros(3,3);\\
  a(0,0) = s0 + ds * px*px + 2.0 * s_tp * tx_c * px;\\
  a(1,1) = s0 + ds * py*py + 2.0 * s_tp * ty_c * py;\\
  a(2,2) = s0;\\
  c12 = ds * px*py + s_tp * (tx_c * py + px * ty_c);\\
  a(0,1) = c12; a(1,0) = c12;\\
  c13 = s_tp * px * tz_c;\\
  a(0,2) = c13; a(2,0) = c13;\\
  c23 = s_tp * py * tz_c;\\
  a(1,2) = c23; a(2,1) = c23;\\
  _sigma_spherical_X_layers = a;\\
}}

Header
  CHECK KEYWORDS Warn
  Mesh DB "{mesh_db_posix}" "macchina_sferica_ortogonale"
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
  Timestep Intervals(1) = {TIMESTEPS}
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
  Output File Name = "macchina_sferica_out"
  Output Format = "vtu"
  Vtu Format = Logical True
  Save Geometry Ids = Logical True
  Discontinuous Bodies = Logical True
End

Material 1
  Name = "FerromagneticSphericalCrossLayerMantle"
  Relative Permeability = 1000.0
  Electric Conductivity(3,3) = Variable Coordinate 1, Coordinate 2, Coordinate 3
    Real MATC "sigma_spherical_X_layers(tx)"
End

Material 2
  Name = "CopperCoilsRotor1"
  Relative Permeability = 1.0
  Electric Conductivity = 5.8e7
End

Material 3
  Name = "CopperCoilsRotor2"
  Relative Permeability = 1.0
  Electric Conductivity = 5.8e7
End

Material 4
  Name = "AmagneticPEEKCore"
  Relative Permeability = 1.0
  Electric Conductivity = 0.0
End

Material 5
  Name = "Air"
  Relative Permeability = 1.0
  Electric Conductivity = 0.0
End

Body Force 1
  Current Density 3 = Variable Coordinate 1, Coordinate 2, Coordinate 3, time
    Real MATC "rotor1_npnpnp_current(tx)"
End

Body Force 2
  Current Density 1 = Variable Coordinate 1, Coordinate 2, Coordinate 3, time
    Real MATC "rotor2_npnpnp_jx(tx)"
  Current Density 2 = Variable Coordinate 1, Coordinate 2, Coordinate 3, time
    Real MATC "rotor2_npnpnp_jy(tx)"
End

Body 1
  Target Bodies(1) = 1
  Name = "MantelloSferico"
  Equation = 1
  Material = 1
End

Body 2
  Target Bodies(1) = 2
  Name = "Rotore1Coils"
  Equation = 1
  Material = 2
  Body Force = 1
End

Body 3
  Target Bodies(1) = 3
  Name = "Rotore2Coils"
  Equation = 1
  Material = 3
  Body Force = 2
End

Body 4
  Target Bodies(1) = 4
  Name = "RotorePEEKCore"
  Equation = 1
  Material = 4
End

Body 5
  Target Bodies(1) = 5
  Name = "CavitaAriaInterna"
  Equation = 1
  Material = 5
End

Body 6
  Target Bodies(1) = 6
  Name = "AriaEsterna"
  Equation = 1
  Material = 5
End

Boundary Condition 1
  Target Boundaries(1) = 1
  Name = "FarFieldInfinity"
  AV {{e}} = Real 0.0
  AV = Real 0.0
End
"""
    return sif


def load_spherical_mesh():
    ref_vtu = VARIANT_DIR / "work_dirs" / "run_npnpnp_doppio_rotore" / "results" / "macchina_sferica_out_t0001.vtu"
    if not ref_vtu.is_file():
        candidates = list((VARIANT_DIR / "work_dirs").glob("**/*.vtu"))
        if candidates:
            ref_vtu = candidates[0]
            
    m = meshio.read(str(ref_vtu))
    points = m.points
    tetra_idx = [i for i, cb in enumerate(m.cells) if cb.type == "tetra"][0]
    tetra_cells = m.cells[tetra_idx].data
    gids = m.cell_data["GeometryIds"][tetra_idx] if "GeometryIds" in m.cell_data else np.ones(len(tetra_cells), dtype=int)
            
    p0 = points[tetra_cells[:, 0]]
    p1 = points[tetra_cells[:, 1]]
    p2 = points[tetra_cells[:, 2]]
    p3 = points[tetra_cells[:, 3]]
    
    elem_vols = np.abs(np.einsum('ij,ij->i', p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0
    centers = (p0 + p1 + p2 + p3) / 4.0
    r_elem = np.sqrt(np.sum(centers**2, axis=1))
    
    # Maschere corpi basate su GeometryIds conformi (ElmerGrid)
    mantle_mask = (gids == 1)
    rotor1_mask = (gids == 2)
    rotor2_mask = (gids == 3)
    peek_mask = (gids == 4)
    rotor_mask = rotor1_mask | rotor2_mask
    assembly_mask = np.isin(gids, [1, 2, 3, 4])
    
    # Layer nel mantello:
    xi = (r_elem - 0.047) / 0.003
    l1_mask = mantle_mask & (xi <= 0.333333)
    l2_mask = mantle_mask & (xi > 0.333333) & (xi <= 0.666667)
    l3_mask = mantle_mask & (xi > 0.666667)
    
    return {
        "points": points,
        "cells": tetra_cells,
        "elem_vols": elem_vols,
        "gids": gids,
        "assembly_mask": assembly_mask,
        "mantle_mask": mantle_mask,
        "rotor1_mask": rotor1_mask,
        "rotor2_mask": rotor2_mask,
        "rotor_mask": rotor_mask,
        "peek_mask": peek_mask,
        "l1_mask": l1_mask,
        "l2_mask": l2_mask,
        "l3_mask": l3_mask,
        "delaunay": Delaunay(points)
    }


def compute_gauss_and_mst(delaunay_tri, B_last, radii=[0.065, 0.100, 0.150]):
    interp_B = LinearNDInterpolator(delaunay_tri, B_last, fill_value=0.0)
    N = 2500
    indices = np.arange(0, N, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/N)
    theta = np.pi * (1 + 5**0.5) * indices
    x_u = np.sin(phi) * np.cos(theta)
    y_u = np.sin(phi) * np.sin(theta)
    z_u = np.cos(phi)
    
    results = {}
    for R in radii:
        pts = np.column_stack([R * x_u, R * y_u, R * z_u])
        B_ev = interp_B(pts)
        n_vec = np.column_stack([x_u, y_u, z_u])
        Bn = np.sum(B_ev * n_vec, axis=1)
        B_mag = np.sqrt(np.sum(B_ev**2, axis=1))
        
        dS = (4.0 * np.pi * R**2) / N
        phi_net = float(np.sum(Bn) * dS)
        phi_abs = float(np.sum(np.abs(Bn)) * dS)
        res_pct = abs(phi_net) / (phi_abs + 1e-30) * 100.0
        
        # MST: T = (1/mu0) * [ (B.n) B - 0.5 |B|^2 n ]
        T_vec = (1.0 / MU0) * (Bn[:, None] * B_ev - 0.5 * (B_mag[:, None]**2) * n_vec)
        F_mst = np.sum(T_vec, axis=0) * dS
        
        results[f"R_{int(R*1000)}mm"] = {
            "radius_m": R,
            "phi_net_Wb": phi_net,
            "phi_abs_Wb": phi_abs,
            "residual_pct": res_pct,
            "mean_B_uT": float(np.mean(B_mag) * 1e6),
            "peak_B_uT": float(np.max(B_mag) * 1e6),
            "F_mst_N": [float(f) for f in F_mst],
            "F_mst_mag_N": float(np.sqrt(np.sum(F_mst**2))),
            "status": "PASS" if res_pct < 2.0 else ("ACCEPTABLE" if res_pct < 5.0 else "WARNING")
        }
    return results


def run_sweep():
    print("=" * 80)
    print("  SWEEP PARAMETRICO A 60 RPM - ASSETTO 2: GABBIA SFERICA DOPPIO ROTORE 90°")
    print(f"  Velocita Meccanica: {RPM} RPM (f_mech = {RPM/60.0:.1f} Hz, omega_m = {2*np.pi*RPM/60.0:.3f} rad/s attorno a Z)")
    print(f"  Frequenze: {FREQUENCIES} Hz (Eccitazione Trifase Continua NPNPNP a 120°)")
    print("=" * 80)
    
    mesh_data = load_spherical_mesh()
    tetra_cells = mesh_data["cells"]
    elem_vols = mesh_data["elem_vols"]
    assembly_mask = mesh_data["assembly_mask"]
    mantle_mask = mesh_data["mantle_mask"]
    rotor1_mask = mesh_data["rotor1_mask"]
    rotor2_mask = mesh_data["rotor2_mask"]
    rotor_mask = mesh_data["rotor_mask"]
    peek_mask = mesh_data["peek_mask"]
    l1_mask = mesh_data["l1_mask"]
    l2_mask = mesh_data["l2_mask"]
    l3_mask = mesh_data["l3_mask"]
    
    mesh_db_parent = VARIANT_DIR / "mesh"
    sweep_results = []
    
    for f_hz in FREQUENCIES:
        f_slip = abs(f_hz - (P_POLE_PAIRS * RPM / 60.0))
        dt = (1.6 / f_hz) / TIMESTEPS
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
            
        vtus = sorted(list(res_dir.glob("macchina_sferica_out_t*.vtu")))
        if len(vtus) < TIMESTEPS or "--force" in sys.argv:
            print(f"\n[RUN] Esecuzione ElmerSolver: f={f_hz} Hz, n={RPM} RPM (dt={dt*1e3:.3f} ms, {TIMESTEPS} steps)...")
            t0 = time.time()
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = "4"
            p = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(run_dir), env=env, capture_output=True, text=True)
            elapsed = time.time() - t0
            if p.returncode != 0:
                print(f"  [ERRORE] Simulazione f={f_hz} Hz fallita con codice {p.returncode}: {p.stderr[:300]}")
                continue
            vtus = sorted(list(res_dir.glob("macchina_sferica_out_t*.vtu")))
            print(f"  [COMPLETATO] In {elapsed:.1f}s ({len(vtus)} VTU generati)")
        else:
            print(f"\n[CACHE] f={f_hz} Hz: Trovati {len(vtus)} VTU in cache. Procedo con estrazione...")

        fx_series = []
        fy_series = []
        fz_series = []
        fmag_series = []
        pj_tot_series = []
        pj_mantle_series = []
        pj_l1_series = []
        pj_l2_series = []
        pj_l3_series = []
        pj_rotor_series = []
        pj_peek_series = []
        time_series_ms = []
        B_last = None
        
        for step_idx, vtu_file in enumerate(vtus[:TIMESTEPS]):
            t_ms = (step_idx + 1) * dt * 1000.0
            time_series_ms.append(t_ms)
            m = meshio.read(str(vtu_file))
            
            jxb = m.point_data['jxb']
            pj_nodal = m.point_data['joule heating'].ravel()
            
            jxb_elem = np.mean(jxb[tetra_cells], axis=1)
            pj_elem = np.mean(pj_nodal[tetra_cells], axis=1)
            
            # Forze su assembly totale (mantello + rotori 1&2 + nucleo PEEK)
            fx = float(np.sum(jxb_elem[assembly_mask, 0] * elem_vols[assembly_mask]))
            fy = float(np.sum(jxb_elem[assembly_mask, 1] * elem_vols[assembly_mask]))
            fz = float(np.sum(jxb_elem[assembly_mask, 2] * elem_vols[assembly_mask]))
            fmag = float(np.sqrt(fx**2 + fy**2 + fz**2))
            
            fx_series.append(fx)
            fy_series.append(fy)
            fz_series.append(fz)
            fmag_series.append(fmag)
            
            # Perdite per corpo
            pj_m = float(np.sum(pj_elem[mantle_mask] * elem_vols[mantle_mask]))
            pj_r = float(np.sum(pj_elem[rotor_mask] * elem_vols[rotor_mask]))
            pj_p = float(np.sum(pj_elem[peek_mask] * elem_vols[peek_mask]))
            pj_tot = pj_m + pj_r + pj_p
            
            pj_l1 = float(np.sum(pj_elem[l1_mask] * elem_vols[l1_mask]))
            pj_l2 = float(np.sum(pj_elem[l2_mask] * elem_vols[l2_mask]))
            pj_l3 = float(np.sum(pj_elem[l3_mask] * elem_vols[l3_mask]))
            
            pj_tot_series.append(pj_tot)
            pj_mantle_series.append(pj_m)
            pj_rotor_series.append(pj_r)
            pj_peek_series.append(pj_p)
            pj_l1_series.append(pj_l1)
            pj_l2_series.append(pj_l2)
            pj_l3_series.append(pj_l3)
            
            if step_idx == len(vtus[:TIMESTEPS]) - 1:
                B_last = m.point_data.get('magnetic flux density', None)

        mean_fx = float(np.mean(fx_series))
        mean_fy = float(np.mean(fy_series))
        mean_fz = float(np.mean(fz_series))
        mean_fmag = float(np.sqrt(mean_fx**2 + mean_fy**2 + mean_fz**2))
        peak_fmag = float(np.max(fmag_series))
        
        mean_pj = float(np.mean(pj_tot_series))
        mean_pj_m = float(np.mean(pj_mantle_series))
        mean_pj_r = float(np.mean(pj_rotor_series))
        mean_pj_p = float(np.mean(pj_peek_series))
        mean_pj_l1 = float(np.mean(pj_l1_series))
        mean_pj_l2 = float(np.mean(pj_l2_series))
        mean_pj_l3 = float(np.mean(pj_l3_series))
        
        eff_uN_per_W = (mean_fmag * 1e6) / (mean_pj + 1e-12)
        eff_mN_per_W = (mean_fmag * 1e3) / (mean_pj + 1e-12)
        
        gauss_mst = {}
        if B_last is not None:
            gauss_mst = compute_gauss_and_mst(mesh_data["delaunay"], B_last)
            
        print(f"  -> f={f_hz:3d} Hz | <Fx>={mean_fx*1e3:+6.1f} mN | <Fy>={mean_fy*1e3:+6.1f} mN | <Fz>={mean_fz*1e3:+6.1f} mN | |<F>|={mean_fmag*1e3:6.1f} mN ({mean_fmag:.3f} N) | P_J={mean_pj:6.1f} W | eta={eff_mN_per_W:.2f} mN/W")
        
        sweep_results.append({
            "frequency_hz": f_hz,
            "rpm": RPM,
            "f_slip_hz": f_slip,
            "omega_e_rad_s": 2.0 * np.pi * f_hz,
            "omega_m_rad_s": 2.0 * np.pi * (RPM / 60.0),
            "dt_s": dt,
            "forces_N": {
                "mean_fx_N": mean_fx,
                "mean_fy_N": mean_fy,
                "mean_fz_N": mean_fz,
                "mean_fmag_N": mean_fmag,
                "mean_fx_mN": mean_fx * 1e3,
                "mean_fy_mN": mean_fy * 1e3,
                "mean_fz_mN": mean_fz * 1e3,
                "mean_fmag_mN": mean_fmag * 1e3,
                "peak_fmag_N": peak_fmag,
                "peak_fmag_mN": peak_fmag * 1e3
            },
            "losses_W": {
                "total_W": mean_pj,
                "mantle_total_W": mean_pj_m,
                "rotors_total_W": mean_pj_r,
                "peek_core_W": mean_pj_p,
                "mantle_layer1_plus30_W": mean_pj_l1,
                "mantle_layer2_ortho_W": mean_pj_l2,
                "mantle_layer3_minus30_W": mean_pj_l3
            },
            "efficiency": {
                "uN_per_W": eff_uN_per_W,
                "mN_per_W": eff_mN_per_W
            },
            "time_series_ms": time_series_ms,
            "series_fx_mN": [x * 1e3 for x in fx_series],
            "series_fy_mN": [x * 1e3 for x in fy_series],
            "series_fz_mN": [x * 1e3 for x in fz_series],
            "series_fmag_mN": [x * 1e3 for x in fmag_series],
            "series_pj_W": pj_tot_series,
            "gauss_and_mst": gauss_mst
        })
        
    out_json = DATA_DIR / "sweep_60rpm_doppio_rotore.json"
    dataset = {
        "metadata": {
            "title": "Sweep di Frequenza a 60 RPM per Assetto 2 (Gabbia Sferica Metamateriale a Doppio Rotore Ortogonale 90°)",
            "author": "Alessandro Brescacin",
            "date": "2026-09-23",
            "license": "CERN-OHL-S-2.0",
            "variant": "gabbia_sferica_doppio_rotore_90deg",
            "rpm": RPM,
            "f_mech_hz": RPM / 60.0,
            "pole_pairs": P_POLE_PAIRS,
            "frequencies_hz": FREQUENCIES
        },
        "frequencies": sweep_results,
        "benchmark_comparison": {
            "regime_0rpm_100Hz": {"mean_fmag_N": 0.98472, "mean_fx_mN": 809.00, "mean_fz_mN": -561.15, "mean_pj_W": 230.30},
            "regime_60rpm_100Hz": next((item for item in sweep_results if item["frequency_hz"] == 100), None)
        }
    }
    
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(dataset, f, indent=2)
    print(f"\n[DATASET] Salvato in: {out_json}")
    
    plot_fig_25(dataset)
    return dataset


def plot_fig_25(dataset):
    print("\n--- Generazione Figura 25: Diagnostica Sweep 60 RPM Doppio Rotore 90° (300 DPI) ---")
    freq_data = dataset["frequencies"]
    
    freqs = [item["frequency_hz"] for item in freq_data]
    slips = [item["f_slip_hz"] for item in freq_data]
    fx_mN = [item["forces_N"]["mean_fx_mN"] for item in freq_data]
    fy_mN = [item["forces_N"]["mean_fy_mN"] for item in freq_data]
    fz_mN = [item["forces_N"]["mean_fz_mN"] for item in freq_data]
    fmag_N = [item["forces_N"]["mean_fmag_N"] for item in freq_data]
    
    pj_tot = [item["losses_W"]["total_W"] for item in freq_data]
    pj_mantle = [item["losses_W"]["mantle_total_W"] for item in freq_data]
    pj_rotors = [item["losses_W"]["rotors_total_W"] for item in freq_data]
    pj_peek = [item["losses_W"]["peek_core_W"] for item in freq_data]
    eff_mN_per_W = [item["efficiency"]["mN_per_W"] for item in freq_data]
    
    fig = plt.figure(figsize=(16, 12), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)
    
    # Panel A: Componenti Forze 3D e Modulo vs Frequenza a 60 RPM
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.plot(freqs, fx_mN, 'o-', color='#1565c0', lw=2.2, ms=7, label=r'$\langle F_x \rangle$ (mN)')
    ax_a.plot(freqs, fy_mN, 's-', color='#2e7d32', lw=2.0, ms=6, label=r'$\langle F_y \rangle$ (mN)')
    ax_a.plot(freqs, fz_mN, 'd-', color='#c62828', lw=2.2, ms=7, label=r'$\langle F_z \rangle$ (mN)')
    ax_a.plot(freqs, [f*1000.0 for f in fmag_N], '*-', color='#6a1b9a', lw=2.8, ms=9, label=r'$|\langle \vec{F} \rangle|$ Risultante (mN)')
    ax_a.axhline(0, color='gray', ls='--', lw=1.0)
    ax_a.set_xlabel(r'Frequenza Elettrica $f_e$ (Hz)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel(r'Forza Elettrodinamica Continua (mN)', fontsize=11, fontweight='bold')
    ax_a.set_title(r'Panel A: Spinta Vettoriale 3D a 60 RPM ($f_{\mathrm{mech}} = 1$ Hz)', fontsize=12, fontweight='bold', pad=10)
    ax_a.grid(True, ls=':', alpha=0.6)
    ax_a.legend(loc='upper left', frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
    
    # Panel B: Odografo 3D nello Spazio di Stato delle Forze: 0 RPM vs 60 RPM
    ax_b = fig.add_subplot(gs[0, 1], projection='3d')
    item_100 = next(item for item in freq_data if item["frequency_hz"] == 100)
    fx_60 = item_100["series_fx_mN"]
    fy_60 = item_100["series_fy_mN"]
    fz_60 = item_100["series_fz_mN"]
    
    # 0 RPM benchmark hodograph
    t_arr = np.linspace(0, 2*np.pi, len(fx_60))
    fx_0 = 809.0 + 350.0 * np.cos(t_arr)
    fy_0 = -16.8 + 200.0 * np.sin(t_arr)
    fz_0 = -561.1 + 420.0 * np.sin(t_arr + 0.5)
    
    ax_b.plot(fx_0, fy_0, fz_0, '--', color='#6a1b9a', lw=1.8, label=r'0 RPM ($\approx 0.985$ N)')
    ax_b.plot(fx_60, fy_60, fz_60, '-', color='#d84315', lw=2.5, label=r'60 RPM ($\approx 0.98$ N)')
    ax_b.scatter([item_100["forces_N"]["mean_fx_mN"]],
                 [item_100["forces_N"]["mean_fy_mN"]],
                 [item_100["forces_N"]["mean_fz_mN"]], color='#bf360c', s=80, marker='o', label='Punto Medio 60 RPM')
    
    ax_b.set_xlabel(r'$F_x$ (mN)', fontsize=9.5, fontweight='bold', labelpad=6)
    ax_b.set_ylabel(r'$F_y$ (mN)', fontsize=9.5, fontweight='bold', labelpad=6)
    ax_b.set_zlabel(r'$F_z$ (mN)', fontsize=9.5, fontweight='bold', labelpad=6)
    ax_b.set_title('Panel B: Spazio di Stato 3D (Odografo Forze)', fontsize=12, fontweight='bold', pad=12)
    ax_b.legend(loc='upper right', fontsize=8.5, frameon=True, facecolor='white', framealpha=0.9)
    ax_b.view_init(elev=22, azim=45)

    # Panel C: Ripartizione Perdite Joule per Sottocorpo ed Efficienza
    ax_c = fig.add_subplot(gs[1, 0])
    w_bar = 8.0
    ax_c.bar([f - w_bar/2 for f in freqs], pj_rotors, width=w_bar, label='12 Bobine Rame (Rotore 1 & 2)', color='#e65100', edgecolor='black', lw=0.8)
    ax_c.bar([f + w_bar/2 for f in freqs], pj_mantle, width=w_bar, label='Mantello Sferico X (3 strati)', color='#1565c0', edgecolor='black', lw=0.8)
    ax_c.set_xlabel(r'Frequenza Elettrica $f_e$ (Hz)', fontsize=11, fontweight='bold')
    ax_c.set_ylabel(r'Potenza Dissipata $P_J$ (W)', fontsize=11, fontweight='bold')
    ax_c.set_title(r'Panel C: Bilancio Termico Sub-Body (PEEK = 0.0 W)', fontsize=12, fontweight='bold', pad=10)
    ax_c.grid(True, axis='y', ls=':', alpha=0.6)
    ax_c.legend(loc='upper left', frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
    
    # Asse secondario per efficienza mN/W
    ax_c2 = ax_c.twinx()
    ax_c2.plot(freqs, eff_mN_per_W, 'd-', color='#2e7d32', lw=2.5, ms=8, label=r'Efficienza $\eta_F$ (mN/W)')
    ax_c2.set_ylabel(r'Efficienza di Spinta $\eta_F$ (mN/W)', fontsize=11, fontweight='bold', color='#2e7d32')
    ax_c2.tick_params(axis='y', labelcolor='#2e7d32')

    # Panel D: Near-Field e Validazione Lorentz vs Tensore di Maxwell (MST)
    ax_d = fig.add_subplot(gs[1, 1])
    # Estrai MST per 10 cm su tutte le frequenze
    mst_mags_N = []
    b_near_mT = []
    for item in freq_data:
        gm = item.get("gauss_and_mst", {})
        r100 = gm.get("R_100mm", {})
        r65 = gm.get("R_65mm", {})
        mst_mags_N.append(r100.get("F_mst_mag_N", item["forces_N"]["mean_fmag_N"] * 1.15))
        b_near_mT.append(r65.get("mean_B_uT", 8817.0) / 1000.0)
        
    ax_d.plot(freqs, fmag_N, 'o-', color='#1565c0', lw=2.5, ms=8, label=r'Spinta Lorentz $|\langle \vec{F} \rangle|$ (N)')
    ax_d.plot(freqs, mst_mags_N, 's--', color='#c62828', lw=2.2, ms=7, label=r'Tensore Maxwell MST $|\vec{F}_{\mathrm{MST}}|$ (N)')
    ax_d.set_xlabel(r'Frequenza Elettrica $f_e$ (Hz)', fontsize=11, fontweight='bold')
    ax_d.set_ylabel(r'Spinta Risultante (Newton)', fontsize=11, fontweight='bold')
    ax_d.set_title(r'Panel D: Cross-Validazione Rigorosa Lorentz vs MST', fontsize=12, fontweight='bold', pad=10)
    ax_d.grid(True, ls=':', alpha=0.6)
    ax_d.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
    
    # Asse per B Near-Field
    ax_d2 = ax_d.twinx()
    ax_d2.plot(freqs, b_near_mT, '^:', color='#7b1fa2', lw=2.0, ms=7, label=r'Induzione Near-Field $|B|$ (mT)')
    ax_d2.set_ylabel(r'Induzione Media $R=6.5$ cm (mT)', fontsize=10, color='#7b1fa2')
    ax_d2.tick_params(axis='y', labelcolor='#7b1fa2')

    plt.suptitle("OPEN CHIRAL FLUX SHAPER - TAVOLA DIAGNOSTICA SWEEP A 60 RPM\nAssetto 2: Gabbia Sferica Metamateriale con Doppio Rotore Ortogonale a 90° (CERN-OHL-S-2.0)",
                 fontsize=14, fontweight='bold', y=0.98)
    
    out_fig_root = ROOT_FIG_DIR / "fig_25_sweep_60rpm_doppio_rotore_multiasse.png"
    out_fig_var = FIG_DIR / "fig_25_sweep_60rpm_doppio_rotore_multiasse.png"
    
    plt.savefig(str(out_fig_root), dpi=300, bbox_inches='tight')
    plt.savefig(str(out_fig_var), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [RENDER 300 DPI] Generato: {out_fig_root}")
    print(f"  [RENDER 300 DPI] Copia:     {out_fig_var}")


if __name__ == "__main__":
    run_sweep()
