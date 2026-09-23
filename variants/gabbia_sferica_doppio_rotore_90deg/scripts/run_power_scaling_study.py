#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico di Scalatura di Potenza, Verifica di Saturazione e Audit Termico Radiativo
Configurazione: Gabbia Sferica Metamateriale a Doppio Rotore Ortogonale a 90°
Eccitazione: NPNPNP Trifase a Terzi Continua (120°) con Quadratura Temporale 90° tra Rotori Z ed X
Frequenza: 100 Hz, 64 timestep (dt = 0.25 ms, T_tot = 16.0 ms, >1.6 cicli completi)

Punti di corrente analizzati:
  1. J0 = 1.0e5 A/m^2 (Baseline: P_J ~ 230 W, F ~ 0.985 N)
  2. J0 = 1.5e5 A/m^2 (Intermedio: P_J ~ 518 W, F ~ 2.22 N)
  3. J0 = 2.0e5 A/m^2 (Raddoppio Corrente: P_J ~ 921 W, F ~ 3.94 N)

Verifiche fisiche:
  - Legge di proporzionalità quadratica F ~ J0^2 e F ~ P_elec
  - Margine di saturazione ferromagnetica B_max vs B_sat (~1.5 T)
  - Dissipazione Joule per sottocorpo (dissezione mantello a 3 strati e schermatura esterna -30°)
  - Bilancio termico radiativo di Stefan-Boltzmann nel vuoto profondo (T_eq e area radiatore ausiliario)
  - Conservazione del flusso di Gauss div(B) = 0 su sfere di Fibonacci a 2500 punti e MST

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
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

# Percorsi
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent.parent.parent
ROOT_FIGURES_DIR = ROOT_DIR / "figures"
ROOT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

VAR_DIR = SCRIPT_DIR.parent
BASE_CONFIG_SIF = VAR_DIR / "config" / "case_doppio_rotore_npnpnp_terzi.sif"
MESH_DIR = VAR_DIR / "mesh"
MESH_NAME = "macchina_sferica_ortogonale"
DATA_DIR = VAR_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON = DATA_DIR / "power_scaling_study.json"

FIG_VAR_DIR = VAR_DIR / "figures"
FIG_VAR_DIR.mkdir(parents=True, exist_ok=True)

ELMER_SOLVER = r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"
MU0 = 4.0 * np.pi * 1e-7
SIGMA_SB = 5.670374419e-8  # Stefan-Boltzmann W / (m^2 K^4)
EMISSIVITY_MANTLE = 0.85
R_EXT_SPHERE = 0.050  # 50 mm
A_RAD_SPHERE = 4.0 * np.pi * (R_EXT_SPHERE**2)  # 0.031416 m^2

TIMESTEPS = 64
DT = 0.00025  # 0.25 ms
F_HZ = 100.0
N_FIBONACCI = 2500
RADII = [0.065, 0.100, 0.150]
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]

# Punti di scalatura
SCALING_CASES = [
    {
        "id": "j10_baseline",
        "label": "J0 = 1.0x10^5 A/m² (1.0x Nominale)",
        "scale_factor": 1.0,
        "j0_val": 1.0e5,
        "work_dir": VAR_DIR / "work_dirs" / "run_npnpnp_doppio_rotore"
    },
    {
        "id": "j15_intermediate",
        "label": "J0 = 1.5x10^5 A/m² (1.5x Intermedio)",
        "scale_factor": 1.5,
        "j0_val": 1.5e5,
        "work_dir": VAR_DIR / "work_dirs" / "power_scaling_j15"
    },
    {
        "id": "j20_doubled",
        "label": "J0 = 2.0x10^5 A/m² (2.0x Raddoppio)",
        "scale_factor": 2.0,
        "j0_val": 2.0e5,
        "work_dir": VAR_DIR / "work_dirs" / "power_scaling_j20"
    }
]


def generate_sif_for_j0(base_sif_path, dest_sif_path, j0_val, mesh_dir, mesh_name):
    """Crea una configurazione Elmer SIF modificata con la corrente j0 specificata."""
    text = base_sif_path.read_text(encoding="utf-8")
    lines = text.splitlines()
    new_lines = []
    
    mesh_db_rel = os.path.relpath(str(mesh_dir), str(dest_sif_path.parent)).replace('\\', '/')
    
    for line in lines:
        if 'Mesh DB' in line:
            new_lines.append(f'  Mesh DB "{mesh_db_rel}" "{mesh_name}"')
        elif 'Results Directory' in line:
            new_lines.append('  Results Directory "results"')
        elif 'j0 = 1.0e+05;' in line:
            new_lines.append(f'  rc = 0.035; rwsq = 2.5e-05; j0 = {j0_val:.4e};\\')
        else:
            new_lines.append(line)
            
    dest_sif_path.write_text("\n".join(new_lines), encoding="utf-8")


def run_simulation_if_needed(case):
    """Esegue ElmerSolver per il caso specificato se i VTU non sono presenti."""
    work_dir = case["work_dir"]
    res_dir = work_dir / "results"
    work_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)
    
    vtus = sorted(list(res_dir.glob("macchina_sferica_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] {case['label']}: Trovati {len(vtus)} VTU in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]
        
    print(f"\n{'='*75}")
    print(f"  [SOLVER] Esecuzione ElmerSolver per {case['label']}")
    print(f"  [DIR] Work Dir: {work_dir}")
    print(f"{'='*75}")
    
    for f in res_dir.glob("macchina_sferica_out_t*.vtu"):
        f.unlink()
        
    sif_path = work_dir / "case.sif"
    generate_sif_for_j0(BASE_CONFIG_SIF, sif_path, case["j0_val"], MESH_DIR, MESH_NAME)
    (work_dir / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")
    
    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    proc = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(work_dir), env=env, capture_output=True, text=True)
    duration = time.time() - t0
    
    if proc.returncode != 0:
        print(f"  [ERRORE] ElmerSolver fallito per {case['label']} (codice {proc.returncode}):")
        print("STDERR:\n", proc.stderr[-1000:])
        print("STDOUT:\n", proc.stdout[-1000:])
        sys.exit(1)
        
    vtus = sorted(list(res_dir.glob("macchina_sferica_out_t*.vtu")))
    print(f"  [SOLVER] Completato in {duration:.1f}s ({len(vtus)} VTU generati).")
    if len(vtus) < TIMESTEPS:
        print(f"  [ERRORE] Generati solo {len(vtus)} VTU su {TIMESTEPS} attesi.")
        sys.exit(1)
    return vtus[:TIMESTEPS]


def generate_fibonacci_sphere(n_points, radius):
    """Genera una distribuzione sferica uniforme di Fibonacci a N punti."""
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * indices / n_points)
    theta = np.pi * (1.0 + 5.0**0.5) * indices
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)
    points = np.column_stack([x, y, z])
    normals = points / radius
    area_weight = 4.0 * np.pi * (radius**2) / n_points
    return points, normals, area_weight


def analyze_fibonacci_spheres(vtus, sample_idx=32):
    """Calcola Gauss e MST su sfere di Fibonacci."""
    vtu_file = vtus[sample_idx]
    m = meshio.read(str(vtu_file))
    pts = m.points
    b_field = m.point_data['magnetic flux density']
    
    delaunay_tri = Delaunay(pts)
    interp_b = LinearNDInterpolator(delaunay_tri, b_field, fill_value=0.0)
    
    results = {}
    for r_val, label in zip(RADII, RADII_LABELS):
        pts_fib, normals, d_area = generate_fibonacci_sphere(N_FIBONACCI, r_val)
        b_samples = interp_b(pts_fib)
        
        b_dot_n = np.einsum('ij,ij->i', b_samples, normals)
        flux_net = float(np.sum(b_dot_n) * d_area)
        flux_abs = float(np.sum(np.abs(b_dot_n)) * d_area)
        residual_pct = float(abs(flux_net) / (flux_abs + 1e-15) * 100.0)
        
        b_mag = np.linalg.norm(b_samples, axis=1)
        
        t_mst_x = (b_samples[:, 0] * b_dot_n - 0.5 * (b_mag**2) * normals[:, 0]) / MU0
        t_mst_y = (b_samples[:, 1] * b_dot_n - 0.5 * (b_mag**2) * normals[:, 1]) / MU0
        t_mst_z = (b_samples[:, 2] * b_dot_n - 0.5 * (b_mag**2) * normals[:, 2]) / MU0
        
        fx_mst = float(np.sum(t_mst_x) * d_area)
        fy_mst = float(np.sum(t_mst_y) * d_area)
        fz_mst = float(np.sum(t_mst_z) * d_area)
        f_mst_mag = float(np.sqrt(fx_mst**2 + fy_mst**2 + fz_mst**2))
        
        results[label] = {
            "radius_m": r_val,
            "gauss_net_weber": flux_net,
            "gauss_abs_weber": flux_abs,
            "gauss_residual_pct": residual_pct,
            "gauss_pass": bool(residual_pct < 2.0),
            "b_mean_uT": float(np.mean(b_mag) * 1e6),
            "b_max_uT": float(np.max(b_mag) * 1e6),
            "mst_fx_N": fx_mst,
            "mst_fy_N": fy_mst,
            "mst_fz_N": fz_mst,
            "mst_mag_N": f_mst_mag
        }
    return results


def process_case(case, vtus):
    """Elabora le forze di Lorentz, perdite Joule e saturazione magnetica per un caso."""
    print(f"\n  [ANALYSIS] Elaborazione dati per {case['label']}...")
    m0 = meshio.read(str(vtus[0]))
    pts = m0.points
    cells = None
    for cb in m0.cells:
        if cb.type == 'tetra':
            cells = cb.data
            break
            
    v0, v1, v2, v3 = pts[cells[:, 0]], pts[cells[:, 1]], pts[cells[:, 2]], pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0
    elem_com = (v0 + v1 + v2 + v3) / 4.0
    x, y, z = elem_com[:, 0], elem_com[:, 1], elem_com[:, 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    rho = np.sqrt(x**2 + y**2)
    rho_yz = np.sqrt(y**2 + z**2)
    
    # Maschere geometriche
    mask_mantle = (r >= 0.0465) & (r <= 0.0505)
    mask_l1 = mask_mantle & (r < 0.0480)
    mask_l2 = mask_mantle & (r >= 0.0480) & (r < 0.0490)
    mask_l3 = mask_mantle & (r >= 0.0490)
    mask_core = r <= 0.0125
    mask_coils1 = (r < 0.0465) & (r > 0.025) & (np.abs(z) <= 0.008) & (rho >= 0.028) & (rho <= 0.042)
    mask_coils2 = (r < 0.0465) & (r > 0.025) & (np.abs(x) <= 0.008) & (rho_yz >= 0.028) & (rho_yz <= 0.042)
    mask_assembly = mask_mantle | mask_coils1 | mask_coils2 | mask_core
    
    # Maschere nodali per induzione B
    pts_r = np.sqrt(pts[:, 0]**2 + pts[:, 1]**2 + pts[:, 2]**2)
    pts_mantle = (pts_r >= 0.0465) & (pts_r <= 0.0505)
    
    masks = {
        "mantle_total": mask_mantle,
        "layer1_plus30": mask_l1,
        "layer2_ortho": mask_l2,
        "layer3_minus30": mask_l3,
        "rotor1_coils": mask_coils1,
        "rotor2_coils": mask_coils2,
        "peek_core": mask_core,
        "total": mask_assembly
    }
    
    forces = {k: {"fx": [], "fy": [], "fz": []} for k in masks}
    joule = {k: [] for k in masks}
    
    b_mantle_max_series = []
    b_mantle_mean_series = []
    b_global_max_series = []
    
    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        jxb = m.point_data['jxb']
        pj = m.point_data['joule heating'].ravel()
        b_pts = m.point_data['magnetic flux density']
        b_mag_pts = np.linalg.norm(b_pts, axis=1)
        
        b_mantle_max_series.append(float(np.max(b_mag_pts[pts_mantle])))
        b_mantle_mean_series.append(float(np.mean(b_mag_pts[pts_mantle])))
        b_global_max_series.append(float(np.max(b_mag_pts)))
        
        jxb_elem = np.mean(jxb[cells], axis=1)
        pj_elem = np.mean(pj[cells], axis=1)
        
        for k, mask in masks.items():
            if np.sum(mask) > 0:
                vols = elem_vols[mask]
                fx_val = float(np.sum(vols * jxb_elem[mask, 0]))
                fy_val = float(np.sum(vols * jxb_elem[mask, 1]))
                fz_val = float(np.sum(vols * jxb_elem[mask, 2]))
                pj_val = float(np.sum(vols * pj_elem[mask]))
            else:
                fx_val, fy_val, fz_val, pj_val = 0.0, 0.0, 0.0, 0.0
                
            forces[k]["fx"].append(fx_val)
            forces[k]["fy"].append(fy_val)
            forces[k]["fz"].append(fz_val)
            joule[k].append(pj_val)
            
    stats = {}
    for k in masks:
        fx = np.array(forces[k]["fx"])
        fy = np.array(forces[k]["fy"])
        fz = np.array(forces[k]["fz"])
        pj = np.array(joule[k])
        stats[k] = {
            "mean_fx_N": float(np.mean(fx)),
            "mean_fy_N": float(np.mean(fy)),
            "mean_fz_N": float(np.mean(fz)),
            "peak_fx_N": float(np.max(fx)),
            "peak_fy_N": float(np.max(fy)),
            "peak_fz_N": float(np.max(fz)),
            "std_fx_N": float(np.std(fx)),
            "std_fy_N": float(np.std(fy)),
            "std_fz_N": float(np.std(fz)),
            "mean_pj_W": float(np.mean(pj)),
            "peak_pj_W": float(np.max(pj)),
            "series_fx_N": fx.tolist(),
            "series_fy_N": fy.tolist(),
            "series_fz_N": fz.tolist(),
            "series_pj_W": pj.tolist()
        }
        
    tot = stats["total"]
    f_mean_vec = np.array([tot["mean_fx_N"], tot["mean_fy_N"], tot["mean_fz_N"]])
    f_mag_mean_N = float(np.linalg.norm(f_mean_vec))
    f_peak_inst_N = float(np.max(np.sqrt(np.array(tot["series_fx_N"])**2 + np.array(tot["series_fy_N"])**2 + np.array(tot["series_fz_N"])**2)))
    
    total_pj_W = tot["mean_pj_W"]
    eta_F_mN_per_W = (f_mag_mean_N * 1000.0) / (total_pj_W + 1e-12)
    
    # Bilancio termico radiativo nel vuoto profondo (Stefan-Boltzmann)
    # T_eq^4 = P_J / (eps * sigma_SB * A_rad)
    t_eq_kelvin = (total_pj_W / (EMISSIVITY_MANTLE * SIGMA_SB * A_RAD_SPHERE))**0.25
    t_eq_celsius = t_eq_kelvin - 273.15
    
    # Area radiante totale richiesta per stabilizzare a T <= 100 °C (373.15 K) e T <= 80 °C (353.15 K)
    t_target_100k = 373.15
    a_req_100c = total_pj_W / (EMISSIVITY_MANTLE * SIGMA_SB * (t_target_100k**4))
    a_aux_100c = max(0.0, a_req_100c - A_RAD_SPHERE)
    
    t_target_80k = 353.15
    a_req_80c = total_pj_W / (EMISSIVITY_MANTLE * SIGMA_SB * (t_target_80k**4))
    a_aux_80c = max(0.0, a_req_80c - A_RAD_SPHERE)
    
    # Margine di saturazione rispetto a B_sat = 1.5 T
    b_sat_ref = 1.50  # Tesla
    b_mantle_peak_abs = float(np.max(b_mantle_max_series))
    b_mantle_mean_abs = float(np.mean(b_mantle_mean_series))
    b_global_peak_abs = float(np.max(b_global_max_series))
    sat_margin_pct = float((1.0 - (b_mantle_peak_abs / b_sat_ref)) * 100.0)
    
    # Fibonacci spheres
    fib_res = analyze_fibonacci_spheres(vtus, sample_idx=32)
    
    case_summary = {
        "id": case["id"],
        "label": case["label"],
        "scale_factor": case["scale_factor"],
        "j0_A_m2": case["j0_val"],
        "thrust": {
            "mean_fx_N": tot["mean_fx_N"],
            "mean_fy_N": tot["mean_fy_N"],
            "mean_fz_N": tot["mean_fz_N"],
            "resultant_mag_N": f_mag_mean_N,
            "peak_instantaneous_N": f_peak_inst_N,
            "mst_mag_midfield_N": fib_res["Mid-Field (R=10.0 cm)"]["mst_mag_N"]
        },
        "power_and_efficiency": {
            "total_joule_W": total_pj_W,
            "rotor1_W": stats["rotor1_coils"]["mean_pj_W"],
            "rotor2_W": stats["rotor2_coils"]["mean_pj_W"],
            "mantle_W": stats["mantle_total"]["mean_pj_W"],
            "mantle_l1_plus30_W": stats["layer1_plus30"]["mean_pj_W"],
            "mantle_l2_ortho_W": stats["layer2_ortho"]["mean_pj_W"],
            "mantle_l3_minus30_W": stats["layer3_minus30"]["mean_pj_W"],
            "peek_core_W": stats["peek_core"]["mean_pj_W"],
            "thrust_efficiency_mN_per_W": eta_F_mN_per_W
        },
        "magnetic_saturation": {
            "b_mantle_peak_T": b_mantle_peak_abs,
            "b_mantle_mean_T": b_mantle_mean_abs,
            "b_global_peak_T": b_global_peak_abs,
            "b_sat_ref_T": b_sat_ref,
            "saturation_margin_pct": sat_margin_pct,
            "linear_regime_pass": bool(b_mantle_peak_abs < b_sat_ref)
        },
        "thermal_vacuum_equilibrium": {
            "t_eq_kelvin": float(t_eq_kelvin),
            "t_eq_celsius": float(t_eq_celsius),
            "mantle_radiating_area_m2": float(A_RAD_SPHERE),
            "mantle_emissivity": float(EMISSIVITY_MANTLE),
            "area_req_100C_m2": float(a_req_100c),
            "area_aux_100C_m2": float(a_aux_100c),
            "area_req_80C_m2": float(a_req_80c),
            "area_aux_80C_m2": float(a_aux_80c)
        },
        "fibonacci_spheres": fib_res,
        "detailed_stats": stats
    }
    
    print(f"    - Spinta Risultante |<F>| : {f_mag_mean_N:.4f} N (Picco: {f_peak_inst_N:.2f} N)")
    print(f"    - Potenza Joule Totale   : {total_pj_W:.2f} W (Efficienza: {eta_F_mN_per_W:.2f} mN/W)")
    print(f"    - Induzione Mantello B    : Peak={b_mantle_peak_abs*1e3:.1f} mT, Mean={b_mantle_mean_abs*1e3:.1f} mT (Margine Sat: {sat_margin_pct:.1f}%)")
    print(f"    - Equilibrio Vuoto T_eq  : {t_eq_kelvin:.1f} K ({t_eq_celsius:.1f} °C) | Area Aux (100°C)={a_aux_100c:.3f} m²")
    print(f"    - Gauss Residuo (Mid)    : {fib_res['Mid-Field (R=10.0 cm)']['gauss_residual_pct']:.3f}% ({'PASS' if fib_res['Mid-Field (R=10.0 cm)']['gauss_pass'] else 'FAIL'})")
    return case_summary


def render_diagnostic_plate_fig26(results, out_fig_paths):
    """Genera la tavola comparativa a 4 pannelli a 300 DPI nativi."""
    print("\n  [PLOTTING] Rendering Tavola Diagnostica a 300 DPI (Fig. 26)...")
    
    plt.rcParams.update({
        'font.size': 10,
        'font.family': 'sans-serif',
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.titlesize': 14,
        'axes.grid': True,
        'grid.alpha': 0.4,
        'grid.linestyle': '--'
    })
    
    fig = plt.figure(figsize=(16, 12), dpi=300)
    gs = gridspec.GridSpec(2, 2, hspace=0.28, wspace=0.25)
    
    scales = [r["scale_factor"] for r in results]
    j0_vals = [r["j0_A_m2"] / 1e5 for r in results]
    f_mag = [r["thrust"]["resultant_mag_N"] for r in results]
    fx_vals = [r["thrust"]["mean_fx_N"] for r in results]
    fz_vals = [r["thrust"]["mean_fz_N"] for r in results]
    f_peak = [r["thrust"]["peak_instantaneous_N"] for r in results]
    f_mst = [r["thrust"]["mst_mag_midfield_N"] for r in results]
    pj_vals = [r["power_and_efficiency"]["total_joule_W"] for r in results]
    eta_vals = [r["power_and_efficiency"]["thrust_efficiency_mN_per_W"] for r in results]
    
    # ----------------------------------------------------
    # Panel A: Spinta Risultante e Componenti vs Scalatura di Corrente e Potenza
    # ----------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    
    # Teoria quadratica ideale basata sul punto 1.0x
    j_cont = np.linspace(0.8, 2.2, 100)
    f_theory = f_mag[0] * (j_cont / scales[0])**2
    
    ax_a.plot(j_cont, f_theory, 'k--', label=r'Teoria Quadratica Ideale $\propto J_0^2$', alpha=0.7, lw=1.8)
    ax_a.plot(scales, f_mag, 'o-', color='#1f77b4', lw=2.5, markersize=8, label=r'Risultante $|\langle\vec{F}\rangle|$ FEM (Lorentz)')
    ax_a.plot(scales, f_mst, 's--', color='#17becf', lw=1.8, markersize=7, label=r'Tensore di Maxwell $|\vec{F}_{\mathrm{MST}}|$')
    ax_a.plot(scales, fx_vals, 'd-', color='#2ca02c', lw=1.8, markersize=6, label=r'$\langle F_x \rangle$ (Spinta Longitudinale)')
    ax_a.plot(scales, fz_vals, 'v-', color='#d62728', lw=1.8, markersize=6, label=r'$\langle F_z \rangle$ (Spinta Verticale)')
    
    for s, f_val, p_val in zip(scales, f_mag, pj_vals):
        ax_a.annotate(f"{f_val:.3f} N\n({p_val:.1f} W)",
                     (s, f_val), textcoords="offset points", xytext=(0, 12),
                     ha='center', fontsize=9, fontweight='bold',
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#1f77b4", alpha=0.85))
                     
    ax_a.set_title("PANEL A: Scalatura Spinta Vettoriale 3D vs Corrente ($J_0$)", fontweight='bold')
    ax_a.set_xlabel(r"Fattore di Scalatura Corrente $J_0 / J_{0,\mathrm{nom}}$ ($1.0 \times 10^5\ \mathrm{A/m}^2$)")
    ax_a.set_ylabel("Spinta Stazionaria [N]")
    ax_a.set_xlim(0.8, 2.2)
    ax_a.set_ylim(-1.5, 4.6)
    ax_a.legend(loc='upper left', framealpha=0.9)
    
    # ----------------------------------------------------
    # Panel B: Verifica di Saturazione Magnetica & Margine B_sat
    # ----------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    
    b_peak_mantle = [r["magnetic_saturation"]["b_mantle_peak_T"] for r in results]
    b_mean_mantle = [r["magnetic_saturation"]["b_mantle_mean_T"] for r in results]
    b_peak_global = [r["magnetic_saturation"]["b_global_peak_T"] for r in results]
    b_sat = results[0]["magnetic_saturation"]["b_sat_ref_T"]
    
    # Curva di magnetizzazione B-H indicativa
    h_norm = np.linspace(0, 3.0, 200)
    # B = Bsat * tanh(mu_r * mu0 * H / Bsat)
    b_sat_curve = b_sat * np.tanh(1.2 * h_norm)
    ax_b.plot(h_norm, b_sat_curve, 'k:', label=r'Curva Saturazione Metamateriale ($B_{\mathrm{sat}} = 1.50\ \mathrm{T}$)', alpha=0.6, lw=2.0)
    
    ax_b.axhline(b_sat, color='#d62728', linestyle='--', lw=2.0, label=r'Soglia Saturazione Fisica $B_{\mathrm{sat}} = 1.5\ \mathrm{T}$')
    ax_b.plot(scales, b_peak_mantle, 'o-', color='#e377c2', lw=2.2, markersize=8, label=r'$|B|_{\mathrm{max}}$ Mantello Sferico (Local Peak)')
    ax_b.plot(scales, b_mean_mantle, 's-', color='#9467bd', lw=2.2, markersize=7, label=r'$\langle |B| \rangle$ Mantello Sferico (Media)')
    ax_b.plot(scales, b_peak_global, '^--', color='#ff7f0e', lw=1.8, markersize=7, label=r'$|B|_{\mathrm{max}}$ Traferro / Interno Bobine')
    
    for s, b_pk, b_mn in zip(scales, b_peak_mantle, b_mean_mantle):
        margin = (1.0 - b_pk / b_sat) * 100.0
        ax_b.annotate(f"Peak: {b_pk*1e3:.1f} mT\nMargine: {margin:.1f}%",
                     (s, b_pk), textcoords="offset points", xytext=(0, 12),
                     ha='center', fontsize=8.5, fontweight='bold',
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#e377c2", alpha=0.85))
                     
    ax_b.set_title("PANEL B: Monitoraggio Induzione e Margine di Saturazione Magnetica", fontweight='bold')
    ax_b.set_xlabel(r"Fattore di Scalatura Corrente $J_0 / J_{0,\mathrm{nom}}$")
    ax_b.set_ylabel("Induzione Magnetica $B$ [T]")
    ax_b.set_xlim(0.8, 2.2)
    ax_b.set_ylim(0.0, 1.8)
    ax_b.legend(loc='center right', framealpha=0.9)
    
    # ----------------------------------------------------
    # Panel C: Ripartizione Dissipazione Joule per Sottocorpo ed Efficienza
    # ----------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    
    bar_width = 0.22
    x_indices = np.arange(len(results))
    
    p_r1 = [r["power_and_efficiency"]["rotor1_W"] for r in results]
    p_r2 = [r["power_and_efficiency"]["rotor2_W"] for r in results]
    p_l1 = [r["power_and_efficiency"]["mantle_l1_plus30_W"] for r in results]
    p_l2 = [r["power_and_efficiency"]["mantle_l2_ortho_W"] for r in results]
    p_l3 = [r["power_and_efficiency"]["mantle_l3_minus30_W"] for r in results]
    p_core = [r["power_and_efficiency"]["peek_core_W"] for r in results]
    
    # Barre impilate
    b1 = ax_c.bar(x_indices, p_r1, bar_width, label='Rotore 1 (6 Bobine Z)', color='#1f77b4', edgecolor='black')
    b2 = ax_c.bar(x_indices, p_r2, bar_width, bottom=p_r1, label='Rotore 2 (6 Bobine X)', color='#aec7e8', edgecolor='black')
    
    cum_r = [r1 + r2 for r1, r2 in zip(p_r1, p_r2)]
    b3 = ax_c.bar(x_indices, p_l1, bar_width, bottom=cum_r, label='Mantello Strato 1 (+30°)', color='#ff7f0e', edgecolor='black')
    
    cum_l1 = [c + l1 for c, l1 in zip(cum_r, p_l1)]
    b4 = ax_c.bar(x_indices, p_l2, bar_width, bottom=cum_l1, label='Mantello Strato 2 (0°)', color='#ffbb78', edgecolor='black')
    
    cum_l2 = [c + l2 for c, l2 in zip(cum_l1, p_l2)]
    b5 = ax_c.bar(x_indices, p_l3, bar_width, bottom=cum_l2, label='Mantello Strato 3 (-30° - Schermato: 0 W)', color='#2ca02c', edgecolor='black')
    
    for i, p_tot, eta in zip(x_indices, pj_vals, eta_vals):
        ax_c.annotate(f"Totale: {p_tot:.1f} W\nη = {eta:.2f} mN/W",
                     (i, p_tot), textcoords="offset points", xytext=(0, 10),
                     ha='center', fontsize=9, fontweight='bold')
                     
    ax_c.set_title("PANEL C: Bilancio Joule per Sottocorpo & Efficienza di Spinta $\eta_F$", fontweight='bold')
    ax_c.set_xticks(x_indices)
    ax_c.set_xticklabels([f"{r['label']}\n($J_0 = {r['j0_A_m2']/1e5:.1f}\\times 10^5$)" for r in results])
    ax_c.set_ylabel("Potenza Dissipata $P_J$ [W]")
    ax_c.set_ylim(0, 1150)
    ax_c.legend(loc='upper left', framealpha=0.9)
    
    # ----------------------------------------------------
    # Panel D: Audit Termico Radiativo nel Vuoto Spaziale & Dimensionamento Radiatori
    # ----------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    
    p_sweep = np.linspace(100, 1200, 100)
    # Stefan-Boltzmann: T_eq = (P / (eps * sigma * A))^(1/4)
    t_eq_mantle_sweep_c = ((p_sweep / (EMISSIVITY_MANTLE * SIGMA_SB * A_RAD_SPHERE))**0.25) - 273.15
    
    ax_d.plot(p_sweep, t_eq_mantle_sweep_c, 'k-', lw=2.2, label=r'Radiazione Solo Mantello Sferico ($A = 314\ \mathrm{cm}^2, \epsilon = 0.85$)')
    ax_d.axhline(100.0, color='#2ca02c', linestyle='--', lw=1.8, label=r'Limite Termico Operativo Tecnopolimero ($100^\circ\mathrm{C}$)')
    ax_d.axhline(80.0, color='#17becf', linestyle=':', lw=1.8, label=r'Target Operativo Ottimale Coils ($80^\circ\mathrm{C}$)')
    
    t_celsius = [r["thermal_vacuum_equilibrium"]["t_eq_celsius"] for r in results]
    a_aux_100 = [r["thermal_vacuum_equilibrium"]["area_aux_100C_m2"] for r in results]
    
    ax_d.scatter(pj_vals, t_celsius, color='#d62728', s=80, zorder=5, label='Punti Simulati FEM')
    for p_val, t_c, a_aux in zip(pj_vals, t_celsius, a_aux_100):
        ax_d.annotate(f"{t_c:.1f} °C\n(Req. Aux Radiator: {a_aux:.2f} m²)",
                     (p_val, t_c), textcoords="offset points", xytext=(0, 12),
                     ha='center', fontsize=8.5, fontweight='bold',
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#d62728", alpha=0.85))
                     
    ax_d.set_title("PANEL D: Audit Radiativo nel Vuoto Profondo & Radiatori Ausiliari", fontweight='bold')
    ax_d.set_xlabel("Potenza Dissipata Totale $P_J$ [W]")
    ax_d.set_ylabel(r"Temperatura di Equilibrio nel Vuoto $T_{\mathrm{eq}}$ [°C]")
    ax_d.set_xlim(100, 1150)
    ax_d.set_ylim(0, 750)
    ax_d.legend(loc='lower right', framealpha=0.9)
    
    # Titolo Generale
    fig.suptitle("Gabbia Sferica Metamateriale a Doppio Rotore 90°: Scalatura di Potenza, Margine di Saturazione e Bilancio Termico",
                 fontsize=14, fontweight='bold', y=0.98)
                 
    for p in out_fig_paths:
        fig.savefig(str(p), dpi=300, bbox_inches='tight')
        print(f"  [SAVED] {p}")
    plt.close(fig)


def main():
    print("=" * 80)
    print("STUDIO ELETTRODINAMICO: SCALATURA DI POTENZA & VERIFICA SATURAZIONE")
    print("GABBIA SFERICA METAMATERIALE A DOPPIO ROTORE ORTOGONALE A 90°")
    print("=" * 80)
    
    all_results = []
    
    for case in SCALING_CASES:
        vtus = run_simulation_if_needed(case)
        summary = process_case(case, vtus)
        all_results.append(summary)
        
    # Salvataggio JSON
    OUT_JSON.write_text(json.dumps(all_results, indent=2), encoding="utf-8")
    print(f"\n  [JSON] Dataset salvato in: {OUT_JSON}")
    
    # Rendering 300 DPI
    fig_paths = [
        ROOT_FIGURES_DIR / "fig_26_scalatura_potenza_saturazione.png",
        FIG_VAR_DIR / "fig_26_scalatura_potenza_saturazione.png"
    ]
    render_diagnostic_plate_fig26(all_results, fig_paths)
    
    # Stampa Riepilogo Esecutivo
    print("\n" + "=" * 80)
    print("RIEPILOGO COMPARATIVO DELLA SCALATURA DI POTENZA")
    print("=" * 80)
    header = f"{'Caso':<22} | {'J0 (A/m^2)':<12} | {'|<F>| (N)':<10} | {'F_peak (N)':<10} | {'P_J (W)':<9} | {'eta (mN/W)':<10} | {'B_pk (mT)':<10} | {'Marg. Bsat':<11} | {'T_eq (C)':<9}"
    print(header)
    print("-" * len(header))
    for r in all_results:
        print(f"{r['label'][:22]:<22} | {r['j0_A_m2']:<12.1e} | {r['thrust']['resultant_mag_N']:<10.3f} | {r['thrust']['peak_instantaneous_N']:<10.2f} | {r['power_and_efficiency']['total_joule_W']:<9.1f} | {r['power_and_efficiency']['thrust_efficiency_mN_per_W']:<10.2f} | {r['magnetic_saturation']['b_mantle_peak_T']*1e3:<10.1f} | {r['magnetic_saturation']['saturation_margin_pct']:<10.1f}% | {r['thermal_vacuum_equilibrium']['t_eq_celsius']:<9.1f}")
    print("=" * 80)
    print("Studio completato con successo!")


if __name__ == "__main__":
    main()
