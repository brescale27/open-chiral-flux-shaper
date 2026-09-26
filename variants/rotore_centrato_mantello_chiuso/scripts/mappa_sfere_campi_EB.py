#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mappatura Sferica Completa a 360° del Campo Magnetico B, Elettrico E e Flusso di Poynting S.
Variante Mantello Chiuso a Barattolo (Top & Bottom Lids a Z = ±H/2).

1. Campionamento su 3 superfici sferiche concentriche di Fibonacci (N = 1200 punti per sfera):
   - Sfera Near-Field (R = 6.5 cm)
   - Sfera Mid-Field (R = 10.0 cm)
   - Sfera Far-Field (R = 15.0 cm)
2. Estrazione per ciascuna sfera:
   - Modulo e componenti B: |B|, B_rad, B_theta, B_phi
   - Modulo e componenti E: |E|, E_rad, E_theta, E_phi
   - Flusso attivo del vettore di Poynting: S_rad = (1/mu0) * (E x B)_rad
   - Verifica di conservazione del flusso di Gauss: integrale chiuso di B · n dA (residuo < 2%)
   - Potenza totale irradiata attraverso ciascuna sfera: P_rad = ∮ S_rad dA
3. Generazione grafici a 300 DPI:
   - fig_03_mappatura_sferica_3d_campo_B_E.png: sfere 3D con falsi colori e vettori freccia 360°
   - fig_04_diagramma_radiazione_poynting_360.png: diagramma polare equatoriale e meridiano con effetto di confinamento dei coperchi
4. Salvataggio dati in data/mappatura_sfere_campi_EB.json

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
from scipy.interpolate import LinearNDInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_DIR = SCRIPT_DIR.parent
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
MU0 = 4.0 * np.pi * 1e-7
TIMESTEPS = 20
DT = 0.0005
N_FIBONACCI = 1200
RADII = [0.065, 0.100, 0.150]  # Metri: 6.5 cm, 10.0 cm, 15.0 cm
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]


def ensure_simulation_data(work_sub):
    """Verifica la presenza dei file VTU o esegue ElmerSolver se necessario."""
    res_dir = work_sub / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} VTU in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]
        
    print(f"  [SOLVER] Esecuzione ElmerSolver per case_plus30deg.sif ({TIMESTEPS} timestep)...")
    sif_src = CONFIG_DIR / "case_plus30deg.sif"
    sif_dest = work_sub / "case.sif"
    
    mesh_db_rel = os.path.relpath(str(MESH_DIR), str(work_sub)).replace('\\', '/')
    lines = sif_src.read_text(encoding="utf-8").splitlines()
    new_lines = []
    for line in lines:
        if 'Mesh DB' in line:
            new_lines.append(f'  Mesh DB "{mesh_db_rel}" "macchina_mantello_chiuso"')
        elif 'Results Directory' in line:
            new_lines.append('  Results Directory "results"')
        else:
            new_lines.append(line)
            
    sif_dest.write_text("\n".join(new_lines), encoding="utf-8")
    (work_sub / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")
    
    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    proc = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(work_sub), env=env, capture_output=True, text=True)
    duration = time.time() - t0
    
    if proc.returncode != 0:
        print(f"  [ERRORE] ElmerSolver fallito (codice {proc.returncode}):")
        print(proc.stderr[-1000:])
        sys.exit(1)
        
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    print(f"  [SOLVER] Calcolo completato in {duration:.1f}s ({len(vtus)} file VTU generati).")
    return vtus[:TIMESTEPS]


def generate_fibonacci_sphere(r, n_points=1200):
    """
    Genera una distribuzione quasi-uniforme di punti su una sfera di raggio r
    utilizzando il reticolo sferico di Fibonacci.
    Restituisce coordinate cartesiane, sferiche, versori di base e pesi di area dA.
    """
    i = np.arange(0, n_points, dtype=float) + 0.5
    phi_polar = np.arccos(1.0 - 2.0 * i / n_points)  # theta polare in [0, pi]
    golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0
    theta_azimuth = (2.0 * np.pi * i / golden_ratio) % (2.0 * np.pi)  # phi azimutale in [0, 2pi]
    
    # Coordinate cartesiane
    x = r * np.sin(phi_polar) * np.cos(theta_azimuth)
    y = r * np.sin(phi_polar) * np.sin(theta_azimuth)
    z = r * np.cos(phi_polar)
    coords = np.column_stack([x, y, z])
    
    # Versori ortonormali (r_hat, theta_hat, phi_hat)
    r_hat = coords / r
    
    sin_p = np.sin(phi_polar)
    cos_p = np.cos(phi_polar)
    sin_t = np.sin(theta_azimuth)
    cos_t = np.cos(theta_azimuth)
    
    theta_hat = np.column_stack([cos_p * cos_t, cos_p * sin_t, -sin_p])
    phi_hat = np.column_stack([-sin_t, cos_t, np.zeros_like(sin_t)])
    
    # Area totale = 4 * pi * r^2, area per elemento dA
    da = (4.0 * np.pi * r**2) / n_points
    
    return {
        "r": r,
        "n_points": n_points,
        "coords": coords,
        "x": x, "y": y, "z": z,
        "theta_polar": phi_polar,
        "phi_azimuth": theta_azimuth,
        "r_hat": r_hat,
        "theta_hat": theta_hat,
        "phi_hat": phi_hat,
        "da": da,
        "total_area": 4.0 * np.pi * r**2
    }


def sample_fields_on_sphere(sphere, vtus, pts_mesh):
    """
    Interpola i campi B(t), E(t) e S(t) sui punti sferici per ogni timestep
    e calcola medie temporali, componenti proiettate, residuo di Gauss e potenza di Poynting.
    """
    coords = sphere["coords"]
    r_hat = sphere["r_hat"]
    theta_hat = sphere["theta_hat"]
    phi_hat = sphere["phi_hat"]
    da = sphere["da"]
    n_pts = sphere["n_points"]
    
    # Array per serie temporali sui punti sferici: shape (TIMESTEPS, n_pts, 3)
    b_time = np.zeros((TIMESTEPS, n_pts, 3))
    e_time = np.zeros((TIMESTEPS, n_pts, 3))
    s_time = np.zeros((TIMESTEPS, n_pts, 3))
    
    for t_idx, vtu_file in enumerate(vtus):
        m = meshio.read(str(vtu_file))
        b_nodal = m.point_data["magnetic flux density"]
        e_nodal = m.point_data["electric field"]
        
        # Costruisci interpolatore lineare 3D
        interp_b = LinearNDInterpolator(pts_mesh, b_nodal, fill_value=0.0)
        interp_e = LinearNDInterpolator(pts_mesh, e_nodal, fill_value=0.0)
        
        b_eval = interp_b(coords)
        e_eval = interp_e(coords)
        
        # Flusso di Poynting istantaneo: S = (1/mu0) * (E x B)
        s_eval = np.cross(e_eval, b_eval) / MU0
        
        b_time[t_idx] = b_eval
        e_time[t_idx] = e_eval
        s_time[t_idx] = s_eval

    # 1. Medie temporali
    b_mean = np.mean(b_time, axis=0)       # shape (n_pts, 3)
    e_mean = np.mean(e_time, axis=0)
    s_mean = np.mean(s_time, axis=0)
    
    # Valori RMS temporali
    b_rms = np.sqrt(np.mean(np.sum(b_time**2, axis=-1), axis=0))  # shape (n_pts,)
    e_rms = np.sqrt(np.mean(np.sum(e_time**2, axis=-1), axis=0))
    s_rms = np.sqrt(np.mean(np.sum(s_time**2, axis=-1), axis=0))
    
    # 2. Decomposizione nelle componenti sferiche (r, theta, phi)
    # Su vettore medio temporale
    b_rad = np.sum(b_mean * r_hat, axis=1)
    b_theta = np.sum(b_mean * theta_hat, axis=1)
    b_phi = np.sum(b_mean * phi_hat, axis=1)
    b_mag = np.linalg.norm(b_mean, axis=1)
    
    e_rad = np.sum(e_mean * r_hat, axis=1)
    e_theta = np.sum(e_mean * theta_hat, axis=1)
    e_phi = np.sum(e_mean * phi_hat, axis=1)
    e_mag = np.linalg.norm(e_mean, axis=1)
    
    s_rad = np.sum(s_mean * r_hat, axis=1)
    s_theta = np.sum(s_mean * theta_hat, axis=1)
    s_phi = np.sum(s_mean * phi_hat, axis=1)
    s_mag = np.linalg.norm(s_mean, axis=1)
    
    # 3. Verifica Teorema di Gauss per B (Solenoidalita)
    # Su ciascun timestep e sul campo medio
    phi_net_steps = []
    phi_abs_steps = []
    rel_res_steps = []
    for t_idx in range(TIMESTEPS):
        bn_t = np.sum(b_time[t_idx] * r_hat, axis=1)
        p_net = float(np.sum(bn_t * da))
        p_abs = float(np.sum(np.abs(bn_t) * da))
        res_t = abs(p_net) / (p_abs + 1e-15)
        phi_net_steps.append(p_net)
        phi_abs_steps.append(p_abs)
        rel_res_steps.append(res_t)
        
    phi_net_mean = float(np.sum(b_rad * da))
    phi_abs_mean = float(np.sum(np.abs(b_rad) * da))
    rel_residual_mean = abs(phi_net_mean) / (phi_abs_mean + 1e-15)
    max_step_residual = float(np.max(rel_res_steps))
    mean_step_residual = float(np.mean(rel_res_steps))
    gauss_pass = bool(rel_residual_mean < 0.02)  # Soglia 2%
    
    # 4. Potenza totale irradiata attraverso la sfera (∮ S · n dA)
    p_rad_steps = []
    for t_idx in range(TIMESTEPS):
        sn_t = np.sum(s_time[t_idx] * r_hat, axis=1)
        p_rad_steps.append(float(np.sum(sn_t * da)))
        
    p_rad_total_W = float(np.sum(s_rad * da))
    p_rad_total_uW = p_rad_total_W * 1e6
    p_rad_peak_uW = float(np.max(p_rad_steps) * 1e6)
    
    return {
        "radius_m": sphere["r"],
        "radius_cm": sphere["r"] * 100.0,
        "area_m2": sphere["total_area"],
        "n_points": n_pts,
        # Gauss Conservation
        "gauss": {
            "phi_net_Tm2": phi_net_mean,
            "phi_abs_Tm2": phi_abs_mean,
            "relative_residual": rel_residual_mean,
            "relative_residual_percent": rel_residual_mean * 100.0,
            "mean_step_residual_percent": mean_step_residual * 100.0,
            "max_step_residual_percent": max_step_residual * 100.0,
            "certified_pass": gauss_pass
        },
        # Poynting Power
        "poynting": {
            "p_rad_total_W": p_rad_total_W,
            "p_rad_total_uW": p_rad_total_uW,
            "p_rad_peak_uW": p_rad_peak_uW,
            "mean_S_rad_uW_m2": float(np.mean(s_rad)) * 1e6,
            "max_S_rad_uW_m2": float(np.max(s_rad)) * 1e6,
            "min_S_rad_uW_m2": float(np.min(s_rad)) * 1e6,
            "mean_S_mag_uW_m2": float(np.mean(s_mag)) * 1e6
        },
        # Magnetic Field Statistics
        "B_field": {
            "mean_mag_uT": float(np.mean(b_mag)) * 1e6,
            "max_mag_uT": float(np.max(b_mag)) * 1e6,
            "rms_mag_uT": float(np.mean(b_rms)) * 1e6,
            "mean_Brad_uT": float(np.mean(b_rad)) * 1e6,
            "mean_abs_Brad_uT": float(np.mean(np.abs(b_rad))) * 1e6,
            "mean_Btheta_uT": float(np.mean(b_theta)) * 1e6,
            "mean_Bphi_uT": float(np.mean(b_phi)) * 1e6
        },
        # Electric Field Statistics
        "E_field": {
            "mean_mag_mV_m": float(np.mean(e_mag)) * 1e3,
            "max_mag_mV_m": float(np.max(e_mag)) * 1e3,
            "rms_mag_mV_m": float(np.mean(e_rms)) * 1e3,
            "mean_Erad_mV_m": float(np.mean(e_rad)) * 1e3,
            "mean_Etheta_mV_m": float(np.mean(e_theta)) * 1e3,
            "mean_Ephi_mV_m": float(np.mean(e_phi)) * 1e3
        },
        # Raw arrays per plotting
        "_arrays": {
            "coords": coords,
            "b_mean": b_mean, "b_mag": b_mag, "b_rad": b_rad, "b_rms": b_rms,
            "e_mean": e_mean, "e_mag": e_mag, "e_rms": e_rms,
            "s_mean": s_mean, "s_mag": s_mag, "s_rad": s_rad
        }
    }


def sample_planar_circles(vtus, pts_mesh, r_val=0.10, n_angles=360):
    """
    Campiona i campi sul piano equatoriale (Z=0, phi in [0, 360])
    e sul piano meridiano (X-Z, theta in [0, 360]) per i diagrammi polari.
    """
    angles = np.linspace(0, 2.0 * np.pi, n_angles, endpoint=False)
    
    # 1. Piano Equatoriale (Z = 0)
    x_eq = r_val * np.cos(angles)
    y_eq = r_val * np.sin(angles)
    z_eq = np.zeros_like(angles)
    pts_eq = np.column_stack([x_eq, y_eq, z_eq])
    
    # 2. Piano Meridiano (Y = 0, X = r*sin(alpha), Z = r*cos(alpha))
    # alpha=0: Polo Nord (+Z), alpha=90: Equatore (+X), alpha=180: Polo Sud (-Z), alpha=270: Equatore (-X)
    x_mer = r_val * np.sin(angles)
    y_mer = np.zeros_like(angles)
    z_mer = r_val * np.cos(angles)
    pts_mer = np.column_stack([x_mer, y_mer, z_mer])
    
    s_eq_time = np.zeros((TIMESTEPS, n_angles, 3))
    s_mer_time = np.zeros((TIMESTEPS, n_angles, 3))
    b_mer_time = np.zeros((TIMESTEPS, n_angles, 3))
    
    for t_idx, vtu_file in enumerate(vtus):
        m = meshio.read(str(vtu_file))
        b_n = m.point_data["magnetic flux density"]
        e_n = m.point_data["electric field"]
        
        ib = LinearNDInterpolator(pts_mesh, b_n, fill_value=0.0)
        ie = LinearNDInterpolator(pts_mesh, e_n, fill_value=0.0)
        
        # Equatoriale
        b_eq = ib(pts_eq)
        e_eq = ie(pts_eq)
        s_eq = np.cross(e_eq, b_eq) / MU0
        s_eq_time[t_idx] = s_eq
        
        # Meridiano
        b_m = ib(pts_mer)
        e_m = ie(pts_mer)
        s_m = np.cross(e_m, b_m) / MU0
        s_mer_time[t_idx] = s_m
        b_mer_time[t_idx] = b_m
        
    s_eq_mean = np.mean(s_eq_time, axis=0)
    s_mer_mean = np.mean(s_mer_time, axis=0)
    b_mer_mean = np.mean(b_mer_time, axis=0)
    
    # Componente radiale equatoriale
    r_hat_eq = np.column_stack([np.cos(angles), np.sin(angles), np.zeros_like(angles)])
    s_rad_eq = np.sum(s_eq_mean * r_hat_eq, axis=1)
    
    # Componente radiale meridiana (r_hat = [sin(alpha), 0, cos(alpha)])
    r_hat_mer = np.column_stack([np.sin(angles), np.zeros_like(angles), np.cos(angles)])
    s_rad_mer = np.sum(s_mer_mean * r_hat_mer, axis=1)
    b_rad_mer = np.sum(b_mer_mean * r_hat_mer, axis=1)
    
    return {
        "angles": angles,
        "s_rad_eq": s_rad_eq,
        "s_rad_mer": s_rad_mer,
        "b_rad_mer": b_rad_mer
    }


def plot_fig03_spherical_3d(spheres_data, out_path):
    """
    fig_03_mappatura_sferica_3d_campo_B_E.png:
    Visualizzazione 3D delle sfere con falsi colori di |B|, |E| e S_rad e frecce orientate a 360°.
    """
    fig = plt.figure(figsize=(19, 6), dpi=300)
    
    # Sfera da visualizzare: Mid-Field (R = 10 cm)
    sp = spheres_data[1]  # Mid-field
    arr = sp["_arrays"]
    coords = arr["coords"]
    x, y, z = coords[:, 0] * 100, coords[:, 1] * 100, coords[:, 2] * 100  # in cm
    
    # Sottocampionamento per frecce 3D (quiver)
    stride = 18
    sub_x, sub_y, sub_z = x[::stride], y[::stride], z[::stride]
    sub_bx = arr["b_mean"][::stride, 0]
    sub_by = arr["b_mean"][::stride, 1]
    sub_bz = arr["b_mean"][::stride, 2]
    b_len = np.sqrt(sub_bx**2 + sub_by**2 + sub_bz**2) + 1e-12
    sub_bx, sub_by, sub_bz = sub_bx / b_len, sub_by / b_len, sub_bz / b_len
    
    sub_ex = arr["e_mean"][::stride, 0]
    sub_ey = arr["e_mean"][::stride, 1]
    sub_ez = arr["e_mean"][::stride, 2]
    e_len = np.sqrt(sub_ex**2 + sub_ey**2 + sub_ez**2) + 1e-12
    sub_ex, sub_ey, sub_ez = sub_ex / e_len, sub_ey / e_len, sub_ez / e_len
    
    sub_sx = arr["s_mean"][::stride, 0]
    sub_sy = arr["s_mean"][::stride, 1]
    sub_sz = arr["s_mean"][::stride, 2]
    s_len = np.sqrt(sub_sx**2 + sub_sy**2 + sub_sz**2) + 1e-12
    sub_sx, sub_sy, sub_sz = sub_sx / s_len, sub_sy / s_len, sub_sz / s_len

    # --- Pannello (a): Campo Magnetico |B| & Frecce 3D ---
    ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    b_vals_uT = arr["b_mag"] * 1e6
    p1 = ax1.scatter(x, y, z, c=b_vals_uT, cmap='viridis', s=16, alpha=0.85, edgecolors='none')
    ax1.quiver(sub_x, sub_y, sub_z, sub_bx, sub_by, sub_bz, length=1.8, normalize=True, color='red', lw=0.9, alpha=0.8)
    cb1 = fig.colorbar(p1, ax=ax1, shrink=0.55, aspect=14, pad=0.08)
    cb1.set_label(r'$|\vec{B}_{\mathrm{mean}}|\;[\mu\mathrm{T}]$', fontsize=11, fontweight='bold')
    ax1.set_title(r'(a) Densità Flusso Magnetico $\vec{B}$ (R=10 cm)' + '\nFrecce rosse: Direzione orientata $\hat{B}$', fontsize=12, fontweight='bold')
    ax1.set_xlabel('X [cm]')
    ax1.set_ylabel('Y [cm]')
    ax1.set_zlabel('Z [cm]')
    ax1.set_box_aspect([1, 1, 1])

    # --- Pannello (b): Campo Elettrico Indotto |E| & Frecce 3D ---
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    e_vals_mV = arr["e_mag"] * 1e3
    p2 = ax2.scatter(x, y, z, c=e_vals_mV, cmap='plasma', s=16, alpha=0.85, edgecolors='none')
    ax2.quiver(sub_x, sub_y, sub_z, sub_ex, sub_ey, sub_ez, length=1.8, normalize=True, color='cyan', lw=0.9, alpha=0.8)
    cb2 = fig.colorbar(p2, ax=ax2, shrink=0.55, aspect=14, pad=0.08)
    cb2.set_label(r'$|\vec{E}_{\mathrm{mean}}|\;[\mathrm{mV/m}]$', fontsize=11, fontweight='bold')
    ax2.set_title(r'(b) Campo Elettrico Indotto $\vec{E}$ (R=10 cm)' + '\nFrecce ciano: Vettore azimutale $\hat{E}$', fontsize=12, fontweight='bold')
    ax2.set_xlabel('X [cm]')
    ax2.set_ylabel('Y [cm]')
    ax2.set_zlabel('Z [cm]')
    ax2.set_box_aspect([1, 1, 1])

    # --- Pannello (c): Flusso Attivo di Poynting S_rad & Frecce 3D ---
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    s_rad_uW = arr["s_rad"] * 1e6
    p3 = ax3.scatter(x, y, z, c=s_rad_uW, cmap='inferno', s=16, alpha=0.85, edgecolors='none')
    ax3.quiver(sub_x, sub_y, sub_z, sub_sx, sub_sy, sub_sz, length=1.8, normalize=True, color='lime', lw=0.9, alpha=0.8)
    cb3 = fig.colorbar(p3, ax=ax3, shrink=0.55, aspect=14, pad=0.08)
    cb3.set_label(r'$S_{\mathrm{rad}}\;[\mu\mathrm{W/m}^2]$', fontsize=11, fontweight='bold')
    ax3.set_title(r'(c) Flusso Attivo di Poynting $\vec{S}$ (R=10 cm)' + '\nFrecce verdi: Flusso energetico proiettato', fontsize=12, fontweight='bold')
    ax3.set_xlabel('X [cm]')
    ax3.set_ylabel('Y [cm]')
    ax3.set_zlabel('Z [cm]')
    ax3.set_box_aspect([1, 1, 1])

    plt.suptitle("Open Chiral Flux Shaper (Mantello Chiuso a Barattolo) — Mappatura Sferica 3D Campi Elettrodinamici",
                 fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(str(out_path), dpi=300)
    plt.close()
    print(f"  [FIGURE] Salvata figura 3D: {out_path}")


def plot_fig04_polar_radiation(spheres_data, planar_dict, out_path):
    """
    fig_04_diagramma_radiazione_poynting_360.png:
    Diagramma polare/radiale a 360° dell'emissione di Poynting sul piano equatoriale (Z=0)
    e sul piano meridiano (Z-R) per dimostrare l'effetto di confinamento dei coperchi.
    """
    fig = plt.figure(figsize=(18, 5.5), dpi=300)
    angles = planar_dict["angles"]
    s_eq_uW = planar_dict["s_rad_eq"] * 1e6
    s_mer_uW = planar_dict["s_rad_mer"] * 1e6
    b_mer_uT = planar_dict["b_rad_mer"] * 1e6
    
    # 1. Pannello Polare Equatoriale (Z = 0)
    ax1 = fig.add_subplot(1, 3, 1, projection='polar')
    # Assicurati valori positivi o traslati per diagramma polare radiale
    s_eq_plot = np.maximum(s_eq_uW, 0.0)
    ax1.plot(angles, s_eq_plot, color='tab:red', lw=2.2, label=r'$S_{\mathrm{rad}}(\phi)$ equatoriale')
    ax1.fill(angles, s_eq_plot, color='tab:red', alpha=0.25)
    ax1.set_theta_zero_location("E")
    ax1.set_title(r'(a) Irradiazione Equatoriale ($Z=0$, $R=10$ cm)' + '\n' + r'$S_{\mathrm{rad}}(\phi)$ [$\mu\mathrm{W/m}^2$] a 360°',
                  fontsize=11, fontweight='bold', pad=15)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='lower right', fontsize=9)

    # 2. Pannello Polare Meridiano (X-Z) - EFFETTO CONFINAMENTO DEI COPERCHI
    ax2 = fig.add_subplot(1, 3, 2, projection='polar')
    # alpha=0: +Z (Top Lid), alpha=pi/2: +X (Equatore), alpha=pi: -Z (Bottom Lid), alpha=3pi/2: -X (Equatore)
    s_mer_plot = np.maximum(s_mer_uW, 0.0)
    ax2.plot(angles, s_mer_plot, color='tab:blue', lw=2.2, label=r'$S_{\mathrm{rad}}(\theta_{\mathrm{mer}})$')
    ax2.fill(angles, s_mer_plot, color='tab:blue', alpha=0.25)
    ax2.set_theta_zero_location("N")  # Polo Nord (+Z) in alto
    ax2.set_theta_direction(-1)       # Senso orario: 0=Top, 90=Equatore (+X), 180=Bottom (-Z), 270=Equatore (-X)
    
    # Annotazioni sui coperchi
    ax2.annotate('Coperchio Superiore\n(Flusso Confinato ~ 0)', xy=(0, 0.05 * np.max(s_mer_plot)),
                 xytext=(np.radians(45), 0.85 * np.max(s_mer_plot)),
                 arrowprops=dict(facecolor='darkred', shrink=0.08, width=1.2, headwidth=6),
                 fontsize=8.5, fontweight='bold', color='darkred')
    ax2.annotate('Lobi Equatoriali\n(Espulsione Mantello)', xy=(np.pi/2, 0.95 * np.max(s_mer_plot)),
                 xytext=(np.radians(70), 1.25 * np.max(s_mer_plot)),
                 arrowprops=dict(facecolor='darkblue', shrink=0.08, width=1.2, headwidth=6),
                 fontsize=8.5, fontweight='bold', color='darkblue')
                 
    ax2.set_title(r'(b) Confinamento Meridiano ($X-Z$, $R=10$ cm)' + '\n' + r'Schermatura assiale dei coperchi ($Z=\pm H/2$)',
                  fontsize=11, fontweight='bold', pad=15)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='lower left', fontsize=9)

    # 3. Pannello (c): Decadimento Radiale di Potenza e Flusso sulle 3 Sfere
    ax3 = fig.add_subplot(1, 3, 3)
    radii_cm = [sp["radius_cm"] for sp in spheres_data]
    p_rad_uW = [sp["poynting"]["p_rad_total_uW"] for sp in spheres_data]
    b_mean_uT = [sp["B_field"]["mean_mag_uT"] for sp in spheres_data]
    
    color1 = 'tab:orange'
    ax3.set_xlabel('Raggio Sfera di Valutazione R [cm]', fontsize=11, fontweight='bold')
    ax3.set_ylabel(r'Potenza Irradiata Netta $P_{\mathrm{rad}}$ [$\mu\mathrm{W}$] (log)', color=color1, fontsize=11, fontweight='bold')
    ax3.set_yscale('log')
    line1 = ax3.plot(radii_cm, p_rad_uW, marker='o', color=color1, lw=2.2, label=r'Potenza Irradiata $P_{\mathrm{rad}}$')
    ax3.tick_params(axis='y', labelcolor=color1)
    ax3.grid(True, which='both', linestyle=':', alpha=0.6)
    
    for r_c, p_val in zip(radii_cm, p_rad_uW):
        if p_val >= 1e6:
            lbl = f"{p_val/1e6:.2f} W"
        elif p_val >= 1e3:
            lbl = f"{p_val/1e3:.2f} mW"
        else:
            lbl = f"{p_val:.1f} $\mu$W"
        ax3.annotate(lbl, (r_c, p_val), textcoords="offset points", xytext=(0, 10), ha='center',
                     fontweight='bold', color=color1, fontsize=9.5)

    ax3_twin = ax3.twinx()
    color2 = 'tab:purple'
    ax3_twin.set_ylabel(r'Induzione Magnetica Media $|\vec{B}|$ [$\mu\mathrm{T}$] (log)', color=color2, fontsize=11, fontweight='bold')
    ax3_twin.set_yscale('log')
    line2 = ax3_twin.plot(radii_cm, b_mean_uT, marker='s', linestyle='--', color=color2, lw=2.0, label=r'Induzione Media $|\vec{B}|$')
    ax3_twin.tick_params(axis='y', labelcolor=color2)
    
    for r_c, b_val in zip(radii_cm, b_mean_uT):
        if b_val >= 1000:
            lbl = f"{b_val/1e3:.2f} mT"
        else:
            lbl = f"{b_val:.1f} $\mu$T"
        ax3_twin.annotate(lbl, (r_c, b_val), textcoords="offset points", xytext=(0, -18), ha='center',
                          fontweight='bold', color=color2, fontsize=9.5)

    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='lower left', fontsize=9.5)
    ax3.set_title(r'(c) Profilo di Confinamento e Decadimento Radiale' + '\n' + r'Near-Field $\to$ Mid-Field $\to$ Far-Field', fontsize=11, fontweight='bold')

    plt.suptitle("Open Chiral Flux Shaper — Diagramma Polare di Irradiazione a 360° ed Effetto Confinamento Coperchi",
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(str(out_path), dpi=300)
    plt.close()
    print(f"  [FIGURE] Salvata figura diagramma di radiazione: {out_path}")


def main():
    print("=" * 80)
    print("MAPPATURA SFERICA COMPLETA A 360° DEI CAMPI B, E E POYNTING (MANTELLO CHIUSO)")
    print("=" * 80)
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    work_sub = WORK_DIR / "run_plus30deg"
    
    # 1. Assicura la presenza dei file VTU della soluzione FEM
    vtus = ensure_simulation_data(work_sub)
    
    # 2. Carica geometria dei nodi della mesh
    m0 = meshio.read(str(vtus[0]))
    pts_mesh = m0.points
    print(f"\n[1/4] Mesh caricata: {len(pts_mesh)} nodi, {len(m0.cells_dict.get('tetra', []))} tetraedri.")
    
    # 3. Campionamento sferico su 3 sfere di Fibonacci
    print(f"\n[2/4] Generazione reticolo sferico di Fibonacci ({N_FIBONACCI} punti per sfera)...")
    spheres_data = []
    for r_val, r_label in zip(RADII, RADII_LABELS):
        print(f"  -> Campionamento {r_label}...")
        sp = generate_fibonacci_sphere(r_val, n_points=N_FIBONACCI)
        sp_res = sample_fields_on_sphere(sp, vtus, pts_mesh)
        spheres_data.append(sp_res)
        
        # Stampa a video dei risultati salienti
        g = sp_res["gauss"]
        p = sp_res["poynting"]
        b = sp_res["B_field"]
        e = sp_res["E_field"]
        pass_str = "PASS" if g["certified_pass"] else "FAIL"
        print(f"     * Teorema di Gauss (|Phi_net|/Phi_abs): {g['relative_residual_percent']:.3f}% [{pass_str}] (Net: {g['phi_net_Tm2']:+.3e} T*m^2)")
        print(f"     * Potenza Irradiata Poynting P_rad: {p['p_rad_total_uW']:.3f} uW (Picco: {p['p_rad_peak_uW']:.3f} uW)")
        print(f"     * Campo Magnetico Medio |B|: {b['mean_mag_uT']:.2f} uT (B_rad: {b['mean_Brad_uT']:+.2f} uT, B_theta: {b['mean_Btheta_uT']:+.2f} uT)")
        print(f"     * Campo Elettrico Medio |E|: {e['mean_mag_mV_m']:.2f} mV/m (E_phi: {e['mean_Ephi_mV_m']:+.2f} mV/m)")

    # 4. Campionamento piani ortogonali per diagrammi polari
    print("\n[3/4] Campionamento piani equatoriale (Z=0) e meridiano (X-Z) a 360°...")
    planar_dict = sample_planar_circles(vtus, pts_mesh, r_val=0.10, n_angles=360)

    # 5. Generazione Grafici 300 DPI
    print("\n[4/4] Generazione delle figure ad alta risoluzione (300 DPI)...")
    fig03_path = FIGURES_DIR / "fig_03_mappatura_sferica_3d_campo_B_E.png"
    fig04_path = FIGURES_DIR / "fig_04_diagramma_radiazione_poynting_360.png"
    plot_fig03_spherical_3d(spheres_data, fig03_path)
    plot_fig04_polar_radiation(spheres_data, planar_dict, fig04_path)

    # 6. Esportazione JSON consolidato
    json_path = DATA_DIR / "mappatura_sfere_campi_EB.json"
    export_dict = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "simulation_variant": "rotore_centrato_mantello_chiuso",
        "description": "Mappatura sferica 360 gradi dei campi B, E e vettore di Poynting su sfere di Fibonacci con verifica di Gauss",
        "frequency_Hz": 100,
        "rpm": 1200,
        "timesteps": TIMESTEPS,
        "dt_s": DT,
        "n_fibonacci_points": N_FIBONACCI,
        "spheres": []
    }
    for sp in spheres_data:
        sp_copy = {k: v for k, v in sp.items() if k != "_arrays"}
        export_dict["spheres"].append(sp_copy)
        
    json_path.write_text(json.dumps(export_dict, indent=2), encoding="utf-8")
    print(f"  [DATA] Dataset consolidato salvato in: {json_path}")
    
    # 7. Pulizia cartelle pesanti di calcolo intermedio per git
    if "--keep-vtus" not in sys.argv:
        print("\n[CLEANUP] Rimozione cartella temporanea VTU per mantenere leggero il repository git...")
        shutil.rmtree(str(WORK_DIR), ignore_errors=True)
        print("  [CLEANUP] Rimossi file pesanti intermedi.")
        
    print("\n" + "=" * 80)
    print("MAPPATURA SFERICA COMPLETATA CON SUCCESSO!")
    print("=" * 80)


if __name__ == "__main__":
    main()
