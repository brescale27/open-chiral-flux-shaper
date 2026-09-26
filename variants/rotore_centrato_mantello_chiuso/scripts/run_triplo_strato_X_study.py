#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico Comparativo Avanzato a 4 Vie:
1. Mantello Chiuso in Alluminio Anisotropo (+30°, mu_r = 1.0, Core Ferro mu_r = 1000)
2. Gabbia Ferromagnetica in Ferro Pieno Isotropo (mu_r = 1000, sigma = 1e7, Core Ferro mu_r = 1000)
3. Rete Stirata Ferromagnetica Monostrato (+30°, mu_r = 1000, Core Ferro mu_r = 1000)
4. Metasuperficie a Triplo Strato Incrociato a "X" (+30°/0°/-30°, mu_r = 1000) con ROTORE AMAGNETICO (PEEK, mu_r = 1.0)

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

# Percorsi di lavoro
SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_DIR = SCRIPT_DIR.parent
CONFIG_DIR = VARIANT_DIR / "config"
DATA_DIR = VARIANT_DIR / "data"
FIGURES_DIR = VARIANT_DIR / "figures"
WORK_DIR = VARIANT_DIR / "work_dirs" / "run_triplo_strato_X_rotore_amagnetico"
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
TIMESTEPS = 64
DT = 0.00025  # 0.25 ms
F_HZ = 100.0
W_RAD = 2.0 * np.pi * F_HZ

N_FIBONACCI = 1200
RADII = [0.065, 0.100, 0.150]
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]


def run_simulation():
    """Configura ed esegue la simulazione Elmer FEM per la variante a triplo strato X con rotore amagnetico."""
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    res_dir = WORK_DIR / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} file VTU in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]
        
    print(f"  [SOLVER] Esecuzione ElmerSolver per triplo strato X + rotore amagnetico ({TIMESTEPS} timestep da {DT*1000:.2f} ms)...")
    sif_src = CONFIG_DIR / "case_triplo_strato_X_rotore_amagnetico.sif"
    sif_dest = WORK_DIR / "case.sif"
    
    mesh_db_rel = os.path.relpath(str(MESH_DIR), str(WORK_DIR)).replace('\\', '/')
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
    (WORK_DIR / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")
    
    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    proc = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(WORK_DIR), env=env, capture_output=True, text=True)
    duration = time.time() - t0
    
    if proc.returncode != 0:
        print(f"  [ERRORE] ElmerSolver fallito (codice {proc.returncode}):")
        print(proc.stderr[-1000:])
        sys.exit(1)
        
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    print(f"  [SOLVER] Simulazione completata con successo in {duration:.1f}s ({len(vtus)} VTU generati).")
    return vtus[:TIMESTEPS]


def extract_forces_and_losses(vtus, mesh_data):
    """Calcola le perdite Joule P_J(t) e la forza di Lorentz F_z(t) sui settori fisici e sui 3 sottostrati."""
    cells = mesh_data["cells"]
    elem_vols = mesh_data["elem_vols"]
    mask_lat = mesh_data["mask_lat"]
    mask_l1 = mesh_data["mask_l1"]
    mask_l2 = mesh_data["mask_l2"]
    mask_l3 = mesh_data["mask_l3"]
    mask_top = mesh_data["mask_top"]
    mask_bot = mesh_data["mask_bot"]
    mask_can = mask_lat | mask_top | mask_bot
    mask_core = mesh_data["mask_core"]
    mask_total = mask_can | mask_core
    
    fz_series = {
        "layer1_plus30": [], "layer2_ortho": [], "layer3_minus30": [],
        "lateral": [], "top_cap": [], "bottom_cap": [], "can_total": [],
        "amagnetic_core": [], "total": []
    }
    pj_series = {
        "layer1_plus30": [], "layer2_ortho": [], "layer3_minus30": [],
        "lateral": [], "top_cap": [], "bottom_cap": [], "can_total": [],
        "amagnetic_core": [], "total": []
    }
    
    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        jxb_z = m.point_data['jxb'][:, 2]
        pj = m.point_data['joule heating'].ravel()
        
        jxb_z_elem = np.mean(jxb_z[cells], axis=1)
        pj_elem = np.mean(pj[cells], axis=1)
        
        for k, mask in [
            ("layer1_plus30", mask_l1), ("layer2_ortho", mask_l2), ("layer3_minus30", mask_l3),
            ("lateral", mask_lat), ("top_cap", mask_top), ("bottom_cap", mask_bot),
            ("can_total", mask_can), ("amagnetic_core", mask_core), ("total", mask_total)
        ]:
            if np.sum(mask) > 0:
                fz_series[k].append(float(np.sum(elem_vols[mask] * jxb_z_elem[mask])))
                pj_series[k].append(float(np.sum(elem_vols[mask] * pj_elem[mask])))
            else:
                fz_series[k].append(0.0)
                pj_series[k].append(0.0)
                
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


def generate_fibonacci_sphere(r, n_points=1200):
    """Genera reticolo quasi-uniforme di Fibonacci su sfera di raggio r."""
    i = np.arange(0, n_points, dtype=float) + 0.5
    phi_polar = np.arccos(1.0 - 2.0 * i / n_points)
    golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0
    theta_azimuth = (2.0 * np.pi * i / golden_ratio) % (2.0 * np.pi)
    
    x = r * np.sin(phi_polar) * np.cos(theta_azimuth)
    y = r * np.sin(phi_polar) * np.sin(theta_azimuth)
    z = r * np.cos(phi_polar)
    coords = np.column_stack([x, y, z])
    r_hat = coords / r
    
    sin_p, cos_p = np.sin(phi_polar), np.cos(phi_polar)
    sin_t, cos_t = np.sin(theta_azimuth), np.cos(theta_azimuth)
    theta_hat = np.column_stack([cos_p * cos_t, cos_p * sin_t, -sin_p])
    phi_hat = np.column_stack([-sin_t, cos_t, np.zeros_like(sin_t)])
    da = (4.0 * np.pi * r**2) / n_points
    
    return {
        "r": r, "n_points": n_points, "coords": coords,
        "x": x, "y": y, "z": z,
        "theta_polar": phi_polar, "phi_azimuth": theta_azimuth,
        "r_hat": r_hat, "theta_hat": theta_hat, "phi_hat": phi_hat,
        "da": da, "total_area": 4.0 * np.pi * r**2
    }


def sample_sphere_data(sphere, vtus, delaunay_tri):
    """Campiona B, E e S su una sfera di Fibonacci lungo i 64 timestep."""
    coords = sphere["coords"]
    r_hat = sphere["r_hat"]
    theta_hat = sphere["theta_hat"]
    phi_hat = sphere["phi_hat"]
    da = sphere["da"]
    n_pts = sphere["n_points"]
    
    b_time = np.zeros((len(vtus), n_pts, 3))
    e_time = np.zeros((len(vtus), n_pts, 3))
    s_time = np.zeros((len(vtus), n_pts, 3))
    
    for t_idx, vtu_file in enumerate(vtus):
        m = meshio.read(str(vtu_file))
        b_nodal = m.point_data["magnetic flux density"]
        e_nodal = m.point_data["electric field"]
        
        interp_b = LinearNDInterpolator(delaunay_tri, b_nodal, fill_value=0.0)
        interp_e = LinearNDInterpolator(delaunay_tri, e_nodal, fill_value=0.0)
        
        b_eval = interp_b(coords)
        e_eval = interp_e(coords)
        s_eval = np.cross(e_eval, b_eval) / MU0
        
        b_time[t_idx] = b_eval
        e_time[t_idx] = e_eval
        s_time[t_idx] = s_eval

    b_mean = np.mean(b_time, axis=0)
    e_mean = np.mean(e_time, axis=0)
    s_mean = np.mean(s_time, axis=0)
    
    b_rad = np.sum(b_mean * r_hat, axis=1)
    b_theta = np.sum(b_mean * theta_hat, axis=1)
    b_phi = np.sum(b_mean * phi_hat, axis=1)
    b_mag = np.linalg.norm(b_mean, axis=1)
    
    e_rad = np.sum(e_mean * r_hat, axis=1)
    e_theta = np.sum(e_mean * theta_hat, axis=1)
    e_phi = np.sum(e_mean * phi_hat, axis=1)
    e_mag = np.linalg.norm(e_mean, axis=1)
    
    s_rad = np.sum(s_mean * r_hat, axis=1)
    s_mag = np.linalg.norm(s_mean, axis=1)
    
    phi_net_mean = float(np.sum(b_rad * da))
    phi_abs_mean = float(np.sum(np.abs(b_rad) * da))
    rel_res = abs(phi_net_mean) / (phi_abs_mean + 1e-15)
    p_rad_total_W = float(np.sum(s_rad * da))
    
    return {
        "sphere": sphere,
        "b_mean": b_mean, "e_mean": e_mean, "s_mean": s_mean,
        "b_mag": b_mag, "b_rad": b_rad, "b_theta": b_theta, "b_phi": b_phi,
        "e_mag": e_mag, "e_rad": e_rad, "e_theta": e_theta, "e_phi": e_phi,
        "s_mag": s_mag, "s_rad": s_rad,
        "gauss_net_Wb": phi_net_mean,
        "gauss_abs_Wb": phi_abs_mean,
        "gauss_residual": rel_res,
        "gauss_pass": bool(rel_res < 0.05),
        "P_rad_W": p_rad_total_W,
        "P_rad_uW": p_rad_total_W * 1e6
    }


def plot_figure_12_four_way_comparison(al_data, fe_solid_data, fe_mesh_data, x_data,
                                       al_spheres, fe_solid_spheres, fe_mesh_spheres, x_spheres):
    """Genera fig_12_confronto_triplo_strato_X_vs_precedenti.png (300 DPI)."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 11), dpi=300)
    time_ms = np.arange(1, TIMESTEPS + 1) * DT * 1000.0
    
    # 1. Dissipazione Termica Joule P_J(t)
    pj_al = al_data["total"]["Pj_series_mW"]
    pj_solid = fe_solid_data["total"]["Pj_series_mW"]
    pj_rf = fe_mesh_data["total"]["Pj_series_mW"]
    pj_x = x_data["total"]["Pj_series_mW"]
    
    m_al = al_data['total']['mean_Pj_mW']
    m_solid = fe_solid_data['total']['mean_Pj_mW']
    m_rf = fe_mesh_data['total']['mean_Pj_mW']
    m_x = x_data['total']['mean_Pj_mW']
    
    ax1.plot(time_ms, pj_al, color='#1f77b4', lw=1.8, label=f"1. Alluminio Aniso 30°: {m_al*1000:.1f} µW")
    ax1.plot(time_ms, pj_solid, color='#d62728', lw=1.8, label=f"2. Ferro Pieno Isotropo: {m_solid*1000:.1f} µW")
    ax1.plot(time_ms, pj_rf, color='#2ca02c', lw=1.8, ls='--', label=f"3. Rete Ferro Monostrato: {m_rf*1000:.1f} µW")
    ax1.plot(time_ms, pj_x, color='#9467bd', lw=2.2, label=f"4. Triplo Strato X + Rotore PEEK: {m_x*1000:.1f} µW")
    
    ax1.set_xlabel("Tempo t [ms]", fontsize=10, fontweight='bold')
    ax1.set_ylabel(r"Dissipazione Joule Totale $P_J(t)$ [mW]", fontsize=10, fontweight='bold')
    ax1.set_title("1. Dissipazione Joule: Abbattimento Parassita con Struttura a 'X'", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=8.2)
    
    # 2. Forza Assiale Lorentz F_z(t)
    fz_al = al_data["total"]["Fz_series_uN"]
    fz_solid = fe_solid_data["total"]["Fz_series_uN"]
    fz_rf = fe_mesh_data["total"]["Fz_series_uN"]
    fz_x = x_data["total"]["Fz_series_uN"]
    
    ax2.plot(time_ms, fz_al, color='#1f77b4', lw=1.8, label=rf"1. Alluminio: $\langle F_z \rangle$ = {al_data['total']['mean_Fz_uN']:+.3f} µN")
    ax2.plot(time_ms, fz_solid, color='#d62728', lw=1.8, label=rf"2. Ferro Pieno: $\langle F_z \rangle$ = {fe_solid_data['total']['mean_Fz_uN']:+.3f} µN")
    ax2.plot(time_ms, fz_rf, color='#2ca02c', lw=1.8, ls='--', label=rf"3. Rete Ferro: $\langle F_z \rangle$ = {fe_mesh_data['total']['mean_Fz_uN']:+.3f} µN")
    ax2.plot(time_ms, fz_x, color='#9467bd', lw=2.2, label=rf"4. Triplo Strato X: $\langle F_z \rangle$ = {x_data['total']['mean_Fz_uN']:+.3f} µN")
    
    ax2.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax2.set_xlabel("Tempo t [ms]", fontsize=10, fontweight='bold')
    ax2.set_ylabel(r"Forza Assiale Lorentz $F_z(t)$ [µN]", fontsize=10, fontweight='bold')
    ax2.set_title("2. Spinta Assiale Ponderomotrice Lorentz F_z(t)", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=8.2)
    
    # 3. Attenuazione Campo Magnetico |B| vs Raggio
    radii_cm = [s["radius_cm"] for s in al_spheres]
    b_al = [s["B_mean_mag_uT"] for s in al_spheres]
    b_solid = [s["B_mean_mag_uT"] for s in fe_solid_spheres]
    b_rf = [s["B_mean_mag_uT"] for s in fe_mesh_spheres]
    b_x = [s["B_mean_mag_uT"] for s in x_spheres]
    
    ax3.plot(radii_cm, b_al, marker='o', lw=1.8, color='#1f77b4', label=r"1. Alluminio ($\mu_r = 1$, Metasuperficie)")
    ax3.plot(radii_cm, b_solid, marker='s', lw=1.8, color='#d62728', label=r"2. Ferro Pieno ($\mu_r = 1000$, Shunt Solido)")
    ax3.plot(radii_cm, b_rf, marker='^', lw=1.8, ls='--', color='#2ca02c', label=r"3. Rete Ferro Monostrato ($\mu_r = 1000$)")
    ax3.plot(radii_cm, b_x, marker='D', lw=2.2, color='#9467bd', label=r"4. Triplo Strato X + Rotore PEEK")
    
    ax3.set_yscale('log')
    ax3.set_xlabel("Raggio Sfera R [cm]", fontsize=10, fontweight='bold')
    ax3.set_ylabel(r"Induzione Magnetica Media $|\vec{B}|$ [µT]", fontsize=10, fontweight='bold')
    ax3.set_title("3. Profilo Radiale Campo Magnetico |B| (Scala Log)", fontsize=11, fontweight='bold')
    ax3.grid(True, which="both", linestyle=":", alpha=0.6)
    ax3.legend(loc="upper right", fontsize=8.2)
    
    # 4. Potenza Irradiata Poynting P_rad vs Raggio
    p_al = [s["P_rad_uW"] for s in al_spheres]
    p_solid = [s["P_rad_uW"] for s in fe_solid_spheres]
    p_rf = [s["P_rad_uW"] for s in fe_mesh_spheres]
    p_x = [s["P_rad_uW"] for s in x_spheres]
    
    x = np.arange(len(radii_cm))
    width = 0.20
    ax4.bar(x - 1.5*width, p_al, width, label='1. Alluminio', color='#1f77b4', alpha=0.85)
    ax4.bar(x - 0.5*width, p_solid, width, label='2. Ferro Pieno', color='#d62728', alpha=0.85)
    ax4.bar(x + 0.5*width, p_rf, width, label='3. Rete Ferro', color='#2ca02c', alpha=0.85)
    ax4.bar(x + 1.5*width, p_x, width, label='4. Triplo Strato X', color='#9467bd', alpha=0.85)
    
    ax4.set_xticks(x)
    ax4.set_xticklabels([f"R={r:.1f} cm" for r in radii_cm], fontsize=9, fontweight='bold')
    ax4.set_ylabel(r"Potenza Radiativa Poynting $P_{\text{rad}}$ [µW]", fontsize=10, fontweight='bold')
    ax4.set_title("4. Potenza Attiva Irradiata Esternamente P_rad", fontsize=11, fontweight='bold')
    ax4.grid(True, linestyle=":", alpha=0.6)
    ax4.legend(loc="upper right", fontsize=8.2)
    
    plt.tight_layout()
    out_png = FIGURES_DIR / "fig_12_confronto_triplo_strato_X_vs_precedenti.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"  [FIGURA] Salvata: {out_png}")


def plot_figure_13_triplo_strato_3d(sphere_data):
    """Genera fig_13_mappatura_3d_triplo_strato_X.png (300 DPI)."""
    fig = plt.figure(figsize=(13, 10), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    
    # Disegna guscio barattolo
    z_cyl = np.linspace(-0.05, 0.05, 30)
    th_cyl = np.linspace(0, 2*np.pi, 50)
    Z_mesh, TH_mesh = np.meshgrid(z_cyl, th_cyl)
    X_cyl = 0.050 * np.cos(TH_mesh)
    Y_cyl = 0.050 * np.sin(TH_mesh)
    ax.plot_surface(X_cyl, Y_cyl, Z_mesh, alpha=0.15, color="#9467bd", edgecolor="none")
    
    r_disc = np.linspace(0, 0.050, 20)
    R_disc, TH_disc = np.meshgrid(r_disc, th_cyl)
    X_disc = R_disc * np.cos(TH_disc)
    Y_disc = R_disc * np.sin(TH_disc)
    ax.plot_surface(X_disc, Y_disc, np.full_like(X_disc, 0.050), alpha=0.18, color="#9467bd")
    ax.plot_surface(X_disc, Y_disc, np.full_like(X_disc, -0.050), alpha=0.18, color="#9467bd")
    
    # Campioni vettoriali sulla sfera Near-Field (R = 6.5 cm)
    s_data = sphere_data[0]
    coords = s_data["sphere"]["coords"]
    b_mean = s_data["b_mean"]
    b_mag = s_data["b_mag"]
    b_rad = s_data["b_rad"]
    
    step = 4
    x_sub = coords[::step, 0]
    y_sub = coords[::step, 1]
    z_sub = coords[::step, 2]
    u_sub = b_mean[::step, 0]
    v_sub = b_mean[::step, 1]
    w_sub = b_mean[::step, 2]
    rad_sub = b_rad[::step]
    mag_sub = b_mag[::step]
    
    scale = 0.015 / (np.max(mag_sub) + 1e-12)
    norm = plt.Normalize(vmin=-np.max(np.abs(rad_sub)), vmax=np.max(np.abs(rad_sub)))
    cmap = plt.cm.coolwarm
    colors = cmap(norm(rad_sub))
    
    ax.quiver(x_sub, y_sub, z_sub, u_sub*scale, v_sub*scale, w_sub*scale,
              colors=colors, arrow_length_ratio=0.35, linewidth=1.1, alpha=0.85)
    ax.scatter(x_sub, y_sub, z_sub, c=rad_sub, cmap='coolwarm', norm=norm, s=12, alpha=0.7)
    
    ax.text(0, 0, 0.080, "Triplo Strato Incrociato a 'X' (+30°/0°/-30°)\nRotore Amagnetico (PEEK, mu_r = 1.0)", color="#4a148c", fontsize=9, fontweight='bold', ha='center')
    ax.text(0.080, 0, 0, "Guida Magnetica X-Chiral", color="#4a148c", fontsize=9, fontweight='bold', ha='left')
    
    ax.set_xlabel("X [m]", fontsize=10, fontweight='bold')
    ax.set_ylabel("Y [m]", fontsize=10, fontweight='bold')
    ax.set_zlabel("Z [m]", fontsize=10, fontweight='bold')
    ax.set_title("Mappatura 3D Campo B su Metasuperficie a Triplo Strato Incrociato ('X')\nAnnullamento Bias Unidirezionali e Focalizzazione Radiale a 360°", fontsize=12, fontweight='bold')
    
    lim = 0.09
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    ax.set_zlim([-lim, lim])
    ax.view_init(elev=22, azim=48)
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label(r"Induzione Radiale $B_{\text{rad}}$ [T]", fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    out_png = FIGURES_DIR / "fig_13_mappatura_3d_triplo_strato_X.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"  [FIGURA] Salvata: {out_png}")


def main():
    print("=" * 85)
    print("STUDIO ELETTRODINAMICO: METASUPERFICIE A TRIPLO STRATO INCROCIATO A 'X'")
    print("Mantello mu_r=1000 (+30°/0°/-30°) con Rotore Amagnetico PEEK (mu_r=1.0)")
    print("=" * 85)
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    
    mesh_msh = MESH_DIR / "macchina_mantello_chiuso.msh"
    if not mesh_msh.is_file():
        print("[ERRORE] Mesh non trovata in mesh/macchina_mantello_chiuso.msh")
        sys.exit(1)
        
    m = meshio.read(str(mesh_msh))
    cells = m.cells_dict['tetra']
    pts = m.points
    v0, v1, v2, v3 = pts[cells[:, 0]], pts[cells[:, 1]], pts[cells[:, 2]], pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0
    elem_com = (v0 + v1 + v2 + v3) / 4.0
    elem_r = np.hypot(elem_com[:, 0], elem_com[:, 1])
    elem_z = elem_com[:, 2]
    
    mask_lat = (elem_r >= 0.0465) & (elem_r <= 0.0505) & (np.abs(elem_z) <= 0.050)
    mask_l1 = mask_lat & (elem_r <= 0.0480)
    mask_l2 = mask_lat & (elem_r > 0.0480) & (elem_r <= 0.0490)
    mask_l3 = mask_lat & (elem_r > 0.0490)
    mask_top = (elem_z > 0.050) & (elem_z <= 0.0535) & (elem_r <= 0.0505)
    mask_bot = (elem_z < -0.050) & (elem_z >= -0.0535) & (elem_r <= 0.0505)
    mask_core = (elem_r < 0.0475) & (np.abs(elem_z) <= 0.008)
    
    mesh_data = {
        "cells": cells, "pts": pts, "elem_vols": elem_vols,
        "mask_lat": mask_lat, "mask_l1": mask_l1, "mask_l2": mask_l2, "mask_l3": mask_l3,
        "mask_top": mask_top, "mask_bot": mask_bot, "mask_core": mask_core
    }
    
    # 1. Esecuzione simulazione
    vtus = run_simulation()
    
    # 2. Estrazione perdite e forze
    print("  [POST-PROC] Estrazione perdite Joule e forza assiale nei settori e sottostrati...")
    x_electro = extract_forces_and_losses(vtus, mesh_data)
    
    # 3. Campionamento sferico di Fibonacci
    print("  [POST-PROC] Campionamento sferico 3D a 360° (N=1200 punti per sfera)...")
    m_vtu0 = meshio.read(str(vtus[0]))
    pts_vtu = m_vtu0.points
    print("    - Costruzione triangolazione Delaunay 3D dei nodi VTU...")
    delaunay_tri = Delaunay(pts_vtu)
    
    x_sphere_raw = []
    x_spheres_json = []
    
    for r_m, r_lbl in zip(RADII, RADII_LABELS):
        print(f"    - Campionamento {r_lbl}...")
        sp = generate_fibonacci_sphere(r_m, n_points=N_FIBONACCI)
        s_data = sample_sphere_data(sp, vtus, delaunay_tri)
        x_sphere_raw.append(s_data)
        
        x_spheres_json.append({
            "label": r_lbl,
            "radius_m": r_m,
            "radius_cm": r_m * 100.0,
            "B_mean_mag_uT": float(np.mean(s_data["b_mag"]) * 1e6),
            "E_mean_mag_mVm": float(np.mean(s_data["e_mag"]) * 1000.0),
            "P_rad_uW": s_data["P_rad_uW"],
            "gauss_net_Wb": s_data["gauss_net_Wb"],
            "gauss_abs_Wb": s_data["gauss_abs_Wb"],
            "gauss_residual_pct": float(s_data["gauss_residual"] * 100.0),
            "gauss_pass": s_data["gauss_pass"]
        })
        
    # 4. Caricamento dati precedenti (Alluminio, Ferro Pieno e Rete Monostrato)
    al_json_path = DATA_DIR / "risultati_handover_terzi.json"
    fe_solid_json_path = DATA_DIR / "risultati_mantello_ferromagnetico.json"
    rf_json_path = DATA_DIR / "risultati_rete_ferromagnetica.json"
    
    with open(al_json_path, "r", encoding="utf-8") as f:
        al_json = json.load(f)
    with open(fe_solid_json_path, "r", encoding="utf-8") as f:
        fe_solid_json = json.load(f)
    with open(rf_json_path, "r", encoding="utf-8") as f:
        rf_json = json.load(f)
        
    al_electro = al_json["electrodynamic_results"]
    al_spheres = al_json["spherical_radiation"]
    fe_solid_electro = fe_solid_json["electrodynamic_results_ferro"]
    fe_solid_spheres = fe_solid_json["spherical_radiation_ferro"]
    rf_electro = rf_json["electrodynamic_results_mesh_ferro"]
    rf_spheres = rf_json["spherical_radiation_mesh_ferro"]
    
    # 5. Generazione figure a 300 DPI
    print("  [PLOTS] Generazione grafici comparativi a 4 vie a 300 DPI...")
    plot_figure_12_four_way_comparison(al_electro, fe_solid_electro, rf_electro, x_electro,
                                       al_spheres, fe_solid_spheres, rf_spheres, x_spheres_json)
    plot_figure_13_triplo_strato_3d(x_sphere_raw)
    
    # 6. Salvataggio dataset consolidato
    out_json_path = DATA_DIR / "risultati_triplo_strato_X.json"
    results_consolidated = {
        "metadata": {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "author": "Alessandro Brescacin",
            "license": "CERN-OHL-S-2.0",
            "description": "Confronto a 4 vie: Alluminio Anisotropo vs Ferro Pieno Isotropo vs Rete Ferromagnetica Monostrato vs Triplo Strato Incrociato a 'X' (+30°/0°/-30°) con Rotore Amagnetico PEEK"
        },
        "material_properties": {
            "mantle": {
                "architecture": "Triplo Strato Incrociato a X",
                "relative_permeability": 1000.0,
                "layer1_inner": "+30 deg (sigma_tz = +3.031e6 S/m)",
                "layer2_mid": "transition orthogonal (sigma_tz = 0.0 S/m)",
                "layer3_outer": "-30 deg (sigma_tz = -3.031e6 S/m)"
            },
            "rotor_core": {
                "material": "PEEK / Tecnopolimero Strutturale Amagnetico",
                "relative_permeability": 1.0,
                "electrical_conductivity": 0.0
            }
        },
        "electrodynamic_results_triplo_strato_X": x_electro,
        "spherical_radiation_triplo_strato_X": x_spheres_json,
        "four_way_comparison": {
            "Pj_mW": {
                "aluminum": al_electro["total"]["mean_Pj_mW"],
                "solid_iron": fe_solid_electro["total"]["mean_Pj_mW"],
                "ferro_mesh_monolayer": rf_electro["total"]["mean_Pj_mW"],
                "triplo_strato_X_peek": x_electro["total"]["mean_Pj_mW"]
            },
            "Fz_mean_uN": {
                "aluminum": al_electro["total"]["mean_Fz_uN"],
                "solid_iron": fe_solid_electro["total"]["mean_Fz_uN"],
                "ferro_mesh_monolayer": rf_electro["total"]["mean_Fz_uN"],
                "triplo_strato_X_peek": x_electro["total"]["mean_Fz_uN"]
            },
            "Fz_peak_uN": {
                "aluminum": al_electro["total"]["peak_Fz_uN"],
                "solid_iron": fe_solid_electro["total"]["peak_Fz_uN"],
                "ferro_mesh_monolayer": rf_electro["total"]["peak_Fz_uN"],
                "triplo_strato_X_peek": x_electro["total"]["peak_Fz_uN"]
            },
            "B_far_field_uT": {
                "aluminum": al_spheres[-1]["B_mean_mag_uT"],
                "solid_iron": fe_solid_spheres[-1]["B_mean_mag_uT"],
                "ferro_mesh_monolayer": rf_spheres[-1]["B_mean_mag_uT"],
                "triplo_strato_X_peek": x_spheres_json[-1]["B_mean_mag_uT"]
            }
        }
    }
    
    with open(out_json_path, "w", encoding="utf-8") as f:
        json.dump(results_consolidated, f, indent=2)
    print(f"  [DATA] Dati salvati con successo in: {out_json_path}")
    
    # 7. Report a video
    fw = results_consolidated["four_way_comparison"]
    print("=" * 85)
    print("SINTESI COMPARATIVA A 4 VIE:")
    print(f"  - Spinta Lorentz media <F_z>:")
    print(f"    * 1. Alluminio Anisotropo:    {fw['Fz_mean_uN']['aluminum']:+.4f} µN")
    print(f"    * 2. Ferro Pieno Isotropo:    {fw['Fz_mean_uN']['solid_iron']:+.4f} µN")
    print(f"    * 3. Rete Ferro Monostrato:   {fw['Fz_mean_uN']['ferro_mesh_monolayer']:+.4f} µN")
    print(f"    * 4. Triplo Strato X + PEEK:  {fw['Fz_mean_uN']['triplo_strato_X_peek']:+.4f} µN")
    print(f"  - Spinta nei 3 sottostrati del mantello X:")
    print(f"    * Strato 1 (+30° interno):    {x_electro['layer1_plus30']['mean_Fz_uN']:+.4f} µN")
    print(f"    * Strato 2 (0° intermedio):   {x_electro['layer2_ortho']['mean_Fz_uN']:+.4f} µN")
    print(f"    * Strato 3 (-30° esterno):    {x_electro['layer3_minus30']['mean_Fz_uN']:+.4f} µN")
    print(f"    * Nucleo Amagnetico (PEEK):   {x_electro['amagnetic_core']['mean_Fz_uN']:+.4f} µN")
    print(f"  - Campo Far-Field |B| (15 cm):")
    print(f"    * 1. Alluminio:               {fw['B_far_field_uT']['aluminum']:.3f} µT")
    print(f"    * 2. Ferro Pieno:             {fw['B_far_field_uT']['solid_iron']:.3f} µT")
    print(f"    * 3. Rete Ferro Monostrato:   {fw['B_far_field_uT']['ferro_mesh_monolayer']:.3f} µT")
    print(f"    * 4. Triplo Strato X + PEEK:  {fw['B_far_field_uT']['triplo_strato_X_peek']:.3f} µT")
    print("=" * 85)


if __name__ == "__main__":
    main()
