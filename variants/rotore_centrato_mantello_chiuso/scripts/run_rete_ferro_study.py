#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico Comparativo a 3 Vie:
1. Mantello Chiuso in Alluminio Anisotropo (EN AW-1050A, mu_r = 1, sigma_theta_z = 3.03e6 S/m)
2. Gabbia Ferromagnetica in Ferro Pieno (Armco Iron, mu_r = 1000, sigma = 1e7 S/m isotropo)
3. Rete Stirata Ferromagnetica Ibrida (Armco Iron, mu_r = 1000, sigma_theta_z = 3.03e6 S/m anisotropo)

1. Esecuzione simulazione Elmer FEM transiente per Rete Ferromagnetica (64 timestep, dt = 0.25 ms).
2. Estrazione perdite Joule P_J(t) e forze assiali F_z(t) sui 4 settori fisici.
3. Campionamento sferico di Fibonacci a 360° su sfere concentriche (6.5 cm, 10.0 cm, 15.0 cm).
4. Generazione deliverable comparativi a 3 vie (300 DPI):
   - fig_10_confronto_tre_materiali_mantello.png
   - fig_11_mappatura_3d_rete_ferromagnetica.png
5. Salvataggio dataset in data/risultati_rete_ferromagnetica.json.

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
WORK_DIR = VARIANT_DIR / "work_dirs" / "run_pulsed_third_handover_rete_ferro"
MESH_DIR = VARIANT_DIR / "mesh"

ELMER_SOLVER = r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"
MU0 = 4.0 * np.pi * 1e-7
TIMESTEPS = 64
DT = 0.00025  # 0.25 ms
F_HZ = 100.0
W_RAD = 2.0 * np.pi * F_HZ

N_FIBONACCI = 1200
RADII = [0.065, 0.100, 0.150]
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]


def run_simulation():
    """Configura ed esegue la simulazione Elmer FEM per la rete ferromagnetica."""
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    res_dir = WORK_DIR / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} file VTU in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]
        
    print(f"  [SOLVER] Esecuzione ElmerSolver per rete ferromagnetica ({TIMESTEPS} timestep da {DT*1000:.2f} ms)...")
    sif_src = CONFIG_DIR / "case_pulsed_third_handover_rete_ferro.sif"
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
    print(f"  [SOLVER] Simulazione rete ferromagnetica completata in {duration:.1f}s ({len(vtus)} VTU generati).")
    return vtus[:TIMESTEPS]


def extract_forces_and_losses(vtus, mesh_data):
    """Calcola le perdite Joule P_J(t) e la forza di Lorentz F_z(t) sui 4 settori."""
    cells = mesh_data["cells"]
    elem_vols = mesh_data["elem_vols"]
    mask_lat = mesh_data["mask_lat"]
    mask_top = mesh_data["mask_top"]
    mask_bot = mesh_data["mask_bot"]
    mask_inner = mesh_data["mask_inner"]
    mask_can = mask_lat | mask_top | mask_bot
    mask_total = mask_can | mask_inner
    
    fz_series = {"lateral": [], "top_cap": [], "bottom_cap": [], "can_total": [], "inner_core": [], "total": []}
    pj_series = {"lateral": [], "top_cap": [], "bottom_cap": [], "can_total": [], "inner_core": [], "total": []}
    
    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        jxb_z = m.point_data['jxb'][:, 2]
        pj = m.point_data['joule heating'].ravel()
        
        jxb_z_elem = np.mean(jxb_z[cells], axis=1)
        pj_elem = np.mean(pj[cells], axis=1)
        
        fz_series["lateral"].append(float(np.sum(elem_vols[mask_lat] * jxb_z_elem[mask_lat])))
        pj_series["lateral"].append(float(np.sum(elem_vols[mask_lat] * pj_elem[mask_lat])))
        
        fz_series["top_cap"].append(float(np.sum(elem_vols[mask_top] * jxb_z_elem[mask_top])))
        pj_series["top_cap"].append(float(np.sum(elem_vols[mask_top] * pj_elem[mask_top])))
        
        fz_series["bottom_cap"].append(float(np.sum(elem_vols[mask_bot] * jxb_z_elem[mask_bot])))
        pj_series["bottom_cap"].append(float(np.sum(elem_vols[mask_bot] * pj_elem[mask_bot])))
        
        fz_series["can_total"].append(float(np.sum(elem_vols[mask_can] * jxb_z_elem[mask_can])))
        pj_series["can_total"].append(float(np.sum(elem_vols[mask_can] * pj_elem[mask_can])))
        
        fz_series["inner_core"].append(float(np.sum(elem_vols[mask_inner] * jxb_z_elem[mask_inner])))
        pj_series["inner_core"].append(float(np.sum(elem_vols[mask_inner] * pj_elem[mask_inner])))
        
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


def plot_figure_10_three_materials(al_data, fe_solid_data, fe_mesh_data, al_spheres, fe_solid_spheres, fe_mesh_spheres):
    """Genera fig_10_confronto_tre_materiali_mantello.png (300 DPI)."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
    
    time_ms = np.arange(1, TIMESTEPS + 1) * DT * 1000.0
    
    # 1. Dissipazione Termica Joule P_J(t)
    pj_al = al_data["total"]["Pj_series_mW"]
    pj_solid = fe_solid_data["total"]["Pj_series_mW"]
    pj_mesh = fe_mesh_data["total"]["Pj_series_mW"]
    
    m_al = al_data['total']['mean_Pj_mW']
    m_solid = fe_solid_data['total']['mean_Pj_mW']
    m_mesh = fe_mesh_data['total']['mean_Pj_mW']
    
    ax1.plot(time_ms, pj_al, color='#1f77b4', lw=2.0, label=f"Alluminio Anisotropo (30°): {m_al*1000:.1f} µW")
    ax1.plot(time_ms, pj_solid, color='#d62728', lw=2.0, label=f"Ferro Pieno Isotropo: {m_solid*1000:.1f} µW")
    ax1.plot(time_ms, pj_mesh, color='#2ca02c', lw=2.2, ls='--', label=f"Rete Ferromagnetica (30°): {m_mesh*1000:.1f} µW")
    
    ax1.set_xlabel("Tempo t [ms]", fontsize=10, fontweight='bold')
    ax1.set_ylabel(r"Dissipazione Joule Totale $P_J(t)$ [mW]", fontsize=10, fontweight='bold')
    ax1.set_title("Confronto Perdite Termiche Joule: Effetto della Maglia Chirale", fontsize=11, fontweight='bold')
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=8.5)
    
    # 2. Forza Assiale Lorentz F_z(t)
    fz_al = al_data["total"]["Fz_series_uN"]
    fz_solid = fe_solid_data["total"]["Fz_series_uN"]
    fz_mesh = fe_mesh_data["total"]["Fz_series_uN"]
    
    ax2.plot(time_ms, fz_al, color='#1f77b4', lw=2.0, label=rf"Alluminio: $\langle F_z \rangle$ = {al_data['total']['mean_Fz_uN']:+.3f} µN")
    ax2.plot(time_ms, fz_solid, color='#d62728', lw=2.0, label=rf"Ferro Pieno: $\langle F_z \rangle$ = {fe_solid_data['total']['mean_Fz_uN']:+.3f} µN")
    ax2.plot(time_ms, fz_mesh, color='#2ca02c', lw=2.2, ls='--', label=rf"Rete Ferro: $\langle F_z \rangle$ = {fe_mesh_data['total']['mean_Fz_uN']:+.3f} µN")
    
    ax2.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax2.set_xlabel("Tempo t [ms]", fontsize=10, fontweight='bold')
    ax2.set_ylabel(r"Forza Assiale Lorentz $F_z(t)$ [µN]", fontsize=10, fontweight='bold')
    ax2.set_title("Forza Assiale Ponderomotrice di Lorentz F_z(t)", fontsize=11, fontweight='bold')
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=8.5)
    
    # 3. Attenuazione Campo Magnetico |B| vs Raggio
    radii_cm = [s["radius_cm"] for s in al_spheres]
    b_al = [s["B_mean_mag_uT"] for s in al_spheres]
    b_solid = [s["B_mean_mag_uT"] for s in fe_solid_spheres]
    b_mesh = [s["B_mean_mag_uT"] for s in fe_mesh_spheres]
    
    ax3.plot(radii_cm, b_al, marker='o', lw=2.0, color='#1f77b4', label=r"Alluminio ($\mu_r = 1$, Metasuperficie)")
    ax3.plot(radii_cm, b_solid, marker='s', lw=2.0, color='#d62728', label=r"Ferro Pieno ($\mu_r = 1000$, Shunt Solido)")
    ax3.plot(radii_cm, b_mesh, marker='^', lw=2.2, ls='--', color='#2ca02c', label=r"Rete Ferromagnetica ($\mu_r = 1000$, Ibrida)")
    
    ax3.set_yscale('log')
    ax3.set_xlabel("Raggio Sfera R [cm]", fontsize=10, fontweight='bold')
    ax3.set_ylabel(r"Induzione Magnetica Media $|\vec{B}|$ [µT]", fontsize=10, fontweight='bold')
    ax3.set_title("Profilo Radiale Induzione Magnetica (Scala Log)", fontsize=11, fontweight='bold')
    ax3.grid(True, which="both", linestyle=":", alpha=0.6)
    ax3.legend(loc="upper right", fontsize=8.5)
    
    # 4. Potenza Irradiata Poynting P_rad vs Raggio
    p_al = [s["P_rad_uW"] for s in al_spheres]
    p_solid = [s["P_rad_uW"] for s in fe_solid_spheres]
    p_mesh = [s["P_rad_uW"] for s in fe_mesh_spheres]
    
    x = np.arange(len(radii_cm))
    width = 0.25
    ax4.bar(x - width, p_al, width, label='Alluminio', color='#1f77b4', alpha=0.85)
    ax4.bar(x, p_solid, width, label='Ferro Pieno', color='#d62728', alpha=0.85)
    ax4.bar(x + width, p_mesh, width, label='Rete Ferro', color='#2ca02c', alpha=0.85)
    
    ax4.set_xticks(x)
    ax4.set_xticklabels([f"R={r:.1f} cm" for r in radii_cm], fontsize=9, fontweight='bold')
    ax4.set_ylabel(r"Potenza Radiativa Poynting $P_{\text{rad}}$ [µW]", fontsize=10, fontweight='bold')
    ax4.set_title("Potenza Attiva Irradiata Esternamente", fontsize=11, fontweight='bold')
    ax4.grid(True, linestyle=":", alpha=0.6)
    ax4.legend(loc="upper right", fontsize=8.5)
    
    plt.tight_layout()
    out_png = FIGURES_DIR / "fig_10_confronto_tre_materiali_mantello.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"  [FIGURA] Salvata: {out_png}")


def plot_figure_11_ferro_mesh_3d(sphere_data):
    """Genera fig_11_mappatura_3d_rete_ferromagnetica.png (300 DPI)."""
    fig = plt.figure(figsize=(13, 10), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    
    # Disegna guscio barattolo
    z_cyl = np.linspace(-0.05, 0.05, 30)
    th_cyl = np.linspace(0, 2*np.pi, 50)
    Z_mesh, TH_mesh = np.meshgrid(z_cyl, th_cyl)
    X_cyl = 0.050 * np.cos(TH_mesh)
    Y_cyl = 0.050 * np.sin(TH_mesh)
    ax.plot_surface(X_cyl, Y_cyl, Z_mesh, alpha=0.15, color="#2ca02c", edgecolor="none")
    
    r_disc = np.linspace(0, 0.050, 20)
    R_disc, TH_disc = np.meshgrid(r_disc, th_cyl)
    X_disc = R_disc * np.cos(TH_disc)
    Y_disc = R_disc * np.sin(TH_disc)
    ax.plot_surface(X_disc, Y_disc, np.full_like(X_disc, 0.050), alpha=0.18, color="#2ca02c")
    ax.plot_surface(X_disc, Y_disc, np.full_like(X_disc, -0.050), alpha=0.18, color="#2ca02c")
    
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
    
    ax.text(0, 0, 0.080, "Rete Ferromagnetica (mu_r = 1000, 30°)\nEspulsione Cuspide + Guida Magnetica", color="#1b5e20", fontsize=9, fontweight='bold', ha='center')
    ax.text(0.080, 0, 0, "Guida Anisotropa", color="#1b5e20", fontsize=9, fontweight='bold', ha='left')
    
    ax.set_xlabel("X [m]", fontsize=10, fontweight='bold')
    ax.set_ylabel("Y [m]", fontsize=10, fontweight='bold')
    ax.set_zlabel("Z [m]", fontsize=10, fontweight='bold')
    ax.set_title("Mappatura 3D Campo B su Rete Stirata Ferromagnetica Anisotropa\nGuida del Flusso ad Alta Permeabilità e Rottura delle Correnti di Foucault", fontsize=12, fontweight='bold')
    
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
    out_png = FIGURES_DIR / "fig_11_mappatura_3d_rete_ferromagnetica.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"  [FIGURA] Salvata: {out_png}")


def main():
    print("=" * 85)
    print("STUDIO ELETTRODINAMICO: RETE STIRATA FERROMAGNETICA ANISOTROPA")
    print("Variante Mantello Chiuso con Rete a 30° ad Alta Permeabilità (mu_r = 1000)")
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
    
    mesh_data = {
        "cells": cells, "pts": pts, "elem_vols": elem_vols,
        "mask_lat": (elem_r >= 0.0465) & (elem_r <= 0.0505) & (np.abs(elem_z) <= 0.050),
        "mask_top": (elem_z > 0.050) & (elem_z <= 0.0535) & (elem_r <= 0.0505),
        "mask_bot": (elem_z < -0.050) & (elem_z >= -0.0535) & (elem_r <= 0.0505),
        "mask_inner": (elem_r < 0.047) & (np.abs(elem_z) <= 0.050)
    }
    
    # 1. Esecuzione simulazione
    vtus = run_simulation()
    
    # 2. Estrazione perdite e forze
    print("  [POST-PROC] Estrazione perdite Joule e forza assiale nei 4 settori...")
    mesh_electro = extract_forces_and_losses(vtus, mesh_data)
    
    # 3. Campionamento sferico di Fibonacci
    print("  [POST-PROC] Campionamento sferico 3D a 360° (N=1200 punti per sfera)...")
    m_vtu0 = meshio.read(str(vtus[0]))
    pts_vtu = m_vtu0.points
    print("    - Costruzione triangolazione Delaunay 3D dei nodi VTU...")
    delaunay_tri = Delaunay(pts_vtu)
    
    mesh_sphere_raw = []
    mesh_spheres_json = []
    
    for r_m, r_lbl in zip(RADII, RADII_LABELS):
        print(f"    - Campionamento {r_lbl}...")
        sp = generate_fibonacci_sphere(r_m, n_points=N_FIBONACCI)
        s_data = sample_sphere_data(sp, vtus, delaunay_tri)
        mesh_sphere_raw.append(s_data)
        
        mesh_spheres_json.append({
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
        
    # 4. Caricamento dati precedenti (Alluminio e Ferro Pieno)
    al_json_path = DATA_DIR / "risultati_handover_terzi.json"
    fe_solid_json_path = DATA_DIR / "risultati_mantello_ferromagnetico.json"
    
    with open(al_json_path, "r", encoding="utf-8") as f:
        al_json = json.load(f)
    with open(fe_solid_json_path, "r", encoding="utf-8") as f:
        fe_solid_json = json.load(f)
        
    al_electro = al_json["electrodynamic_results"]
    al_spheres = al_json["spherical_radiation"]
    fe_solid_electro = fe_solid_json["electrodynamic_results_ferro"]
    fe_solid_spheres = fe_solid_json["spherical_radiation_ferro"]
    
    # 5. Generazione figure a 300 DPI
    print("  [PLOTS] Generazione grafici comparativi a 3 vie a 300 DPI...")
    plot_figure_10_three_materials(al_electro, fe_solid_electro, mesh_electro, al_spheres, fe_solid_spheres, mesh_spheres_json)
    plot_figure_11_ferro_mesh_3d(mesh_sphere_raw)
    
    # 6. Salvataggio dataset consolidato
    out_json_path = DATA_DIR / "risultati_rete_ferromagnetica.json"
    results_consolidated = {
        "metadata": {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "author": "Alessandro Brescacin",
            "license": "CERN-OHL-S-2.0",
            "description": "Confronto a 3 vie: Alluminio Anisotropo vs Ferro Pieno Isotropo vs Rete Stirata Ferromagnetica Anisotropa a 30 gradi"
        },
        "material_properties": {
            "aluminum_mesh": {
                "relative_permeability": 1.0,
                "conductivity_tensor": "anisotropic 30 deg"
            },
            "solid_iron": {
                "relative_permeability": 1000.0,
                "conductivity": "isotropic 1.0e7 S/m"
            },
            "ferromagnetic_mesh": {
                "relative_permeability": 1000.0,
                "conductivity_tensor": "anisotropic 30 deg (sigma_tz = 3.031e6 S/m)"
            }
        },
        "electrodynamic_results_mesh_ferro": mesh_electro,
        "spherical_radiation_mesh_ferro": mesh_spheres_json,
        "three_way_comparison": {
            "Pj_mW": {
                "aluminum": al_electro["total"]["mean_Pj_mW"],
                "solid_iron": fe_solid_electro["total"]["mean_Pj_mW"],
                "ferro_mesh": mesh_electro["total"]["mean_Pj_mW"],
                "saving_vs_solid_pct": (fe_solid_electro["total"]["mean_Pj_mW"] - mesh_electro["total"]["mean_Pj_mW"]) / fe_solid_electro["total"]["mean_Pj_mW"] * 100.0
            },
            "Fz_mean_uN": {
                "aluminum": al_electro["total"]["mean_Fz_uN"],
                "solid_iron": fe_solid_electro["total"]["mean_Fz_uN"],
                "ferro_mesh": mesh_electro["total"]["mean_Fz_uN"]
            },
            "Fz_peak_uN": {
                "aluminum": al_electro["total"]["peak_Fz_uN"],
                "solid_iron": fe_solid_electro["total"]["peak_Fz_uN"],
                "ferro_mesh": mesh_electro["total"]["peak_Fz_uN"]
            },
            "B_far_field_uT": {
                "aluminum": al_spheres[-1]["B_mean_mag_uT"],
                "solid_iron": fe_solid_spheres[-1]["B_mean_mag_uT"],
                "ferro_mesh": mesh_spheres_json[-1]["B_mean_mag_uT"]
            }
        }
    }
    
    with open(out_json_path, "w", encoding="utf-8") as f:
        json.dump(results_consolidated, f, indent=2)
    print(f"  [DATA] Dati salvati con successo in: {out_json_path}")
    
    # 7. Report a video
    tw = results_consolidated["three_way_comparison"]
    print("=" * 85)
    print("SINTESI COMPARATIVA A 3 VIE:")
    print(f"  - Perdite Joule medie P_J:")
    print(f"    * Alluminio Anisotropo:    {tw['Pj_mW']['aluminum']:.4f} mW ({tw['Pj_mW']['aluminum']*1000:.1f} µW)")
    print(f"    * Ferro Pieno Isotropo:    {tw['Pj_mW']['solid_iron']:.4f} mW ({tw['Pj_mW']['solid_iron']*1000:.1f} µW)")
    print(f"    * Rete Ferromagnetica:     {tw['Pj_mW']['ferro_mesh']:.4f} mW ({tw['Pj_mW']['ferro_mesh']*1000:.1f} µW)")
    print(f"    -> Risparmio Rete vs Ferro Pieno: {tw['Pj_mW']['saving_vs_solid_pct']:.2f}% di perdite abbattute!")
    print(f"  - Spinta Lorentz media <F_z>:")
    print(f"    * Alluminio:               {tw['Fz_mean_uN']['aluminum']:+.4f} µN")
    print(f"    * Ferro Pieno:             {tw['Fz_mean_uN']['solid_iron']:+.4f} µN")
    print(f"    * Rete Ferromagnetica:     {tw['Fz_mean_uN']['ferro_mesh']:+.4f} µN")
    print(f"  - Campo Far-Field |B| (15 cm):")
    print(f"    * Alluminio:               {tw['B_far_field_uT']['aluminum']:.3f} µT")
    print(f"    * Ferro Pieno:             {tw['B_far_field_uT']['solid_iron']:.3f} µT")
    print(f"    * Rete Ferromagnetica:     {tw['B_far_field_uT']['ferro_mesh']:.3f} µT")
    print("=" * 85)


if __name__ == "__main__":
    main()
