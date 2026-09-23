#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico e Diagnostica 3D:
Sequenza ad Handover a Terzi (33.3% / 66.7% / 100%) a 3 Coppie Diametrali (180°).
Variante Mantello Chiuso a Barattolo (Top & Bottom Lids a Z = ±H/2).

1. Esecuzione simulazione Elmer FEM transiente (64 timestep, dt = 0.25 ms, T_tot = 16 ms).
2. Estrazione forme d'onda temporali di corrente, perdite Joule P_J(t) e forza di Lorentz F_z(t).
3. Mappatura sferica di Fibonacci a 360° su 3 raggi (6.5 cm, 10.0 cm, 15.0 cm).
4. Verifica legge di Gauss (solenoidalità) e calcolo potenza irradiata di Poynting P_rad.
5. Generazione deliverable grafici a 300 DPI:
   - fig_05_forme_onda_terzi_handover.png
   - fig_06_mappatura_3d_cuspide_divergente.png
   - fig_07_diagrammi_radiazione_handover.png
6. Salvataggio dataset in data/risultati_handover_terzi.json.

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

# Configurazione percorsi
SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_DIR = SCRIPT_DIR.parent
CONFIG_DIR = VARIANT_DIR / "config"
DATA_DIR = VARIANT_DIR / "data"
FIGURES_DIR = VARIANT_DIR / "figures"
WORK_DIR = VARIANT_DIR / "work_dirs" / "run_pulsed_third_handover"
MESH_DIR = VARIANT_DIR / "mesh"

ELMER_SOLVER = r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"
MU0 = 4.0 * np.pi * 1e-7
TIMESTEPS = 64
DT = 0.00025  # 0.25 ms
F_HZ = 100.0
W_RAD = 2.0 * np.pi * F_HZ
J0 = 1.0e5    # A/m^2

# Parametri handover a terzi
TH1 = np.arcsin(1.0 / 3.0)
TH2 = np.pi - np.arcsin(2.0 / 3.0)
TH3 = np.pi / 2.0
DELTA_TH = TH1 + TH2 + TH3
T_CYC = DELTA_TH / W_RAD
TAU = np.pi / W_RAD
TB0 = TH1 / W_RAD
TC0 = (TH1 + TH2) / W_RAD

N_FIBONACCI = 1200
RADII = [0.065, 0.100, 0.150]
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]


def get_analytical_current(t, t_offset):
    """Calcola la densità di corrente normalizzata per un canale con trigger offset."""
    if t < t_offset:
        return 0.0
    t_shift = t - t_offset
    t_rel = t_shift - np.floor(t_shift / T_CYC) * T_CYC
    if t_rel <= TAU:
        return np.sin(W_RAD * t_rel)
    else:
        return 0.0


def run_simulation():
    """Configura ed esegue la simulazione Elmer FEM transiente."""
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    res_dir = WORK_DIR / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    
    vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} file VTU in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]
        
    print(f"  [SOLVER] Preparazione ed esecuzione ElmerSolver ({TIMESTEPS} timestep da {DT*1000:.2f} ms)...")
    sif_src = CONFIG_DIR / "case_pulsed_third_handover.sif"
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
    
    # Verifica Gauss
    phi_net_mean = float(np.sum(b_rad * da))
    phi_abs_mean = float(np.sum(np.abs(b_rad) * da))
    rel_res = abs(phi_net_mean) / (phi_abs_mean + 1e-15)
    
    # Potenza irradiata totale
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
        "gauss_pass": bool(rel_res < 0.02),
        "P_rad_W": p_rad_total_W,
        "P_rad_uW": p_rad_total_W * 1e6
    }


def plot_figure_05_waveforms(time_ms, j_series, electro_res):
    """Genera fig_05_forme_onda_terzi_handover.png (300 DPI)."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), dpi=300, sharex=True)
    
    # 1. Grafico Correnti
    ax1.plot(time_ms, j_series["j_A"] * 100, color="#1f77b4", lw=2.2, label=r"Coppia A (Bobina 1: 0°, Bobina 4: 180°)")
    ax1.plot(time_ms, j_series["j_B"] * 100, color="#ff7f0e", lw=2.2, label=r"Coppia B (Bobina 2: 60°, Bobina 5: 240°)")
    ax1.plot(time_ms, j_series["j_C"] * 100, color="#2ca02c", lw=2.2, label=r"Coppia C (Bobina 3: 120°, Bobina 6: 300°)")
    
    # Linee di soglia orizzontali
    ax1.axhline(33.333, color="#ff7f0e", ls="--", lw=1.2, alpha=0.85, label="Soglia Trigger Coppia B: 33.3% (1/3)")
    ax1.axhline(66.666, color="#2ca02c", ls="--", lw=1.2, alpha=0.85, label="Soglia Trigger Coppia C: 66.7% (2/3)")
    ax1.axhline(100.0, color="#1f77b4", ls=":", lw=1.2, alpha=0.85, label="Soglia Chiusura Ciclo: 100% (Picco)")
    
    # Linee di transizione verticali e annotazioni
    ax1.axvline(TB0 * 1000, color="#ff7f0e", ls="-.", lw=1.0, alpha=0.7)
    ax1.axvline(TC0 * 1000, color="#2ca02c", ls="-.", lw=1.0, alpha=0.7)
    ax1.axvline(T_CYC * 1000, color="#1f77b4", ls="-.", lw=1.0, alpha=0.7)
    
    ax1.scatter([TB0 * 1000], [33.333], color="#ff7f0e", s=60, zorder=5)
    ax1.scatter([TC0 * 1000], [66.666], color="#2ca02c", s=60, zorder=5)
    ax1.scatter([T_CYC * 1000], [100.0], color="#2ca02c", s=60, zorder=5)
    ax1.scatter([T_CYC * 1000], [0.0], color="#1f77b4", s=60, zorder=5)
    
    ax1.annotate(f"Trigger Coppia B (33.3% salita)\nt = {TB0*1000:.2f} ms",
                 xy=(TB0 * 1000, 33.333), xytext=(TB0 * 1000 + 0.35, 45),
                 arrowprops=dict(arrowstyle="->", color="#ff7f0e", lw=1.2),
                 fontsize=9, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ff7f0e", alpha=0.9))
                 
    ax1.annotate(f"Trigger Coppia C (66.7% discesa)\nt = {TC0*1000:.2f} ms",
                 xy=(TC0 * 1000, 66.666), xytext=(TC0 * 1000 - 2.8, 80),
                 arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1.2),
                 fontsize=9, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#2ca02c", alpha=0.9))

    ax1.annotate(f"Chiusura Ciclo su A (100% picco C)\nt = {T_CYC*1000:.2f} ms",
                 xy=(T_CYC * 1000, 100.0), xytext=(T_CYC * 1000 + 0.35, 92),
                 arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1.2),
                 fontsize=9, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#1f77b4", alpha=0.9))

    ax1.set_ylabel(r"Intensità di Corrente $I_k / I_{\max}$ [%]", fontsize=11, fontweight='bold')
    ax1.set_title(r"Sequenza ad Handover Asimmetrico a Terzi (33.3% / 66.7% / 100%) - 3 Coppie Diametrali 180°", fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
    ax1.set_ylim(-5, 115)
    
    # 2. Grafico Perdite Joule e Forza Assiale Fz
    ax2_r = ax2.twinx()
    
    p_pj, = ax2.plot(time_ms, electro_res["total"]["Pj_series_mW"], color="#d62728", lw=2.0, label="Perdite Joule Istantanee P_J(t)")
    p_fz, = ax2_r.plot(time_ms, electro_res["total"]["Fz_series_uN"], color="#9467bd", lw=2.0, ls="-", label="Forza Assiale Lorentz F_z(t)")
    
    mean_pj = electro_res["total"]["mean_Pj_mW"]
    mean_fz = electro_res["total"]["mean_Fz_uN"]
    
    ax2.axhline(mean_pj, color="#d62728", ls=":", lw=1.2, alpha=0.8, label=f"P_J medio = {mean_pj:.3f} mW ({mean_pj*1000:.1f} µW)")
    ax2_r.axhline(mean_fz, color="#9467bd", ls=":", lw=1.2, alpha=0.8, label=rf"$\langle F_z \rangle$ = {mean_fz:.3f} µN")
    
    ax2.set_xlabel(r"Tempo $t$ [ms]", fontsize=11, fontweight='bold')
    ax2.set_ylabel(r"Potenza Dissipata $P_J$ [mW]", fontsize=11, fontweight='bold', color="#d62728")
    ax2_r.set_ylabel(r"Forza di Lorentz $F_z$ [µN]", fontsize=11, fontweight='bold', color="#9467bd")
    
    ax2.tick_params(axis='y', labelcolor="#d62728")
    ax2_r.tick_params(axis='y', labelcolor="#9467bd")
    ax2.grid(True, linestyle=":", alpha=0.6)
    
    lines = [p_pj, p_fz]
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc="upper right", fontsize=9, framealpha=0.95)
    
    ax2.set_xlim(0, time_ms[-1])
    plt.tight_layout()
    
    out_png = FIGURES_DIR / "fig_05_forme_onda_terzi_handover.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"  [FIGURA] Salvata: {out_png}")


def plot_figure_06_cusp_field_3d(sphere_results):
    """Genera fig_06_mappatura_3d_cuspide_divergente.png (300 DPI)."""
    fig = plt.figure(figsize=(13, 10), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    
    # 1. Disegna silhouette trasparente del barattolo (R = 50 mm, H = 100 mm)
    z_cyl = np.linspace(-0.05, 0.05, 30)
    th_cyl = np.linspace(0, 2*np.pi, 50)
    Z_mesh, TH_mesh = np.meshgrid(z_cyl, th_cyl)
    X_cyl = 0.050 * np.cos(TH_mesh)
    Y_cyl = 0.050 * np.sin(TH_mesh)
    ax.plot_surface(X_cyl, Y_cyl, Z_mesh, alpha=0.10, color="gray", edgecolor="none")
    
    # Coperchi del barattolo
    r_disc = np.linspace(0, 0.050, 20)
    R_disc, TH_disc = np.meshgrid(r_disc, th_cyl)
    X_disc = R_disc * np.cos(TH_disc)
    Y_disc = R_disc * np.sin(TH_disc)
    ax.plot_surface(X_disc, Y_disc, np.full_like(X_disc, 0.050), alpha=0.12, color="gray")
    ax.plot_surface(X_disc, Y_disc, np.full_like(X_disc, -0.050), alpha=0.12, color="gray")
    
    # 2. Campioni vettoriali sulla sfera Near-Field (R = 6.5 cm)
    s_data = sphere_results[0]  # Near-Field
    coords = s_data["sphere"]["coords"]
    b_mean = s_data["b_mean"]
    b_mag = s_data["b_mag"]
    b_rad = s_data["b_rad"]
    
    # Seleziona sottocampione diradato per leggibilità 3D (ogni 4 punti)
    step = 4
    x_sub = coords[::step, 0]
    y_sub = coords[::step, 1]
    z_sub = coords[::step, 2]
    u_sub = b_mean[::step, 0]
    v_sub = b_mean[::step, 1]
    w_sub = b_mean[::step, 2]
    mag_sub = b_mag[::step]
    rad_sub = b_rad[::step]
    
    # Normalizza lunghezza frecce per visualizzazione chiara
    scale = 0.015 / (np.max(mag_sub) + 1e-12)
    
    # Colora per componente radiale (divergente / uscente vs entrante)
    norm = plt.Normalize(vmin=-np.max(np.abs(rad_sub)), vmax=np.max(np.abs(rad_sub)))
    cmap = plt.cm.coolwarm
    colors = cmap(norm(rad_sub))
    
    q = ax.quiver(x_sub, y_sub, z_sub, u_sub*scale, v_sub*scale, w_sub*scale,
                  colors=colors, arrow_length_ratio=0.35, linewidth=1.1, alpha=0.85)
    
    # Aggiungi punti di campionamento
    ax.scatter(x_sub, y_sub, z_sub, c=rad_sub, cmap='coolwarm', norm=norm, s=12, alpha=0.7)
    
    # Annotazioni frecce di espulsione divergente
    ax.text(0, 0, 0.085, "↑ Flusso Uscente Superiore (+Z)\n(Espulsione Coperchio Sup.)", color="#b2182b", fontsize=9, fontweight='bold', ha='center')
    ax.text(0, 0, -0.095, "↓ Flusso Uscente Inferiore (-Z)\n(Espulsione Coperchio Inf.)", color="#b2182b", fontsize=9, fontweight='bold', ha='center')
    ax.text(0.080, 0, 0, "→ Espulsione Radiale Equatoriale", color="#b2182b", fontsize=9, fontweight='bold', ha='left')
    
    ax.set_xlabel("X [m]", fontsize=10, fontweight='bold')
    ax.set_ylabel("Y [m]", fontsize=10, fontweight='bold')
    ax.set_zlabel("Z [m]", fontsize=10, fontweight='bold')
    ax.set_title("Mappatura 3D Campo Magnetico B e Topologia a Cuspide Divergente (360°)\nSequenza Pulsata ad Handover a Terzi su Mantello Chiuso", fontsize=12, fontweight='bold')
    
    # Imposta limiti assi isotropi
    lim = 0.09
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    ax.set_zlim([-lim, lim])
    ax.view_init(elev=22, azim=48)
    
    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label(r"Induzione Radiale $B_{\text{rad}}$ [T] (Rosso = Uscente / Divergente)", fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    out_png = FIGURES_DIR / "fig_06_mappatura_3d_cuspide_divergente.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"  [FIGURA] Salvata: {out_png}")


def plot_figure_07_radiation_patterns(sphere_results):
    """Genera fig_07_diagrammi_radiazione_handover.png (300 DPI)."""
    fig, (ax_eq, ax_mer) = plt.subplots(1, 2, figsize=(14, 6.5), dpi=300, subplot_kw=dict(polar=True))
    
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    
    # 1. Diagramma Polare Equatoriale (Z = 0, variazione lungo azimut phi)
    for idx, s_res in enumerate(sphere_results):
        sp = s_res["sphere"]
        theta_p = sp["theta_polar"]
        phi_a = sp["phi_azimuth"]
        s_rad = s_res["s_rad"]
        
        # Filtra fascia equatoriale (|theta_p - pi/2| < 0.25 rad)
        mask_eq = np.abs(theta_p - np.pi/2) < 0.25
        phi_sort = np.argsort(phi_a[mask_eq])
        phi_eq = phi_a[mask_eq][phi_sort]
        s_eq = s_rad[mask_eq][phi_sort]
        
        # Chiudi cerchio
        phi_eq = np.append(phi_eq, phi_eq[0] + 2*np.pi)
        s_eq = np.append(s_eq, s_eq[0])
        
        # Scala mW/m^2
        ax_eq.plot(phi_eq, np.abs(s_eq) * 1000, color=colors[idx], lw=2.0, label=RADII_LABELS[idx])
        
    ax_eq.set_title("Diagramma Polare Equatoriale (Piano XY, Z=0)\nFlusso Radiale di Poynting |S_rad| [mW/m²]", fontsize=11, fontweight='bold', pad=15)
    ax_eq.set_theta_zero_location('E')
    ax_eq.grid(True, linestyle=":", alpha=0.6)
    ax_eq.legend(loc="lower right", bbox_to_anchor=(1.15, -0.15), fontsize=8.5)
    
    # 2. Diagramma Polare Meridiano (variazione lungo angolo polare theta)
    for idx, s_res in enumerate(sphere_results):
        sp = s_res["sphere"]
        theta_p = sp["theta_polar"]
        phi_a = sp["phi_azimuth"]
        s_rad = s_res["s_rad"]
        
        # Filtra sezione meridiana (phi vicino a 0 o pi)
        mask_mer1 = (phi_a < 0.3) | (phi_a > (2*np.pi - 0.3))
        mask_mer2 = (np.abs(phi_a - np.pi) < 0.3)
        
        th_mer1 = theta_p[mask_mer1]
        s_mer1 = s_rad[mask_mer1]
        th_mer2 = 2*np.pi - theta_p[mask_mer2]
        s_mer2 = s_rad[mask_mer2]
        
        th_all = np.concatenate([th_mer1, th_mer2])
        s_all = np.concatenate([s_mer1, s_mer2])
        
        idx_sort = np.argsort(th_all)
        th_sorted = th_all[idx_sort]
        s_sorted = s_all[idx_sort]
        
        th_sorted = np.append(th_sorted, th_sorted[0] + 2*np.pi)
        s_sorted = np.append(s_sorted, s_sorted[0])
        
        ax_mer.plot(th_sorted, np.abs(s_sorted) * 1000, color=colors[idx], lw=2.0, label=RADII_LABELS[idx])
        
    ax_mer.set_title("Diagramma Polare Meridiano (Piano XZ, Y=0)\nConfinamento Assiale Lobi da Coperchi [mW/m²]", fontsize=11, fontweight='bold', pad=15)
    ax_mer.set_theta_zero_location('N')  # 0 deg = +Z
    ax_mer.set_theta_direction(-1)       # Senso orario
    ax_mer.grid(True, linestyle=":", alpha=0.6)
    ax_mer.legend(loc="lower right", bbox_to_anchor=(1.15, -0.15), fontsize=8.5)
    
    plt.tight_layout()
    out_png = FIGURES_DIR / "fig_07_diagrammi_radiazione_handover.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"  [FIGURA] Salvata: {out_png}")


def main():
    print("=" * 85)
    print("STUDIO ELETTRODINAMICO: SEQUENZA AD HANDOVER A TERZI A 3 COPPIE DIAMETRALI")
    print("Variante Mantello Chiuso a Barattolo con Coperchi in Rete Stirata a Z = ±H/2")
    print("=" * 85)
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Caricamento mesh
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
    
    # 2. Esecuzione simulazione
    vtus = run_simulation()
    
    # 3. Estrazione perdite e forze
    print("  [POST-PROC] Estrazione perdite Joule e forza assiale nei 4 settori...")
    electro_res = extract_forces_and_losses(vtus, mesh_data)
    
    # 4. Forme d'onda di corrente
    time_pts = np.arange(1, TIMESTEPS + 1) * DT
    time_ms = time_pts * 1000.0
    
    j_A = np.array([get_analytical_current(t, 0.0) for t in time_pts])
    j_B = np.array([get_analytical_current(t, TB0) for t in time_pts])
    j_C = np.array([get_analytical_current(t, TC0) for t in time_pts])
    
    j_series = {
        "time_ms": time_ms.tolist(),
        "j_A": j_A, "j_B": j_B, "j_C": j_C
    }
    
    # 5. Campionamento sferico di Fibonacci
    print("  [POST-PROC] Campionamento sferico 3D a 360° (N=1200 punti per sfera)...")
    m_vtu0 = meshio.read(str(vtus[0]))
    pts_vtu = m_vtu0.points
    print("    - Costruzione triangolazione Delaunay 3D dei nodi VTU...")
    delaunay_tri = Delaunay(pts_vtu)
    sphere_results = []
    spheres_data_json = []
    
    for r_m, r_lbl in zip(RADII, RADII_LABELS):
        print(f"    - Campionamento {r_lbl}...")
        sp = generate_fibonacci_sphere(r_m, n_points=N_FIBONACCI)
        s_data = sample_sphere_data(sp, vtus, delaunay_tri)
        sphere_results.append(s_data)
        
        spheres_data_json.append({
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
        
    # 6. Generazione Figure a 300 DPI
    print("  [PLOTS] Generazione grafici a 300 DPI...")
    plot_figure_05_waveforms(time_ms, j_series, electro_res)
    plot_figure_06_cusp_field_3d(sphere_results)
    plot_figure_07_radiation_patterns(sphere_results)
    
    # 7. Salvataggio dati consolidati JSON
    output_json = DATA_DIR / "risultati_handover_terzi.json"
    results_consolidated = {
        "metadata": {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "author": "Alessandro Brescacin",
            "license": "CERN-OHL-S-2.0",
            "description": "Simulazione transiente sequenza ad handover a terzi (33.3% / 66.7% / 100%) con diodi contrapposti su mantello chiuso"
        },
        "handover_parameters": {
            "f_base_Hz": F_HZ,
            "omega_rad_s": W_RAD,
            "threshold_1_pct": 33.333,
            "threshold_2_pct": 66.666,
            "threshold_3_pct": 100.0,
            "T_cycle_ms": float(T_CYC * 1000.0),
            "t_B0_ms": float(TB0 * 1000.0),
            "t_C0_ms": float(TC0 * 1000.0),
            "pulse_duration_tau_ms": float(TAU * 1000.0),
            "timesteps": TIMESTEPS,
            "dt_ms": float(DT * 1000.0),
            "total_time_ms": float(TIMESTEPS * DT * 1000.0),
            "cycles_completed": float(TIMESTEPS * DT / T_CYC)
        },
        "electrodynamic_results": electro_res,
        "spherical_radiation": spheres_data_json
    }
    
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(results_consolidated, f, indent=2)
    print(f"  [DATA] Dati salvati con successo in: {output_json}")
    
    print("=" * 85)
    print("RISULTATI SINTETICI:")
    print(f"  - Perdite Joule medie P_J: {electro_res['total']['mean_Pj_mW']:.3f} mW ({electro_res['total']['mean_Pj_mW']*1000:.1f} µW)")
    print(f"  - Forza Assiale Lorentz <F_z>: {electro_res['total']['mean_Fz_uN']:+.3f} µN (Peak: {electro_res['total']['peak_Fz_uN']:+.3f} µN)")
    for s in spheres_data_json:
        print(f"  - Sfera {s['radius_cm']} cm: |B| = {s['B_mean_mag_uT']:.1f} µT, P_rad = {s['P_rad_uW']:.1f} µW, Gauss Residuo = {s['gauss_residual_pct']:.2f}% (PASS={s['gauss_pass']})")
    print("=" * 85)


if __name__ == "__main__":
    main()
