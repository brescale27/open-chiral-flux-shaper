#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline di Calcolo, Simulazione Elmer FEM e Post-Processing 3D per:
Gabbia Sferica Metamateriale con Doppio Rotore Ortogonale a 90° (Multi-Axis Vector Shaper).

Funzionalità:
1. Esecuzione simulazione transiente Elmer FEM su 64 timestep (dt = 0.25 ms, f = 100 Hz, T_tot = 16.0 ms)
2. Estrazione serie temporale delle forze vettoriali di Lorentz F_x(t), F_y(t), F_z(t) e perdite Joule P_J(t)
   scomposte sui 6 corpi fisici e sui 3 sottostrati del mantello (+30° / 0° / -30°)
3. Calcolo indipendente del Tensore degli Sforzi di Maxwell (MST) su superficie sferica di controllo in aria
4. Campionamento sferico di Fibonacci ad alta densità (N = 2500 punti) su 3 sfere concentriche (R = 6.5, 10.0, 15.0 cm):
   - Verifica della solenoidalità di Gauss (div B = 0, residuo < 2.0%)
   - Calcolo del flusso attivo tridimensionale del vettore di Poynting S = (E x B) / mu0
5. Generazione deliverable diagnostici ad alta risoluzione (300 DPI):
   - fig_16_mappatura_3D_doppio_rotore_sferico.png
   - fig_17_matrice_forze_ortogonali.png
6. Esportazione dataset in data/risultati_gabbia_sferica_ortogonale.json

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
ROOT_FIGURES_DIR = VARIANT_DIR.parent.parent / "figures"
WORK_DIR = VARIANT_DIR / "work_dirs" / "run_gabbia_sferica"
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

N_FIBONACCI = 2500
RADII = [0.065, 0.100, 0.150]
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]


def run_simulation():
    """Configura ed esegue la simulazione Elmer FEM per la variante sferica a doppio rotore ortogonale."""
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    res_dir = WORK_DIR / "results"
    res_dir.mkdir(parents=True, exist_ok=True)

    vtus = sorted(list(res_dir.glob("macchina_sferica_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} file VTU in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]

    print(f"  [SOLVER] Esecuzione ElmerSolver per gabbia sferica ({TIMESTEPS} timestep da {DT*1000:.2f} ms)...")
    sif_src = CONFIG_DIR / "case_gabbia_sferica_doppio_rotore.sif"
    sif_dest = WORK_DIR / "case.sif"

    mesh_db_rel = os.path.relpath(str(MESH_DIR), str(WORK_DIR)).replace('\\', '/')
    lines = sif_src.read_text(encoding="utf-8").splitlines()
    new_lines = []
    for line in lines:
        if 'Mesh DB' in line:
            new_lines.append(f'  Mesh DB "{mesh_db_rel}" "macchina_sferica_ortogonale"')
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
        print(proc.stdout[-1500:])
        sys.exit(1)

    vtus = sorted(list(res_dir.glob("macchina_sferica_out_t*.vtu")))
    print(f"  [SOLVER] Simulazione completata con successo in {duration:.1f}s ({len(vtus)} VTU generati).")
    return vtus[:TIMESTEPS]


def load_mesh_structure(sample_vtu):
    """Carica la mesh ed estrae i filtri geometrici per ciascun corpo fisico e per i 3 sottostrati del mantello."""
    m = meshio.read(str(sample_vtu))
    pts = m.points
    cells = None
    for cb in m.cells:
        if cb.type == 'tetra':
            cells = cb.data
            break

    centroids = np.mean(pts[cells], axis=1)
    x = centroids[:, 0]
    y = centroids[:, 1]
    z = centroids[:, 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    rho = np.sqrt(x**2 + y**2)

    # Calcolo volumi tetraedri
    v0 = pts[cells[:, 0]]
    v1 = pts[cells[:, 1]]
    v2 = pts[cells[:, 2]]
    v3 = pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0

    # Maschere per corpi fisici conformi
    # Mantello: guscio sferico tra 47 mm e 50 mm
    mask_mantle = (r >= 0.0465) & (r <= 0.0505)
    # 3 sottostrati radiali del mantello (1 mm ciascuno)
    mask_l1 = mask_mantle & (r < 0.0480)
    mask_l2 = mask_mantle & (r >= 0.0480) & (r < 0.0490)
    mask_l3 = mask_mantle & (r >= 0.0490)

    # Nucleo centrale PEEK: sfera r <= 12 mm
    mask_core = r <= 0.0125

    # Rotore 1 (Coils Z): raggio coil ~35 mm, vicine al piano Z=0
    mask_coils1 = (r < 0.0465) & (r > 0.025) & (np.abs(z) <= 0.008) & (rho >= 0.028) & (rho <= 0.042)
    
    # Rotore 2 (Coils X): raggio coil ~35 mm nel piano Y-Z (vicine a X=0)
    rho_yz = np.sqrt(y**2 + z**2)
    mask_coils2 = (r < 0.0465) & (r > 0.025) & (np.abs(x) <= 0.008) & (rho_yz >= 0.028) & (rho_yz <= 0.042)

    # Cavità interna restante
    mask_cavity = (r < 0.0465) & (~mask_core) & (~mask_coils1) & (~mask_coils2)
    # Aria esterna
    mask_air_ext = r > 0.0505

    print(f"  [MESH] Tetraedri: {len(cells)}, Nodi: {len(pts)}")
    print(f"    - Mantello Sferico : {np.sum(mask_mantle)} elem (Vol: {np.sum(elem_vols[mask_mantle])*1e6:.2f} cm3)")
    print(f"      * Strato 1 (+30°): {np.sum(mask_l1)} elem")
    print(f"      * Strato 2 (0°)  : {np.sum(mask_l2)} elem")
    print(f"      * Strato 3 (-30°): {np.sum(mask_l3)} elem")
    print(f"    - Rotore 1 (Coils Z): {np.sum(mask_coils1)} elem")
    print(f"    - Rotore 2 (Coils X): {np.sum(mask_coils2)} elem")
    print(f"    - Nucleo PEEK       : {np.sum(mask_core)} elem (Vol: {np.sum(elem_vols[mask_core])*1e6:.2f} cm3)")
    print(f"    - Cavita Aria       : {np.sum(mask_cavity)} elem")
    print(f"    - Aria Esterna      : {np.sum(mask_air_ext)} elem")

    return {
        "points": pts,
        "cells": cells,
        "centroids": centroids,
        "elem_vols": elem_vols,
        "mask_mantle": mask_mantle,
        "mask_l1": mask_l1,
        "mask_l2": mask_l2,
        "mask_l3": mask_l3,
        "mask_coils1": mask_coils1,
        "mask_coils2": mask_coils2,
        "mask_core": mask_core,
        "mask_cavity": mask_cavity,
        "mask_air_ext": mask_air_ext,
    }


def extract_forces_and_losses(vtus, mesh_data):
    """Calcola le componenti vettoriali della forza di Lorentz (Fx, Fy, Fz) e perdite Joule P_J su ciascun settore."""
    cells = mesh_data["cells"]
    elem_vols = mesh_data["elem_vols"]
    masks = {
        "mantle": mesh_data["mask_mantle"],
        "layer1_plus30": mesh_data["mask_l1"],
        "layer2_ortho": mesh_data["mask_l2"],
        "layer3_minus30": mesh_data["mask_l3"],
        "rotor1_coils": mesh_data["mask_coils1"],
        "rotor2_coils": mesh_data["mask_coils2"],
        "peek_core": mesh_data["mask_core"],
        "cavity_air": mesh_data["mask_cavity"],
        "assembly_total": mesh_data["mask_mantle"] | mesh_data["mask_coils1"] | mesh_data["mask_coils2"] | mesh_data["mask_core"]
    }

    forces = {k: {"fx": [], "fy": [], "fz": []} for k in masks}
    joule = {k: [] for k in masks}

    for idx, vtu_file in enumerate(vtus):
        m = meshio.read(str(vtu_file))
        jxb = m.point_data['jxb']  # shape (N, 3)
        pj = m.point_data['joule heating'].ravel()

        jxb_elem = np.mean(jxb[cells], axis=1)  # shape (n_elem, 3)
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

    return forces, joule


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


def compute_spherical_sampling_and_mst(vtus, sample_idx=32):
    """Esegue il campionamento di Fibonacci (N=2500) e calcola Gauss, Poynting e Tensore di Maxwell MST."""
    vtu_file = vtus[sample_idx]
    m = meshio.read(str(vtu_file))
    pts = m.points
    b_field = m.point_data['magnetic flux density']

    # Calcolo fisico del campo elettrico indotto da Faraday: E_ind = -dA/dt
    if sample_idx > 0:
        m_prev = meshio.read(str(vtus[sample_idx - 1]))
        e_field = -(m.point_data['magnetic vector potential'] - m_prev.point_data['magnetic vector potential']) / DT
    else:
        e_field = -(m.point_data['magnetic vector potential']) / DT

    interp_b = LinearNDInterpolator(pts, b_field, fill_value=0.0)
    interp_e = LinearNDInterpolator(pts, e_field, fill_value=0.0)

    results = {}
    for r_val, label in zip(RADII, RADII_LABELS):
        pts_fib, normals, d_area = generate_fibonacci_sphere(N_FIBONACCI, r_val)
        b_samples = interp_b(pts_fib)
        e_samples = interp_e(pts_fib)

        # Gauss Solenoidality: phi = B . n
        b_dot_n = np.einsum('ij,ij->i', b_samples, normals)
        flux_net = float(np.sum(b_dot_n) * d_area)
        flux_abs = float(np.sum(np.abs(b_dot_n)) * d_area)
        residual_pct = float(abs(flux_net) / flux_abs * 100.0) if flux_abs > 0 else 0.0

        # Poynting Vector: S = (E x B) / mu0
        s_vec = np.cross(e_samples, b_samples) / MU0
        s_rad = np.einsum('ij,ij->i', s_vec, normals)
        p_rad_watts = float(np.sum(s_rad) * d_area)

        b_mag = np.linalg.norm(b_samples, axis=1)
        e_mag = np.linalg.norm(e_samples, axis=1)

        # Maxwell Stress Tensor on this sphere:
        # T_ij = (1/mu0) * [ B_i B_j - 0.5 * delta_ij * |B|^2 ]
        # Traction: t_i = T_ij n_j = (1/mu0) * [ B_i (B . n) - 0.5 * |B|^2 n_i ]
        t_mst_x = (b_samples[:, 0] * b_dot_n - 0.5 * (b_mag**2) * normals[:, 0]) / MU0
        t_mst_y = (b_samples[:, 1] * b_dot_n - 0.5 * (b_mag**2) * normals[:, 1]) / MU0
        t_mst_z = (b_samples[:, 2] * b_dot_n - 0.5 * (b_mag**2) * normals[:, 2]) / MU0

        fx_mst = float(np.sum(t_mst_x) * d_area)
        fy_mst = float(np.sum(t_mst_y) * d_area)
        fz_mst = float(np.sum(t_mst_z) * d_area)

        results[label] = {
            "radius_m": r_val,
            "gauss_net_weber": flux_net,
            "gauss_abs_weber": flux_abs,
            "gauss_residual_pct": residual_pct,
            "gauss_pass": bool(residual_pct < 2.0),
            "b_mean_uT": float(np.mean(b_mag) * 1e6),
            "b_max_uT": float(np.max(b_mag) * 1e6),
            "e_mean_mV_m": float(np.mean(e_mag) * 1e3),
            "poynting_rad_mW": float(p_rad_watts * 1e3),
            "mst_fx_uN": float(fx_mst * 1e6),
            "mst_fy_uN": float(fy_mst * 1e6),
            "mst_fz_uN": float(fz_mst * 1e6),
            "b_samples": b_samples,
            "normals": normals,
            "pts": pts_fib
        }

    return results


def plot_diagnostics(forces, joule, sampling_results, mesh_data):
    """Genera le due figure ad alta risoluzione (300 DPI): fig_16 e fig_17."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    ROOT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    t_ms = np.arange(TIMESTEPS) * DT * 1000.0

    # -------------------------------------------------------------------------
    # FIGURA 16: Mappatura 3D Doppio Rotore Sferico e Cupola di Flusso
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(16, 8), dpi=300)
    fig.patch.set_facecolor('#0d1117')

    # Subplot 1: Vista 3D Architettura e Bobine Ortogonali
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.set_facecolor('#0d1117')

    # Wireframe sfera mantello
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 15)
    r_can = 50.0  # mm
    xs = r_can * np.outer(np.cos(u), np.sin(v))
    ys = r_can * np.outer(np.sin(u), np.sin(v))
    zs = r_can * np.outer(np.ones(np.size(u)), np.cos(v))
    ax1.plot_wireframe(xs, ys, zs, color='#58a6ff', alpha=0.12, linewidth=0.5)

    # Bobine Rotore 1 (Asse Z, Piano Z=0, R=35 mm) -> Cyan/Blu
    for deg in [0, 60, 120, 180, 240, 300]:
        rad = np.radians(deg)
        xc = 35.0 * np.cos(rad)
        yc = 35.0 * np.sin(rad)
        zc = np.linspace(-7, 7, 10)
        theta_c = np.linspace(0, 2*np.pi, 12)
        X_c = xc + 4.0 * np.cos(theta_c)[:, None]
        Y_c = yc + 4.0 * np.sin(theta_c)[:, None]
        Z_c = np.repeat(zc[None, :], 12, axis=0)
        ax1.plot_surface(X_c, Y_c, Z_c, color='#388bfd', alpha=0.85, shade=True)

    # Bobine Rotore 2 (Asse X, Piano X=0, R=35 mm, sfasate 30°) -> Arancio/Rosso
    for deg in [30, 90, 150, 210, 270, 330]:
        rad = np.radians(deg)
        yc = 35.0 * np.cos(rad)
        zc = 35.0 * np.sin(rad)
        xc = np.linspace(-7, 7, 10)
        theta_c = np.linspace(0, 2*np.pi, 12)
        Y_c = yc + 4.0 * np.cos(theta_c)[:, None]
        Z_c = zc + 4.0 * np.sin(theta_c)[:, None]
        X_c = np.repeat(xc[None, :], 12, axis=0)
        ax1.plot_surface(X_c, Y_c, Z_c, color='#f0883e', alpha=0.85, shade=True)

    # Nucleo Centrale PEEK (Verde)
    rc = 12.0
    xc_core = rc * np.outer(np.cos(u), np.sin(v))
    yc_core = rc * np.outer(np.sin(u), np.sin(v))
    zc_core = rc * np.outer(np.ones(np.size(u)), np.cos(v))
    ax1.plot_surface(xc_core, yc_core, zc_core, color='#3fb950', alpha=0.6, shade=True)

    ax1.set_xlim(-60, 60)
    ax1.set_ylim(-60, 60)
    ax1.set_zlim(-60, 60)
    ax1.set_xlabel('X [mm]', color='#c9d1d9', fontsize=10)
    ax1.set_ylabel('Y [mm]', color='#c9d1d9', fontsize=10)
    ax1.set_zlabel('Z [mm]', color='#c9d1d9', fontsize=10)
    ax1.tick_params(colors='#8b949e', labelsize=8)
    ax1.set_title("Architettura CAD 3D: Gabbia Sferica e Rotori a 90°\nBlu: Rotore 1 (Asse Z) | Arancio: Rotore 2 (Asse X) | Verde: PEEK Core", color='#f0f6fc', fontsize=11, fontweight='bold')

    # Subplot 2: Campo Magnetico 3D a Cupola e Fibonacci Sampling a R=6.5 cm
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.set_facecolor('#0d1117')

    res_nf = sampling_results["Near-Field (R=6.5 cm)"]
    pts_nf = res_nf["pts"] * 1000.0  # mm
    b_nf = res_nf["b_samples"]
    b_mag_nf = np.linalg.norm(b_nf, axis=1) * 1e6  # uT

    sc = ax2.scatter(pts_nf[:, 0], pts_nf[:, 1], pts_nf[:, 2], c=b_mag_nf, cmap='plasma', s=8, alpha=0.8)
    cbar = fig.colorbar(sc, ax=ax2, shrink=0.6, pad=0.1)
    cbar.set_label('|B| [µT] a R = 65 mm', color='#c9d1d9', fontsize=10)
    cbar.ax.yaxis.set_tick_params(color='#8b949e')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='#c9d1d9')

    # Campione vettori B
    skip = 60
    b_unit = b_nf / (np.linalg.norm(b_nf, axis=1, keepdims=True) + 1e-12) * 12.0
    ax2.quiver(pts_nf[::skip, 0], pts_nf[::skip, 1], pts_nf[::skip, 2],
               b_unit[::skip, 0], b_unit[::skip, 1], b_unit[::skip, 2],
               color='#7ee787', length=1.0, normalize=False, linewidth=1.0, alpha=0.9)

    ax2.set_xlim(-80, 80)
    ax2.set_ylim(-80, 80)
    ax2.set_zlim(-80, 80)
    ax2.set_xlabel('X [mm]', color='#c9d1d9', fontsize=10)
    ax2.set_ylabel('Y [mm]', color='#c9d1d9', fontsize=10)
    ax2.set_zlabel('Z [mm]', color='#c9d1d9', fontsize=10)
    ax2.tick_params(colors='#8b949e', labelsize=8)
    ax2.set_title(f"Mappatura 3D Cupola di Flusso (Fibonacci N=2500, R=6.5 cm)\nMean |B| = {res_nf['b_mean_uT']:.1f} µT, Peak = {res_nf['b_max_uT']:.1f} µT", color='#f0f6fc', fontsize=11, fontweight='bold')

    plt.tight_layout()
    fig16_path = FIGURES_DIR / "fig_16_mappatura_3D_doppio_rotore_sferico.png"
    plt.savefig(fig16_path, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    shutil.copy(fig16_path, ROOT_FIGURES_DIR / "fig_16_mappatura_3D_doppio_rotore_sferico.png")
    plt.close()
    print(f"  [FIGURE] Salvata {fig16_path}")

    # -------------------------------------------------------------------------
    # FIGURA 17: Matrice Forze Vettoriali sui 3 Assi, Odografo e MST
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(15, 11), dpi=300)
    fig.patch.set_facecolor('#0d1117')

    fx_tot = np.array(forces["assembly_total"]["fx"]) * 1e6
    fy_tot = np.array(forces["assembly_total"]["fy"]) * 1e6
    fz_tot = np.array(forces["assembly_total"]["fz"]) * 1e6

    # Pannello 1: Forme d'onda delle Forze sui Tre Assi Fx(t), Fy(t), Fz(t)
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_facecolor('#161b22')
    ax1.plot(t_ms, fx_tot, color='#f0883e', linewidth=2.0, label='Fx (Asse Trasversale X)')
    ax1.plot(t_ms, fy_tot, color='#58a6ff', linewidth=1.5, linestyle='--', label='Fy (Asse Intermedio Y)')
    ax1.plot(t_ms, fz_tot, color='#7ee787', linewidth=2.0, label='Fz (Asse Verticale Z)')
    ax1.axhline(0, color='#8b949e', linestyle=':', alpha=0.6)
    ax1.set_xlabel('Tempo [ms]', color='#c9d1d9', fontsize=10)
    ax1.set_ylabel('Forza di Lorentz [µN]', color='#c9d1d9', fontsize=10)
    ax1.set_title('(a) Dinamica Forze Vettoriali sui 3 Assi', color='#f0f6fc', fontsize=11, fontweight='bold')
    ax1.grid(True, color='#30363d', alpha=0.5)
    ax1.tick_params(colors='#8b949e')
    ax1.legend(loc='upper right', facecolor='#0d1117', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=9)

    # Pannello 2: Odografo Vettoriale 3D (Fx, Fy, Fz)
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    ax2.set_facecolor('#161b22')
    ax2.plot(fx_tot, fy_tot, fz_tot, color='#d2a8ff', linewidth=2.0)
    ax2.scatter(fx_tot[0], fy_tot[0], fz_tot[0], color='#7ee787', s=50, label='Start (t=0)')
    ax2.scatter(np.mean(fx_tot), np.mean(fy_tot), np.mean(fz_tot), color='#ff7b72', marker='*', s=150, label='Mean Offset')
    ax2.set_xlabel('Fx [µN]', color='#c9d1d9', fontsize=9)
    ax2.set_ylabel('Fy [µN]', color='#c9d1d9', fontsize=9)
    ax2.set_zlabel('Fz [µN]', color='#c9d1d9', fontsize=9)
    ax2.tick_params(colors='#8b949e', labelsize=8)
    ax2.set_title('(b) Odografo Spazio-Forze 3D (Traiettoria Vettoriale)', color='#f0f6fc', fontsize=11, fontweight='bold')
    ax2.legend(loc='upper right', facecolor='#0d1117', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=8)

    # Pannello 3: Matrice di Confronto Forze per Settore e Tensore di Maxwell MST
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_facecolor('#161b22')
    labels = ['Mantello\nSferico', 'Rotore 1\n(Coils Z)', 'Rotore 2\n(Coils X)', 'PEEK\nCore', 'Totale\nLorentz', 'MST\n(Mid-Field)']
    fx_means = [np.mean(forces["mantle"]["fx"])*1e6, np.mean(forces["rotor1_coils"]["fx"])*1e6, np.mean(forces["rotor2_coils"]["fx"])*1e6, np.mean(forces["peek_core"]["fx"])*1e6, np.mean(fx_tot), sampling_results["Mid-Field (R=10.0 cm)"]["mst_fx_uN"]]
    fy_means = [np.mean(forces["mantle"]["fy"])*1e6, np.mean(forces["rotor1_coils"]["fy"])*1e6, np.mean(forces["rotor2_coils"]["fy"])*1e6, np.mean(forces["peek_core"]["fy"])*1e6, np.mean(fy_tot), sampling_results["Mid-Field (R=10.0 cm)"]["mst_fy_uN"]]
    fz_means = [np.mean(forces["mantle"]["fz"])*1e6, np.mean(forces["rotor1_coils"]["fz"])*1e6, np.mean(forces["rotor2_coils"]["fz"])*1e6, np.mean(forces["peek_core"]["fz"])*1e6, np.mean(fz_tot), sampling_results["Mid-Field (R=10.0 cm)"]["mst_fz_uN"]]

    x_idx = np.arange(len(labels))
    w_bar = 0.25
    ax3.bar(x_idx - w_bar, fx_means, w_bar, label='⟨Fx⟩ [µN]', color='#f0883e', alpha=0.9)
    ax3.bar(x_idx, fy_means, w_bar, label='⟨Fy⟩ [µN]', color='#58a6ff', alpha=0.9)
    ax3.bar(x_idx + w_bar, fz_means, w_bar, label='⟨Fz⟩ [µN]', color='#7ee787', alpha=0.9)
    ax3.set_xticks(x_idx)
    ax3.set_xticklabels(labels, color='#c9d1d9', fontsize=9)
    ax3.set_ylabel('Forza Media [µN]', color='#c9d1d9', fontsize=10)
    ax3.set_title('(c) Ripartizione Elettrodinamica per Settore e MST', color='#f0f6fc', fontsize=11, fontweight='bold')
    ax3.grid(True, color='#30363d', alpha=0.5, axis='y')
    ax3.tick_params(colors='#8b949e')
    ax3.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=9)

    # Pannello 4: Bilancio Termico Joule P_J(t) e Potenza Irradiata Poynting
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.set_facecolor('#161b22')
    pj_mantle = np.array(joule["mantle"]) * 1e6
    pj_r1 = np.array(joule["rotor1_coils"]) * 1e6
    pj_r2 = np.array(joule["rotor2_coils"]) * 1e6
    pj_tot = np.array(joule["assembly_total"]) * 1e6

    ax4.plot(t_ms, pj_tot, color='#ff7b72', linewidth=2.0, label=f'Totale Macchina (⟨P_J⟩ = {np.mean(pj_tot):.2f} µW)')
    ax4.plot(t_ms, pj_mantle, color='#d2a8ff', linewidth=1.5, linestyle='--', label=f'Mantello Sferico (⟨P_J⟩ = {np.mean(pj_mantle):.2f} µW)')
    ax4.plot(t_ms, pj_r1, color='#58a6ff', linewidth=1.2, linestyle=':', label='Bobine Rotore 1')
    ax4.plot(t_ms, pj_r2, color='#f0883e', linewidth=1.2, linestyle=':', label='Bobine Rotore 2')

    ax4.set_xlabel('Tempo [ms]', color='#c9d1d9', fontsize=10)
    ax4.set_ylabel('Perdite Joule [µW]', color='#c9d1d9', fontsize=10)
    p_rad = sampling_results["Far-Field (R=15.0 cm)"]["poynting_rad_mW"]
    ax4.set_title(f"(d) Perdite Joule e Irradiazione (P_rad a 15 cm = {p_rad:.3f} mW)", color='#f0f6fc', fontsize=11, fontweight='bold')
    ax4.grid(True, color='#30363d', alpha=0.5)
    ax4.tick_params(colors='#8b949e')
    ax4.legend(facecolor='#0d1117', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=9)

    plt.tight_layout()
    fig17_path = FIGURES_DIR / "fig_17_matrice_forze_ortogonali.png"
    plt.savefig(fig17_path, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    shutil.copy(fig17_path, ROOT_FIGURES_DIR / "fig_17_matrice_forze_ortogonali.png")
    plt.close()
    print(f"  [FIGURE] Salvata {fig17_path}")


def main():
    print("=" * 80)
    print("STUDIO ELETTRODINAMICO 3D: GABBIA SFERICA METAMATERIALE A DOPPIO ROTORE ORTOGONALE")
    print("=" * 80)

    # 1. Esecuzione o verifica simulazione FEM
    vtus = run_simulation()

    # 2. Caricamento mesh e filtri geometrici
    print("\n[ANALISI] Caricamento struttura geometrica e volumi...")
    mesh_data = load_mesh_structure(vtus[0])

    # 3. Estrazione serie temporali forze vettoriali e perdite Joule
    print("\n[ANALISI] Estrazione forze vettoriali di Lorentz (Fx, Fy, Fz) e perdite Joule su 64 timestep...")
    forces, joule = extract_forces_and_losses(vtus, mesh_data)

    # 4. Campionamento sferico di Fibonacci ed MST
    print(f"\n[ANALISI] Campionamento sferico di Fibonacci (N={N_FIBONACCI} punti) su 3 sfere...")
    sampling_results = compute_spherical_sampling_and_mst(vtus, sample_idx=32)

    # Stampa sintesi dei risultati
    print("\n" + "=" * 80)
    print("RISULTATI CAMPIONAMENTO SFERICO ED ELETTRODINAMICA MULTI-ASSE:")
    print("=" * 80)
    for label, data in sampling_results.items():
        print(f"  * {label}:")
        print(f"    - Modulo B medio: {data['b_mean_uT']:.2f} uT (Picco: {data['b_max_uT']:.2f} uT)")
        print(f"    - Campo E medio : {data['e_mean_mV_m']:.2f} mV/m")
        print(f"    - Gauss Residual: {data['gauss_residual_pct']:.3f}% ({'PASS' if data['gauss_pass'] else 'FAIL'})")
        print(f"    - Poynting Rad  : {data['poynting_rad_mW']:.3f} mW")
        print(f"    - MST Forces    : Fx = {data['mst_fx_uN']:+.2f} uN, Fy = {data['mst_fy_uN']:+.2f} uN, Fz = {data['mst_fz_uN']:+.2f} uN")

    fx_mean = np.mean(forces["assembly_total"]["fx"]) * 1e6
    fy_mean = np.mean(forces["assembly_total"]["fy"]) * 1e6
    fz_mean = np.mean(forces["assembly_total"]["fz"]) * 1e6
    pj_mean = np.mean(joule["assembly_total"]) * 1e6

    print("\nBILANCIO MEDIO TOTALE MACCHINA:")
    print(f"  - Forza Lorentz [Fx]: {fx_mean:+.4f} uN (Picco: {np.max(np.abs(forces['assembly_total']['fx']))*1e6:.2f} uN)")
    print(f"  - Forza Lorentz [Fy]: {fy_mean:+.4f} uN (Picco: {np.max(np.abs(forces['assembly_total']['fy']))*1e6:.2f} uN)")
    print(f"  - Forza Lorentz [Fz]: {fz_mean:+.4f} uN (Picco: {np.max(np.abs(forces['assembly_total']['fz']))*1e6:.2f} uN)")
    print(f"  - Dissipazione Joule [P_J]: {pj_mean:.4f} uW (Mantello: {np.mean(joule['mantle'])*1e6:.4f} uW)")

    # 5. Generazione grafici
    print("\n[GRAFICA] Generazione deliverable 300 DPI...")
    plot_diagnostics(forces, joule, sampling_results, mesh_data)

    # 6. Esportazione JSON strutturato
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    json_path = DATA_DIR / "risultati_gabbia_sferica_ortogonale.json"
    
    clean_sampling = {}
    for k, v in sampling_results.items():
        clean_sampling[k] = {
            "radius_m": v["radius_m"],
            "gauss_net_weber": v["gauss_net_weber"],
            "gauss_abs_weber": v["gauss_abs_weber"],
            "gauss_residual_pct": v["gauss_residual_pct"],
            "gauss_pass": v["gauss_pass"],
            "b_mean_uT": v["b_mean_uT"],
            "b_max_uT": v["b_max_uT"],
            "e_mean_mV_m": v["e_mean_mV_m"],
            "poynting_rad_mW": v["poynting_rad_mW"],
            "mst_fx_uN": v["mst_fx_uN"],
            "mst_fy_uN": v["mst_fy_uN"],
            "mst_fz_uN": v["mst_fz_uN"],
        }

    dataset = {
        "variant": "gabbia_sferica_doppio_rotore_90deg",
        "description": "Gabbia Sferica Metamateriale a Triplo Strato X (+30/0/-30) con Doppio Rotore Ortogonale a 90 gradi",
        "parameters": {
            "r_ext_m": 0.050,
            "r_int_m": 0.047,
            "thickness_m": 0.003,
            "r_coil_m": 0.035,
            "r_core_m": 0.012,
            "mu_r_mantle": 1000.0,
            "mu_r_core": 1.0,
            "frequency_hz": F_HZ,
            "dt_s": DT,
            "timesteps": TIMESTEPS,
            "t_cycle_s": 0.00687947,
            "n_fibonacci": N_FIBONACCI
        },
        "lorentz_forces_mean_uN": {
            "assembly_total": {"fx": fx_mean, "fy": fy_mean, "fz": fz_mean},
            "mantle": {
                "fx": float(np.mean(forces["mantle"]["fx"]) * 1e6),
                "fy": float(np.mean(forces["mantle"]["fy"]) * 1e6),
                "fz": float(np.mean(forces["mantle"]["fz"]) * 1e6),
            },
            "layer1_plus30": {
                "fx": float(np.mean(forces["layer1_plus30"]["fx"]) * 1e6),
                "fy": float(np.mean(forces["layer1_plus30"]["fy"]) * 1e6),
                "fz": float(np.mean(forces["layer1_plus30"]["fz"]) * 1e6),
            },
            "layer2_ortho": {
                "fx": float(np.mean(forces["layer2_ortho"]["fx"]) * 1e6),
                "fy": float(np.mean(forces["layer2_ortho"]["fy"]) * 1e6),
                "fz": float(np.mean(forces["layer2_ortho"]["fz"]) * 1e6),
            },
            "layer3_minus30": {
                "fx": float(np.mean(forces["layer3_minus30"]["fx"]) * 1e6),
                "fy": float(np.mean(forces["layer3_minus30"]["fy"]) * 1e6),
                "fz": float(np.mean(forces["layer3_minus30"]["fz"]) * 1e6),
            },
            "rotor1_coils": {
                "fx": float(np.mean(forces["rotor1_coils"]["fx"]) * 1e6),
                "fy": float(np.mean(forces["rotor1_coils"]["fy"]) * 1e6),
                "fz": float(np.mean(forces["rotor1_coils"]["fz"]) * 1e6),
            },
            "rotor2_coils": {
                "fx": float(np.mean(forces["rotor2_coils"]["fx"]) * 1e6),
                "fy": float(np.mean(forces["rotor2_coils"]["fy"]) * 1e6),
                "fz": float(np.mean(forces["rotor2_coils"]["fz"]) * 1e6),
            },
            "peek_core": {
                "fx": float(np.mean(forces["peek_core"]["fx"]) * 1e6),
                "fy": float(np.mean(forces["peek_core"]["fy"]) * 1e6),
                "fz": float(np.mean(forces["peek_core"]["fz"]) * 1e6),
            }
        },
        "joule_losses_mean_uW": {
            "assembly_total": pj_mean,
            "mantle": float(np.mean(joule["mantle"]) * 1e6),
            "layer1_plus30": float(np.mean(joule["layer1_plus30"]) * 1e6),
            "layer2_ortho": float(np.mean(joule["layer2_ortho"]) * 1e6),
            "layer3_minus30": float(np.mean(joule["layer3_minus30"]) * 1e6),
            "rotor1_coils": float(np.mean(joule["rotor1_coils"]) * 1e6),
            "rotor2_coils": float(np.mean(joule["rotor2_coils"]) * 1e6),
            "peek_core": float(np.mean(joule["peek_core"]) * 1e6)
        },
        "spherical_sampling": clean_sampling
    }

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(dataset, f, indent=2)
    print(f"  [DATA] Dataset salvato in {json_path}")
    print("Studio completato con successo.")


if __name__ == "__main__":
    main()
