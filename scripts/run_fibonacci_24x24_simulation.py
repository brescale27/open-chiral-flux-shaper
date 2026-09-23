#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico e Termico su Elmer FEM 3D:
Gabbia Sferica Ibrida con Architettura a 24 Gruppi (Mappa di Fibonacci & Radice Numerica).
Potenza Calibrata: 100 W per Solenoide (Potenza Totale Array Solenoidi: 2.40 kW).

Caratteristiche della simulazione:
- 24 gruppi di solenoidi distribuiti sull'equatore a Delta_theta = 15°
- Pilotaggio secondo la radice numerica dei primi 24 numeri di Fibonacci (periodo Pisano mod 9):
    v_k = [1, 1, 2, 3, 5, 8, 4, 3, 7, 1, 8, 9, 8, 8, 7, 6, 4, 1, 5, 6, 2, 8, 1, 9]
    phi_k = (v_k / 9) * 2*pi
  con auto-bilanciamento per coniugazione di fase diametrale phi_{k+12} = -phi_k
- Mantello metamateriale sferico a triplo strato X ferromagnetico (mu_r = 1000.0, spessore 3 mm)
- Nucleo amagnetico sferico centrale in PEEK (R_core = 12 mm)
- Risoluzione transiente su 64 timestep (f = 100 Hz, dt = 0.25 ms, T_tot = 16.0 ms)
- Analisi forze di Lorentz F_x, F_y, F_z, bilancio energetico (100 W/bobina, 2.40 kW array), saturazione B_max, Gauss e MST
- Generazione Tavola Diagnostica 300 DPI (Figura 27)
- Generazione Video Animato ad Alta Risoluzione della dinamica vettoriale di forze e campi

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
from PIL import Image

# Percorsi principali
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
ROOT_FIGURES_DIR = ROOT_DIR / "figures"
ROOT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

VAR_DIR = ROOT_DIR / "variants" / "gabbia_sferica_fibonacci_24x24"
CONFIG_DIR = VAR_DIR / "config"
BASE_SIF = CONFIG_DIR / "case_fibonacci_24x24_100w.sif"
MESH_DIR = VAR_DIR / "mesh"
MESH_NAME = "macchina_fibonacci_24x24"
WORK_DIR = VAR_DIR / "work_dirs" / "run_fibonacci_24x24"
RES_DIR = WORK_DIR / "results"
DATA_DIR = VAR_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON = DATA_DIR / "fibonacci_24x24_100w.json"

FIG_VAR_DIR = VAR_DIR / "figures"
FIG_VAR_DIR.mkdir(parents=True, exist_ok=True)
FIG_27_NAME = "fig_27_architettura_fibonacci_24x24_100w.png"
VIDEO_NAME = "video_dinamica_fibonacci_24x24.gif"

ELMER_SOLVER = r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"
MU0 = 4.0 * np.pi * 1e-7
SIGMA_SB = 5.670374419e-8  # Stefan-Boltzmann W / (m^2 K^4)
EMISSIVITY_MANTLE = 0.85
R_EXT_SPHERE = 0.050  # 50 mm
A_RAD_SPHERE = 4.0 * np.pi * (R_EXT_SPHERE**2)  # 0.031416 m^2

TIMESTEPS = 64
DT = 0.00025  # 0.25 ms
F_HZ = 100.0
T_PERIOD = 0.010  # 10 ms (1 ciclo)
N_FIBONACCI = 2500
RADII = [0.065, 0.100, 0.150]
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]

# Sequenza Pisano periodo 24 mod 9
FIB_NUMS = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368]
PISANO_ROOTS = [1, 1, 2, 3, 5, 8, 4, 3, 7, 1, 8, 9, 8, 8, 7, 6, 4, 1, 5, 6, 2, 8, 1, 9]
PHASES_RAD = [(r / 9.0) * 2.0 * np.pi for r in PISANO_ROOTS]
PHASES_DEG = [np.degrees(p) % 360.0 for p in PHASES_RAD]


def prepare_and_run_solver():
    """Configura l'ambiente ed esegue ElmerSolver se necessario."""
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)

    vtus = sorted(list(RES_DIR.glob("macchina_fibonacci_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} VTU in {RES_DIR}. Skip solver.")
        return vtus[:TIMESTEPS]

    print("=" * 80)
    print("  [SOLVER] Esecuzione Elmer FEM 3D: Gabbia Sferica Fibonacci 24x24 (100 W/bobina)")
    print(f"  [DIR] Work Dir: {WORK_DIR}")
    print("=" * 80)

    for f in RES_DIR.glob("macchina_fibonacci_out_t*.vtu"):
        f.unlink()

    sif_text = BASE_SIF.read_text(encoding="utf-8")
    mesh_rel = os.path.relpath(str(MESH_DIR), str(WORK_DIR)).replace('\\', '/')
    sif_lines = sif_text.splitlines()
    new_lines = []
    for line in sif_lines:
        if 'Mesh DB' in line:
            new_lines.append(f'  Mesh DB "{mesh_rel}" "{MESH_NAME}"')
        elif 'Results Directory' in line:
            new_lines.append('  Results Directory "results"')
        else:
            new_lines.append(line)

    dest_sif = WORK_DIR / "case.sif"
    dest_sif.write_text("\n".join(new_lines), encoding="utf-8")
    (WORK_DIR / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")

    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    proc = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(WORK_DIR), env=env, capture_output=True, text=True)
    dur = time.time() - t0

    if proc.returncode != 0:
        print(f"  [ERRORE] ElmerSolver terminato con errore {proc.returncode}:")
        print("STDERR tail:\n", proc.stderr[-1000:])
        print("STDOUT tail:\n", proc.stdout[-1000:])
        sys.exit(1)

    vtus = sorted(list(RES_DIR.glob("macchina_fibonacci_out_t*.vtu")))
    print(f"  [SOLVER] Concluso con successo in {dur:.1f}s ({len(vtus)} VTU generati).")
    if len(vtus) < TIMESTEPS:
        print(f"  [ERRORE] Generati {len(vtus)} VTU su {TIMESTEPS} attesi.")
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
    """Calcola la divergenza di Gauss e il Tensore di Maxwell MST su sfere concentriche."""
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
            "mst_fx_uN": fx_mst * 1e6,
            "mst_fy_uN": fy_mst * 1e6,
            "mst_fz_uN": fz_mst * 1e6,
            "mst_mag_uN": f_mst_mag * 1e6
        }
    return results


def process_simulation_data(vtus):
    """Elabora le grandezze fisiche dai 64 timestep VTU."""
    print("  [ANALYSIS] Lettura geometria e calcolo volumi tetraedrici con GeometryIds...")
    m0 = meshio.read(str(vtus[0]))
    geom_ids = m0.cell_data['GeometryIds'][0]
    pts = m0.points
    cells = m0.cells[0].data

    v0, v1, v2, v3 = pts[cells[:, 0]], pts[cells[:, 1]], pts[cells[:, 2]], pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0
    elem_com = (v0 + v1 + v2 + v3) / 4.0
    x, y, z = elem_com[:, 0], elem_com[:, 1], elem_com[:, 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y, x)
    theta = np.where(theta < 0.0, theta + 2.0 * np.pi, theta)

    # Maschere fisiche basate direttamente su GeometryIds
    mask_mantle = (geom_ids == 1)
    mask_l1 = mask_mantle & (r < 0.0480)
    mask_l2 = mask_mantle & (r >= 0.0480) & (r < 0.0490)
    mask_l3 = mask_mantle & (r >= 0.0490)
    mask_coils = (geom_ids == 2)
    mask_core = (geom_ids == 3)
    mask_assembly = mask_mantle | mask_coils | mask_core

    # Maschere per i 24 singoli settori di bobina
    d_th = 2.0 * np.pi / 24.0
    coil_masks = []
    for k in range(24):
        ang_k = k * d_th
        d_ang = np.abs((theta - ang_k + np.pi) % (2.0 * np.pi) - np.pi)
        mask_k = mask_coils & (d_ang <= (d_th / 2.0))
        coil_masks.append(mask_k)

    # Nodi per induzione sul mantello
    pts_r = np.sqrt(pts[:, 0]**2 + pts[:, 1]**2 + pts[:, 2]**2)
    pts_mantle = (pts_r >= 0.0465) & (pts_r <= 0.0505)

    masks = {
        "mantle_total": mask_mantle,
        "layer1_plus30": mask_l1,
        "layer2_ortho": mask_l2,
        "layer3_minus30": mask_l3,
        "coils_total": mask_coils,
        "peek_core": mask_core,
        "total": mask_assembly
    }

    forces = {k: {"fx": [], "fy": [], "fz": []} for k in masks}
    joule_eddy = {k: [] for k in masks}

    b_mantle_max_series = []
    b_mantle_mean_series = []
    b_global_max_series = []

    print(f"  [ANALYSIS] Elaborazione di {len(vtus)} timestep VTU...")
    for idx, vtu_file in enumerate(vtus):
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
            joule_eddy[k].append(pj_val)

        if (idx + 1) % 16 == 0 or idx == len(vtus) - 1:
            f_tot_uN = np.array([forces['total']['fx'][-1], forces['total']['fy'][-1], forces['total']['fz'][-1]]) * 1e6
            print(f"    - Timestep {idx+1:2d}/{len(vtus)}: F_tot = ({f_tot_uN[0]:.2f}, {f_tot_uN[1]:.2f}, {f_tot_uN[2]:.2f}) uN, P_eddy_mantle = {joule_eddy['mantle_total'][-1]:.3f} W")

    # Medie integrali nel tempo (ciclo stazionario)
    half_idx = TIMESTEPS // 2
    time_s = [i * DT for i in range(TIMESTEPS)]

    total_fx = np.array(forces["total"]["fx"])
    total_fy = np.array(forces["total"]["fy"])
    total_fz = np.array(forces["total"]["fz"])
    total_fmag = np.sqrt(total_fx**2 + total_fy**2 + total_fz**2)

    mean_fx = float(np.mean(total_fx[half_idx:]))
    mean_fy = float(np.mean(total_fy[half_idx:]))
    mean_fz = float(np.mean(total_fz[half_idx:]))
    mean_fmag = float(np.sqrt(mean_fx**2 + mean_fy**2 + mean_fz**2))
    peak_fmag = float(np.max(total_fmag))

    # Potenza attiva imposta calibrata: 100 W per solenoide (totale array = 2.40 kW)
    p_elec_per_coil = 100.0  # W
    p_elec_array_tot = 24.0 * p_elec_per_coil  # 2400.0 W = 2.40 kW
    per_coil_powers = [p_elec_per_coil for _ in range(24)]

    # Perdite eddy dal risolutore
    mean_eddy_mantle = float(np.mean(joule_eddy["mantle_total"][half_idx:]))
    mean_eddy_l1 = float(np.mean(joule_eddy["layer1_plus30"][half_idx:]))
    mean_eddy_l2 = float(np.mean(joule_eddy["layer2_ortho"][half_idx:]))
    mean_eddy_l3 = float(np.mean(joule_eddy["layer3_minus30"][half_idx:]))
    mean_eddy_core = float(np.mean(joule_eddy["peek_core"][half_idx:]))

    total_system_power_W = p_elec_array_tot + mean_eddy_mantle

    # Efficienza di spinta microN / W
    thrust_eff_uN_per_W = float(mean_fmag * 1e6 / total_system_power_W)

    # Margine di saturazione
    peak_b_mantle = float(np.max(b_mantle_max_series))
    mean_b_mantle = float(np.mean(b_mantle_mean_series[half_idx:]))
    peak_b_global = float(np.max(b_global_max_series))
    b_sat_ref = 1.5  # Tesla
    sat_margin_pct = float((1.0 - peak_b_mantle / b_sat_ref) * 100.0)

    # Analisi termica di radiazione
    t_eq_k = float((total_system_power_W / (EMISSIVITY_MANTLE * SIGMA_SB * A_RAD_SPHERE))**0.25)
    t_eq_c = float(t_eq_k - 273.15)
    target_t_k = 350.0  # 76.85 °C
    area_req_350k = float(total_system_power_W / (EMISSIVITY_MANTLE * SIGMA_SB * (target_t_k**4 - 3.0**4)))
    area_aux_350k = float(max(0.0, area_req_350k - A_RAD_SPHERE))

    # Sfera di Fibonacci per conservazione Gauss & MST
    print("  [ANALYSIS] Calcolo sfere di Fibonacci (Gauss div(B)=0 e MST)...")
    fib_res = analyze_fibonacci_spheres(vtus, sample_idx=32)

    data_summary = {
        "architecture": {
            "name": "Gabbia Sferica Ibrida a 24 Gruppi (Mappa di Fibonacci)",
            "n_groups": 24,
            "delta_theta_deg": 15.0,
            "hierarchical_nodes": 576,
            "pisano_period": 24,
            "pisano_mod": 9,
            "roots": PISANO_ROOTS,
            "phases_deg": PHASES_DEG,
            "calibrated_power_per_coil_W": p_elec_per_coil,
            "target_total_power_coils_kW": p_elec_array_tot / 1000.0
        },
        "simulation_parameters": {
            "frequency_hz": F_HZ,
            "timesteps": TIMESTEPS,
            "timestep_size_s": DT,
            "total_time_s": TIMESTEPS * DT,
            "mantle_mu_r": 1000.0,
            "j0_nominal_A_m2": 2.4405e5
        },
        "thrust_and_forces": {
            "mean_fx_N": mean_fx,
            "mean_fy_N": mean_fy,
            "mean_fz_N": mean_fz,
            "resultant_mag_N": mean_fmag,
            "peak_instantaneous_N": peak_fmag,
            "mean_fx_uN": mean_fx * 1e6,
            "mean_fy_uN": mean_fy * 1e6,
            "mean_fz_uN": mean_fz * 1e6,
            "resultant_mag_uN": mean_fmag * 1e6,
            "peak_instantaneous_uN": peak_fmag * 1e6,
            "time_series_fx_N": forces["total"]["fx"],
            "time_series_fy_N": forces["total"]["fy"],
            "time_series_fz_N": forces["total"]["fz"],
            "time_series_fmag_N": total_fmag.tolist()
        },
        "power_and_efficiency": {
            "total_system_power_W": total_system_power_W,
            "coils_input_power_kW": p_elec_array_tot / 1000.0,
            "average_per_coil_W": p_elec_per_coil,
            "per_coil_powers_W": per_coil_powers,
            "mantle_eddy_total_W": mean_eddy_mantle,
            "mantle_l1_plus30_W": mean_eddy_l1,
            "mantle_l2_ortho_W": mean_eddy_l2,
            "mantle_l3_minus30_W": mean_eddy_l3,
            "peek_core_eddy_W": mean_eddy_core,
            "thrust_efficiency_uN_per_W": thrust_eff_uN_per_W
        },
        "magnetic_saturation": {
            "b_mantle_peak_T": peak_b_mantle,
            "b_mantle_mean_T": mean_b_mantle,
            "b_global_peak_T": peak_b_global,
            "b_sat_ref_T": b_sat_ref,
            "saturation_margin_pct": sat_margin_pct,
            "linear_regime_pass": bool(peak_b_mantle < b_sat_ref)
        },
        "thermal_vacuum_equilibrium": {
            "t_eq_kelvin": t_eq_k,
            "t_eq_celsius": t_eq_c,
            "mantle_radiating_area_m2": A_RAD_SPHERE,
            "mantle_emissivity": EMISSIVITY_MANTLE,
            "area_req_350k_m2": area_req_350k,
            "area_aux_350k_m2": area_aux_350k
        },
        "fibonacci_spheres": fib_res
    }

    # Salvataggio JSON
    OUT_JSON.write_text(json.dumps(data_summary, indent=2), encoding="utf-8")
    print(f"  [DATA] File JSON salvato in: {OUT_JSON}")

    return data_summary, time_s, forces, joule_eddy, per_coil_powers, (b_mantle_max_series, b_global_max_series)


def generate_figure_27(data, time_s, forces, joule_eddy, per_coil_powers):
    """Genera la Tavola Diagnostica a 300 DPI per l'Architettura Fibonacci 24x24 (Figura 27)."""
    print("=" * 80)
    print("  [PLOTTING] Generazione Figura 27: Diagnostica Architettura Fibonacci 24x24 (300 DPI)")
    print("=" * 80)

    plt.style.use('default')
    fig = plt.figure(figsize=(18, 14), dpi=300)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.26)

    # ----------------------------------------------------
    # PANNELLO A: Mappa dei 24 Settori di Fibonacci e Coniugazione di Fase
    # ----------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0], projection='polar')

    angles = np.linspace(0, 2.0 * np.pi, 24, endpoint=False)
    roots = data["architecture"]["roots"]

    norm_roots = np.array(roots) / 9.0
    colors = plt.cm.plasma(norm_roots)

    bars = ax_a.bar(angles, np.array(roots), width=2.0*np.pi/24.0*0.88, bottom=0.0,
                    color=colors, edgecolor='black', linewidth=1.2, alpha=0.85)

    for k in range(24):
        ang = angles[k]
        r_val = roots[k]
        ax_a.text(ang, r_val + 0.65, f"B{k+1}\n({r_val})", ha='center', va='center',
                  fontsize=8, fontweight='bold', color='#1a1a2e')

    for k in range(12):
        opp_k = k + 12
        ax_a.plot([angles[k], angles[opp_k]], [roots[k], roots[opp_k]],
                  color='#2b580c', linestyle=':', linewidth=0.8, alpha=0.6)

    ax_a.set_title(r"PANEL A: Mappa Settoriale Pisano mod 9 (24 Gruppi, $\Delta\theta = 15^\circ$)" + "\n" +
                   r"Coniugazione di Fase Diametrale: $\phi_{k+12} = -\phi_k$ (Bilanciamento Reattivo)",
                   fontsize=11, fontweight='bold', pad=15)
    ax_a.set_rticks([3, 6, 9])
    ax_a.set_rlabel_position(45)
    ax_a.set_yticklabels(['3 (60°)', '6 (120°)', '9 (0°/360°)'], fontsize=8)
    ax_a.grid(True, linestyle='--', alpha=0.5)

    # ----------------------------------------------------
    # PANNELLO B: Dinamica delle Forze di Lorentz (Waveforms & Odografo 3D)
    # ----------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])

    t_ms = np.array(time_s) * 1000.0
    fx_uN = np.array(forces["total"]["fx"]) * 1e6
    fy_uN = np.array(forces["total"]["fy"]) * 1e6
    fz_uN = np.array(forces["total"]["fz"]) * 1e6
    fmag_uN = np.sqrt(fx_uN**2 + fy_uN**2 + fz_uN**2)

    ax_b.plot(t_ms, fx_uN, color='#d90429', linewidth=2.0, label=r"$F_x(t)$ (Trasversale X)")
    ax_b.plot(t_ms, fy_uN, color='#2b9348', linewidth=2.0, label=r"$F_y(t)$ (Laterale Y)")
    ax_b.plot(t_ms, fz_uN, color='#0077b6', linewidth=2.0, label=r"$F_z(t)$ (Assiale Z)")
    ax_b.plot(t_ms, fmag_uN, color='#2b2d42', linewidth=2.4, linestyle='--', label=r"$|\mathbf{F}(t)|$ Risultante")

    mean_f_uN = data["thrust_and_forces"]["resultant_mag_uN"]
    peak_f_uN = data["thrust_and_forces"]["peak_instantaneous_uN"]
    ax_b.axhline(mean_f_uN, color='#6c757d', linestyle=':', linewidth=1.5,
                 label=f"Media Stazionaria = {mean_f_uN:.2f} uN")

    ax_b.set_title(f"PANEL B: Waveform delle Forze di Lorentz Transienti\nSpinta Risultante Media: |F| = {mean_f_uN:.2f} uN (Picco: {peak_f_uN:.2f} uN)",
                   fontsize=11, fontweight='bold')
    ax_b.set_xlabel("Tempo [ms]", fontsize=10, fontweight='bold')
    ax_b.set_ylabel(r"Forza Elettrodinamica [$\mu$N]", fontsize=10, fontweight='bold')
    ax_b.grid(True, linestyle='--', alpha=0.6)
    ax_b.legend(loc='upper right', fontsize=8.5, framealpha=0.95)
    ax_b.set_xlim(0, 16.0)

    # ----------------------------------------------------
    # PANNELLO C: Dissipazione Joule Calibrata (100 W/Bobina) e Ripartizione Mantello
    # ----------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])

    coil_indices = np.arange(1, 25)
    tot_coils_kw = data["power_and_efficiency"]["coils_input_power_kW"]

    bars_c = ax_c.bar(coil_indices, per_coil_powers, color='#028090', edgecolor='#004e64',
                      linewidth=1.0, alpha=0.85, label='Potenza Calibrata per Bobina')
    ax_c.axhline(100.0, color='#e63946', linestyle='--', linewidth=2.0,
                 label='Target Nominale: 100.0 W/bobina (Totale = 2.40 kW)')

    ax_c.set_title(f"PANEL C: Bilancio di Potenza per Solenoide (N=24) e Dissipazione Mantello\nPotenza Totale Solenoidi: P_array = {tot_coils_kw:.2f} kW (Calibrazione: 100.0 W/bobina)",
                   fontsize=11, fontweight='bold')
    ax_c.set_xlabel("Indice Bobina Fibonacci [1 - 24]", fontsize=10, fontweight='bold')
    ax_c.set_ylabel("Potenza Elettrica [W]", fontsize=10, fontweight='bold')
    ax_c.set_xticks(range(1, 25))
    ax_c.set_ylim(0, 130)
    ax_c.grid(True, linestyle='--', alpha=0.6)

    # Tabella inset con breakdown strati mantello
    pj_m_tot = data["power_and_efficiency"]["mantle_eddy_total_W"]
    pj_l1 = data["power_and_efficiency"]["mantle_l1_plus30_W"]
    pj_l2 = data["power_and_efficiency"]["mantle_l2_ortho_W"]
    pj_l3 = data["power_and_efficiency"]["mantle_l3_minus30_W"]

    table_data = [
        ["Array 24 Bobine", f"{tot_coils_kw:.2f} kW"],
        ["Mantello Totale", f"{pj_m_tot:.3f} W"],
        ["Strato 1 (+30°)", f"{pj_l1:.3f} W"],
        ["Strato 2 (0°)", f"{pj_l2:.3f} W"],
        ["Strato 3 (-30°)", f"{pj_l3:.4f} W"],
        ["Nucleo PEEK", f"{data['power_and_efficiency']['peek_core_eddy_W']:.6f} W"]
    ]
    tb = ax_c.table(cellText=table_data, colLabels=["Componente", "Potenza / Eddy"],
                    loc='upper left', bbox=[0.05, 0.45, 0.42, 0.45])
    tb.auto_set_font_size(False)
    tb.set_fontsize(8)
    ax_c.legend(loc='lower right', fontsize=8.5, framealpha=0.95)

    # ----------------------------------------------------
    # PANNELLO D: Verifica Gauss div(B)=0, Tensore MST e Margine Saturazione
    # ----------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])

    fib_spheres = data["fibonacci_spheres"]
    r_labels = list(fib_spheres.keys())
    residuals = [fib_spheres[k]["gauss_residual_pct"] for k in r_labels]
    b_peaks = [fib_spheres[k]["b_max_uT"] / 1000.0 for k in r_labels]  # mT

    x_pos = np.arange(len(r_labels))
    width = 0.35

    ax_d1 = ax_d
    ax_d2 = ax_d.twinx()

    b1 = ax_d1.bar(x_pos - width/2, residuals, width, label='Residuo Flusso Gauss [%]',
                   color='#3a0ca3', alpha=0.85, edgecolor='black')
    b2 = ax_d2.bar(x_pos + width/2, b_peaks, width, label=r'Picco Induzione $|\mathbf{B}|$ [mT]',
                   color='#f72585', alpha=0.85, edgecolor='black')

    ax_d1.axhline(2.0, color='#d90429', linestyle=':', linewidth=1.8, label='Soglia Solenoidalita (< 2%)')

    ax_d1.set_title(r"PANEL D: Solenoidalita $\nabla \cdot \mathbf{B} = 0$ e Induzione alle Sfere di Fibonacci" + "\n" +
                    f"Margine Saturazione Mantello: {data['magnetic_saturation']['saturation_margin_pct']:.1f}% (B_max = {data['magnetic_saturation']['b_mantle_peak_T']*1000.0:.2f} mT < 1500 mT)",
                    fontsize=11, fontweight='bold')
    ax_d1.set_xticks(x_pos)
    ax_d1.set_xticklabels(["Near-Field\n(R = 6.5 cm)", "Mid-Field\n(R = 10.0 cm)", "Far-Field\n(R = 15.0 cm)"],
                          fontsize=9, fontweight='bold')
    ax_d1.set_ylabel("Residuo di Gauss [%]", color='#3a0ca3', fontsize=10, fontweight='bold')
    ax_d2.set_ylabel("Induzione Massima [mT]", color='#f72585', fontsize=10, fontweight='bold')
    ax_d1.set_ylim(0, max(residuals) * 1.35)
    ax_d2.set_ylim(0, max(b_peaks) * 1.35)
    ax_d1.grid(True, linestyle='--', alpha=0.5)

    lines1, labels1 = ax_d1.get_legend_handles_labels()
    lines2, labels2 = ax_d2.get_legend_handles_labels()
    ax_d1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=8.5, framealpha=0.95)

    plt.suptitle("ARCHITETTURA FIBONACCI 24x24 A RADICE NUMERICA CON MANTELLO METAMATERIALE A TRIPLO STRATO X\nCampagna Elmer FEM 3D Transiente (64 Timestep, f = 100 Hz, 100 W/Bobina, Potenza Totale Array = 2.40 kW)",
                 fontsize=14, fontweight='bold', y=0.98)

    out_root = ROOT_FIGURES_DIR / FIG_27_NAME
    out_var = FIG_VAR_DIR / FIG_27_NAME
    fig.savefig(out_root, dpi=300, bbox_inches='tight')
    fig.savefig(out_var, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  [SAVED] Figura 27 salvata in:")
    print(f"    - {out_root}")
    print(f"    - {out_var}")


def generate_animated_video(vtus, data, time_s, forces):
    """
    Genera un video animato (GIF ad alta risoluzione) che evidenzia:
    1. Campo magnetico 2D |B|(x,y) nel piano equatoriale con vortice chirale e le 24 bobine.
    2. Odografo 3D dinamico della forza di Lorentz con vettore istantaneo e traiettoria orbitale.
    3. Waveform dinamica temporale con cursore sincronizzato.
    """
    print("=" * 80)
    print("  [VIDEO] Generazione Video Dinamico 3D di Forze e Campi Elettrodinamici")
    print("=" * 80)

    grid_res = 120
    lim = 0.065  # 65 mm (cattura mantello e bobine)
    xi = np.linspace(-lim, lim, grid_res)
    yi = np.linspace(-lim, lim, grid_res)
    XI, YI = np.meshgrid(xi, yi)

    m0 = meshio.read(str(vtus[0]))
    pts = m0.points
    tri = Delaunay(pts)

    t_ms = np.array(time_s) * 1000.0
    fx_uN = np.array(forces["total"]["fx"]) * 1e6
    fy_uN = np.array(forces["total"]["fy"]) * 1e6
    fz_uN = np.array(forces["total"]["fz"]) * 1e6
    fmag_uN = np.sqrt(fx_uN**2 + fy_uN**2 + fz_uN**2)

    temp_frame_dir = WORK_DIR / "video_frames"
    temp_frame_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []

    rc = 0.035
    r_mantle_ext = 0.050
    r_mantle_int = 0.047
    r_core = 0.012
    coil_angles = np.linspace(0, 2.0 * np.pi, 24, endpoint=False)
    coil_xs = rc * np.cos(coil_angles)
    coil_ys = rc * np.sin(coil_angles)

    print(f"  [VIDEO] Rendering di {len(vtus)} frame sincronizzati...")

    max_f = max(np.max(np.abs(fx_uN)), np.max(np.abs(fy_uN)), np.max(np.abs(fz_uN))) * 1.25

    for idx, vtu_file in enumerate(vtus):
        m = meshio.read(str(vtu_file))
        b_pts = m.point_data['magnetic flux density']
        interp_b = LinearNDInterpolator(tri, b_pts, fill_value=0.0)

        eval_pts = np.column_stack([XI.ravel(), YI.ravel(), np.zeros_like(XI.ravel())])
        b_grid = interp_b(eval_pts)
        b_grid_x = b_grid[:, 0].reshape(XI.shape)
        b_grid_y = b_grid[:, 1].reshape(XI.shape)
        b_mag_grid = np.linalg.norm(b_grid, axis=1).reshape(XI.shape) * 1000.0  # in mT

        fig = plt.figure(figsize=(15, 6.5), dpi=140)
        gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1.25, 1.0], hspace=0.35, wspace=0.25)

        # Subplot 1: Campo magnetico B(x,y) equatoriale con vettori e 24 bobine
        ax_field = fig.add_subplot(gs[:, 0])
        cf = ax_field.contourf(XI * 1000.0, YI * 1000.0, b_mag_grid, levels=40, cmap='inferno', vmin=0, vmax=5.0)
        cbar = fig.colorbar(cf, ax=ax_field, fraction=0.046, pad=0.04)
        cbar.set_label(r"$|\mathbf{B}|$ [mT]", fontsize=9, fontweight='bold')

        skip = 5
        ax_field.quiver(XI[::skip, ::skip] * 1000.0, YI[::skip, ::skip] * 1000.0,
                        b_grid_x[::skip, ::skip], b_grid_y[::skip, ::skip],
                        color='white', alpha=0.55, scale=0.15)

        theta_c = np.linspace(0, 2.0 * np.pi, 200)
        ax_field.plot(r_mantle_ext * 1000.0 * np.cos(theta_c), r_mantle_ext * 1000.0 * np.sin(theta_c),
                      color='cyan', linewidth=1.2, linestyle='--', label='Mantello (R=50 mm)')
        ax_field.plot(r_mantle_int * 1000.0 * np.cos(theta_c), r_mantle_int * 1000.0 * np.sin(theta_c),
                      color='cyan', linewidth=0.8, linestyle=':')
        ax_field.plot(r_core * 1000.0 * np.cos(theta_c), r_core * 1000.0 * np.sin(theta_c),
                      color='gray', linewidth=1.0, label='PEEK Core (R=12 mm)')

        w_val = 2.0 * np.pi * F_HZ
        t_curr = time_s[idx]
        current_phases = [np.sin(w_val * t_curr + p) for p in PHASES_RAD]
        for ck in range(24):
            c_val = current_phases[ck]
            c_color = '#ef233c' if c_val > 0 else '#0466c8'
            ax_field.plot(coil_xs[ck] * 1000.0, coil_ys[ck] * 1000.0, 'o',
                          color=c_color, markersize=5.0 + 3.0 * abs(c_val),
                          markeredgecolor='black', markeredgewidth=0.8)

        ax_field.set_title(f"CAMPO MAGNETICO EQUATORIALE (Z=0) - t = {t_ms[idx]:.2f} ms\nArray 24 Bobine Fibonacci (Pisano mod 9, 100 W/bobina)",
                           fontsize=10, fontweight='bold')
        ax_field.set_xlabel("X [mm]", fontsize=9, fontweight='bold')
        ax_field.set_ylabel("Y [mm]", fontsize=9, fontweight='bold')
        ax_field.set_xlim(-lim * 1000.0, lim * 1000.0)
        ax_field.set_ylim(-lim * 1000.0, lim * 1000.0)
        ax_field.set_aspect('equal')

        # Subplot 2: Odografo 3D della Forza di Lorentz
        ax_3d = fig.add_subplot(gs[0, 1], projection='3d')
        ax_3d.plot(fx_uN[:idx+1], fy_uN[:idx+1], fz_uN[:idx+1], color='#3a0ca3', linewidth=1.8, alpha=0.85)
        ax_3d.quiver(0, 0, 0, fx_uN[idx], fy_uN[idx], fz_uN[idx], color='#d90429', linewidth=2.5,
                     arrow_length_ratio=0.25)
        ax_3d.scatter([fx_uN[idx]], [fy_uN[idx]], [fz_uN[idx]], color='#d90429', s=45, edgecolor='black')

        ax_3d.set_title(f"ODOGRAFO FORZA 3D: F(t) = ({fx_uN[idx]:.1f}, {fy_uN[idx]:.1f}, {fz_uN[idx]:.1f}) uN",
                        fontsize=9, fontweight='bold')
        ax_3d.set_xlabel(r"$F_x$ [$\mu$N]", fontsize=8)
        ax_3d.set_ylabel(r"$F_y$ [$\mu$N]", fontsize=8)
        ax_3d.set_zlabel(r"$F_z$ [$\mu$N]", fontsize=8)
        ax_3d.set_xlim(-max_f, max_f)
        ax_3d.set_ylim(-max_f, max_f)
        ax_3d.set_zlim(-max_f, max_f)
        ax_3d.view_init(elev=25, azim=45 + idx * (360.0 / TIMESTEPS))

        # Subplot 3: Waveform Temporale con Cursore Dinamico
        ax_time = fig.add_subplot(gs[1, 1])
        ax_time.plot(t_ms, fx_uN, color='#d90429', linewidth=1.2, alpha=0.7, label=r"$F_x$")
        ax_time.plot(t_ms, fy_uN, color='#2b9348', linewidth=1.2, alpha=0.7, label=r"$F_y$")
        ax_time.plot(t_ms, fz_uN, color='#0077b6', linewidth=1.2, alpha=0.7, label=r"$F_z$")
        ax_time.plot(t_ms, fmag_uN, color='#2b2d42', linewidth=1.5, linestyle='--', label=r"$|\mathbf{F}|$")

        ax_time.axvline(t_ms[idx], color='#ffb703', linewidth=2.0)
        ax_time.scatter([t_ms[idx]], [fmag_uN[idx]], color='#d90429', s=35, zorder=5)

        ax_time.set_title(f"WAVEFORM FORZA ELETTRODINAMICA (f = 100 Hz, |F|_inst = {fmag_uN[idx]:.1f} uN)",
                          fontsize=9, fontweight='bold')
        ax_time.set_xlabel("Tempo [ms]", fontsize=8, fontweight='bold')
        ax_time.set_ylabel(r"Forza [$\mu$N]", fontsize=8, fontweight='bold')
        ax_time.set_xlim(0, 16.0)
        ax_time.set_ylim(-max_f * 1.05, max_f * 1.05)
        ax_time.grid(True, linestyle='--', alpha=0.5)
        ax_time.legend(loc='upper right', fontsize=7, ncol=4, framealpha=0.9)

        frame_file = temp_frame_dir / f"frame_{idx:03d}.png"
        fig.savefig(frame_file, bbox_inches='tight')
        plt.close(fig)
        frame_paths.append(frame_file)

        if (idx + 1) % 16 == 0 or idx == len(vtus) - 1:
            print(f"    - Frame {idx+1:2d}/{len(vtus)} renderizzato.")

    print("  [VIDEO] Assemblaggio sequenza animata con Pillow...")
    images = [Image.open(f) for f in frame_paths]

    out_video_root = ROOT_FIGURES_DIR / VIDEO_NAME
    out_video_var = FIG_VAR_DIR / VIDEO_NAME

    images[0].save(
        out_video_root,
        save_all=True,
        append_images=images[1:],
        duration=80,
        loop=0,
        optimize=True
    )
    shutil.copy(out_video_root, out_video_var)

    print(f"  [SAVED] Video animato salvato in:")
    print(f"    - {out_video_root}")
    print(f"    - {out_video_var}")

    for f in frame_paths:
        f.unlink()
    temp_frame_dir.rmdir()
    print("  [CLEANUP] Frame temporanei rimossi con successo.")


def main():
    print("=" * 80)
    print("CAMPAIGN RUNNER: GABBIA SFERICA FIBONACCI 24x24 (100 W/BOBINA)")
    print("=" * 80)

    # 1. Esecuzione simulatore (salta se già presenti i 64 VTU)
    vtus = prepare_and_run_solver()

    # 2. Elaborazione dati fisici
    data_summary, time_s, forces, joule_eddy, per_coil_powers, b_series = process_simulation_data(vtus)

    # 3. Generazione Figura 27 a 300 DPI
    generate_figure_27(data_summary, time_s, forces, joule_eddy, per_coil_powers)

    # 4. Generazione Video Animato
    generate_animated_video(vtus, data_summary, time_s, forces)

    print("\n" + "=" * 80)
    print("TUTTE LE OPERAZIONI COMPLETATE CON SUCCESSO!")
    print(f"  - Dataset JSON: {OUT_JSON}")
    print(f"  - Figura 27 (300 DPI): {ROOT_FIGURES_DIR / FIG_27_NAME}")
    print(f"  - Video Animato: {ROOT_FIGURES_DIR / VIDEO_NAME}")
    print("=" * 80)


if __name__ == "__main__":
    main()
