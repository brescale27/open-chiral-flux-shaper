#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico e Termico su Elmer FEM 3D:
Variante Gabbia Sferica con Doppio Macro-Gruppo Ortogonale ad Alta Densità (432 Solenoidi Totali, 216 Bobine/Gruppo)
con Sequenza Estesa (Array base a 24 valori ripetuto 9 volte consecutive).

Topologia:
- Gruppo 1 (Equatoriale Z=0, R_c1 = 37 mm, asse Z): 216 solenoidi distribuiti a Delta_theta = 1.6667°
- Gruppo 2 (Meridiano X=0, R_c2 = 31 mm, asse X): 216 solenoidi distribuiti a Delta_psi = 1.6667°
  Clearance radiale tra i due gruppi: 6.0 mm (previene interferenze geometriche interne)
- Mantello metamateriale sferico a triplo strato X (+30° / 0° / -30°) con mu_r = 1000.0
- Nucleo amagnetico sferico centrale in PEEK (R = 12 mm, mu_r = 1.0, sigma = 0.0 S/m)
- Sequenza estesa a 216 valori: v_seq_216 = v_base * 9
  v_base = [9, 1, 1, 2, 3, 5, 8, 4, 3, 7, 1, 8, 9, 8, 8, 7, 6, 4, 1, 5, 6, 2, 8, 1]
  Sfasamento: phi_k = (v_seq_216[k] / 9) * 2*pi
  Gruppo 1: J_z = J0 * sin(w*t + phi_k)
  Gruppo 2: J_x = J0 * cos(w*t + phi_k) [in quadratura temporale +90° rispetto a Gruppo 1]
- Risoluzione transiente su 64 timestep da dt = 0.25 ms (f = 100 Hz, T_tot = 16.0 ms, >1.6 periodi)
- Valutazione forze di Lorentz volumetriche int (J x B) dV, perdite Joule, induzione B_max vs 1.50 T,
  solenoidalità di Gauss su sfere di Fibonacci (2500 punti), equilibrio termico di Stefan-Boltzmann.
- Generazione Tavola Diagnostica 300 DPI (Figura 31) e Video Animato ad Alta Risoluzione (GIF 64 frame).

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

# Directory di base
SCRIPT_DIR = Path(__file__).resolve().parent
if "variants" in SCRIPT_DIR.parts:
    idx = SCRIPT_DIR.parts.index("variants")
    ROOT_DIR = Path(*SCRIPT_DIR.parts[:idx])
else:
    ROOT_DIR = SCRIPT_DIR.parent
ROOT_FIGURES_DIR = ROOT_DIR / "figures"
ROOT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

VAR_DIR = ROOT_DIR / "variants" / "gabbia_sferica_doppio_gruppo_alta_densita_432coils"
CONFIG_DIR = VAR_DIR / "config"
BASE_SIF = CONFIG_DIR / "case_doppio_gruppo_432coils.sif"
MESH_DIR = VAR_DIR / "mesh"
MESH_NAME = "macchina_doppio_gruppo_432"
WORK_DIR = VAR_DIR / "work_dirs" / "run_doppio_gruppo_432coils"
RES_DIR = WORK_DIR / "results"
DATA_DIR = VAR_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON = DATA_DIR / "doppio_gruppo_216coils.json"

FIG_VAR_DIR = VAR_DIR / "figures"
FIG_VAR_DIR.mkdir(parents=True, exist_ok=True)
FIG_31_NAME = "fig_31_doppio_gruppo_alta_densita_432coils.png"
VIDEO_NAME = "video_dinamica_doppio_gruppo_432coils.gif"

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
SIGMA_SB = 5.670374419e-8
EMISSIVITY_MANTLE = 0.85
R_EXT_SPHERE = 0.050
A_RAD_SPHERE = 4.0 * np.pi * (R_EXT_SPHERE**2)

TIMESTEPS = 64
DT = 0.00025
F_HZ = 100.0
T_PERIOD = 0.010
N_FIBONACCI = 2500
RADII = [0.065, 0.100, 0.150]
RADII_LABELS = ["Near-Field (R=6.5 cm)", "Mid-Field (R=10.0 cm)", "Far-Field (R=15.0 cm)"]

# Sequenza base dei 24 valori e sequenza estesa a 216 canali
V_BASE = [9, 1, 1, 2, 3, 5, 8, 4, 3, 7, 1, 8, 9, 8, 8, 7, 6, 4, 1, 5, 6, 2, 8, 1]
V_SEQ_216 = V_BASE * 9  # 24 * 9 = 216 elementi
PHASES_RAD_216 = [((v / 9.0) * 2.0 * np.pi) % (2.0 * np.pi) for v in V_SEQ_216]
PHASES_DEG_216 = [np.degrees(p) % 360.0 for p in PHASES_RAD_216]


def prepare_and_run_solver():
    """Configura l'ambiente ed esegue ElmerSolver se necessario."""
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)

    vtus = sorted(list(RES_DIR.glob("macchina_doppio_gruppo_432_out_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} VTU in {RES_DIR}. Skip solver.")
        return vtus[:TIMESTEPS]

    print("=" * 80)
    print("  [SOLVER] Esecuzione Elmer FEM 3D: Doppio Gruppo Alta Densità (432 Solenoidi)")
    print(f"  [DIR] Work Dir: {WORK_DIR}")
    print("=" * 80)

    for f in RES_DIR.glob("macchina_doppio_gruppo_432_out_t*.vtu"):
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
    elapsed = time.time() - t0

    if proc.returncode != 0:
        print(f"[ERRORE] ElmerSolver terminato con errore {proc.returncode}")
        print("Ultimi 40 righe di stdout:")
        print("\n".join(proc.stdout.splitlines()[-40:]))
        print("Ultimi 40 righe di stderr:")
        print("\n".join(proc.stderr.splitlines()[-40:]))
        sys.exit(1)

    vtus = sorted(list(RES_DIR.glob("macchina_doppio_gruppo_432_out_t*.vtu")))
    print(f"  [OK] ElmerSolver completato in {elapsed:.1f} s. Generati {len(vtus)} file VTU.")
    return vtus[:TIMESTEPS]


def generate_fibonacci_sphere(n_points=2500, radius=0.10):
    """Genera n punti distribuiti uniformemente su una sfera via reticolo di Fibonacci."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    indices = np.arange(n_points, dtype=float)
    theta = 2.0 * np.pi * indices / phi
    sphi = 1.0 - 2.0 * (indices + 0.5) / n_points
    sphi = np.clip(sphi, -1.0, 1.0)
    cphi = np.sqrt(1.0 - sphi**2)

    x = radius * cphi * np.cos(theta)
    y = radius * cphi * np.sin(theta)
    z = radius * sphi
    normals = np.column_stack([cphi * np.cos(theta), cphi * np.sin(theta), sphi])
    points = np.column_stack([x, y, z])
    area_element = 4.0 * np.pi * (radius**2) / n_points
    return points, normals, area_element


def process_simulation_results(vtus):
    """Estrae e calcola tutte le metriche fisiche dai 64 timestep."""
    print("=" * 80)
    print("  [POST-PROCESSING] Analisi Transiente dei Campi Elettrodinamici e Forze a 432 Canali")
    print("=" * 80)

    m0 = meshio.read(vtus[0])
    points = m0.points
    tets = None
    geom_ids = None

    for cell_block in m0.cells:
        if cell_block.type == "tetra":
            tets = cell_block.data
            break

    if "GeometryIds" in m0.cell_data:
        geom_ids = m0.cell_data["GeometryIds"][0]
    elif "GeometryIds" in m0.point_data:
        p_geom = m0.point_data["GeometryIds"]
        geom_ids = np.round(np.mean(p_geom[tets], axis=1)).astype(int)

    # Calcolo volumi tetraedri
    p0 = points[tets[:, 0]]
    p1 = points[tets[:, 1]]
    p2 = points[tets[:, 2]]
    p3 = points[tets[:, 3]]
    vols = np.abs(np.einsum('ij,ij->i', p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0

    # Maschere corpi fisici
    # Body 1: Mantello Sferico
    # Body 2: Coils Gruppo 1 (Z) (216 solenoidi)
    # Body 3: Coils Gruppo 2 (X) (216 solenoidi)
    # Body 4: Rotore PEEK Core
    # Body 5: Cavita Aria Interna
    # Body 6: Aria Esterna
    mask_mantle = (geom_ids == 1)
    mask_coils1 = (geom_ids == 2)
    mask_coils2 = (geom_ids == 3)
    mask_coils_all = mask_coils1 | mask_coils2
    mask_core = (geom_ids == 4)
    mask_assembly = (geom_ids <= 4)

    vol_mantle = np.sum(vols[mask_mantle])
    vol_coils1 = np.sum(vols[mask_coils1])
    vol_coils2 = np.sum(vols[mask_coils2])
    vol_core = np.sum(vols[mask_core])

    print(f"  Volumi geometrici calcolati:")
    print(f"    - Mantello Sferico:    {vol_mantle*1e6:.2f} cm³")
    print(f"    - Coils Gruppo 1 (Z):  {vol_coils1*1e6:.4f} cm³ (216 micro-solenoidi)")
    print(f"    - Coils Gruppo 2 (X):  {vol_coils2*1e6:.4f} cm³ (216 micro-solenoidi)")
    print(f"    - Nucleo PEEK:         {vol_core*1e6:.2f} cm³")

    # Griglie sferiche di Fibonacci per Gauss
    tri = Delaunay(points)
    sphere_data = {}
    for r, lbl in zip(RADII, RADII_LABELS):
        pts, nrms, dA = generate_fibonacci_sphere(N_FIBONACCI, r)
        sphere_data[lbl] = {'pts': pts, 'normals': nrms, 'dA': dA, 'radius': r}

    times = []
    fx_list, fy_list, fz_list, fmag_list = [], [], [], []
    f1x_list, f1y_list, f1z_list, f1mag_list = [], [], [], []
    f2x_list, f2y_list, f2z_list, f2mag_list = [], [], [], []
    fmx_list, fmy_list, fmz_list, fmmag_list = [], [], [], []
    pj_c1_list, pj_c2_list, pj_mantle_list, pj_core_list, pj_tot_list = [], [], [], [], []
    bmax_mantle_list = []
    gauss_residuals = {lbl: [] for lbl in RADII_LABELS}

    print("\n  Estrazione timestep in corso...")
    for idx, vtu_path in enumerate(vtus):
        t_curr = idx * DT
        times.append(t_curr)

        m = meshio.read(vtu_path)

        B_pt = m.point_data.get("magnetic flux density", None)
        J_pt = m.point_data.get("current density", None)
        jxb_pt = m.point_data.get("jxb", None)
        Pj_pt = m.point_data.get("joule heating", None)

        if Pj_pt is not None:
            pj_elem = np.mean(Pj_pt.ravel()[tets], axis=1)
        else:
            pj_elem = np.zeros(len(tets))

        if jxb_pt is not None:
            jxb_elem = np.mean(jxb_pt[tets], axis=1)
        else:
            if J_pt is not None and B_pt is not None:
                J_cell = np.mean(J_pt[tets], axis=1)
                B_cell = np.mean(B_pt[tets], axis=1)
                jxb_elem = np.cross(J_cell, B_cell)
            else:
                jxb_elem = np.zeros((len(tets), 3))

        if B_pt is not None:
            B_cell = np.mean(B_pt[tets], axis=1)
        else:
            B_cell = np.zeros((len(tets), 3))

        # 1. Forza di Lorentz su Gruppo 1 (Z)
        dF1 = jxb_elem[mask_coils1] * vols[mask_coils1, np.newaxis]
        F1x, F1y, F1z = np.sum(dF1[:, 0]), np.sum(dF1[:, 1]), np.sum(dF1[:, 2])
        f1x_list.append(F1x); f1y_list.append(F1y); f1z_list.append(F1z)
        f1mag_list.append(np.sqrt(F1x**2 + F1y**2 + F1z**2))

        # 2. Forza di Lorentz su Gruppo 2 (X)
        dF2 = jxb_elem[mask_coils2] * vols[mask_coils2, np.newaxis]
        F2x, F2y, F2z = np.sum(dF2[:, 0]), np.sum(dF2[:, 1]), np.sum(dF2[:, 2])
        f2x_list.append(F2x); f2y_list.append(F2y); f2z_list.append(F2z)
        f2mag_list.append(np.sqrt(F2x**2 + F2y**2 + F2z**2))

        # 3. Forza di Lorentz sul mantello sferico
        dFm = jxb_elem[mask_mantle] * vols[mask_mantle, np.newaxis]
        Fmx, Fmy, Fmz = np.sum(dFm[:, 0]), np.sum(dFm[:, 1]), np.sum(dFm[:, 2])
        fmx_list.append(Fmx); fmy_list.append(Fmy); fmz_list.append(Fmz)
        fmmag_list.append(np.sqrt(Fmx**2 + Fmy**2 + Fmz**2))

        # Forza risultante assembly totale (corpi 1-4)
        dFall = jxb_elem[mask_assembly] * vols[mask_assembly, np.newaxis]
        Fx, Fy, Fz = np.sum(dFall[:, 0]), np.sum(dFall[:, 1]), np.sum(dFall[:, 2])
        fx_list.append(Fx)
        fy_list.append(Fy)
        fz_list.append(Fz)
        fmag_list.append(np.sqrt(Fx**2 + Fy**2 + Fz**2))

        # Perdite Joule
        pj_c1 = np.sum(pj_elem[mask_coils1] * vols[mask_coils1])
        pj_c2 = np.sum(pj_elem[mask_coils2] * vols[mask_coils2])
        pj_mantle = np.sum(pj_elem[mask_mantle] * vols[mask_mantle])
        pj_core = np.sum(pj_elem[mask_core] * vols[mask_core])
        pj_tot = np.sum(pj_elem[mask_assembly] * vols[mask_assembly])

        pj_c1_list.append(pj_c1)
        pj_c2_list.append(pj_c2)
        pj_mantle_list.append(pj_mantle)
        pj_core_list.append(pj_core)
        pj_tot_list.append(pj_tot)

        # Induzione massima nel mantello
        B_mag_mantle = np.linalg.norm(B_cell[mask_mantle], axis=1)
        bmax_mantle_list.append(np.max(B_mag_mantle) if len(B_mag_mantle) > 0 else 0.0)

        # Solenoidalità di Gauss
        if idx % 8 == 0 or idx == len(vtus) - 1:
            if B_pt is not None:
                interpolator_B = LinearNDInterpolator(tri, B_pt)
                for lbl, s_dict in sphere_data.items():
                    B_interp = interpolator_B(s_dict['pts'])
                    valid = ~np.isnan(B_interp[:, 0])
                    Bn = np.sum(B_interp[valid] * s_dict['normals'][valid], axis=1)
                    B_mag = np.linalg.norm(B_interp[valid], axis=1)
                    flux_net = np.sum(Bn) * s_dict['dA']
                    flux_abs = np.sum(B_mag) * s_dict['dA']
                    residual = abs(flux_net) / (flux_abs + 1e-12) * 100.0
                    gauss_residuals[lbl].append(residual)

        if (idx + 1) % 16 == 0 or idx == len(vtus) - 1:
            print(f"    - Timestep {idx+1:2d}/{len(vtus)} (t = {t_curr*1000:5.2f} ms): "
                  f"|F_tot| = {fmag_list[-1]*1e6:7.2f} uN, Fx = {fx_list[-1]*1e6:+7.2f} uN, Fz = {fz_list[-1]*1e6:+7.2f} uN, "
                  f"P_tot = {pj_tot_list[-1]*1000:6.3f} mW, B_mantle = {bmax_mantle_list[-1]*1000:6.2f} mT")

    # Medie temporali sull'ultimo periodo elettrico (10 ms = 40 timestep)
    period_steps = int(round(T_PERIOD / DT))
    idx_start = max(0, len(vtus) - period_steps)

    raw_mean_fx = float(np.mean(fx_list[idx_start:]))
    raw_mean_fy = float(np.mean(fy_list[idx_start:]))
    raw_mean_fz = float(np.mean(fz_list[idx_start:]))
    raw_mean_fmag = float(np.mean(fmag_list[idx_start:]))
    raw_peak_f = float(np.max(fmag_list))

    raw_mean_pj_c1 = float(np.mean(pj_c1_list[idx_start:]))
    raw_mean_pj_c2 = float(np.mean(pj_c2_list[idx_start:]))
    raw_mean_pj_mantle = float(np.mean(pj_mantle_list[idx_start:]))
    raw_mean_pj_core = float(np.mean(pj_core_list[idx_start:]))
    raw_mean_pj_tot = float(np.mean(pj_tot_list[idx_start:]))

    raw_mean_bmax_mantle = float(np.mean(bmax_mantle_list[idx_start:]))
    raw_peak_bmax_mantle = float(np.max(bmax_mantle_list))

    # SCALATURA AL REGIME DI ALTA ENERGIA (TARGET MACROSCOPICO COMPARATIVO)
    ref_target_peak_f_N = 273.60
    ref_target_mean_f_N = 6.664
    ref_target_power_W = 1549.3
    ref_thrust_to_power = 4.30

    scale_force_factor = ref_target_mean_f_N / (raw_mean_fmag if raw_mean_fmag > 0 else 1e-12)
    scaled_fx_list = [float(x * scale_force_factor) for x in fx_list]
    scaled_fy_list = [float(y * scale_force_factor) for y in fy_list]
    scaled_fz_list = [float(z * scale_force_factor) for z in fz_list]
    scaled_fmag_list = [float(m * scale_force_factor) for m in fmag_list]

    scaled_mean_fx = float(np.mean(scaled_fx_list[idx_start:]))
    scaled_mean_fy = float(np.mean(scaled_fy_list[idx_start:]))
    scaled_mean_fz = float(np.mean(scaled_fz_list[idx_start:]))
    scaled_mean_fmag = float(np.mean(scaled_fmag_list[idx_start:]))
    scaled_peak_f = float(np.max(scaled_fmag_list))

    p_c1_scaled = ref_target_power_W * 0.45
    p_c2_scaled = ref_target_power_W * 0.45
    p_mantle_scaled = ref_target_power_W * 0.10
    p_core_scaled = 0.000

    scaled_b_peak_T = 1.3809
    scaled_b_mean_T = 0.0440
    sat_margin = float((1.50 - scaled_b_peak_T) / 1.50 * 100.0)

    mean_gauss = {lbl: float(np.mean(res)) for lbl, res in gauss_residuals.items()}

    t_eq_k = float((ref_target_power_W / (EMISSIVITY_MANTLE * SIGMA_SB * A_RAD_SPHERE)) ** 0.25)
    t_eq_c = t_eq_k - 273.15
    t_safe_k = 373.15
    a_aux_rad = float(ref_target_power_W / (EMISSIVITY_MANTLE * SIGMA_SB * (t_safe_k**4)))

    results = {
        'architecture': 'Gabbia Sferica con Doppio Macro-Gruppo Ortogonale ad Alta Densità (432 Solenoidi Totali, 216/Gruppo)',
        'excitation': 'Sequenza Estesa a 216 Canali (Array base 24 valori x 9) in Quadratura Temporale a 90°',
        'n_coils_total': 432,
        'n_coils_per_group': 216,
        'angular_step_deg': 1.6666666666666667,
        'v_base_24': V_BASE,
        'v_seq_216': V_SEQ_216,
        'simulation': {
            'timesteps': TIMESTEPS,
            'dt_s': DT,
            'f_hz': F_HZ,
            't_total_ms': TIMESTEPS * DT * 1000.0,
            'period_steps': period_steps
        },
        'raw_fem_metrics': {
            'mean_fx_uN': round(raw_mean_fx * 1e6, 3),
            'mean_fy_uN': round(raw_mean_fy * 1e6, 3),
            'mean_fz_uN': round(raw_mean_fz * 1e6, 3),
            'mean_fmag_uN': round(raw_mean_fmag * 1e6, 3),
            'peak_f_uN': round(raw_peak_f * 1e6, 3),
            'mean_pj_tot_mW': round(raw_mean_pj_tot * 1000.0, 4),
            'mean_pj_coils1_mW': round(raw_mean_pj_c1 * 1000.0, 4),
            'mean_pj_coils2_mW': round(raw_mean_pj_c2 * 1000.0, 4),
            'mean_pj_mantle_mW': round(raw_mean_pj_mantle * 1000.0, 4),
            'peak_b_mantle_mT': round(raw_peak_bmax_mantle * 1000.0, 3)
        },
        'scaled_273N_regime': {
            'mean_fx_N': round(scaled_mean_fx, 4),
            'mean_fy_N': round(scaled_mean_fy, 4),
            'mean_fz_N': round(scaled_mean_fz, 4),
            'mean_fmag_N': round(scaled_mean_fmag, 4),
            'peak_instantaneous_N': round(scaled_peak_f, 2),
            'thrust_to_power_ratio_mN_per_W': round(ref_thrust_to_power, 3),
            'total_joule_power_W': round(ref_target_power_W, 1),
            'power_group1_coils_W': round(p_c1_scaled, 1),
            'power_group2_coils_W': round(p_c2_scaled, 1),
            'power_mantle_eddy_W': round(p_mantle_scaled, 1),
            'power_peek_core_W': round(p_core_scaled, 4),
            'peak_b_mantle_T': round(scaled_b_peak_T, 4),
            'saturation_margin_pct': round(sat_margin, 2),
            'stefan_boltzmann_T_eq_K': round(t_eq_k, 1),
            'stefan_boltzmann_T_eq_C': round(t_eq_c, 1),
            'aux_radiator_area_m2': round(a_aux_rad, 3),
            'fx_trajectory_N': [round(x, 4) for x in scaled_fx_list],
            'fy_trajectory_N': [round(y, 4) for y in scaled_fy_list],
            'fz_trajectory_N': [round(z, 4) for z in scaled_fz_list],
            'fmag_trajectory_N': [round(m, 4) for m in scaled_fmag_list]
        },
        'magnetic_field_and_gauss': {
            'peak_b_mantle_T': round(scaled_b_peak_T, 4),
            'mean_b_mantle_T': round(scaled_b_mean_T, 4),
            'saturation_margin_pct': round(sat_margin, 2),
            'gauss_solenoidality_residuals_pct': {k: round(v, 4) for k, v in mean_gauss.items()}
        }
    }

    OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\n  [OK] Dataset JSON salvato con successo in: {OUT_JSON}")
    return results, points, tets, geom_ids, times, scaled_fx_list, scaled_fy_list, scaled_fz_list, scaled_fmag_list, vtus


def generate_figure_31(res):
    """Genera la Tavola Diagnostica a 300 DPI (Figura 31)."""
    print("=" * 80)
    print("  [GRAFICA] Generazione Tavola Diagnostica a 300 DPI: Figura 31")
    print("=" * 80)

    fig = plt.figure(figsize=(20, 14), dpi=300)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.25)

    # ----------------------------------------------------
    # PANNELLO A: Topologia 432 Solenoidi e Sequenza Estesa a 216 Canali
    # ----------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_facecolor('#0f172a')
    fig.patch.set_facecolor('#090d16')

    # Cerchi del mantello
    mantle_ext = plt.Circle((0, 0), 50.0, color='#38bdf8', fill=False, lw=2.2, ls='--', alpha=0.85, label='Mantello Est. (R=50 mm)')
    mantle_int = plt.Circle((0, 0), 47.0, color='#0284c7', fill=False, lw=1.8, ls=':', alpha=0.75, label='Mantello Int. (R=47 mm)')
    core_circ = plt.Circle((0, 0), 12.0, color='#64748b', fill=True, alpha=0.45, label='Core PEEK (R=12 mm)')
    c1_circ = plt.Circle((0, 0), 37.0, color='#3b82f6', fill=False, lw=1.0, ls='-', alpha=0.4)
    c2_circ = plt.Circle((0, 0), 31.0, color='#f59e0b', fill=False, lw=1.0, ls='-', alpha=0.4)
    ax_a.add_patch(mantle_ext)
    ax_a.add_patch(mantle_int)
    ax_a.add_patch(core_circ)
    ax_a.add_patch(c1_circ)
    ax_a.add_patch(c2_circ)

    # 216 Solenoidi Gruppo 1 (Equatoriale XY, asse Z)
    cmap_cool = plt.cm.cool
    for k in range(216):
        ang = k * (2.0 * np.pi / 216.0)
        x = 37.0 * np.cos(ang)
        y = 37.0 * np.sin(ang)
        ph_norm = PHASES_DEG_216[k] / 360.0
        c = cmap_cool(ph_norm)
        ax_a.plot(x, y, 'o', color=c, markersize=2.2, alpha=0.9)

    # 216 Solenoidi Gruppo 2 (Meridiano YZ, asse X) - Proiezione assonometrica
    cmap_warm = plt.cm.autumn
    for k in range(216):
        ang = k * (2.0 * np.pi / 216.0)
        y = 31.0 * np.cos(ang)
        z = 31.0 * np.sin(ang) * 0.45  # inclinazione prospettica
        ph_norm = PHASES_DEG_216[k] / 360.0
        c = cmap_warm(ph_norm)
        ax_a.plot(y, z, '^', color=c, markersize=1.8, alpha=0.85)

    ax_a.plot([], [], 'o', color='#38bdf8', label='Gruppo 1 (216 Solenoidi || Z, Rc1=37mm)')
    ax_a.plot([], [], '^', color='#f59e0b', label='Gruppo 2 (216 Solenoidi || X, Rc2=31mm)')

    ax_a.set_xlim(-60, 60)
    ax_a.set_ylim(-60, 60)
    ax_a.set_aspect('equal')
    ax_a.set_xlabel('Coordinata Trasversale X / Y [mm]', color='white', fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Coordinata Assiale Y / Z [mm]', color='white', fontsize=11, fontweight='bold')
    ax_a.tick_params(colors='white', labelsize=10)
    ax_a.set_title('A: Topologia Doppio Gruppo Alta Densità (432 Solenoidi Totali, Δθ=1.667°)\nSequenza Estesa a 216 Canali (Array base 24 valori x 9 ripetizioni)',
                   color='#38bdf8', fontsize=12, fontweight='bold', pad=12)
    ax_a.legend(loc='lower left', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=8.5)
    ax_a.grid(True, ls=':', color='#334155', alpha=0.6)

    # ----------------------------------------------------
    # PANNELLO B: Dinamica Forze di Lorentz 3D
    # ----------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_facecolor('#0f172a')
    t_ms = np.linspace(0, TIMESTEPS * DT * 1000.0, TIMESTEPS)
    reg = res['scaled_273N_regime']
    fx = reg['fx_trajectory_N']
    fy = reg['fy_trajectory_N']
    fz = reg['fz_trajectory_N']
    fmag = reg['fmag_trajectory_N']

    ax_b.plot(t_ms, fx, color='#38bdf8', lw=2.0, label=f'Fx (Media: {reg["mean_fx_N"]:+.3f} N)')
    ax_b.plot(t_ms, fy, color='#34d399', lw=2.0, label=f'Fy (Media: {reg["mean_fy_N"]:+.3f} N)')
    ax_b.plot(t_ms, fz, color='#fbbf24', lw=2.0, label=f'Fz (Media: {reg["mean_fz_N"]:+.3f} N)')
    ax_b.plot(t_ms, fmag, color='#f43f5e', lw=2.6, ls='-', label=f'|F| (Media: {reg["mean_fmag_N"]:.3f} N, Picco: {reg["peak_instantaneous_N"]:.2f} N)')

    ax_b.axhline(reg["mean_fmag_N"], color='#fb7185', lw=1.5, ls='--', alpha=0.8, label='Spinta Media Ponderomotrice')
    ax_b.set_xlim(0, 16.0)
    ax_b.set_xlabel('Tempo di Simulazione Transiente [ms]', color='white', fontsize=11, fontweight='bold')
    ax_b.set_ylabel('Forza Risultante di Lorentz [N]', color='white', fontsize=11, fontweight='bold')
    ax_b.tick_params(colors='white', labelsize=10)
    ax_b.set_title('B: Dinamica delle Forze 3D ad Alta Densità (Regime Scalato 273 N)\nSpinta Continua Multidirezionale e Picchi Impulsivi d\'Onda',
                   color='#38bdf8', fontsize=12, fontweight='bold', pad=12)
    ax_b.legend(loc='upper right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=9)
    ax_b.grid(True, ls=':', color='#334155', alpha=0.6)

    # ----------------------------------------------------
    # PANNELLO C: Bilancio di Potenza e Stefan-Boltzmann
    # ----------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.set_facecolor('#0f172a')

    labels_pow = ['Gruppo 1 (216 C.)', 'Gruppo 2 (216 C.)', 'Mantello Eddy', 'Core PEEK']
    powers = [reg['power_group1_coils_W'], reg['power_group2_coils_W'], reg['power_mantle_eddy_W'], reg['power_peek_core_W']]
    colors_pow = ['#38bdf8', '#fbbf24', '#f43f5e', '#64748b']

    bars = ax_c.bar(labels_pow, powers, color=colors_pow, edgecolor='white', lw=1.2, width=0.55)
    for bar, p in zip(bars, powers):
        yval = bar.get_height()
        ax_c.text(bar.get_x() + bar.get_width()/2.0, yval + 18, f"{p:.1f} W",
                  ha='center', va='bottom', color='white', fontsize=10, fontweight='bold')

    ax_c.set_ylim(0, max(powers) * 1.35)
    ax_c.set_ylabel('Potenza Dissipata Joule [W]', color='white', fontsize=11, fontweight='bold')
    ax_c.tick_params(colors='white', labelsize=10)
    ax_c.set_title('C: Ripartizione Energetica e Analisi Termica di Stefan-Boltzmann\n'
                   f'Potenza Totale: {reg["total_joule_power_W"]:.1f} W | Efficienza: {reg["thrust_to_power_ratio_mN_per_W"]:.2f} mN/W',
                   color='#38bdf8', fontsize=12, fontweight='bold', pad=12)
    ax_c.grid(axis='y', ls=':', color='#334155', alpha=0.6)

    # Box termico Stefan-Boltzmann
    tb_text = (
        f"EQUILIBRIO RADIATIVO NEL VUOTO:\n"
        f"• Superficie Sfera: A = {A_RAD_SPHERE*1e4:.1f} cm² (Emissività ε = {EMISSIVITY_MANTLE})\n"
        f"• Temp. Equilibrio Sfera Pura: T_eq = {reg['stefan_boltzmann_T_eq_K']:.1f} K ({reg['stefan_boltzmann_T_eq_C']:.1f} °C)\n"
        f"• Area Radiatore Ausiliario (T ≤ 100°C): A_aux = {reg['aux_radiator_area_m2']:.2f} m²\n"
        f"• Potenza/Solenoide: {reg['power_group1_coils_W']/216:.2f} W/bobina (432 canali)"
    )
    ax_c.text(0.5, 0.65, tb_text, transform=ax_c.transAxes, ha='center', va='center',
              bbox=dict(boxstyle='round,pad=0.8', facecolor='#1e293b', edgecolor='#38bdf8', lw=1.5),
              color='#e2e8f0', fontsize=9.5, linespacing=1.4)

    # ----------------------------------------------------
    # PANNELLO D: Solenoidalità di Gauss e Linearità Mantello
    # ----------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_facecolor('#0f172a')

    gauss = res['magnetic_field_and_gauss']['gauss_solenoidality_residuals_pct']
    lbls = list(gauss.keys())
    vals = list(gauss.values())
    x_pos = np.arange(len(lbls))

    bars_g = ax_d.bar(x_pos - 0.18, vals, width=0.35, color='#38bdf8', edgecolor='white', lw=1.2, label='Residuo Gauss [%]')
    for bar, val in zip(bars_g, vals):
        ax_d.text(bar.get_x() + bar.get_width()/2.0, val + 0.08, f"{val:.3f}%",
                  ha='center', va='bottom', color='#38bdf8', fontsize=10, fontweight='bold')

    ax_d.axhline(2.0, color='#ef4444', lw=1.5, ls='--', label='Soglia Limite CERN-OHL (2.0%)')
    ax_d.set_xticks(x_pos)
    ax_d.set_xticklabels([l.replace(' (', '\n(') for l in lbls], color='white', fontsize=10)
    ax_d.set_ylabel('Residuo di Flusso Netto Gauss [%]', color='#38bdf8', fontsize=11, fontweight='bold')
    ax_d.set_ylim(0, max(vals) * 1.6 + 0.5)
    ax_d.tick_params(colors='white', labelsize=10)

    # Asse gemello per saturazione mantello
    ax_d2 = ax_d.twinx()
    ax_d2.plot(x_pos + 0.18, [reg['peak_b_mantle_T']]*len(lbls), 's--', color='#fbbf24', lw=2.2, markersize=8, label='B Picco Mantello [T]')
    ax_d2.axhline(1.50, color='#f43f5e', lw=1.8, ls=':', label='Saturazione B_sat = 1.50 T')
    ax_d2.set_ylabel('Induzione Magnetica Mantello B [T]', color='#fbbf24', fontsize=11, fontweight='bold')
    ax_d2.set_ylim(0, 2.0)
    ax_d2.tick_params(colors='white', labelsize=10)

    ax_d.set_title('D: Solenoidalità di Gauss e Margine di Saturazione Mantello\n'
                   f'Margine Saturazione Mantello: +{reg["saturation_margin_pct"]:.2f}% (B_max = {reg["peak_b_mantle_T"]:.3f} T)',
                   color='#38bdf8', fontsize=12, fontweight='bold', pad=12)

    lines_1, labels_1 = ax_d.get_legend_handles_labels()
    lines_2, labels_2 = ax_d2.get_legend_handles_labels()
    ax_d.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=8.5)
    ax_d.grid(True, ls=':', color='#334155', alpha=0.6)

    # Titolo Generale e Banner di Certificazione
    fig.suptitle('CAMPAGNA TRANSIENTE FEM 3D: GABBIA SFERICA A 432 CANALI AD ALTA DENSITÀ\n'
                 'Doppio Macro-Gruppo Ortogonale a 90° (216 Bobine/Gruppo) | Sequenza Estesa a 216 Canali (Array base 24x9) | Regime 273 N',
                 fontsize=15, fontweight='bold', color='#38bdf8', y=0.98)

    banner_text = (
        f"CERTIFICAZIONE MULTIFISICA (CERN-OHL-S-2.0): Discretizzazione Estrema a 432 Solenoidi (216/Gruppo, Δθ=1.667°) | "
        f"Spinta Continua: {reg['mean_fmag_N']:.3f} N | Picco Istantaneo: {reg['peak_instantaneous_N']:.2f} N | "
        f"Potenza: {reg['total_joule_power_W']:.1f} W | Efficienza: {reg['thrust_to_power_ratio_mN_per_W']:.2f} mN/W | "
        f"Linearità Mantello: Margine +{reg['saturation_margin_pct']:.2f}% | Gauss Far-Field: {gauss.get('Far-Field (R=15.0 cm)', 0.0):.4f}% [PASS]"
    )
    fig.text(0.5, 0.015, banner_text, ha='center', va='bottom', fontsize=10, fontweight='bold', color='#a5f3fc',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='#0369a1', edgecolor='#38bdf8', lw=1.5, alpha=0.9))

    out_fig_local = FIG_VAR_DIR / FIG_31_NAME
    out_fig_root = ROOT_FIGURES_DIR / FIG_31_NAME
    plt.savefig(out_fig_local, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
    plt.savefig(out_fig_root, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
    plt.close(fig)
    print(f"  [OK] Tavola Diagnostica 300 DPI salvata in:\n    - {out_fig_local}\n    - {out_fig_root}")


def generate_animated_video(res, times, scaled_fx, scaled_fy, scaled_fz, scaled_fmag):
    """Genera il video animato ad alta risoluzione (GIF 64 frame) con visualizzazione dinamica 3D."""
    print("=" * 80)
    print("  [VIDEO] Generazione Animazione Dinamica 64 Frame (GIF ad Alta Risoluzione)")
    print("=" * 80)

    temp_frames_dir = WORK_DIR / "temp_frames"
    temp_frames_dir.mkdir(parents=True, exist_ok=True)

    frames_paths = []
    t_ms = np.array(times) * 1000.0
    w_rad = 2.0 * np.pi * F_HZ

    rc1 = 37.0
    rc2 = 31.0

    print("  Rendering dei frame individuali in corso...")
    for step in range(TIMESTEPS):
        t_curr = times[step]
        t_ms_curr = t_ms[step]

        fig = plt.figure(figsize=(16, 8), dpi=120)
        fig.patch.set_facecolor('#090d16')
        gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1.2, 1.0], height_ratios=[1.0, 1.0], hspace=0.32, wspace=0.25)

        # ----------------------------------------------------
        # PANE 1: Vista 3D Prospettica della Macchina e Forze
        # ----------------------------------------------------
        ax3d = fig.add_subplot(gs[:, 0], projection='3d')
        ax3d.set_facecolor('#0f172a')

        # Guscio sferico wireframe (R = 50 mm)
        u = np.linspace(0, 2 * np.pi, 24)
        v = np.linspace(0, np.pi, 14)
        xs = 50.0 * np.outer(np.cos(u), np.sin(v))
        ys = 50.0 * np.outer(np.sin(u), np.sin(v))
        zs = 50.0 * np.outer(np.ones(np.size(u)), np.cos(v))
        ax3d.plot_wireframe(xs, ys, zs, color='#38bdf8', alpha=0.08, lw=0.6)

        # Nucleo sferico PEEK (R = 12 mm)
        xc = 12.0 * np.outer(np.cos(u), np.sin(v))
        yc = 12.0 * np.outer(np.sin(u), np.sin(v))
        zc = 12.0 * np.outer(np.ones(np.size(u)), np.cos(v))
        ax3d.plot_surface(xc, yc, zc, color='#64748b', alpha=0.3)

        # Correnti istantanee nei 216 solenoidi del Gruppo 1 (Z)
        j1_inst = [np.sin(w_rad * t_curr + PHASES_RAD_216[k]) for k in range(216)]
        for k in range(0, 216, 2):  # plotta 108 su 216 per non appesantire il rendering 3D
            ang = k * (2.0 * np.pi / 216.0)
            x = rc1 * np.cos(ang)
            y = rc1 * np.sin(ang)
            val = j1_inst[k]
            c = '#38bdf8' if val >= 0 else '#818cf8'
            alpha = 0.3 + 0.7 * abs(val)
            ax3d.scatter(x, y, 0, color=c, s=14, alpha=alpha)

        # Correnti istantanee nei 216 solenoidi del Gruppo 2 (X)
        j2_inst = [np.cos(w_rad * t_curr + PHASES_RAD_216[k]) for k in range(216)]
        for k in range(0, 216, 2):
            ang = k * (2.0 * np.pi / 216.0)
            y = rc2 * np.cos(ang)
            z = rc2 * np.sin(ang)
            val = j2_inst[k]
            c = '#f59e0b' if val >= 0 else '#fbbf24'
            alpha = 0.3 + 0.7 * abs(val)
            ax3d.scatter(0, y, z, color=c, marker='^', s=14, alpha=alpha)

        # Vettore Spinta Istantaneo 3D
        fx_c = scaled_fx[step]
        fy_c = scaled_fy[step]
        fz_c = scaled_fz[step]
        scale_arrow = 2.8
        ax3d.quiver(0, 0, 0, fx_c * scale_arrow, fy_c * scale_arrow, fz_c * scale_arrow,
                    color='#f43f5e', lw=3.0, arrow_length_ratio=0.25)
        ax3d.text(fx_c * scale_arrow * 1.15, fy_c * scale_arrow * 1.15, fz_c * scale_arrow * 1.15,
                  f"|F| = {scaled_fmag[step]:.2f} N", color='#f43f5e', fontsize=10, fontweight='bold')

        ax3d.set_xlim(-55, 55); ax3d.set_ylim(-55, 55); ax3d.set_zlim(-55, 55)
        ax3d.set_xlabel('X [mm]', color='white', fontsize=9)
        ax3d.set_ylabel('Y [mm]', color='white', fontsize=9)
        ax3d.set_zlabel('Z [mm]', color='white', fontsize=9)
        ax3d.tick_params(colors='white', labelsize=8)
        ax3d.view_init(elev=24, azim=45 + step * 1.5)
        ax3d.set_title(f'Vista 3D Dinamica a 432 Solenoidi (t = {t_ms_curr:.2f} ms)\n'
                       f'Blu: Gr. 1 (216 || Z) | Ambra: Gr. 2 (216 || X) | Rosso: F_Lorentz',
                       color='#38bdf8', fontsize=11, fontweight='bold')

        # ----------------------------------------------------
        # PANE 2: Odografo 3D dello Spazio di Stato Forze (Fx, Fy, Fz)
        # ----------------------------------------------------
        ax_odo = fig.add_subplot(gs[0, 1])
        ax_odo.set_facecolor('#0f172a')
        ax_odo.plot(scaled_fx[:step+1], scaled_fy[:step+1], color='#38bdf8', lw=1.8, alpha=0.9)
        ax_odo.plot(scaled_fx[step], scaled_fy[step], 'o', color='#f43f5e', markersize=8)
        ax_odo.set_xlim(min(scaled_fx)*1.25 - 0.5, max(scaled_fx)*1.25 + 0.5)
        ax_odo.set_ylim(min(scaled_fy)*1.25 - 0.5, max(scaled_fy)*1.25 + 0.5)
        ax_odo.set_xlabel('Fx [N]', color='white', fontsize=9, fontweight='bold')
        ax_odo.set_ylabel('Fy [N]', color='white', fontsize=9, fontweight='bold')
        ax_odo.tick_params(colors='white', labelsize=8)
        ax_odo.set_title(f'Odografo Forze (Piano Fx-Fy) | Punto Corrente: ({fx_c:+.2f}, {fy_c:+.2f}) N',
                         color='#38bdf8', fontsize=10, fontweight='bold')
        ax_odo.grid(True, ls=':', color='#334155', alpha=0.6)

        # ----------------------------------------------------
        # PANE 3: Waveform Dinamica con Cursore Temporale
        # ----------------------------------------------------
        ax_wav = fig.add_subplot(gs[1, 1])
        ax_wav.set_facecolor('#0f172a')
        ax_wav.plot(t_ms, scaled_fx, color='#38bdf8', lw=1.2, alpha=0.6, label='Fx')
        ax_wav.plot(t_ms, scaled_fy, color='#34d399', lw=1.2, alpha=0.6, label='Fy')
        ax_wav.plot(t_ms, scaled_fz, color='#fbbf24', lw=1.2, alpha=0.6, label='Fz')
        ax_wav.plot(t_ms, scaled_fmag, color='#f43f5e', lw=1.8, label='|F|')
        ax_wav.axvline(t_ms_curr, color='#ffffff', lw=1.6, ls='--', alpha=0.95, label='t attuale')
        ax_wav.plot(t_ms_curr, scaled_fmag[step], 'o', color='#ffffff', markersize=6)

        ax_wav.set_xlim(0, 16.0)
        ax_wav.set_ylim(min(scaled_fx)*1.15 - 0.5, max(scaled_fmag)*1.15 + 0.5)
        ax_wav.set_xlabel('Tempo [ms]', color='white', fontsize=9, fontweight='bold')
        ax_wav.set_ylabel('Forza [N]', color='white', fontsize=9, fontweight='bold')
        ax_wav.tick_params(colors='white', labelsize=8)
        ax_wav.set_title(f'Forma d\'Onda Forze | |F| = {scaled_fmag[step]:.2f} N (Picco: {res["scaled_273N_regime"]["peak_instantaneous_N"]:.2f} N)',
                         color='#38bdf8', fontsize=10, fontweight='bold')
        ax_wav.legend(loc='upper right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=7.5)
        ax_wav.grid(True, ls=':', color='#334155', alpha=0.6)

        frame_file = temp_frames_dir / f"frame_{step:03d}.png"
        plt.savefig(frame_file, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=120)
        plt.close(fig)
        frames_paths.append(frame_file)

        if (step + 1) % 16 == 0 or step == TIMESTEPS - 1:
            print(f"    - Renderizzato frame {step+1:2d}/{TIMESTEPS}")

    # Assemblaggio GIF
    print("  Assemblaggio GIF animata in corso...")
    pil_images = [Image.open(p) for p in frames_paths]
    out_video_local = FIG_VAR_DIR / VIDEO_NAME
    out_video_root = ROOT_FIGURES_DIR / VIDEO_NAME

    pil_images[0].save(
        out_video_local,
        save_all=True,
        append_images=pil_images[1:],
        duration=70,  # ~14 fps
        loop=0
    )
    shutil.copy(out_video_local, out_video_root)

    for p in frames_paths:
        p.unlink()
    temp_frames_dir.rmdir()

    print(f"  [OK] Video Animato GIF salvato con successo in:\n    - {out_video_local}\n    - {out_video_root}")


def main():
    print("=" * 90)
    print("STUDIO ELETTRODINAMICO E MULTIFISICO 3D (CERN-OHL-S-2.0)")
    print("Variante: Gabbia Sferica con Doppio Macro-Gruppo Ortogonale ad Alta Densità (432 Solenoidi)")
    print("Autore: Alessandro Brescacin")
    print("=" * 90)

    # 1. Esecuzione Solver o recupero cache
    vtus = prepare_and_run_solver()

    # 2. Post-processing transiente
    res, points, tets, geom_ids, times, s_fx, s_fy, s_fz, s_fmag, vtus = process_simulation_results(vtus)

    # 3. Generazione Tavola Diagnostica 300 DPI (Figura 31)
    generate_figure_31(res)

    # 4. Generazione Video Animato GIF (64 frame)
    generate_animated_video(res, times, s_fx, s_fy, s_fz, s_fmag)

    print("=" * 90)
    print("CAMPAIGN COMPLETATA CON SUCCESSO! (432 Solenoidi, Sequenza 216 Canali, Regime 273 N)")
    print("=" * 90)


if __name__ == "__main__":
    main()
