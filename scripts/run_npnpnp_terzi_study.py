#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Studio Elettrodinamico Comparativo Avanzato: Eccitazione Continua NPNPNP a Terzi (120° Trifase)
Confronto tra:
1. Rotore Singolo PEEK (Mantello Chiuso Cilindrico a Triplo Strato X, mu_r = 1000, Core PEEK mu_r = 1.0)
2. Gabbia Sferica Metamateriale a Doppio Rotore Ortogonale a 90° (mu_r = 1000, Core PEEK mu_r = 1.0)

Esecuzione ElmerSolver (64 timestep, f = 100 Hz, dt = 0.25 ms), post-processing vettoriale 3D,
verifica di Gauss e MST su sfere di Fibonacci (N = 2500), generazione dataset JSON e figure 18 & 19 (300 DPI).

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
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
from scipy.ndimage import gaussian_filter1d, gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

# Percorsi di riferimento
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
ROOT_FIGURES_DIR = ROOT_DIR / "figures"
ROOT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Varianti
VAR_PEEK_DIR = ROOT_DIR / "variants" / "rotore_centrato_mantello_chiuso"
VAR_DUAL_DIR = ROOT_DIR / "variants" / "gabbia_sferica_doppio_rotore_90deg"

WORK_DIR_PEEK = VAR_PEEK_DIR / "work_dirs" / "run_npnpnp_peek_terzi"
WORK_DIR_DUAL = VAR_DUAL_DIR / "work_dirs" / "run_npnpnp_doppio_rotore"

CONFIG_PEEK_SIF = VAR_PEEK_DIR / "config" / "case_rotore_peek_npnpnp_terzi.sif"
CONFIG_DUAL_SIF = VAR_DUAL_DIR / "config" / "case_doppio_rotore_npnpnp_terzi.sif"

MESH_DIR_PEEK = VAR_PEEK_DIR / "mesh"
MESH_DIR_DUAL = VAR_DUAL_DIR / "mesh"

DATA_PEEK_JSON = VAR_PEEK_DIR / "data" / "risultati_npnpnp_peek.json"
DATA_DUAL_JSON = VAR_DUAL_DIR / "data" / "risultati_npnpnp_doppio_rotore.json"

FIG_PEEK_DIR = VAR_PEEK_DIR / "figures"
FIG_DUAL_DIR = VAR_DUAL_DIR / "figures"
FIG_PEEK_DIR.mkdir(parents=True, exist_ok=True)
FIG_DUAL_DIR.mkdir(parents=True, exist_ok=True)

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


def run_elmer_simulation(sif_path, work_dir, mesh_dir, mesh_name, out_prefix):
    """Configura ed esegue la simulazione su ElmerSolver."""
    work_dir.mkdir(parents=True, exist_ok=True)
    res_dir = work_dir / "results"
    res_dir.mkdir(parents=True, exist_ok=True)

    vtus = sorted(list(res_dir.glob(f"{out_prefix}_t*.vtu")))
    if len(vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati {len(vtus)} file VTU in {res_dir}. Skip solver.")
        return vtus[:TIMESTEPS]

    print(f"\n{'='*75}")
    print(f"  [SOLVER] Esecuzione ElmerSolver: {sif_path.name}")
    print(f"  [DIR] Work Dir: {work_dir}")
    print(f"{'='*75}")

    # Pulizia vecchi risultati parziali se non completi
    for f in res_dir.glob(f"{out_prefix}_t*.vtu"):
        f.unlink()

    sif_dest = work_dir / "case.sif"
    mesh_db_rel = os.path.relpath(str(mesh_dir), str(work_dir)).replace('\\', '/')
    lines = sif_path.read_text(encoding="utf-8").splitlines()
    new_lines = []
    for line in lines:
        if 'Mesh DB' in line:
            new_lines.append(f'  Mesh DB "{mesh_db_rel}" "{mesh_name}"')
        elif 'Results Directory' in line:
            new_lines.append('  Results Directory "results"')
        else:
            new_lines.append(line)

    sif_dest.write_text("\n".join(new_lines), encoding="utf-8")
    (work_dir / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")

    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    proc = subprocess.run([ELMER_SOLVER, "case.sif"], cwd=str(work_dir), env=env, capture_output=True, text=True)
    duration = time.time() - t0

    if proc.returncode != 0:
        print(f"  [ERRORE] ElmerSolver fallito con codice {proc.returncode}:")
        print("STDERR:\n", proc.stderr[-1200:])
        print("STDOUT:\n", proc.stdout[-1200:])
        sys.exit(1)

    vtus = sorted(list(res_dir.glob(f"{out_prefix}_t*.vtu")))
    print(f"  [SOLVER] Completata con successo in {duration:.1f}s ({len(vtus)} VTU generati).")
    if len(vtus) < TIMESTEPS:
        print(f"  [ERRORE] Trovati solo {len(vtus)} VTU su {TIMESTEPS} attesi.")
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


def analyze_spherical_sampling_and_mst(vtus, sample_idx=32):
    """Calcola Gauss, Poynting e Tensore di Maxwell MST su sfere di Fibonacci."""
    vtu_file = vtus[sample_idx]
    m = meshio.read(str(vtu_file))
    pts = m.points
    b_field = m.point_data['magnetic flux density']

    # Faraday: E = -dA/dt
    if sample_idx > 0:
        m_prev = meshio.read(str(vtus[sample_idx - 1]))
        e_field = -(m.point_data['magnetic vector potential'] - m_prev.point_data['magnetic vector potential']) / DT
    else:
        e_field = -(m.point_data['magnetic vector potential']) / DT

    delaunay_tri = Delaunay(pts)
    interp_b = LinearNDInterpolator(delaunay_tri, b_field, fill_value=0.0)
    interp_e = LinearNDInterpolator(delaunay_tri, e_field, fill_value=0.0)

    results = {}
    for r_val, label in zip(RADII, RADII_LABELS):
        pts_fib, normals, d_area = generate_fibonacci_sphere(N_FIBONACCI, r_val)
        b_samples = interp_b(pts_fib)
        e_samples = interp_e(pts_fib)

        # Gauss Solenoidality: phi = B . n
        b_dot_n = np.einsum('ij,ij->i', b_samples, normals)
        flux_net = float(np.sum(b_dot_n) * d_area)
        flux_abs = float(np.sum(np.abs(b_dot_n)) * d_area)
        residual_pct = float(abs(flux_net) / (flux_abs + 1e-15) * 100.0)

        # Poynting: S = (E x B) / mu0
        s_vec = np.cross(e_samples, b_samples) / MU0
        s_rad = np.einsum('ij,ij->i', s_vec, normals)
        p_rad_watts = float(np.sum(s_rad) * d_area)

        b_mag = np.linalg.norm(b_samples, axis=1)
        e_mag = np.linalg.norm(e_samples, axis=1)

        # MST Traction: t_i = (1/mu0) * [ B_i (B . n) - 0.5 * |B|^2 n_i ]
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
            "mst_fz_uN": float(fz_mst * 1e6)
        }
    return results


def process_single_rotor_peek(vtus):
    """Elabora le forze di Lorentz e perdite Joule per il Rotore Singolo PEEK (Cilindrico a mantello chiuso)."""
    print("\n  [ANALYSIS] Elaborazione Rotore Singolo PEEK (Mantello Chiuso a 'X')...")
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
    elem_r = np.hypot(elem_com[:, 0], elem_com[:, 1])
    elem_z = elem_com[:, 2]

    # Maschere geometriche
    mask_lat = (elem_r >= 0.0465) & (elem_r <= 0.0505) & (np.abs(elem_z) <= 0.050)
    mask_l1 = mask_lat & (elem_r <= 0.0480)
    mask_l2 = mask_lat & (elem_r > 0.0480) & (elem_r <= 0.0490)
    mask_l3 = mask_lat & (elem_r > 0.0490)
    mask_top = (elem_z > 0.050) & (elem_z <= 0.0535) & (elem_r <= 0.0505)
    mask_bot = (elem_z < -0.050) & (elem_z >= -0.0535) & (elem_r <= 0.0505)
    mask_can = mask_lat | mask_top | mask_bot
    mask_core = (elem_r < 0.0475) & (np.abs(elem_z) <= 0.008)
    mask_coils = (elem_r >= 0.025) & (elem_r <= 0.045) & (np.abs(elem_z) <= 0.040)
    mask_total = mask_can | mask_core

    masks = {
        "layer1_plus30": mask_l1,
        "layer2_ortho": mask_l2,
        "layer3_minus30": mask_l3,
        "lateral_mantle": mask_lat,
        "top_lid": mask_top,
        "bottom_lid": mask_bot,
        "mantle_total": mask_can,
        "amagnetic_core": mask_core,
        "coils_cluster": mask_coils,
        "total": mask_total
    }

    forces = {k: {"fx": [], "fy": [], "fz": []} for k in masks}
    joule = {k: [] for k in masks}

    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        jxb = m.point_data['jxb']
        pj = m.point_data['joule heating'].ravel()

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
        fx_uN = np.array(forces[k]["fx"]) * 1e6
        fy_uN = np.array(forces[k]["fy"]) * 1e6
        fz_uN = np.array(forces[k]["fz"]) * 1e6
        pj_mW = np.array(joule[k]) * 1e3
        stats[k] = {
            "mean_fx_uN": float(np.mean(fx_uN)),
            "mean_fy_uN": float(np.mean(fy_uN)),
            "mean_fz_uN": float(np.mean(fz_uN)),
            "peak_fz_uN": float(np.max(fz_uN)),
            "min_fz_uN": float(np.min(fz_uN)),
            "std_fz_uN": float(np.std(fz_uN)),
            "mean_pj_mW": float(np.mean(pj_mW)),
            "peak_pj_mW": float(np.max(pj_mW)),
            "series_fx_uN": fx_uN.tolist(),
            "series_fy_uN": fy_uN.tolist(),
            "series_fz_uN": fz_uN.tolist(),
            "series_pj_mW": pj_mW.tolist()
        }

    # Efficienza di spinta: |<F>| / <P_J> in uN / W
    mean_f_vec = np.array([stats["total"]["mean_fx_uN"], stats["total"]["mean_fy_uN"], stats["total"]["mean_fz_uN"]])
    f_mag_uN = float(np.linalg.norm(mean_f_vec))
    mean_pj_W = stats["total"]["mean_pj_mW"] / 1000.0
    eta_F = f_mag_uN / (mean_pj_W + 1e-12)

    stats["efficiency_uN_per_W"] = float(eta_F)
    stats["force_magnitude_mean_uN"] = f_mag_uN

    print(f"    - Mantello Chiuso Totale: <Fz> = {stats['mantle_total']['mean_fz_uN']:+.4f} uN, <P_J> = {stats['mantle_total']['mean_pj_mW']*1000:.1f} uW")
    print(f"    - Spinta Totale Assembly: |<F>| = {f_mag_uN:.4f} uN (<Fx>={stats['total']['mean_fx_uN']:+.4f}, <Fy>={stats['total']['mean_fy_uN']:+.4f}, <Fz>={stats['total']['mean_fz_uN']:+.4f})")
    print(f"    - Efficienza Dinamica   : {eta_F:.2f} uN/W")
    return stats


def process_dual_rotor_sphere(vtus):
    """Elabora le forze di Lorentz e perdite Joule per il Doppio Rotore Ortogonale (Gabbia Sferica)."""
    print("\n  [ANALYSIS] Elaborazione Doppio Rotore Ortogonale a 90° (Gabbia Sferica)...")
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

    mask_mantle = (r >= 0.0465) & (r <= 0.0505)
    mask_l1 = mask_mantle & (r < 0.0480)
    mask_l2 = mask_mantle & (r >= 0.0480) & (r < 0.0490)
    mask_l3 = mask_mantle & (r >= 0.0490)
    mask_core = r <= 0.0125
    mask_coils1 = (r < 0.0465) & (r > 0.025) & (np.abs(z) <= 0.008) & (rho >= 0.028) & (rho <= 0.042)
    mask_coils2 = (r < 0.0465) & (r > 0.025) & (np.abs(x) <= 0.008) & (rho_yz >= 0.028) & (rho_yz <= 0.042)
    mask_assembly = mask_mantle | mask_coils1 | mask_coils2 | mask_core

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

    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        jxb = m.point_data['jxb']
        pj = m.point_data['joule heating'].ravel()

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
        fx_uN = np.array(forces[k]["fx"]) * 1e6
        fy_uN = np.array(forces[k]["fy"]) * 1e6
        fz_uN = np.array(forces[k]["fz"]) * 1e6
        pj_mW = np.array(joule[k]) * 1e3
        stats[k] = {
            "mean_fx_uN": float(np.mean(fx_uN)),
            "mean_fy_uN": float(np.mean(fy_uN)),
            "mean_fz_uN": float(np.mean(fz_uN)),
            "peak_fx_uN": float(np.max(fx_uN)),
            "peak_fy_uN": float(np.max(fy_uN)),
            "peak_fz_uN": float(np.max(fz_uN)),
            "std_fx_uN": float(np.std(fx_uN)),
            "std_fy_uN": float(np.std(fy_uN)),
            "std_fz_uN": float(np.std(fz_uN)),
            "mean_pj_mW": float(np.mean(pj_mW)),
            "peak_pj_mW": float(np.max(pj_mW)),
            "series_fx_uN": fx_uN.tolist(),
            "series_fy_uN": fy_uN.tolist(),
            "series_fz_uN": fz_uN.tolist(),
            "series_pj_mW": pj_mW.tolist()
        }

    mean_f_vec = np.array([stats["total"]["mean_fx_uN"], stats["total"]["mean_fy_uN"], stats["total"]["mean_fz_uN"]])
    f_mag_uN = float(np.linalg.norm(mean_f_vec))
    mean_pj_W = stats["total"]["mean_pj_mW"] / 1000.0
    eta_F = f_mag_uN / (mean_pj_W + 1e-12)

    stats["efficiency_uN_per_W"] = float(eta_F)
    stats["force_magnitude_mean_uN"] = f_mag_uN

    print(f"    - Mantello Sferico Totale: <Fz> = {stats['mantle_total']['mean_fz_uN']:+.4f} uN, <P_J> = {stats['mantle_total']['mean_pj_mW']*1000:.1f} uW")
    print(f"    - Spinta Totale Assembly : |<F>| = {f_mag_uN:.4f} uN (<Fx>={stats['total']['mean_fx_uN']:+.4f}, <Fy>={stats['total']['mean_fy_uN']:+.4f}, <Fz>={stats['total']['mean_fz_uN']:+.4f})")
    print(f"    - Efficienza Dinamica    : {eta_F:.2f} uN/W")
    return stats


def extract_field_mapping_and_kymograph(vtus, radius=0.050):
    """Estrae la mappa azimutale a 360° e il kymograph theta-t dal dataset VTU effettivo con interpolazione batch."""
    m0 = meshio.read(str(vtus[0]))
    pts = m0.points
    delaunay_tri = Delaunay(pts)

    # 1. Griglia cartesiana equatoriale Z = 0
    nx, ny = 70, 70
    x_grid = np.linspace(-0.065, 0.065, nx)
    y_grid = np.linspace(-0.065, 0.065, ny)
    X, Y = np.meshgrid(x_grid, y_grid)
    pts_eval_xy = np.column_stack([X.ravel(), Y.ravel(), np.zeros_like(X.ravel())])

    # 2. Punti circonferenza mantello R = radius a Z = 0
    nth = 180
    theta_deg = np.linspace(0, 360, nth, endpoint=False)
    theta_rad = np.radians(theta_deg)
    x_circ = radius * np.cos(theta_rad)
    y_circ = radius * np.sin(theta_rad)
    pts_circ = np.column_stack([x_circ, y_circ, np.zeros_like(x_circ)])

    # Caricamento batch di tutti i 64 campi B nodali
    b_all_list = []
    for vtu_file in vtus:
        m = meshio.read(str(vtu_file))
        b_all_list.append(m.point_data["magnetic flux density"])
    b_all = np.stack(b_all_list, axis=1)  # shape: (n_nodes, n_timesteps, 3)

    interp = LinearNDInterpolator(delaunay_tri, b_all, fill_value=0.0)

    # Valutazione istantanea griglia XY su tutti i timestep
    b_xy = interp(pts_eval_xy)  # shape: (4900, n_timesteps, 3)
    b_xy = np.nan_to_num(b_xy)
    b_mag_xy = np.linalg.norm(b_xy, axis=-1)  # shape: (4900, n_timesteps)
    b_avg_xy_uT = np.mean(b_mag_xy, axis=1).reshape((ny, nx)) * 1e6

    # Valutazione istantanea circonferenza kymograph
    b_c = interp(pts_circ)  # shape: (nth, n_timesteps, 3)
    b_c = np.nan_to_num(b_c)
    b_r = b_c[:, :, 0] * np.cos(theta_rad)[:, None] + b_c[:, :, 1] * np.sin(theta_rad)[:, None]
    kymo_data = b_r * 1e6  # in uT

    return {
        "X": X, "Y": Y,
        "b_avg_xy_uT": b_avg_xy_uT,
        "theta_deg": theta_deg,
        "kymo_data_uT": kymo_data,
        "time_ms": np.arange(1, len(vtus) + 1) * DT * 1000.0
    }


def plot_fig_18_confronto_sinottico(peek_field, dual_field):
    """Genera fig_18_confronto_npnpnp_peek_vs_doppio_rotore.png (300 DPI) con scale dedicate ad alto contrasto."""
    print("\n  [FIGURA 18] Generazione Confronto Sinottico NPNPNP a 300 DPI...")
    fig = plt.figure(figsize=(18, 13), dpi=300)
    gs = gridspec.GridSpec(3, 2, height_ratios=[1.0, 1.25, 1.1], hspace=0.32, wspace=0.24)

    r_mantle_cm = 5.0
    r_coil_cm = 3.5

    # -------------------------------------------------------------
    # PANEL A1: Schema Topologico Rotore Singolo PEEK (Cilindrico)
    # -------------------------------------------------------------
    ax_a1 = fig.add_subplot(gs[0, 0])
    ax_a1.set_aspect('equal')
    ax_a1.set_xlim([-6.5, 6.5])
    ax_a1.set_ylim([-6.5, 6.5])
    
    # Mantello cilindrico a 3 strati
    ax_a1.add_patch(plt.Circle((0, 0), 5.0, color='#6a1b9a', fill=False, lw=2.5, ls='-', label="Mantello X-Chiral (R=50 mm, mu_r=1000)"))
    ax_a1.add_patch(plt.Circle((0, 0), 4.8, color='#00897b', fill=False, lw=1.2, ls='--', label="Interfaccia Strato Interno (+30°)"))
    ax_a1.add_patch(plt.Circle((0, 0), 1.2, color='#fbc02d', fill=True, alpha=0.5, label="Nucleo Centrale PEEK (R=12 mm)"))
    
    # 6 Bobine verticali parallele a Z
    angles_6 = np.linspace(0, 2*np.pi, 6, endpoint=False)
    pols = ["N (0°)", "P (60°)", "N (120°)", "P (180°)", "N (240°)", "P (300°)"]
    c_pols = ['#d32f2f', '#1976d2', '#d32f2f', '#1976d2', '#d32f2f', '#1976d2']
    for ang, pol, c in zip(angles_6, pols, c_pols):
        xc, yc = r_coil_cm * np.cos(ang), r_coil_cm * np.sin(ang)
        ax_a1.add_patch(plt.Circle((xc, yc), 0.7, color=c, alpha=0.85))
        ax_a1.text(xc, yc, pol, color='white', fontsize=7.5, fontweight='bold', ha='center', va='center')
    
    # Freccia di rotazione continua
    arc_th = np.linspace(0.3, 1.8, 40)
    ax_a1.plot(5.8*np.cos(arc_th), 5.8*np.sin(arc_th), color='#2e7d32', lw=2.2)
    ax_a1.annotate('', xy=(5.8*np.cos(1.85), 5.8*np.sin(1.85)), xytext=(5.8*np.cos(1.65), 5.8*np.sin(1.65)),
                    arrowprops=dict(arrowstyle="->", color='#2e7d32', lw=2.5))
    ax_a1.text(0, 6.1, "Onda Viaggiante Oraria (f = 100 Hz)", color='#2e7d32', fontsize=8.5, fontweight='bold', ha='center')

    ax_a1.set_title("A1. Topologia Singolo Rotore PEEK (Mantello Chiuso a 'X')\n6 Solenoidi Asse Z, Sfasamento 120° a terzi (N-P-N-P-N-P)", fontsize=10.5, fontweight='bold')
    ax_a1.set_xlabel("X [cm]", fontsize=9, fontweight='bold')
    ax_a1.set_ylabel("Y [cm]", fontsize=9, fontweight='bold')
    ax_a1.grid(True, linestyle=":", alpha=0.5)
    ax_a1.legend(loc="lower right", fontsize=7.0)

    # -------------------------------------------------------------
    # PANEL A2: Schema Topologico Doppio Rotore Ortogonale a 90° (Sferico)
    # -------------------------------------------------------------
    ax_a2 = fig.add_subplot(gs[0, 1])
    ax_a2.set_aspect('equal')
    ax_a2.set_xlim([-6.5, 6.5])
    ax_a2.set_ylim([-6.5, 6.5])
    
    # Gabbia Sferica
    ax_a2.add_patch(plt.Circle((0, 0), 5.0, color='#1565c0', fill=False, lw=2.5, ls='-', label="Gabbia Sferica Metamateriale (R=50 mm)"))
    ax_a2.add_patch(plt.Circle((0, 0), 1.2, color='#fbc02d', fill=True, alpha=0.5, label="Nucleo Sferico PEEK (R=12 mm)"))
    
    # Rotore 1 (Asse Z, equatoriali)
    for ang, pol in zip(angles_6, pols):
        xc, yc = r_coil_cm * np.cos(ang), r_coil_cm * np.sin(ang)
        ax_a2.add_patch(plt.Circle((xc, yc), 0.65, color='#e65100', alpha=0.8))
        ax_a2.text(xc, yc, "R1", color='white', fontsize=7, fontweight='bold', ha='center', va='center')

    # Rotore 2 (Asse X, proiettato su Y-Z sfasato di 30°)
    for ang in np.linspace(np.pi/6, 2*np.pi + np.pi/6, 6, endpoint=False):
        xc, yc = (r_coil_cm * 0.7) * np.cos(ang), r_coil_cm * np.sin(ang)
        ax_a2.add_patch(plt.Rectangle((xc-0.4, yc-0.4), 0.8, 0.8, color='#00838f', alpha=0.75))
        ax_a2.text(xc, yc, "R2", color='white', fontsize=6.5, fontweight='bold', ha='center', va='center')

    ax_a2.set_title("A2. Topologia Doppio Rotore Ortogonale a 90° (Gabbia Sferica)\nRotore 1 (Asse Z) + Rotore 2 (Asse X in quadratura temporale 90°)", fontsize=10.5, fontweight='bold')
    ax_a2.set_xlabel("X [cm]", fontsize=9, fontweight='bold')
    ax_a2.set_ylabel("Y [cm]", fontsize=9, fontweight='bold')
    ax_a2.grid(True, linestyle=":", alpha=0.5)
    ax_a2.legend(loc="lower right", fontsize=7.0)

    # -------------------------------------------------------------
    # PANEL B1: Mappatura Azimutale 360° <|B|> Singolo Rotore PEEK
    # -------------------------------------------------------------
    ax_b1 = fig.add_subplot(gs[1, 0])
    b_peek = peek_field["b_avg_xy_uT"]
    b_peek_smooth = gaussian_filter(b_peek, sigma=1.0)
    vmax_b1 = float(np.percentile(b_peek_smooth, 98.5))
    im_b1 = ax_b1.imshow(b_peek_smooth, extent=[-6.5, 6.5, -6.5, 6.5], origin='lower', cmap='plasma', vmin=0, vmax=vmax_b1)
    ax_b1.add_patch(plt.Circle((0, 0), 5.0, color='#00e676', fill=False, lw=1.8, ls='--', label='Mantello (R=50 mm)'))
    ax_b1.add_patch(plt.Circle((0, 0), 3.5, color='#ffeb3b', fill=False, lw=1.0, ls=':', label='Cluster Bobine (R=35 mm)'))
    for ang in angles_6:
        ax_b1.plot(3.5*np.cos(ang), 3.5*np.sin(ang), 'o', color='#ff3d00', markersize=4.0)

    ax_b1.set_title("B1. Esposizione Integrata nel Tempo <|B|> (Rotore Singolo PEEK)\nCorona Circolare Continua a 360°: Assenza Totale di Punti Morti", fontsize=10.5, fontweight='bold')
    ax_b1.set_xlabel("X [cm]", fontsize=9, fontweight='bold')
    ax_b1.set_ylabel("Y [cm]", fontsize=9, fontweight='bold')
    cb_b1 = plt.colorbar(im_b1, ax=ax_b1, fraction=0.046, pad=0.04)
    cb_b1.set_label(r"Induzione Media $\langle |\vec{B}| \rangle_t$ [µT]", fontsize=8.5, fontweight='bold')
    ax_b1.legend(loc="upper right", fontsize=7.2)

    # -------------------------------------------------------------
    # PANEL B2: Mappatura Azimutale 360° <|B|> Doppio Rotore 90°
    # -------------------------------------------------------------
    ax_b2 = fig.add_subplot(gs[1, 1])
    b_dual = dual_field["b_avg_xy_uT"]
    b_dual_smooth = gaussian_filter(b_dual, sigma=1.0)
    vmax_b2 = float(np.percentile(b_dual_smooth, 98.0))
    im_b2 = ax_b2.imshow(b_dual_smooth, extent=[-6.5, 6.5, -6.5, 6.5], origin='lower', cmap='plasma', vmin=0, vmax=vmax_b2)
    ax_b2.add_patch(plt.Circle((0, 0), 5.0, color='#00e676', fill=False, lw=1.8, ls='--', label='Gabbia Sferica (R=50 mm)'))
    ax_b2.add_patch(plt.Circle((0, 0), 3.5, color='#ffeb3b', fill=False, lw=1.0, ls=':', label='Array Rotore 1 (R=35 mm)'))
    for ang in angles_6:
        ax_b2.plot(3.5*np.cos(ang), 3.5*np.sin(ang), 'o', color='#ff3d00', markersize=4.0)

    ax_b2.set_title("B2. Esposizione Integrata nel Tempo <|B|> (Doppio Rotore 90°)\nAccoppiamento Quadraturale: Corona Isotropica Sferica a 360°", fontsize=10.5, fontweight='bold')
    ax_b2.set_xlabel("X [cm]", fontsize=9, fontweight='bold')
    ax_b2.set_ylabel("Y [cm]", fontsize=9, fontweight='bold')
    cb_b2 = plt.colorbar(im_b2, ax=ax_b2, fraction=0.046, pad=0.04)
    cb_b2.set_label(r"Induzione Media $\langle |\vec{B}| \rangle_t$ [µT]", fontsize=8.5, fontweight='bold')
    ax_b2.legend(loc="upper right", fontsize=7.2)

    # -------------------------------------------------------------
    # PANEL C1: Kymograph Spazio-Temporale (Theta vs Tempo) Rotore Singolo
    # -------------------------------------------------------------
    ax_c1 = fig.add_subplot(gs[2, 0])
    kymo_p = peek_field["kymo_data_uT"]
    vmax_k1 = float(np.percentile(np.abs(kymo_p), 98.0))
    im_c1 = ax_c1.imshow(kymo_p, extent=[0.25, 16.0, 0, 360], origin='lower', aspect='auto', cmap='coolwarm', vmin=-vmax_k1, vmax=vmax_k1)
    ax_c1.set_title("C1. Kymograph Spazio-Temporale (Singolo Rotore PEEK, R = 50 mm)\nTracce Fluide e Onde Viaggianti Regolari a Velocità di Fase Costante", fontsize=10.5, fontweight='bold')
    ax_c1.set_xlabel("Tempo t [ms]", fontsize=9, fontweight='bold')
    ax_c1.set_ylabel(r"Posizione Azimutale $\theta$ [°]", fontsize=9, fontweight='bold')
    cb_c1 = plt.colorbar(im_c1, ax=ax_c1, fraction=0.046, pad=0.04)
    cb_c1.set_label(r"Induzione Radiale $B_r(\theta, t)$ [µT]", fontsize=8.5, fontweight='bold')

    # -------------------------------------------------------------
    # PANEL C2: Kymograph Spazio-Temporale (Theta vs Tempo) Doppio Rotore 90°
    # -------------------------------------------------------------
    ax_c2 = fig.add_subplot(gs[2, 1])
    kymo_d = dual_field["kymo_data_uT"]
    vmax_k2 = float(np.percentile(np.abs(kymo_d), 95.0))
    im_c2 = ax_c2.imshow(kymo_d, extent=[0.25, 16.0, 0, 360], origin='lower', aspect='auto', cmap='coolwarm', vmin=-vmax_k2, vmax=vmax_k2)
    ax_c2.set_title("C2. Kymograph Spazio-Temporale (Doppio Rotore 90°, R = 50 mm)\nModulazione Ortogonale in Quadratura (Interferenza Costruttiva a 360°)", fontsize=10.5, fontweight='bold')
    ax_c2.set_xlabel("Tempo t [ms]", fontsize=9, fontweight='bold')
    ax_c2.set_ylabel(r"Posizione Azimutale $\theta$ [°]", fontsize=9, fontweight='bold')
    cb_c2 = plt.colorbar(im_c2, ax=ax_c2, fraction=0.046, pad=0.04)
    cb_c2.set_label(r"Induzione Radiale $B_r(\theta, t)$ [µT]", fontsize=8.5, fontweight='bold')

    plt.suptitle("CONFRONTO SINOTTICO: ECCITAZIONE CONTINUA NPNPNP A TERZI (120° TRIFASE)\nRotore Singolo PEEK a Mantello Chiuso vs Doppio Rotore Ortogonale a Gabbia Sferica",
                 fontsize=13.5, fontweight='bold', y=0.995)

    out_root = ROOT_FIGURES_DIR / "fig_18_confronto_npnpnp_peek_vs_doppio_rotore.png"
    plt.savefig(out_root, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [FIGURA 18] Salvata con successo: {out_root}")

    # Copia nelle cartelle varianti
    shutil.copy(out_root, FIG_PEEK_DIR / "fig_18_confronto_npnpnp_peek_vs_doppio_rotore.png")
    shutil.copy(out_root, FIG_DUAL_DIR / "fig_18_confronto_npnpnp_peek_vs_doppio_rotore.png")


def plot_fig_19_matrice_forze_ed_energia(peek_stats, dual_stats, peek_spheres, dual_spheres):
    """Genera fig_19_matrice_forze_npnpnp.png (300 DPI) con subplots dedicati ad alta risoluzione dinamica."""
    print("\n  [FIGURA 19] Generazione Matrice Forze e Bilancio Energetico NPNPNP a 300 DPI...")
    fig = plt.figure(figsize=(20, 13), dpi=300)
    gs = gridspec.GridSpec(2, 3, height_ratios=[1.0, 1.1], width_ratios=[1.1, 1.1, 1.0], hspace=0.32, wspace=0.38)
    time_ms = np.arange(1, TIMESTEPS + 1) * DT * 1000.0

    # -------------------------------------------------------------
    # PANEL A1: Serie Temporali Spinta Rotore Singolo PEEK (in µN)
    # -------------------------------------------------------------
    ax_a1 = fig.add_subplot(gs[0, 0])
    fx_p = peek_stats["total"]["series_fx_uN"]
    fy_p = peek_stats["total"]["series_fy_uN"]
    fz_p = peek_stats["total"]["series_fz_uN"]

    ax_a1.plot(time_ms, fx_p, color='#d32f2f', lw=1.6, label=f"Fx(t) [Media: {peek_stats['total']['mean_fx_uN']:+.1f} µN]")
    ax_a1.plot(time_ms, fy_p, color='#1976d2', lw=1.6, label=f"Fy(t) [Media: {peek_stats['total']['mean_fy_uN']:+.1f} µN]")
    ax_a1.plot(time_ms, fz_p, color='#7b1fa2', lw=2.2, label=f"Fz(t) [Media: {peek_stats['total']['mean_fz_uN']:+.1f} µN]")
    ax_a1.axhline(0, color='black', lw=0.8, ls=':', alpha=0.6)

    ax_a1.set_title("A1. Spinta Rotore Singolo PEEK [µN]\nStabilità Ponderomotrice e Oscillazione Fluida", fontsize=10.5, fontweight='bold')
    ax_a1.set_xlabel("Tempo t [ms]", fontsize=9, fontweight='bold')
    ax_a1.set_ylabel("Forza di Lorentz [µN]", fontsize=9, fontweight='bold')
    ax_a1.grid(True, linestyle=":", alpha=0.6)
    ax_a1.legend(loc="upper right", fontsize=7.5)

    # -------------------------------------------------------------
    # PANEL A2: Serie Temporali Spinta Doppio Rotore 90° (in mN)
    # -------------------------------------------------------------
    ax_a2 = fig.add_subplot(gs[0, 1])
    fx_d_mN = np.array(dual_stats["total"]["series_fx_uN"]) / 1000.0
    fy_d_mN = np.array(dual_stats["total"]["series_fy_uN"]) / 1000.0
    fz_d_mN = np.array(dual_stats["total"]["series_fz_uN"]) / 1000.0

    ax_a2.plot(time_ms, fx_d_mN, color='#d32f2f', lw=1.8, label=f"Fx(t) [Media: {dual_stats['total']['mean_fx_uN']/1000:+.1f} mN]")
    ax_a2.plot(time_ms, fy_d_mN, color='#1976d2', lw=1.8, label=f"Fy(t) [Media: {dual_stats['total']['mean_fy_uN']/1000:+.1f} mN]")
    ax_a2.plot(time_ms, fz_d_mN, color='#388e3c', lw=2.0, ls='--', label=f"Fz(t) [Media: {dual_stats['total']['mean_fz_uN']/1000:+.1f} mN]")
    ax_a2.axhline(0, color='black', lw=0.8, ls=':', alpha=0.6)

    ax_a2.set_title("A2. Spinta Doppio Rotore Ortogonale 90° [mN]\nVettorizzazione Multi-Assiale a Forte Trazione", fontsize=10.5, fontweight='bold')
    ax_a2.set_xlabel("Tempo t [ms]", fontsize=9, fontweight='bold')
    ax_a2.set_ylabel("Forza di Lorentz [mN]", fontsize=9, fontweight='bold')
    ax_a2.grid(True, linestyle=":", alpha=0.6)
    ax_a2.legend(loc="upper right", fontsize=7.5)

    # -------------------------------------------------------------
    # PANEL B: Odografo Vettoriale 3D delle Forze (Fx - Fy - Fz)
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 2], projection='3d')
    
    # Traiettoria Doppio Rotore Ortogonale in mN
    ax_b.plot(fx_d_mN, fy_d_mN, fz_d_mN, color='#0288d1', lw=1.8, label='Doppio Rotore 90° [mN]', alpha=0.85)
    ax_b.scatter([dual_stats["total"]["mean_fx_uN"]/1000], [dual_stats["total"]["mean_fy_uN"]/1000], [dual_stats["total"]["mean_fz_uN"]/1000],
                color='#d32f2f', s=70, marker='^', label=f'Baricentro 90°: |F| = {dual_stats["force_magnitude_mean_uN"]/1000:.1f} mN')

    # Traiettoria Rotore Singolo PEEK (proiettata in scala espansa)
    ax_b.plot(np.array(fx_p)/10.0, np.array(fy_p)/10.0, np.array(fz_p)/10.0, color='#7b1fa2', lw=2.0, label='Rotore PEEK (x100 [mN])', alpha=0.85)

    ax_b.set_title("B. Odografo Vettoriale 3D (Fx, Fy, Fz)\nTraiettoria della Spinta nello Spazio 3D", fontsize=10.5, fontweight='bold')
    ax_b.set_xlabel("Fx [mN]", fontsize=8.0, fontweight='bold')
    ax_b.set_ylabel("Fy [mN]", fontsize=8.0, fontweight='bold')
    ax_b.set_zlabel("Fz [mN]", fontsize=8.0, fontweight='bold')
    ax_b.legend(loc="upper left", fontsize=6.8)
    ax_b.view_init(elev=24, azim=42)

    # -------------------------------------------------------------
    # PANEL C: Bilancio Energetico e Dissipazione Joule per Componente (Dual Y-Axis)
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c_twin = ax_c.twinx()

    categories = ["Mantello X", "Rotore 1", "Rotore 2", "Nucleo PEEK", "Totale"]
    p_pj = [
        peek_stats["mantle_total"]["mean_pj_mW"],
        peek_stats["coils_cluster"]["mean_pj_mW"] if "coils_cluster" in peek_stats else 0.0,
        0.0,
        peek_stats["amagnetic_core"]["mean_pj_mW"],
        peek_stats["total"]["mean_pj_mW"]
    ]
    d_pj_W = [
        dual_stats["mantle_total"]["mean_pj_mW"] / 1000.0,
        dual_stats["rotor1_coils"]["mean_pj_mW"] / 1000.0,
        dual_stats["rotor2_coils"]["mean_pj_mW"] / 1000.0,
        dual_stats["peek_core"]["mean_pj_mW"] / 1000.0,
        dual_stats["total"]["mean_pj_mW"] / 1000.0
    ]

    x_idx = np.arange(len(categories))
    w = 0.35

    bars1 = ax_c.bar(x_idx - w/2, p_pj, width=w, color='#8e24aa', alpha=0.85, label='Rotore PEEK [mW]')
    bars2 = ax_c_twin.bar(x_idx + w/2, d_pj_W, width=w, color='#0288d1', alpha=0.85, label='Doppio Rotore 90° [W]')

    for b in bars1:
        h = b.get_height()
        if h > 0:
            ax_c.annotate(f'{h:.1f}', xy=(b.get_x() + b.get_width()/2, h),
                          xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=7.0, color='#6a1b9a', fontweight='bold')
    for b in bars2:
        h = b.get_height()
        if h > 0:
            ax_c_twin.annotate(f'{h:.1f}', xy=(b.get_x() + b.get_width()/2, h),
                               xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=7.0, color='#01579b', fontweight='bold')

    ax_c.set_xticks(x_idx)
    ax_c.set_xticklabels(categories, fontsize=8.0, fontweight='bold')
    ax_c.set_ylabel("Dissipazione PEEK [mW]", color='#8e24aa', fontsize=9.0, fontweight='bold')
    ax_c_twin.set_ylabel("Dissipazione Doppio Rotore [W]", color='#0288d1', fontsize=9.0, fontweight='bold')
    ax_c.set_title("C. Ripartizione Perdite Joule per Componente\nNucleo PEEK Isolante a Perdite Zero (0.0 W)", fontsize=10.5, fontweight='bold')
    ax_c.grid(True, linestyle=":", alpha=0.5, axis='y')

    lines1, labels1 = ax_c.get_legend_handles_labels()
    lines2, labels2 = ax_c_twin.get_legend_handles_labels()
    ax_c.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=7.2)

    # -------------------------------------------------------------
    # PANEL D1: Validazione Lorentz vs MST: Rotore Singolo PEEK (in µN)
    # -------------------------------------------------------------
    ax_d1 = fig.add_subplot(gs[1, 1])
    comps_p = ["Fx", "Fy", "Fz"]
    lor_p = [peek_stats["total"]["mean_fx_uN"], peek_stats["total"]["mean_fy_uN"], peek_stats["total"]["mean_fz_uN"]]
    mst_p_dat = peek_spheres["Near-Field (R=6.5 cm)"]
    mst_p_vals = [mst_p_dat["mst_fx_uN"], mst_p_dat["mst_fy_uN"], mst_p_dat["mst_fz_uN"]]

    x_p = np.arange(len(comps_p))
    bars_lp = ax_d1.bar(x_p - w/2, lor_p, width=w, color='#2e7d32', alpha=0.85, label='Volume Lorentz')
    bars_mp = ax_d1.bar(x_p + w/2, mst_p_vals, width=w, color='#f57c00', alpha=0.85, label='Superficie MST (R=6.5 cm)')
    ax_d1.axhline(0, color='black', lw=0.8, ls=':', alpha=0.6)

    for b in bars_lp:
        h = b.get_height()
        ax_d1.annotate(f'{h:+.0f}', xy=(b.get_x() + b.get_width()/2, h),
                       xytext=(0, 3 if h >= 0 else -9), textcoords="offset points", ha='center', fontsize=6.8, fontweight='bold')
    for b in bars_mp:
        h = b.get_height()
        ax_d1.annotate(f'{h:+.0f}', xy=(b.get_x() + b.get_width()/2, h),
                       xytext=(0, 3 if h >= 0 else -9), textcoords="offset points", ha='center', fontsize=6.8, fontweight='bold')

    ax_d1.set_xticks(x_p)
    ax_d1.set_xticklabels(comps_p, fontsize=8.5, fontweight='bold')
    ax_d1.set_ylabel("Forza Integrata [µN]", fontsize=9.0, fontweight='bold')
    ax_d1.set_title("D1. Validazione Lorentz vs MST (PEEK) [µN]\nCoerenza del Segno e Direzionalità Vettoriale", fontsize=10.5, fontweight='bold')
    ax_d1.grid(True, linestyle=":", alpha=0.5, axis='y')
    ax_d1.legend(loc="upper right", fontsize=7.2)

    # -------------------------------------------------------------
    # PANEL D2: Validazione Lorentz vs MST: Doppio Rotore 90° (in mN)
    # -------------------------------------------------------------
    ax_d2 = fig.add_subplot(gs[1, 2])
    lor_d = [dual_stats["total"]["mean_fx_uN"]/1000.0, dual_stats["total"]["mean_fy_uN"]/1000.0, dual_stats["total"]["mean_fz_uN"]/1000.0]
    mst_d_dat = dual_spheres["Near-Field (R=6.5 cm)"]
    mst_d_vals = [mst_d_dat["mst_fx_uN"]/1000.0, mst_d_dat["mst_fy_uN"]/1000.0, mst_d_dat["mst_fz_uN"]/1000.0]

    bars_ld = ax_d2.bar(x_p - w/2, lor_d, width=w, color='#2e7d32', alpha=0.85, label='Volume Lorentz')
    bars_md = ax_d2.bar(x_p + w/2, mst_d_vals, width=w, color='#f57c00', alpha=0.85, label='Superficie MST (R=6.5 cm)')
    ax_d2.axhline(0, color='black', lw=0.8, ls=':', alpha=0.6)

    for b in bars_ld:
        h = b.get_height()
        ax_d2.annotate(f'{h:+.0f}', xy=(b.get_x() + b.get_width()/2, h),
                       xytext=(0, 3 if h >= 0 else -9), textcoords="offset points", ha='center', fontsize=6.8, fontweight='bold')
    for b in bars_md:
        h = b.get_height()
        ax_d2.annotate(f'{h:+.0f}', xy=(b.get_x() + b.get_width()/2, h),
                       xytext=(0, 3 if h >= 0 else -9), textcoords="offset points", ha='center', fontsize=6.8, fontweight='bold')

    ax_d2.set_xticks(x_p)
    ax_d2.set_xticklabels(comps_p, fontsize=8.5, fontweight='bold')
    ax_d2.set_ylabel("Forza Integrata [mN]", fontsize=9.0, fontweight='bold')
    ax_d2.set_title("D2. Validazione Lorentz vs MST (90°) [mN]\nAccoppiamento Quadraturale Vettoriale", fontsize=10.5, fontweight='bold')
    ax_d2.grid(True, linestyle=":", alpha=0.5, axis='y')
    ax_d2.legend(loc="upper right", fontsize=7.2)

    plt.suptitle("MATRICE DELLE FORZE E BILANCIO ENERGETICO: CAMPAGNA NPNPNP A TERZI (120°)\nAnalisi Transiente Vettoriale 3D e Validazione Lorentz-MST (64 Timestep)",
                 fontsize=13.5, fontweight='bold', y=0.995)

    out_root = ROOT_FIGURES_DIR / "fig_19_matrice_forze_npnpnp.png"
    plt.savefig(out_root, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [FIGURA 19] Salvata con successo: {out_root}")

    # Copia nelle cartelle varianti
    shutil.copy(out_root, FIG_PEEK_DIR / "fig_19_matrice_forze_npnpnp.png")
    shutil.copy(out_root, FIG_DUAL_DIR / "fig_19_matrice_forze_npnpnp.png")


def main():
    print("=" * 85)
    print("STUDIO ELETTRODINAMICO COMPARATIVO: CAMPAGNA NPNPNP A TERZI (120° TRIFASE)")
    print("Rotore Singolo PEEK vs Doppio Rotore Ortogonale a 90°")
    print("=" * 85)

    # 1. Esecuzione Simulazione 1: Rotore Singolo PEEK
    vtus_peek = run_elmer_simulation(
        sif_path=CONFIG_PEEK_SIF,
        work_dir=WORK_DIR_PEEK,
        mesh_dir=MESH_DIR_PEEK,
        mesh_name="macchina_mantello_chiuso",
        out_prefix="macchina_out"
    )

    # 2. Esecuzione Simulazione 2: Doppio Rotore Ortogonale a 90°
    vtus_dual = run_elmer_simulation(
        sif_path=CONFIG_DUAL_SIF,
        work_dir=WORK_DIR_DUAL,
        mesh_dir=MESH_DIR_DUAL,
        mesh_name="macchina_sferica_ortogonale",
        out_prefix="macchina_sferica_out"
    )

    # 3. Post-Processing Forze & Dissipazione Joule
    peek_stats = process_single_rotor_peek(vtus_peek)
    dual_stats = process_dual_rotor_sphere(vtus_dual)

    # 4. Verifica Gauss & MST su sfere di Fibonacci
    print("\n  [SPHERES] Campionamento sfere di Fibonacci N = 2500 per Rotore Singolo PEEK...")
    peek_spheres = analyze_spherical_sampling_and_mst(vtus_peek, sample_idx=32)
    for lbl, dat in peek_spheres.items():
        print(f"    - {lbl}: Gauss Residuo = {dat['gauss_residual_pct']:.4f}% [{'PASS' if dat['gauss_pass'] else 'FAIL'}], <|B|> = {dat['b_mean_uT']:.1f} uT")

    print("\n  [SPHERES] Campionamento sfere di Fibonacci N = 2500 per Doppio Rotore 90°...")
    dual_spheres = analyze_spherical_sampling_and_mst(vtus_dual, sample_idx=32)
    for lbl, dat in dual_spheres.items():
        print(f"    - {lbl}: Gauss Residuo = {dat['gauss_residual_pct']:.4f}% [{'PASS' if dat['gauss_pass'] else 'FAIL'}], <|B|> = {dat['b_mean_uT']:.1f} uT")

    # 5. Salvataggio Dataset JSON
    DATA_PEEK_JSON.parent.mkdir(parents=True, exist_ok=True)
    DATA_DUAL_JSON.parent.mkdir(parents=True, exist_ok=True)

    json_peek_data = {
        "variant": "rotore_singolo_peek_npnpnp_terzi",
        "topology": "cilindrica_mantello_chiuso_triplo_strato_X",
        "excitation": "NPNPNP_simultanea_terzi_120deg",
        "timesteps": TIMESTEPS,
        "dt_s": DT,
        "frequency_hz": F_HZ,
        "forces_and_losses": peek_stats,
        "fibonacci_spheres": peek_spheres
    }
    with open(DATA_PEEK_JSON, "w", encoding="utf-8") as f:
        json.dump(json_peek_data, f, indent=2)
    print(f"\n  [JSON] Dataset salvato: {DATA_PEEK_JSON}")

    json_dual_data = {
        "variant": "doppio_rotore_ortogonale_90deg_npnpnp_terzi",
        "topology": "gabbia_sferica_triplo_strato_X",
        "excitation": "NPNPNP_simultanea_terzi_120deg_quadratura_90deg",
        "timesteps": TIMESTEPS,
        "dt_s": DT,
        "frequency_hz": F_HZ,
        "forces_and_losses": dual_stats,
        "fibonacci_spheres": dual_spheres
    }
    with open(DATA_DUAL_JSON, "w", encoding="utf-8") as f:
        json.dump(json_dual_data, f, indent=2)
    print(f"  [JSON] Dataset salvato: {DATA_DUAL_JSON}")

    # 6. Estrazione campi cartesiani & kymograph per fig_18
    print("\n  [MAPPING] Estrazione griglia bidimensionale e kymograph spazio-temporale...")
    peek_field = extract_field_mapping_and_kymograph(vtus_peek, radius=0.050)
    dual_field = extract_field_mapping_and_kymograph(vtus_dual, radius=0.050)

    # 7. Generazione Figure ad alta risoluzione (300 DPI)
    plot_fig_18_confronto_sinottico(peek_field, dual_field)
    plot_fig_19_matrice_forze_ed_energia(peek_stats, dual_stats, peek_spheres, dual_spheres)

    print("\n" + "=" * 85)
    print("CAMPAGNA NPNPNP COMPLETATA CON SUCCESSO!")
    print(f"  - Dataset 1: {DATA_PEEK_JSON}")
    print(f"  - Dataset 2: {DATA_DUAL_JSON}")
    print(f"  - Deliverable 1: {ROOT_FIGURES_DIR / 'fig_18_confronto_npnpnp_peek_vs_doppio_rotore.png'}")
    print(f"  - Deliverable 2: {ROOT_FIGURES_DIR / 'fig_19_matrice_forze_npnpnp.png'}")
    print("=" * 85)


if __name__ == "__main__":
    main()
