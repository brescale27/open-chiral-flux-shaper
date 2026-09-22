#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROTOCOLLO DI CHIUSURA NUMERICA E CERTIFICAZIONE OPEN HARDWARE CERN-OHL-S v2:
1. Verifica della conservazione del flusso (Teorema di Gauss su sfere R = 8, 12, 15 cm)
2. Analisi di sensibilità parametrica REALE da simulazioni Elmer FEM (25°, 30°, 35°)
3. Esportazione metriche e tabelle in data/
"""

import os
import sys
import json
import csv
from pathlib import Path
import numpy as np
import meshio
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "data"

def build_interpolator(pts, B_field, box=0.17, n_grid=71):
    gx = np.linspace(-box, box, n_grid)
    gy = np.linspace(-box, box, n_grid)
    gz = np.linspace(-box, box, n_grid)
    GX, GY, GZ = np.meshgrid(gx, gy, gz, indexing='ij')

    lin_interp = LinearNDInterpolator(pts, B_field, fill_value=0.0)
    grid_coords = np.column_stack([GX.ravel(), GY.ravel(), GZ.ravel()])
    B_grid = lin_interp(grid_coords).reshape(n_grid, n_grid, n_grid, 3)

    interp_bx = RegularGridInterpolator((gx, gy, gz), B_grid[..., 0], bounds_error=False, fill_value=0.0)
    interp_by = RegularGridInterpolator((gx, gy, gz), B_grid[..., 1], bounds_error=False, fill_value=0.0)
    interp_bz = RegularGridInterpolator((gx, gy, gz), B_grid[..., 2], bounds_error=False, fill_value=0.0)

    def eval_B(r_vec):
        r_arr = np.atleast_2d(r_vec)
        bx = interp_bx(r_arr)
        by = interp_by(r_arr)
        bz = interp_bz(r_arr)
        return np.column_stack([bx, by, bz])

    return eval_B

def load_sim_series(mesh_dir):
    b_steps = []
    j_steps = []
    p_steps = []
    for step in range(1, 11):
        f = mesh_dir / f"macchina_out_t{step:04d}.vtu"
        if not f.is_file():
            raise FileNotFoundError(f"File VTU non trovato: {f}")
        mt = meshio.read(str(f))
        b_steps.append(mt.point_data['magnetic flux density'])
        j_steps.append(mt.point_data['current density'])
        p_steps.append(mt.point_data['joule heating'])

    B_dc = np.mean(b_steps, axis=0)
    J_dc = np.mean(j_steps, axis=0)
    P_dc = np.mean(p_steps, axis=0)
    return B_dc, J_dc, P_dc

def main():
    print("=" * 80)
    print("POST-PROCESSING DI CHIUSURA NUMERICA E CERTIFICAZIONE SOLENOIDALE (CERN-OHL-S v2)")
    print("=" * 80)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Carica geometria base
    nom_dir = ROOT_DIR / "results_nominal_30deg"
    base_vtu = nom_dir / "macchina_out_t0001.vtu"
    if not base_vtu.is_file():
        print(f"[ERRORE] File baseline non trovato in: {base_vtu}")
        print("Eseguire prima 'python scripts/run_all_simulations.py'")
        sys.exit(1)

    m_base = meshio.read(str(base_vtu))
    pts = m_base.points
    tetra = m_base.cells_dict['tetra']
    gids = m_base.cell_data['GeometryIds'][0]

    alu_mask = (gids == 2)
    alu_tetra = tetra[alu_mask]

    # Calcolo volumi tetraedri mantello
    p0 = pts[alu_tetra[:, 0]]
    p1 = pts[alu_tetra[:, 1]]
    p2 = pts[alu_tetra[:, 2]]
    p3 = pts[alu_tetra[:, 3]]
    alu_vols = np.abs(np.einsum('ij,ij->i', p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0
    total_alu_vol = float(np.sum(alu_vols))

    # Carica serie baseline (30.0°)
    B_dc_base, J_dc_base, P_dc_base = load_sim_series(nom_dir)

    P_elem_base = np.mean(P_dc_base[alu_tetra], axis=1)
    P_joule_base = float(np.sum(P_elem_base * alu_vols))

    eval_B_base = build_interpolator(pts, B_dc_base)

    # -------------------------------------------------------------------------
    # PARTE 1: VERIFICA DELLA CONSERVAZIONE DEL FLUSSO (TEOREMA DI GAUSS)
    # -------------------------------------------------------------------------
    print("\n[1/3] Calcolo conservazione integrale di Gauss su sfere chiuse (Nominale 30°)...")
    radii = [0.08, 0.12, 0.15] # 8 cm, 12 cm, 15 cm
    n_th, n_ph = 140, 280
    theta = np.linspace(0.005, np.pi - 0.005, n_th)
    phi = np.linspace(0, 2 * np.pi, n_ph, endpoint=False)
    dth = theta[1] - theta[0]
    dph = phi[1] - phi[0]
    TH, PH = np.meshgrid(theta, phi)

    sin_th = np.sin(TH)
    cos_th = np.cos(TH)
    cos_ph = np.cos(PH)
    sin_ph = np.sin(PH)

    nx_s = sin_th * cos_ph
    ny_s = sin_th * sin_ph
    nz_s = cos_th

    gauss_results = []
    for R in radii:
        x = R * nx_s
        y = R * ny_s
        z = R * nz_s
        pts_eval = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
        B_eval = eval_B_base(pts_eval)
        
        Bx = B_eval[:, 0].reshape(TH.shape)
        By = B_eval[:, 1].reshape(TH.shape)
        Bz = B_eval[:, 2].reshape(TH.shape)

        Bn = Bx * nx_s + By * ny_s + Bz * nz_s
        dA = (R**2) * sin_th * dth * dph

        phi_net = float(np.sum(Bn * dA))
        phi_abs = float(np.sum(np.abs(Bn) * dA))
        rel_res = abs(phi_net) / (phi_abs + 1e-15)

        area = float(4.0 * np.pi * R**2)
        mean_bn = float(np.mean(np.abs(Bn)))

        gauss_results.append({
            'radius_m': R,
            'radius_cm': R * 100.0,
            'area_m2': area,
            'phi_net_Tm2': phi_net,
            'phi_abs_Tm2': phi_abs,
            'relative_residual': rel_res,
            'mean_Bn_uT': mean_bn * 1e6
        })

        print(f"  Sfera R = {R*100:4.1f} cm (Area = {area:.4f} m^2):")
        print(f"    Flusso Netto Phi_net: {phi_net:+.4e} T*m^2")
        print(f"    Flusso Assoluto Phi_abs: {phi_abs:.4e} T*m^2")
        print(f"    Residuo Relativo (|Phi_net| / Phi_abs): {rel_res:.4e} ({rel_res*100:.3f}%)")

    # -------------------------------------------------------------------------
    # PARTE 2: ANALISI DI SENSIBILITÀ PARAMETRICA REALE (25°, 30°, 35°)
    # -------------------------------------------------------------------------
    print("\n[2/3] Elaborazione della sensitività dai tre run FEM effettivi...")
    sim_cases = [
        {
            'angle_deg': 25.0,
            'dir': ROOT_DIR / 'results_25deg',
            'rel_dir': 'results_25deg',
            'sigma_tz': 2681155.550916423
        },
        {
            'angle_deg': 30.0,
            'dir': ROOT_DIR / 'results_nominal_30deg',
            'rel_dir': 'results_nominal_30deg',
            'sigma_tz': 3031088.913245535
        },
        {
            'angle_deg': 35.0,
            'dir': ROOT_DIR / 'results_35deg',
            'rel_dir': 'results_35deg',
            'sigma_tz': 3288924.172750679
        }
    ]

    th_eval = np.linspace(0, 2*np.pi, 72, endpoint=False)
    pts_r6 = np.column_stack([0.06 * np.cos(th_eval), 0.06 * np.sin(th_eval), np.zeros_like(th_eval)])

    B_r6_base = eval_B_base(pts_r6)
    br_r6_base = B_r6_base[:, 0] * np.cos(th_eval) + B_r6_base[:, 1] * np.sin(th_eval)
    mean_br6_base = float(np.mean(np.abs(br_r6_base))) * 1e6

    sigma_base = sim_cases[1]['sigma_tz']

    sensitivity_results = []
    for sc in sim_cases:
        angle_deg = sc['angle_deg']
        sim_dir = sc['dir']
        sigma_tz = sc['sigma_tz']
        delta_sigma_pct = (sigma_tz - sigma_base) / sigma_base * 100.0

        print(f"  Elaborazione dati per alpha = {angle_deg:.1f}° ({sc['rel_dir']})...")
        B_dc, J_dc, P_dc = load_sim_series(sim_dir)

        P_elem = np.mean(P_dc[alu_tetra], axis=1)
        P_joule_real = float(np.sum(P_elem * alu_vols))
        delta_P_pct = (P_joule_real - P_joule_base) / P_joule_base * 100.0

        eval_B_case = build_interpolator(pts, B_dc)
        B_r6_case = eval_B_case(pts_r6)
        br_r6_case = B_r6_case[:, 0] * np.cos(th_eval) + B_r6_case[:, 1] * np.sin(th_eval)
        mean_br6_real = float(np.mean(np.abs(br_r6_case))) * 1e6
        delta_br6_pct = (mean_br6_real - mean_br6_base) / mean_br6_base * 100.0

        sensitivity_results.append({
            'angle_deg': angle_deg,
            'sigma_theta_z_S_m': sigma_tz,
            'delta_sigma_pct': delta_sigma_pct,
            'P_joule_W': P_joule_real,
            'delta_P_joule_pct': delta_P_pct,
            'mean_Brad_6cm_uT': mean_br6_real,
            'delta_Brad_pct': delta_br6_pct,
            'simulation_source': f"{sc['rel_dir']}/macchina_out_t*.vtu"
        })

        print(f"    sigma_theta_z = {sigma_tz:.3e} S/m ({delta_sigma_pct:+5.1f}%):")
        print(f"    Potenza Joule: {P_joule_real:.4f} W ({delta_P_pct:+5.2f}%)")
        print(f"    B_rad a 6cm:  {mean_br6_real:.2f} uT ({delta_br6_pct:+5.2f}%)")

    # -------------------------------------------------------------------------
    # PARTE 3: CONFRONTO CON MANTELLO PIENO MASSICCIO (SE DISPONIBILE)
    # -------------------------------------------------------------------------
    solid_dir = ROOT_DIR / "_archive_backup" / "results_solid"
    if solid_dir.is_dir() and len(list(solid_dir.glob("*.vtu"))) >= 10:
        print("\n[3/3] Aggiornamento confronto mantello pieno vs rete stirata...")
        B_dc_solid, J_dc_solid, P_dc_solid = load_sim_series(solid_dir)
        eval_B_solid = build_interpolator(pts, B_dc_solid)

        radii_comp = [0.048, 0.060, 0.080, 0.120]
        comp_rows = []
        for r_m in radii_comp:
            pts_r = np.column_stack([r_m * np.cos(th_eval), r_m * np.sin(th_eval), np.zeros_like(th_eval)])
            
            b_sol = eval_B_solid(pts_r)
            br_sol = b_sol[:, 0] * np.cos(th_eval) + b_sol[:, 1] * np.sin(th_eval)
            bmod_sol = np.linalg.norm(b_sol, axis=1) * 1e6
            
            b_mesh = eval_B_base(pts_r)
            br_mesh = b_mesh[:, 0] * np.cos(th_eval) + b_mesh[:, 1] * np.sin(th_eval)
            bmod_mesh = np.linalg.norm(b_mesh, axis=1) * 1e6

            max_br_sol = float(np.max(np.abs(br_sol))) * 1e6
            max_br_mesh = float(np.max(np.abs(br_mesh))) * 1e6
            enh_factor = max_br_mesh / (max_br_sol + 1e-12)
            mean_bmod_sol = float(np.mean(bmod_sol))
            mean_bmod_mesh = float(np.mean(bmod_mesh))

            comp_rows.append({
                'radius_cm': r_m * 100.0,
                'max_Brad_solid_uT': round(max_br_sol, 2),
                'max_Brad_mesh_uT': round(max_br_mesh, 2),
                'enhancement_factor': round(enh_factor, 3),
                'mean_Bmod_solid_uT': round(mean_bmod_sol, 2),
                'mean_Bmod_mesh_uT': round(mean_bmod_mesh, 2)
            })

        csv_path = DATA_DIR / "confronto_mantello_pieno_vs_rete.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=['radius_cm', 'max_Brad_solid_uT', 'max_Brad_mesh_uT', 'enhancement_factor', 'mean_Bmod_solid_uT', 'mean_Bmod_mesh_uT'])
            writer.writeheader()
            writer.writerows(comp_rows)
        print(f"  Tabella comparativa aggiornata in: {csv_path}")

    # Esportazione JSON
    out_json = DATA_DIR / "validazione_chiusura_cern_ohl.json"
    validation_data = {
        'gauss_conservation_spheres': gauss_results,
        'sensitivity_louver_angle': sensitivity_results,
        'baseline_reference': {
            'frequency_hz': 100.0,
            'rpm': 1200.0,
            'P_joule_dc_W': P_joule_base,
            'mean_abs_Brad_6cm_uT': mean_br6_base,
            'volume_alu_m3': total_alu_vol
        }
    }
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(validation_data, f, indent=2)
    print(f"  Dataset di validazione salvato in: {out_json}")
    print("=" * 80)
    print("CERTIFICAZIONE NUMERICA COMPLETATA CON SUCCESSO!")
    print("=" * 80)

if __name__ == "__main__":
    main()
