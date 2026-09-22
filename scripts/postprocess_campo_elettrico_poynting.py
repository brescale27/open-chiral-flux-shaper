#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modulo Multifisico Avanzato: Post-processing Sfasamento, E, Poynting, Harvesting e AM (CERN-OHL-S v2)
1. Analisi comparativa dei 4 regimi (A, B, C, D): B_rad, Ripple Index, perdite Joule P_J.
2. Mappatura campo elettrico E(r, t) e flusso vettoriale di Poynting S(r, t) = (1/mu0) (E x B).
3. Accoppiamento a distanza e virtual harvesting su 3 sonde virtuali (radiale, assiale, capacitiva).
4. Caratterizzazione del respiro dinamico da modulazione AM a 10 Hz.
5. Esportazione data/sweep_sfasamento_risultati.json e generazione figures/04, 05, 06 a 300 DPI.
"""

import os
import sys
import json
from pathlib import Path
import numpy as np
import meshio
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "data"
FIG_DIR = ROOT_DIR / "figures"
MU0 = 4.0 * np.pi * 1e-7
EPS0 = 8.8541878128e-12

def load_sim_series(dir_path, n_steps=10):
    b_steps = []
    e_steps = []
    j_steps = []
    p_steps = []
    for step in range(1, n_steps + 1):
        f = dir_path / f"macchina_out_t{step:04d}.vtu"
        if not f.is_file():
            raise FileNotFoundError(f"File non trovato: {f}")
        m = meshio.read(str(f))
        b_steps.append(m.point_data['magnetic flux density'])
        e_steps.append(m.point_data['electric field'])
        j_steps.append(m.point_data['current density'])
        p_steps.append(m.point_data['joule heating'])
    return np.array(b_steps), np.array(e_steps), np.array(j_steps), np.array(p_steps)

def build_interpolator_3d(pts, field_data, box_xy=0.18, box_z=0.15, n_grid=51):
    gx = np.linspace(-box_xy, box_xy, n_grid)
    gy = np.linspace(-box_xy, box_xy, n_grid)
    gz = np.linspace(-box_z, box_z, n_grid)
    GX, GY, GZ = np.meshgrid(gx, gy, gz, indexing='ij')

    lin_interp = LinearNDInterpolator(pts, field_data, fill_value=0.0)
    grid_coords = np.column_stack([GX.ravel(), GY.ravel(), GZ.ravel()])
    grid_vals = lin_interp(grid_coords).reshape(n_grid, n_grid, n_grid, 3)

    ix = RegularGridInterpolator((gx, gy, gz), grid_vals[..., 0], bounds_error=False, fill_value=0.0)
    iy = RegularGridInterpolator((gx, gy, gz), grid_vals[..., 1], bounds_error=False, fill_value=0.0)
    iz = RegularGridInterpolator((gx, gy, gz), grid_vals[..., 2], bounds_error=False, fill_value=0.0)

    def eval_vec(r_vec):
        r_arr = np.atleast_2d(r_vec)
        vx = ix(r_arr)
        vy = iy(r_arr)
        vz = iz(r_arr)
        return np.column_stack([vx, vy, vz])

    return eval_vec

def main():
    print("=" * 80)
    print("POST-PROCESSING MULTIFISICO AVANZATO: SFASAMENTO, POYNTING, HARVESTING E AM")
    print("=" * 80)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    nom_vtu = ROOT_DIR / "results_regime_A" / "macchina_out_t0001.vtu"
    if not nom_vtu.is_file():
        print(f"[ERRORE] File baseline non trovato: {nom_vtu}")
        sys.exit(1)

    m_base = meshio.read(str(nom_vtu))
    pts = m_base.points
    tetra = m_base.cells_dict['tetra']
    gids = m_base.cell_data['GeometryIds'][0]

    alu_mask = (gids == 2)
    alu_tetra = tetra[alu_mask]
    p0 = pts[alu_tetra[:, 0]]
    p1 = pts[alu_tetra[:, 1]]
    p2 = pts[alu_tetra[:, 2]]
    p3 = pts[alu_tetra[:, 3]]
    alu_vols = np.abs(np.einsum('ij,ij->i', p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0
    tot_vol = np.sum(alu_vols)

    regimes = [
        ("Regime A (Sincrono 0°)", "results_regime_A", "tab:blue"),
        ("Regime B (Co-rotante 60°)", "results_regime_B", "tab:red"),
        ("Regime C (Quadratura 90°/180°)", "results_regime_C", "tab:green"),
        ("Regime D (Contro-rotante 60°)", "results_regime_D", "tab:purple")
    ]

    dt_step = 0.001 # 1.0 ms
    n_steps = 10
    time_arr = np.arange(1, n_steps + 1) * dt_step * 1000.0 # ms

    n_theta = 72
    theta_arr = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)

    results_dict = {"regimes": {}, "virtual_harvesting": {}, "poynting_flux": {}, "am_modulation": {}}
    data_all_regimes = {}

    for name, folder, col in regimes:
        dir_path = ROOT_DIR / folder
        print(f"\n[Elaborazione] {name} ({folder})...")
        B_all, E_all, J_all, P_all = load_sim_series(dir_path, n_steps=10)
        data_all_regimes[name] = {"B": B_all, "E": E_all, "J": J_all, "P": P_all, "col": col}

        # Calcolo perdite Joule nel mantello
        P_joule_t = []
        for step in range(n_steps):
            p_elem = np.mean(P_all[step][alu_tetra], axis=1)
            P_joule_t.append(float(np.sum(p_elem * alu_vols)))
        P_joule_t = np.array(P_joule_t)
        P_joule_mean = float(np.mean(P_joule_t))

        B_mean_t = np.mean(B_all, axis=0)
        eval_B_dc = build_interpolator_3d(pts, B_mean_t)

        # Campionamento temporale a R = 6 cm
        brad_t_series = []
        for step in range(n_steps):
            eval_B_step = build_interpolator_3d(pts, B_all[step])
            pts_6cm = np.column_stack([0.06 * np.cos(theta_arr), 0.06 * np.sin(theta_arr), np.zeros(n_theta)])
            B_eval = eval_B_step(pts_6cm)
            brad = B_eval[:, 0] * np.cos(theta_arr) + B_eval[:, 1] * np.sin(theta_arr)
            brad_t_series.append(np.mean(np.abs(brad)))
        brad_t_series = np.array(brad_t_series)

        brad_mean_6cm = float(np.mean(brad_t_series))
        b_max = float(np.max(brad_t_series))
        b_min = float(np.min(brad_t_series))
        ripple_idx = float((b_max - b_min) / (brad_mean_6cm + 1e-12))

        pts_8cm = np.column_stack([0.08 * np.cos(theta_arr), 0.08 * np.sin(theta_arr), np.zeros(n_theta)])
        B_8cm = eval_B_dc(pts_8cm)
        brad_8cm = B_8cm[:, 0] * np.cos(theta_arr) + B_8cm[:, 1] * np.sin(theta_arr)
        brad_mean_8cm = float(np.mean(np.abs(brad_8cm)))

        results_dict["regimes"][name] = {
            "mean_Brad_6cm_uT": brad_mean_6cm * 1e6,
            "max_Brad_6cm_uT": b_max * 1e6,
            "min_Brad_6cm_uT": b_min * 1e6,
            "ripple_index": ripple_idx,
            "mean_Brad_8cm_uT": brad_mean_8cm * 1e6,
            "mean_P_joule_W": P_joule_mean,
            "P_joule_t_W": P_joule_t.tolist(),
            "Brad_t_series_uT": (brad_t_series * 1e6).tolist()
        }
        print(f"  -> <|Brad|> 6cm: {brad_mean_6cm*1e6:.2f} uT | Ripple: {ripple_idx*100:.2f}% | P_joule: {P_joule_mean:.3f} W")

    # Poynting & E-field per Regime B
    print("\n[Mappatura Vettore di Poynting ed E-field (Regime B)]...")
    opt_B_all = data_all_regimes["Regime B (Co-rotante 60°)"]["B"]
    opt_E_all = data_all_regimes["Regime B (Co-rotante 60°)"]["E"]

    S_all = np.cross(opt_E_all, opt_B_all) / MU0
    S_mean_nodal = np.mean(S_all, axis=0)
    E_mean_nodal = np.mean(opt_E_all, axis=0)
    B_mean_nodal = np.mean(opt_B_all, axis=0)

    eval_S = build_interpolator_3d(pts, S_mean_nodal)
    eval_E = build_interpolator_3d(pts, E_mean_nodal)
    eval_B_dc_opt = build_interpolator_3d(pts, B_mean_nodal)

    # Flusso di Poynting attraverso cilindro di raggio 12 cm, H = 20 cm
    r_cyl = 0.12
    h_cyl = 0.10
    nz_cyl = 25
    nth_cyl = 36
    z_cyl_arr = np.linspace(-h_cyl, h_cyl, nz_cyl)
    th_cyl_arr = np.linspace(0, 2*np.pi, nth_cyl, endpoint=False)
    TH_c, Z_c = np.meshgrid(th_cyl_arr, z_cyl_arr)
    xc = r_cyl * np.cos(TH_c).ravel()
    yc = r_cyl * np.sin(TH_c).ravel()
    zc = Z_c.ravel()
    pts_cyl = np.column_stack([xc, yc, zc])
    S_cyl_eval = eval_S(pts_cyl)
    nr_x = np.cos(TH_c).ravel()
    nr_y = np.sin(TH_c).ravel()
    S_rad = S_cyl_eval[:, 0] * nr_x + S_cyl_eval[:, 1] * nr_y
    dA_cyl = (2.0 * np.pi * r_cyl / nth_cyl) * (2.0 * h_cyl / (nz_cyl - 1))
    P_rad_out = float(np.sum(S_rad) * dA_cyl)

    results_dict["poynting_flux"] = {
        "cylinder_radius_m": r_cyl,
        "cylinder_height_m": 2.0 * h_cyl,
        "net_outward_Poynting_power_mW": P_rad_out * 1000.0,
        "mean_S_magnitude_W_m2": float(np.mean(np.linalg.norm(S_cyl_eval, axis=1)))
    }
    print(f"  -> Potenza netta uscente di Poynting a R=12cm: {P_rad_out*1000.0:.3f} mW")

    # Virtual Probes & Remote Harvesting
    print("\n[Virtual Probes & Remote Harvesting]...")
    r_coil = 0.025
    area_coil = np.pi * r_coil * r_coil
    N_turns = 100
    R_load_opt = 1.34 # Ohm

    phi_sonda1 = []
    phi_sonda2 = []
    er_sonda3 = []

    for step in range(n_steps):
        eval_B_s = build_interpolator_3d(pts, opt_B_all[step])
        eval_E_s = build_interpolator_3d(pts, opt_E_all[step])

        b_s1 = eval_B_s(np.array([[0.10, 0.0, 0.0]]))
        phi_sonda1.append(float(b_s1[0, 0] * area_coil))

        b_s2 = eval_B_s(np.array([[0.15, 0.0, 0.10]]))
        phi_sonda2.append(float(b_s2[0, 2] * area_coil))

        e_s3 = eval_E_s(np.array([[0.12, 0.0, 0.0]]))
        er_sonda3.append(float(e_s3[0, 0]))

    phi_sonda1 = np.array(phi_sonda1)
    phi_sonda2 = np.array(phi_sonda2)
    er_sonda3 = np.array(er_sonda3)

    dphi1_dt = np.gradient(phi_sonda1, dt_step)
    V_ind1 = -N_turns * dphi1_dt
    V_ind1_rms = float(np.sqrt(np.mean(V_ind1**2)))
    P_harvest1 = float((V_ind1_rms**2) / (4.0 * R_load_opt))

    dphi2_dt = np.gradient(phi_sonda2, dt_step)
    V_ind2 = -N_turns * dphi2_dt
    V_ind2_rms = float(np.sqrt(np.mean(V_ind2**2)))
    P_harvest2 = float((V_ind2_rms**2) / (4.0 * R_load_opt))

    area_plate = 0.005 # 50 cm^2
    der_dt = np.gradient(er_sonda3, dt_step)
    Id_s3 = EPS0 * der_dt * area_plate
    Id_rms_uA = float(np.sqrt(np.mean(Id_s3**2)) * 1e6)
    Id_peak_uA = float(np.max(np.abs(Id_s3)) * 1e6)

    # Profilo di decadimento con R in [5.5, 22] cm
    r_sweep = np.linspace(0.055, 0.22, 18)
    v_decay_series = []
    s_decay_series = []
    for r_val in r_sweep:
        b_pt = eval_B_dc_opt(np.array([[r_val, 0.0, 0.0]]))
        s_pt = eval_S(np.array([[r_val, 0.0, 0.0]]))
        v_est = N_turns * (2.0 * np.pi * 100.0) * np.linalg.norm(b_pt) * area_coil / np.sqrt(2.0)
        v_decay_series.append(float(v_est))
        s_decay_series.append(float(np.linalg.norm(s_pt)))

    results_dict["virtual_harvesting"] = {
        "probe1_radial": {
            "radius_cm": 10.0,
            "z_cm": 0.0,
            "turns": N_turns,
            "diameter_cm": 5.0,
            "V_ind_rms_V": V_ind1_rms,
            "V_ind_peak_V": float(np.max(np.abs(V_ind1))),
            "P_matched_mW": P_harvest1 * 1000.0
        },
        "probe2_axial": {
            "radius_cm": 15.0,
            "z_cm": 10.0,
            "turns": N_turns,
            "diameter_cm": 5.0,
            "V_ind_rms_V": V_ind2_rms,
            "V_ind_peak_V": float(np.max(np.abs(V_ind2))),
            "P_matched_mW": P_harvest2 * 1000.0
        },
        "probe3_capacitive": {
            "radius_cm": 12.0,
            "area_cm2": 50.0,
            "Id_rms_uA": Id_rms_uA,
            "Id_peak_uA": Id_peak_uA
        },
        "decay_curve": {
            "r_cm": (r_sweep * 100.0).tolist(),
            "V_rms_V": v_decay_series,
            "S_W_m2": s_decay_series
        }
    }
    print(f"  -> Sonda 1 (Radiale 10cm): V_rms = {V_ind1_rms:.3f} V | Potenza = {P_harvest1*1000.0:.2f} mW")
    print(f"  -> Sonda 2 (Assiale 15cm, Z=10cm): V_rms = {V_ind2_rms:.3f} V | Potenza = {P_harvest2*1000.0:.2f} mW")
    print(f"  -> Sonda 3 (Capacitiva 12cm): I_disp RMS = {Id_rms_uA:.2f} uA | Peak = {Id_peak_uA:.2f} uA")

    # Modulazione AM a 10 Hz
    print("\n[Analisi Modulazione AM (10 Hz)]...")
    am_dir = ROOT_DIR / "results_am_modulation"
    if am_dir.exists() and len(list(am_dir.glob("*.vtu"))) >= 10:
        B_am, E_am, _, _ = load_sim_series(am_dir, n_steps=10)
        t_am = np.arange(1, 11) * 0.005 * 1000.0 # ms

        r_line = np.linspace(0.045, 0.18, 35)
        bell_envelope_t = []
        for step in range(10):
            eval_b_am = build_interpolator_3d(pts, B_am[step])
            pts_line = np.column_stack([r_line, np.zeros_like(r_line), np.zeros_like(r_line)])
            b_vals = np.linalg.norm(eval_b_am(pts_line), axis=1) * 1e6
            idx_thresh = np.where(b_vals >= 20.0)[0]
            if len(idx_thresh) > 0:
                bell_envelope_t.append(float(r_line[idx_thresh[-1]] * 100.0))
            else:
                bell_envelope_t.append(float(r_line[0] * 100.0))

        results_dict["am_modulation"] = {
            "f_carrier_hz": 100.0,
            "f_mod_hz": 10.0,
            "m_depth": 0.5,
            "r_bell_max_cm": float(np.max(bell_envelope_t)),
            "r_bell_min_cm": float(np.min(bell_envelope_t)),
            "breathing_delta_cm": float(np.max(bell_envelope_t) - np.min(bell_envelope_t)),
            "envelope_timeseries_cm": bell_envelope_t
        }
        print(f"  -> Respiro Inviluppo Campana: Max = {np.max(bell_envelope_t):.2f} cm | Min = {np.min(bell_envelope_t):.2f} cm | Delta = {np.max(bell_envelope_t)-np.min(bell_envelope_t):.2f} cm")

    json_path = DATA_DIR / "sweep_sfasamento_risultati.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(results_dict, f, indent=2)
    print(f"\n[OK] Dataset JSON salvato in: {json_path}")

    # =========================================================================
    # GENERAZIONE FIGURE DIAGNOSTICHE UFFICIALI (300 DPI)
    # =========================================================================
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14
    })

    # FIGURA 04: Sweep Sfasamento Confronto
    print("\n[Generazione Grafico 04] figures/04_sweep_sfasamento_confronto.png...")
    fig = plt.figure(figsize=(14, 11), dpi=300)
    ax1 = fig.add_subplot(2, 2, 1)
    ax_pol = fig.add_subplot(2, 2, 2, projection='polar')
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    # Subplot (0,0): Forme d'onda temporali di B_radiale
    for name, fld, col in regimes:
        ax1.plot(time_arr, results_dict["regimes"][name]["Brad_t_series_uT"], label=name, color=col, lw=2.2, marker='o', markersize=4)
    ax1.set_title("Emissione Radiale Media nel Tempo <|Brad|> (R = 6 cm)")
    ax1.set_xlabel("Tempo [ms]")
    ax1.set_ylabel("<|Brad|> [uT]")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right")

    # Subplot (0,1): Profilo Polare 360° DC di Brad a R = 6 cm
    for name, fld, col in regimes:
        eval_b_reg = build_interpolator_3d(pts, np.mean(data_all_regimes[name]["B"], axis=0))
        pts_6 = np.column_stack([0.06 * np.cos(theta_arr), 0.06 * np.sin(theta_arr), np.zeros(n_theta)])
        b_eval = eval_b_reg(pts_6)
        br = np.abs(b_eval[:, 0] * np.cos(theta_arr) + b_eval[:, 1] * np.sin(theta_arr)) * 1e6
        ax_pol.plot(theta_arr, br, label=name, color=col, lw=2.0)
    ax_pol.set_title("Profilo Polare 360° |Brad| a R = 6 cm", va='bottom')
    ax_pol.legend(loc="lower right", bbox_to_anchor=(1.35, -0.1))

    # Subplot (1,0): Istogramma Ripple Index
    names_short = ["Regime A\n(Sincrono)", "Regime B\n(Co-rot. 60°)", "Regime C\n(Quadratura)", "Regime D\n(Contro-rot.)"]
    ripples = [results_dict["regimes"][r[0]]["ripple_index"] * 100 for r in regimes]
    bars = ax3.bar(names_short, ripples, color=["tab:blue", "tab:red", "tab:green", "tab:purple"], width=0.55, edgecolor="black")
    for bar in bars:
        h = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., h + 0.8, f"{h:.1f}%", ha='center', va='bottom', fontweight='bold')
    ax3.set_title("Indice di Ripple di Emissione [(Bmax - Bmin) / Bmean]")
    ax3.set_ylabel("Ripple [%]")
    ax3.set_ylim(0, max(ripples) * 1.25)
    ax3.grid(True, linestyle=":", alpha=0.6, axis='y')

    # Subplot (1,1): Potenza Joule Dissipata nel Mantello
    joules = [results_dict["regimes"][r[0]]["mean_P_joule_W"] for r in regimes]
    bars = ax4.bar(names_short, joules, color=["tab:blue", "tab:red", "tab:green", "tab:purple"], width=0.55, edgecolor="black")
    for bar in bars:
        h = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., h + 0.05, f"{h:.3f} W", ha='center', va='bottom', fontweight='bold')
    ax4.set_title("Dissipazione Termica Media nel Mantello P_J [W]")
    ax4.set_ylabel("Potenza Joule [W]")
    ax4.set_ylim(0, max(joules) * 1.25)
    ax4.grid(True, linestyle=":", alpha=0.6, axis='y')

    fig.suptitle("CARATTERIZZAZIONE DELLO SFASAMENTO SPAZIO-TEMPORALE (REGIMI A, B, C, D)", fontsize=15, fontweight='bold')
    fig.tight_layout()
    fig04_path = FIG_DIR / "04_sweep_sfasamento_confronto.png"
    fig.savefig(str(fig04_path), dpi=300)
    plt.close(fig)
    print(f"  [OK] Salvato: {fig04_path}")

    # FIGURA 05: Vettore di Poynting e Campo Elettrico
    print("\n[Generazione Grafico 05] figures/05_vettore_poynting_e_campo_elettrico.png...")
    fig, axs = plt.subplots(1, 2, figsize=(16, 7), dpi=300)

    ax = axs[0]
    grid_x = np.linspace(-0.16, 0.16, 31)
    grid_z = np.linspace(-0.14, 0.14, 27)
    GX, GZ = np.meshgrid(grid_x, grid_z)
    pts_xz = np.column_stack([GX.ravel(), np.zeros_like(GX.ravel()), GZ.ravel()])
    S_xz = eval_S(pts_xz)
    Sx = S_xz[:, 0].reshape(GX.shape)
    Sz = S_xz[:, 2].reshape(GX.shape)
    S_mag = np.sqrt(Sx**2 + Sz**2)

    im = ax.pcolormesh(GX*100, GZ*100, S_mag, shading='gouraud', cmap='inferno')
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Densità di Potenza Poynting <|S|> [W/m²]")

    Sx_norm = Sx / (S_mag + 1e-6)
    Sz_norm = Sz / (S_mag + 1e-6)
    ax.quiver(GX[::2, ::2]*100, GZ[::2, ::2]*100, Sx_norm[::2, ::2], Sz_norm[::2, ::2], color='white', alpha=0.85, scale=25)

    ax.plot([-5, -5, 5, 5, -5], [-5, 5, 5, -5, -5], 'r--', lw=2, label="Mantello Rete Stirata")
    ax.set_title("Flusso di Poynting Medio <S> (Sezione Meridiana XZ, Y=0)")
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Z [cm]")
    ax.set_xlim(-16, 16)
    ax.set_ylim(-14, 14)
    ax.legend(loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.4)

    ax = axs[1]
    grid_xy = np.linspace(-0.12, 0.12, 31)
    GX_y, GY_x = np.meshgrid(grid_xy, grid_xy)
    pts_xy = np.column_stack([GX_y.ravel(), GY_x.ravel(), np.zeros_like(GX_y.ravel())])
    E_xy = eval_E(pts_xy)
    Ex = E_xy[:, 0].reshape(GX_y.shape)
    Ey = E_xy[:, 1].reshape(GX_y.shape)
    E_mag = np.sqrt(Ex**2 + Ey**2)

    im2 = ax.pcolormesh(GX_y*100, GY_x*100, E_mag, shading='gouraud', cmap='viridis')
    cb2 = fig.colorbar(im2, ax=ax)
    cb2.set_label("Modulo Campo Elettrico Indotto <|E|> [V/m]")

    Ex_n = Ex / (E_mag + 1e-6)
    Ey_n = Ey / (E_mag + 1e-6)
    ax.quiver(GX_y[::2, ::2]*100, GY_x[::2, ::2]*100, Ex_n[::2, ::2], Ey_n[::2, ::2], color='white', alpha=0.85, scale=25)

    mant_circle = plt.Circle((0, 0), 5.0, color='r', fill=False, linestyle='--', lw=2, label="Mantello R=5cm")
    ax.add_patch(mant_circle)
    ax.set_title("Vortice di Campo Elettrico <E> (Piano Equatoriale XY, Z=0)")
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")
    ax.set_xlim(-12, 12)
    ax.set_ylim(-12, 12)
    ax.legend(loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.4)

    fig.suptitle("MAPPATURA MULTIFISICA: CAMPO ELETTRICO E TRASPORTO DI POTENZA POYNTING", fontsize=15, fontweight='bold')
    fig.tight_layout()
    fig05_path = FIG_DIR / "05_vettore_poynting_e_campo_elettrico.png"
    fig.savefig(str(fig05_path), dpi=300)
    plt.close(fig)
    print(f"  [OK] Salvato: {fig05_path}")

    # FIGURA 06: Accoppiamento a Distanza e Virtual Harvesting
    print("\n[Generazione Grafico 06] figures/06_accoppiamento_distanza_harvesting.png...")
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), dpi=300)

    ax = axs[0]
    ax.plot(time_arr, V_ind1 * 1e3, 'b-o', lw=2.2, label=f"Sonda 1 Radiale (R=10cm, Z=0)\nVrms = {V_ind1_rms*1e3:.2f} mV")
    ax.plot(time_arr, V_ind2 * 1e3, 'm-s', lw=2.2, label=f"Sonda 2 Assiale (R=15cm, Z=10cm)\nVrms = {V_ind2_rms*1e3:.2f} mV")
    ax.set_title("Tensione a Vuoto Indotta Vind(t) (N = 100 spire)")
    ax.set_xlabel("Tempo [ms]")
    ax.set_ylabel("f.e.m. Indotta Vind [mV]")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right")

    ax = axs[1]
    ax2 = ax.twinx()
    l1 = ax.plot(results_dict["virtual_harvesting"]["decay_curve"]["r_cm"],
                 results_dict["virtual_harvesting"]["decay_curve"]["V_rms_V"],
                 'tab:red', lw=2.5, marker='d', label="Tensione RMS Indotta V_rms [V]")
    l2 = ax2.plot(results_dict["virtual_harvesting"]["decay_curve"]["r_cm"],
                  results_dict["virtual_harvesting"]["decay_curve"]["S_W_m2"],
                  'tab:cyan', lw=2.0, marker='^', linestyle='--', label="Densità di Potenza Poynting <|S|> [W/m²]")
    ax.set_title("Profilo di Accoppiamento vs Distanza Radiale R")
    ax.set_xlabel("Distanza Radiale R [cm]")
    ax.set_ylabel("Tensione Indotta Vrms [V]", color='tab:red')
    ax2.set_ylabel("Poynting <|S|> [W/m²]", color='tab:cyan')
    ax.grid(True, linestyle=":", alpha=0.6)
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc="upper right")

    ax = axs[2]
    if "am_modulation" in results_dict and "envelope_timeseries_cm" in results_dict["am_modulation"]:
        t_am = np.arange(1, 11) * 5.0 # ms
        ax.plot(t_am, results_dict["am_modulation"]["envelope_timeseries_cm"], 'crimson', lw=2.5, marker='o', label="Raggio Inviluppo Campana (|B|>=20uT)")
        ax.axhline(results_dict["am_modulation"]["r_bell_max_cm"], color='navy', linestyle=':', label=f"R_max = {results_dict['am_modulation']['r_bell_max_cm']:.1f} cm")
        ax.axhline(results_dict["am_modulation"]["r_bell_min_cm"], color='green', linestyle=':', label=f"R_min = {results_dict['am_modulation']['r_bell_min_cm']:.1f} cm")
        ax.set_title("Respiro Dinamico dell'Inviluppo (Modulazione AM 10 Hz)")
        ax.set_xlabel("Tempo [ms]")
        ax.set_ylabel("Raggio Inviluppo R [cm]")
        ax.legend(loc="upper right")
    else:
        ax.plot(time_arr, Id_s3 * 1e6, 'tab:purple', lw=2.2, label=f"Corrente Id (50 cm² a 12cm)\nPeak={Id_peak_uA:.1f} uA")
        ax.set_title("Corrente di Spostamento Sonda Capacitiva")
        ax.set_xlabel("Tempo [ms]")
        ax.set_ylabel("Id [uA]")
        ax.legend(loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.6)

    fig.suptitle("VIRTUAL HARVESTING, ACCOPPIAMENTO A DISTANZA E RESPIRO DINAMICO AM", fontsize=15, fontweight='bold')
    fig.tight_layout()
    fig06_path = FIG_DIR / "06_accoppiamento_distanza_harvesting.png"
    fig.savefig(str(fig06_path), dpi=300)
    plt.close(fig)
    print(f"  [OK] Salvato: {fig06_path}")

    print("\n" + "="*80)
    print("POST-PROCESSING MULTIFISICO AVANZATO COMPLETATO CON SUCCESSO!")
    print("="*80)

if __name__ == "__main__":
    main()
