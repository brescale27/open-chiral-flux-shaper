#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-processing Scientifico: Variante Rotore Centrato a Z=0
6 Bobine a Polarità Alternate Specchiate (N-S-N-S-N-S) con Pilotaggio a Semionde Pulsate (CERN-OHL-S v2)
1. Verifica Solenoidalità del campo B su sfere di controllo (8, 12, 15 cm)
2. Calcolo integrale della forza ponderomotrice assiale netta F_z(t) su 40 timestep (2 periodi)
3. Estrazione profili radiali multi-quota B_rad(R) e mappa azimutale B_rad(360°)
4. Bilancio Vettore di Poynting S = (1/mu0) (E x B) e perdite Joule
5. Generazione Figure Ufficiali a 300 DPI in figures/:
   - fig_01_topologia_poli_specchiati_3d.png
   - fig_02_forme_onda_semionda_e_profilo_radiale.png
   - fig_03_confronto_forza_lift_Fz_impulsi.png
6. Esportazione dataset comparativo in data/confronto_semionda_specchiata.json
"""

import os
import sys
import json
import math
from pathlib import Path
import numpy as np
import meshio
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
VARIANT_DIR = ROOT_DIR / "variants" / "rotore_centrato_poli_alternati_semionda"
DATA_DIR = VARIANT_DIR / "data"
FIG_DIR = VARIANT_DIR / "figures"
RESULTS_DIR = VARIANT_DIR / "results"
BENCHMARK_JSON = ROOT_DIR / "variants" / "rotore_centrato_z0" / "data" / "confronto_variante_centrata.json"

MU0 = 4.0 * np.pi * 1e-7
EPS0 = 8.8541878128e-12

def check_solenoidality(delaunay_tri, B_last):
    print("\n--- 1. Verifica di Solenoidalita del Campo Magnetico (Gauss) ---")
    interp_B = LinearNDInterpolator(delaunay_tri, B_last, fill_value=0.0)
    
    # Campionamento sferico uniforme (spirale di Fibonacci)
    N = 2500
    indices = np.arange(0, N, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/N)
    theta = np.pi * (1 + 5**0.5) * indices
    x_u = np.sin(phi) * np.cos(theta)
    y_u = np.sin(phi) * np.sin(theta)
    z_u = np.cos(phi)
    
    radii = [0.08, 0.12, 0.15]
    solenoid_results = {}
    
    for R in radii:
        pts_sphere = np.column_stack([R * x_u, R * y_u, R * z_u])
        B_eval = interp_B(pts_sphere)
        Bn = np.sum(B_eval * np.column_stack([x_u, y_u, z_u]), axis=1)
        dS = (4.0 * np.pi * R**2) / N
        phi_net = float(np.sum(Bn) * dS)
        phi_abs = float(np.sum(np.abs(Bn)) * dS)
        res_rel = abs(phi_net) / (phi_abs + 1e-30) * 100.0
        
        solenoid_results[f"sphere_{int(R*100)}cm"] = {
            "radius_m": R,
            "area_m2": float(4.0 * np.pi * R**2),
            "phi_net_Wb": phi_net,
            "phi_abs_Wb": phi_abs,
            "residual_pct": res_rel,
            "status": "PASS" if res_rel < 3.0 else "WARNING"
        }
        print(f"  Sfera R = {R*100:4.1f} cm | Phi_net = {phi_net:+10.4e} Wb | Phi_abs = {phi_abs:10.4e} Wb | Residuo = {res_rel:6.3f}% [{solenoid_results[f'sphere_{int(R*100)}cm']['status']}]")
        
    return solenoid_results

def compute_streamlines(interp_B):
    """Calcola le linee di flusso magnetiche tridimensionali con RK45 evidenziando i dipoli alternati."""
    def b_field_norm(t, pos):
        p = np.array([pos])
        b_val = interp_B(p)[0]
        norm = np.linalg.norm(b_val)
        if norm < 1e-12 or np.isnan(norm):
            return [0.0, 0.0, 0.0]
        return (b_val / norm).tolist()
    
    streamlines = []
    # Semi di tracciamento attorno alle 6 bobine e all'equatore
    seeds = []
    
    # Semi posti sulle 6 bobine (R = 3.5 cm) a quote Z = 0, +2 cm, -2 cm
    for k in range(6):
        th = k * (2 * np.pi / 6)
        r_c = 0.035
        for z_c in [-0.02, 0.0, 0.02]:
            seeds.append([r_c * np.cos(th), r_c * np.sin(th), z_c])
            # Semi leggermente spostati verso l'esterno del mantello
            seeds.append([(r_c + 0.012) * np.cos(th), (r_c + 0.012) * np.sin(th), z_c])
            
    # Semi equatoriali nel traferro esterno
    for th in np.linspace(0, 2*np.pi, 12, endpoint=False):
        seeds.append([0.055 * np.cos(th), 0.055 * np.sin(th), 0.0])
        
    print(f"  Integrazione RK45 su {len(seeds)} linee di flusso 3D...")
    for s in seeds:
        # Integrazione in avanti
        sol_fwd = solve_ivp(b_field_norm, (0, 0.15), s, method='RK45', max_step=0.003)
        # Integrazione all'indietro
        sol_bwd = solve_ivp(lambda t, p: [-v for v in b_field_norm(t, p)], (0, 0.15), s, method='RK45', max_step=0.003)
        
        # Concatena percorso
        pts_bwd = sol_bwd.y[:, ::-1]
        pts_fwd = sol_fwd.y[:, 1:]
        pts_full = np.hstack([pts_bwd, pts_fwd])
        
        # Filtra linee valide
        if pts_full.shape[1] > 8:
            streamlines.append(pts_full)
            
    print(f"  Linee di flusso tracciate con successo: {len(streamlines)}")
    return streamlines

def plot_fig_01_topologia_3d(streamlines):
    fig = plt.figure(figsize=(11, 9), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    
    # Cilindro mantello di riferimento (R = 4.85 cm, H = 10 cm)
    z_cyl = np.linspace(-0.05, 0.05, 30)
    theta_cyl = np.linspace(0, 2*np.pi, 40)
    theta_grid, z_grid = np.meshgrid(theta_cyl, z_cyl)
    x_cyl = 0.0485 * np.cos(theta_grid)
    y_cyl = 0.0485 * np.sin(theta_grid)
    ax.plot_surface(x_cyl * 100, y_cyl * 100, z_grid * 100, color='gray', alpha=0.15, edgecolor='none')
    
    # Nucleo ferromagnetico equatoriale (R = 3.5 cm, t = 1.5 cm)
    z_core = np.linspace(-0.0075, 0.0075, 10)
    theta_core, z_core_g = np.meshgrid(theta_cyl, z_core)
    x_core = 0.035 * np.cos(theta_core)
    y_core = 0.035 * np.sin(theta_core)
    ax.plot_surface(x_core * 100, y_core * 100, z_core_g * 100, color='darkred', alpha=0.35, edgecolor='none')
    
    # Plot delle 6 posizioni delle bobine con marcatura di polarita N-S alternata
    th_coils = np.linspace(0, 2*np.pi, 6, endpoint=False)
    for k, th in enumerate(th_coils):
        x_c = 3.5 * np.cos(th)
        y_c = 3.5 * np.sin(th)
        pol_str = "N (+Z)" if k % 2 == 0 else "S (-Z)"
        col_pol = 'crimson' if k % 2 == 0 else 'royalblue'
        ax.scatter([x_c], [y_c], [0], color=col_pol, s=90, edgecolors='black', zorder=10)
        ax.text(x_c * 1.15, y_c * 1.15, 0, f"B{k+1}\n{pol_str}", color=col_pol, fontsize=8, fontweight='bold', ha='center')
        
    # Traccia le linee di flusso raggruppate per orientamento
    for line in streamlines:
        x_pts = line[0] * 100
        y_pts = line[1] * 100
        z_pts = line[2] * 100
        
        # Filtra all'interno della bounding box
        valid = (np.abs(x_pts) <= 15.0) & (np.abs(y_pts) <= 15.0) & (np.abs(z_pts) <= 12.0)
        if np.sum(valid) < 5:
            continue
        xv, yv, zv = x_pts[valid], y_pts[valid], z_pts[valid]
        
        # Colore in base alla quota media e all'avvolgimento
        mean_z = np.mean(zv)
        if mean_z > 1.5:
            color = 'navy'
            alpha = 0.65
        elif mean_z < -1.5:
            color = 'darkorange'
            alpha = 0.65
        else:
            color = 'mediumseagreen'
            alpha = 0.75
            
        ax.plot(xv, yv, zv, color=color, linewidth=1.2, alpha=alpha)
        
    ax.set_title("Open Chiral Flux Shaper - Variante Polarità Alternate N-S-N-S-N-S (Z=0)\nTopologia Elettrodinamica a Coppie Specchiate e Lobi di Espulsione Radiale", fontsize=11, fontweight='bold', pad=15)
    ax.set_xlabel("X [cm]", fontsize=9, labelpad=5)
    ax.set_ylabel("Y [cm]", fontsize=9, labelpad=5)
    ax.set_zlabel("Z [cm]", fontsize=9, labelpad=5)
    ax.set_xlim(-12, 12)
    ax.set_ylim(-12, 12)
    ax.set_zlim(-10, 10)
    ax.view_init(elev=28, azim=45)
    
    # Legenda personalizzata
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color='crimson', marker='o', linestyle='', markersize=7, label='Polo Nord (+Z, Bobine 1, 3, 5)'),
        Line2D([0], [0], color='royalblue', marker='o', linestyle='', markersize=7, label='Polo Sud (-Z, Bobine 2, 4, 6)'),
        Line2D([0], [0], color='navy', lw=2, label='Linee Flusso Lobo Superiore (+Z)'),
        Line2D([0], [0], color='darkorange', lw=2, label='Linee Flusso Lobo Inferiore (-Z)'),
        Line2D([0], [0], color='mediumseagreen', lw=2, label='Chiusura Dipolare & Espulsione Equatoriale (Z=0)'),
        Line2D([0], [0], color='gray', lw=4, alpha=0.3, label='Mantello Rete Stirata Anisotropa 30°')
    ]
    ax.legend(handles=custom_lines, loc='upper right', fontsize=8, framealpha=0.9)
    
    plt.tight_layout()
    out_path = FIG_DIR / "fig_01_topologia_poli_specchiati_3d.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Figura salvata: {out_path}")

def plot_fig_02_forme_onda_e_profilo_radiale(t_arr, j_channels, th_deg, brad_circ_cycle1, brad_circ_cycle2, r_cm, brad_radial):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=300)
    
    # 1. Forme d'onda delle correnti nelle 6 bobine (2 periodi: 0-20 ms)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    for k in range(6):
        pol_tag = "(+Z, N)" if k % 2 == 0 else "(-Z, S)"
        ax1.plot(t_arr, j_channels[k] / 1e5, label=f"Bobina {k+1} {pol_tag}", color=colors[k], lw=1.8)
        
    ax1.axvline(10.0, color='black', linestyle='--', alpha=0.6, label='Inizio 2° Periodo (T=10ms)')
    ax1.set_title("Forme d'Onda delle Correnti Rotoriche: Treno a Semionde Pulsate Sfasate 60°\nPolarità Geometriche Alternate N-S-N-S-N-S", fontsize=10, fontweight='bold')
    ax1.set_xlabel("Tempo [ms]", fontsize=9)
    ax1.set_ylabel(r"Densità di Corrente $J_z$ [$10^5$ A/m$^2$]", fontsize=9)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='lower left', fontsize=7.5, ncol=2)
    ax1.set_xlim(0, 20.0)
    
    # 2. Profilo circonferenziale di induzione radiale B_rad(theta) a R=6cm, Z=0
    ax2.plot(th_deg, brad_circ_cycle1, label=r"Ciclo 1 (t = 5.0 ms)", color='darkorange', lw=2.0, linestyle='--')
    ax2.plot(th_deg, brad_circ_cycle2, label=r"Ciclo 2 a Regime (t = 15.0 ms)", color='navy', lw=2.2)
    
    # Marcatori poli bobine a 60°
    for k in range(6):
        ang = k * 60
        pol_tag = "N" if k % 2 == 0 else "S"
        col_pol = 'crimson' if k % 2 == 0 else 'royalblue'
        ax2.axvline(ang, color=col_pol, linestyle=':', alpha=0.5)
        ax2.text(ang, ax2.get_ylim()[1]*0.85 if hasattr(ax2, 'get_ylim') else 50, f"B{k+1}\n({pol_tag})", 
                 color=col_pol, fontsize=7.5, fontweight='bold', ha='center')
        
    ax2.set_title(r"Profilo Azimutale di Induzione Radiale $|B_{\mathrm{rad}}(\theta)|$ all'Equatore ($R=6$ cm, $Z=0$)" + "\nEvidenza dei 6 Lobi Dipolari Formati dalle Coppie Specchiate", fontsize=10, fontweight='bold')
    ax2.set_xlabel("Angolo Azimutale [gradi]", fontsize=9)
    ax2.set_ylabel(r"Induzione Radiale $|B_{\mathrm{rad}}|$ [$\mu$T]", fontsize=9)
    ax2.set_xlim(0, 360)
    ax2.set_xticks(np.arange(0, 361, 60))
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='upper right', fontsize=8.5)
    
    plt.tight_layout()
    out_path = FIG_DIR / "fig_02_forme_onda_semionda_e_profilo_radiale.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Figura salvata: {out_path}")

def plot_fig_03_confronto_lift(t_ms, fz_tot, fz_mantle, fz_rotor, fz_bench_mean, fz_bench_series):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), dpi=300, gridspec_kw={'height_ratios': [2, 1]})
    
    # 1. Confronto dinamico della forza di Lorentz assiale F_z(t)
    ax1.plot(t_ms, fz_tot * 1e3, label=r"$F_z(t)$ Totale (Poli Alternati Semionde)", color='navy', lw=2.2)
    ax1.plot(t_ms, fz_mantle * 1e3, label=r"$F_z(t)$ Mantello Rete Stirata (60%)", color='mediumseagreen', lw=1.6, linestyle='--')
    ax1.plot(t_ms, fz_rotor * 1e3, label=r"$F_z(t)$ Rotore e Nucleo", color='darkorange', lw=1.4, linestyle=':')
    
    # Medie
    mean_all = float(np.mean(fz_tot))
    mean_cycle2 = float(np.mean(fz_tot[20:]))  # Secondo periodo (regime)
    
    ax1.axhline(mean_all * 1e3, color='blue', linestyle='-.', lw=1.5, 
                label=rf"Media Globale: $\langle F_z \rangle = {mean_all*1e3:+.4f}$ mN ({mean_all*1e6:+.2f} $\mu$N)")
    ax1.axhline(mean_cycle2 * 1e3, color='crimson', linestyle='-', lw=1.8, 
                label=rf"Media Regime Ciclo 2: $\langle F_z \rangle = {mean_cycle2*1e3:+.4f}$ mN ({mean_cycle2*1e6:+.2f} $\mu$N)")
    ax1.axhline(fz_bench_mean, color='gray', linestyle='--', lw=1.8, 
                label=rf"Benchmark Continuo (Regime B): $\langle F_z \rangle = {fz_bench_mean:+.4f}$ mN ({fz_bench_mean*1e3:+.2f} $\mu$N)")
    
    ax1.axvline(10.0, color='black', linestyle=':', alpha=0.7, label='Separazione Ciclo 1 / Ciclo 2')
    
    ax1.set_title("Evoluzione Temporale della Forza Assiale di Lorentz $F_z(t)$ (Lift Ponderomotore)\nConfronto: Pilotaggio Pulsato a Semionde vs Regime Continuo a Sinusoide Piena", fontsize=11, fontweight='bold')
    ax1.set_ylabel("Forza Assiale $F_z$ [mN]", fontsize=9)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=8, ncol=2)
    ax1.set_xlim(0, 20.0)
    
    # 2. Zoom comparativo sul secondo periodo (10-20 ms) con confronto diretto con Benchmark
    t_c2 = t_ms[20:] - 10.0
    ax2.plot(t_c2, fz_tot[20:] * 1e3, label=r"Poli Alternati Semionde (Ciclo 2 a Regime)", color='navy', lw=2.0)
    if len(fz_bench_series) == 20:
        t_bench = np.linspace(0.5, 10.0, 20)
        ax2.plot(t_bench, np.array(fz_bench_series), label=r"Regime B Continuo (1 Ciclo)", color='gray', lw=1.8, linestyle='--')
        
    ax2.axhline(mean_cycle2 * 1e3, color='crimson', linestyle='-', lw=1.5)
    ax2.axhline(fz_bench_mean, color='gray', linestyle='--', lw=1.5)
    
    ax2.set_title("Confronto Dettagliato a Regime di 1 Ciclo (10 ms)", fontsize=10, fontweight='bold')
    ax2.set_xlabel("Tempo Ciclo a Regime [ms]", fontsize=9)
    ax2.set_ylabel("Forza Assiale $F_z$ [mN]", fontsize=9)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='upper right', fontsize=8)
    ax2.set_xlim(0, 10.0)
    
    plt.tight_layout()
    out_path = FIG_DIR / "fig_03_confronto_forza_lift_Fz_impulsi.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Figura salvata: {out_path}")

def main():
    print("=" * 80)
    print("Post-processing Scientifico: Variante Rotore Centrato con Poli Alternate Semionde")
    print("=" * 80)
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    
    # Verifica presenza dei file VTU
    vtu_files = sorted(list(RESULTS_DIR.glob("macchina_out_t*.vtu")))
    if len(vtu_files) < 40:
        print(f"[ERRORE] Risultati incompleti in {RESULTS_DIR}: trovati {len(vtu_files)}/40 VTU.")
        return 1
        
    print(f"Trovati tutti i {len(vtu_files)} file VTU in {RESULTS_DIR}.")
    
    # Carica la mesh di base
    m0 = meshio.read(str(vtu_files[0]))
    pts = m0.points
    cells = m0.cells_dict['tetra']
    n_elems = len(cells)
    
    v0 = pts[cells[:, 0]]
    v1 = pts[cells[:, 1]]
    v2 = pts[cells[:, 2]]
    v3 = pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0
    elem_com = (v0 + v1 + v2 + v3) / 4.0
    elem_r = np.hypot(elem_com[:, 0], elem_com[:, 1])
    elem_z = np.abs(elem_com[:, 2])
    
    active_mask = (elem_r <= 0.051) & (elem_z <= 0.051)
    mantle_mask = (elem_r >= 0.046) & (elem_r <= 0.051) & (elem_z <= 0.051)
    rotor_mask = (elem_r < 0.046) & (elem_z <= 0.051)
    
    delaunay_tri = Delaunay(pts)
    
    # Vettori per serie temporale
    time_ms = np.linspace(0.5, 20.0, 40)
    fz_total = []
    fz_mantle = []
    fz_rotor = []
    pj_total = []
    
    # Cilindro per flusso di Poynting (R=12cm, H=20cm)
    R_cyl = 0.12
    H_cyl = 0.20
    n_th = 36
    n_z = 20
    th_grid = np.linspace(0, 2*np.pi, n_th, endpoint=False)
    z_grid = np.linspace(-H_cyl/2, H_cyl/2, n_z)
    TH, ZG = np.meshgrid(th_grid, z_grid)
    pts_cyl = np.column_stack([R_cyl * np.cos(TH).ravel(), R_cyl * np.sin(TH).ravel(), ZG.ravel()])
    nx_cyl = np.cos(TH).ravel()
    ny_cyl = np.sin(TH).ravel()
    dA_cyl = (2*np.pi*R_cyl / n_th) * (H_cyl / n_z)
    
    poynting_series = []
    s_mag_series = []
    
    # Campionatori specifici
    interp_B_cycle1 = None
    interp_B_cycle2 = None
    interp_B_last = None
    
    print("\n--- 2. Elaborazione dei 40 Timestep: Forza di Lorentz, Poynting e Campi ---")
    for step_idx, f in enumerate(vtu_files, 1):
        mt = meshio.read(str(f))
        jxb = mt.point_data['jxb']
        jxb_z = jxb[:, 2]
        pj = mt.point_data['joule heating'].ravel()
        B = mt.point_data['magnetic flux density']
        E = mt.point_data['electric field']
        
        # Integrazione ponderomotrice
        jxb_z_elem = np.mean(jxb_z[cells], axis=1)
        pj_elem = np.mean(pj[cells], axis=1)
        
        fz_tot = float(np.sum(elem_vols[active_mask] * jxb_z_elem[active_mask]))
        fz_m = float(np.sum(elem_vols[mantle_mask] * jxb_z_elem[mantle_mask]))
        fz_r = float(np.sum(elem_vols[rotor_mask] * jxb_z_elem[rotor_mask]))
        pj_tot = float(np.sum(elem_vols[active_mask] * pj_elem[active_mask]))
        
        fz_total.append(fz_tot)
        fz_mantle.append(fz_m)
        fz_rotor.append(fz_r)
        pj_total.append(pj_tot)
        
        # Flusso di Poynting sul cilindro
        interp_B_t = LinearNDInterpolator(delaunay_tri, B, fill_value=0.0)
        interp_E_t = LinearNDInterpolator(delaunay_tri, E, fill_value=0.0)
        
        if step_idx == 10:
            interp_B_cycle1 = interp_B_t
        if step_idx == 30:
            interp_B_cycle2 = interp_B_t
        if step_idx == 40:
            interp_B_last = interp_B_t
            
        B_c = interp_B_t(pts_cyl)
        E_c = interp_E_t(pts_cyl)
        S_c = np.cross(E_c, B_c) / MU0
        Sn_c = S_c[:, 0] * nx_cyl + S_c[:, 1] * ny_cyl
        poynting_series.append(float(np.sum(Sn_c) * dA_cyl))
        s_mag_series.append(float(np.mean(np.linalg.norm(S_c, axis=1))))
        
        if step_idx % 10 == 0:
            print(f"  Step {step_idx}/40 completato | Fz = {fz_tot*1e3:+.4f} mN | P_out = {poynting_series[-1]*1e3:+.3f} mW")
            
    fz_total = np.array(fz_total)
    fz_mantle = np.array(fz_mantle)
    fz_rotor = np.array(fz_rotor)
    pj_total = np.array(pj_total)
    poynting_series = np.array(poynting_series)
    
    # 3. Solenoidalità del Flusso
    solenoid_data = check_solenoidality(delaunay_tri, mt.point_data['magnetic flux density'])
    
    # 4. Profili di Induzione Radiale
    print("\n--- 3. Estrazione Profili Radiali Multi-Quota e Profilo Azimutale ---")
    r_sweep = np.linspace(0.048, 0.25, 60)
    theta_sweep = np.linspace(0, 2*np.pi, 72, endpoint=False)
    
    # Profilo circonferenziale B_rad(theta) a R=6cm, Z=0 per ciclo 1 (t=5ms) e ciclo 2 (t=15ms)
    pts_circ = np.column_stack([0.06 * np.cos(theta_sweep), 0.06 * np.sin(theta_sweep), np.zeros(72)])
    B_circ_c1 = interp_B_cycle1(pts_circ)
    Brad_circ_c1 = np.abs(B_circ_c1[:, 0] * np.cos(theta_sweep) + B_circ_c1[:, 1] * np.sin(theta_sweep)) * 1e6
    B_circ_c2 = interp_B_cycle2(pts_circ)
    Brad_circ_c2 = np.abs(B_circ_c2[:, 0] * np.cos(theta_sweep) + B_circ_c2[:, 1] * np.sin(theta_sweep)) * 1e6
    
    # Profili radiali alle 3 quote Z = 0, +5 cm, -5 cm per il ciclo 2 a regime
    profiles = {}
    for z_target, label in [(0.0, "Z_0cm"), (0.05, "Z_plus5cm"), (-0.05, "Z_minus5cm")]:
        b_rad_curve = []
        for r_val in r_sweep:
            x_sw = r_val * np.cos(theta_sweep)
            y_sw = r_val * np.sin(theta_sweep)
            z_sw = np.full_like(x_sw, z_target)
            pts_sw = np.column_stack([x_sw, y_sw, z_sw])
            B_sw = interp_B_cycle2(pts_sw)
            Brad_sw = B_sw[:, 0] * np.cos(theta_sweep) + B_sw[:, 1] * np.sin(theta_sweep)
            b_rad_curve.append(float(np.mean(np.abs(Brad_sw))))
        profiles[label] = {
            "r_cm": (r_sweep * 100.0).tolist(),
            "Brad_uT": (np.array(b_rad_curve) * 1e6).tolist(),
            "Brad_at_6cm_uT": float(b_rad_curve[4] * 1e6),
            "Brad_at_10cm_uT": float(b_rad_curve[18] * 1e6),
            "max_Brad_uT": float(np.max(b_rad_curve) * 1e6)
        }
        print(f"  {label:12s} | B_rad(6cm) = {profiles[label]['Brad_at_6cm_uT']:8.2f} uT | B_rad(10cm) = {profiles[label]['Brad_at_10cm_uT']:6.2f} uT | Max = {profiles[label]['max_Brad_uT']:8.2f} uT")
        
    # Forme d'onda correnti nei 6 canali
    w_e = 2 * np.pi * 100.0
    t_dense = np.linspace(0, 0.02, 400)
    j_channels = []
    for k in range(6):
        phi_k = k * (np.pi / 3.0)
        s_k = 1.0 if k % 2 == 0 else -1.0
        v_k = np.sin(w_e * t_dense - phi_k)
        j_k = s_k * 1.0e5 * 0.5 * (v_k + np.abs(v_k))
        j_channels.append(j_k)
        
    # Carica dati di benchmark continuo se presenti
    fz_bench_mean = 0.004669
    fz_bench_series = []
    if BENCHMARK_JSON.is_file():
        with open(BENCHMARK_JSON, "r", encoding="utf-8") as f:
            bdata = json.load(f)
            fz_bench_mean = bdata["electrodynamic_lift_fz"].get("mean_Fz_total_mN", 0.004669)
            fz_bench_series = bdata["electrodynamic_lift_fz"].get("Fz_time_series_mN", [])
            
    # 5. Rendering Figure Diagnostiche a 300 DPI
    print("\n--- 4. Rendering Figure Diagnostiche a 300 DPI ---")
    streamlines = compute_streamlines(interp_B_cycle2)
    plot_fig_01_topologia_3d(streamlines)
    plot_fig_02_forme_onda_e_profilo_radiale(t_dense * 1e3, j_channels, np.degrees(theta_sweep), Brad_circ_c1, Brad_circ_c2, r_sweep * 100, profiles["Z_0cm"]["Brad_uT"])
    plot_fig_03_confronto_lift(time_ms, fz_total, fz_mantle, fz_rotor, fz_bench_mean, fz_bench_series)
    
    # 6. Statistiche di Forza e Confronto
    mean_fz_tot = float(np.mean(fz_total))
    mean_fz_cycle2 = float(np.mean(fz_total[20:]))
    max_fz = float(np.max(fz_total))
    min_fz = float(np.min(fz_total))
    ripple_fz = float((max_fz - min_fz) / (abs(mean_fz_cycle2) + 1e-12))
    
    mean_pj_tot = float(np.mean(pj_total))
    mean_pj_cycle2 = float(np.mean(pj_total[20:]))
    mean_poynting_out = float(np.mean(poynting_series))
    mean_poynting_c2 = float(np.mean(poynting_series[20:]))
    
    print("\n" + "=" * 80)
    print("SINTESI RISULTATI VARIANTE POLI ALTERNATI SEMIONDE (N-S-N-S-N-S):")
    print(f"  Lift Fz Medio Globale:        {mean_fz_tot*1e3:+.4f} mN ({mean_fz_tot*1e6:+.2f} uN)")
    print(f"  Lift Fz Medio Ciclo 2 Regime: {mean_fz_cycle2*1e3:+.4f} mN ({mean_fz_cycle2*1e6:+.2f} uN)")
    print(f"  Benchmark Regime B Continuo:  {fz_bench_mean:+.4f} mN ({fz_bench_mean*1e3:+.2f} uN)")
    print(f"  Picco Fz:                     {max_fz*1e3:+.4f} mN | Min: {min_fz*1e3:+.4f} mN")
    print(f"  Brad all'Equatore (6 cm):     {profiles['Z_0cm']['Brad_at_6cm_uT']:.2f} uT")
    print(f"  Dissipazione Joule Media:     {mean_pj_cycle2:.4f} W")
    print(f"  Potenza Poynting Uscente:     {mean_poynting_c2*1e3:+.4f} mW")
    print("=" * 80)
    
    # Esportazione JSON
    out_json = {
        "metadata": {
            "variant": "rotore_centrato_poli_alternati_semionda",
            "description": "Variante rotore centrato Z=0 con 6 bobine a polarità alternate specchiate N-S-N-S-N-S e pilotaggio a semionde pulsate sfasate di 60 gradi",
            "license": "CERN-OHL-S-2.0",
            "author": "Alessandro Brescacin",
            "date": "2026-09-22",
            "timesteps": 40,
            "dt_s": 0.0005,
            "total_time_ms": 20.0
        },
        "solenoidality_gauss": solenoid_data,
        "electrodynamic_lift_fz": {
            "mean_Fz_global_mN": mean_fz_tot * 1e3,
            "mean_Fz_global_uN": mean_fz_tot * 1e6,
            "mean_Fz_regime_cycle2_mN": mean_fz_cycle2 * 1e3,
            "mean_Fz_regime_cycle2_uN": mean_fz_cycle2 * 1e6,
            "peak_Fz_mN": max_fz * 1e3,
            "min_Fz_mN": min_fz * 1e3,
            "ripple_ratio": ripple_fz,
            "mean_Fz_mantle_regime_mN": float(np.mean(fz_mantle[20:])) * 1e3,
            "mean_Fz_rotor_regime_mN": float(np.mean(fz_rotor[20:])) * 1e3,
            "benchmark_continuous_regime_B_mN": fz_bench_mean,
            "benchmark_continuous_regime_B_uN": fz_bench_mean * 1e3,
            "time_series_ms": time_ms.tolist(),
            "fz_total_series_mN": (fz_total * 1e3).tolist(),
            "fz_mantle_series_mN": (fz_mantle * 1e3).tolist(),
            "fz_rotor_series_mN": (fz_rotor * 1e3).tolist()
        },
        "radial_induction_profiles": profiles,
        "circumferential_profile_R6cm": {
            "theta_deg": np.degrees(theta_sweep).tolist(),
            "Brad_cycle1_t5ms_uT": Brad_circ_c1.tolist(),
            "Brad_cycle2_t15ms_uT": Brad_circ_c2.tolist()
        },
        "energy_and_losses": {
            "mean_joule_heating_regime_W": mean_pj_cycle2,
            "mean_joule_heating_global_W": mean_pj_tot,
            "net_poynting_power_out_regime_mW": mean_poynting_c2 * 1e3,
            "net_poynting_power_out_global_mW": mean_poynting_out * 1e3,
            "mean_poynting_magnitude_W_m2": float(np.mean(s_mag_series[20:])),
            "control_cylinder_R_cm": 12.0,
            "control_cylinder_H_cm": 20.0
        },
        "full_scale_1e7_projection": {
            "scale_factor_current": 100.0,
            "scale_factor_B": 100.0,
            "scale_factor_Force": 10000.0,
            "projected_Fz_regime_N": mean_fz_cycle2 * 10000.0,
            "projected_Brad_at_6cm_mT": profiles["Z_0cm"]["Brad_at_6cm_uT"] * 1e-6 * 100.0 * 1e3,
            "projected_Poynting_W": mean_poynting_c2 * 1e-3 * 10000.0
        }
    }
    
    json_path = DATA_DIR / "confronto_semionda_specchiata.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out_json, f, indent=2)
    print(f"  Dataset salvato: {json_path}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
