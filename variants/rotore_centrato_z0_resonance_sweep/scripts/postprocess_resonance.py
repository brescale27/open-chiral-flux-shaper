#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-processing Scientifico: Sweep Parametrico 2D (Frequenza vs RPM)
Mappatura Superficie di Risonanza e Massimizzazione del Lift di Lorentz (CERN-OHL-S v2)

1. Analisi matrice 5x5: identificazione punto di risonanza ottimale (f_opt, RPM_opt)
2. Calcolo Figure of Merit / Efficienza elettrodinamica eta_F = <Fz> / P_J [uN / W]
3. Fit della curva di dispersione di Kloss in funzione dello scorrimento f_slip
4. Riscalamento quadratico a piena scala di rame J_0 = 1.0e7 A/m^2 (x 10^4)
5. Generazione Figure Ufficiali a 300 DPI in figures/:
   - fig_01_superficie_risonanza_lift_2d.png
   - fig_02_curva_dispersione_vs_slip.png
   - fig_03_confronto_forma_onda_ottimo_vs_nominale.png
6. Esportazione dataset strutturato in data/sweep_risonanza_matrice.json
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from scipy.optimize import curve_fit

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
SWEEP_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0_resonance_sweep"
DATA_DIR = SWEEP_DIR / "data"
FIG_DIR = SWEEP_DIR / "figures"
RAW_JSON = DATA_DIR / "sweep_risonanza_raw.json"

FREQUENCIES = [50, 100, 200, 400, 800]
RPMS = [0, 600, 1200, 2400, 4800]

def kloss_func(f_s, F_max, f_crit):
    """Modello dispersivo di induzione a campana di Kloss."""
    ratio = (f_s + 1e-6) / (f_crit + 1e-6)
    return F_max * (2.0 / (ratio + 1.0 / ratio))

def plot_fig_01_superficie_2d(freqs, rpms, fz_grid, opt_point, nom_point):
    fig = plt.figure(figsize=(15, 6), dpi=300)
    
    # Subplot 1: Mappa di livello 2D (Contour Plot)
    ax1 = fig.add_subplot(1, 2, 1)
    
    # Griglia interpolata densa
    f_log = np.log10(freqs)
    f_dense_log = np.linspace(f_log[0], f_log[-1], 200)
    rpm_dense = np.linspace(rpms[0], rpms[-1], 200)
    FLOG, RPMG = np.meshgrid(f_dense_log, rpm_dense)
    
    # Punti griglia originale
    pts_orig = []
    vals_orig = []
    for i, f in enumerate(freqs):
        for j, n in enumerate(rpms):
            pts_orig.append([np.log10(f), n])
            vals_orig.append(fz_grid[i, j])
            
    FZ_dense = griddata(pts_orig, vals_orig, (FLOG, RPMG), method='cubic')
    
    cf = ax1.contourf(10**FLOG, RPMG, FZ_dense, levels=30, cmap='viridis')
    cbar = plt.colorbar(cf, ax=ax1)
    cbar.set_label(r"Spinta Media $\langle F_z \rangle$ [$\mu$N]", fontsize=9)
    
    cs = ax1.contour(10**FLOG, RPMG, FZ_dense, levels=12, colors='white', linewidths=0.7, alpha=0.6)
    ax1.clabel(cs, inline=True, fontsize=7.5, fmt='%.1f')
    
    # Evidenzia punto nominale e punto ottimale
    f_nom, n_nom, fz_nom = nom_point
    f_opt, n_opt, fz_opt = opt_point
    
    ax1.scatter([f_nom], [n_nom], color='cyan', edgecolors='black', s=120, zorder=10, 
                label=rf"Nominale (100 Hz, 1200 RPM): ${fz_nom:+.2f}\,\mu$N")
    ax1.scatter([f_opt], [n_opt], color='red', marker='*', edgecolors='black', s=250, zorder=10, 
                label=rf"Risonanza Ottima ({f_opt:.0f} Hz, {n_opt:.0f} RPM): ${fz_opt:+.2f}\,\mu$N (+{((fz_opt/fz_nom)-1)*100:.1f}%)")
    
    ax1.set_xscale('log')
    ax1.set_xlabel("Frequenza Elettrica $f$ [Hz]", fontsize=9)
    ax1.set_ylabel("Velocità di Rotazione Meccanica $n$ [RPM]", fontsize=9)
    ax1.set_title("Mappa 2D di Risonanza del Lift di Lorentz $\langle F_z \rangle (f, n)$\n(Simmetria Biconica Equatoriale $Z=0$, Regime B Co-rotante)", fontsize=10, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.5, which='both')
    ax1.legend(loc='upper left', fontsize=8, framealpha=0.9)
    
    # Subplot 2: Superficie 3D
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax2.plot_surface(FLOG, RPMG, FZ_dense, cmap='viridis', edgecolor='none', alpha=0.85)
    
    ax2.scatter([np.log10(f_nom)], [n_nom], [fz_nom], color='cyan', edgecolors='black', s=80, zorder=10)
    ax2.scatter([np.log10(f_opt)], [n_opt], [fz_opt], color='red', marker='*', edgecolors='black', s=180, zorder=10)
    
    ax2.set_xlabel(r"$\log_{10}(f)$ [Hz]", fontsize=8, labelpad=5)
    ax2.set_ylabel("RPM", fontsize=8, labelpad=5)
    ax2.set_zlabel(r"$\langle F_z \rangle$ [$\mu$N]", fontsize=8, labelpad=5)
    ax2.set_title("Superficie Tridimensionale di Spinta Assiale", fontsize=10, fontweight='bold', pad=10)
    ax2.view_init(elev=28, azim=-125)
    
    plt.tight_layout()
    out_path = FIG_DIR / "fig_01_superficie_risonanza_lift_2d.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Figura salvata: {out_path}")

def plot_fig_02_dispersione_slip(slip_data, freqs):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=300)
    
    # Subplot 1: Fz vs f_slip
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    f_slip_all = []
    fz_all = []
    
    for idx, f in enumerate(freqs):
        f_pts = [p for p in slip_data if p['frequency_Hz'] == f]
        f_pts.sort(key=lambda x: x['f_slip_Hz'])
        s_vals = [p['f_slip_Hz'] for p in f_pts]
        fz_vals = [p['mean_Fz_uN'] for p in f_pts]
        
        f_slip_all.extend(s_vals)
        fz_all.extend(fz_vals)
        
        ax1.plot(s_vals, fz_vals, marker='o', lw=1.6, color=colors[idx], label=f"f = {f} Hz")
        
    # Fit di Kloss su tutti i punti o punti dominanti
    f_slip_arr = np.array(f_slip_all)
    fz_arr = np.array(fz_all)
    
    try:
        popt, _ = curve_fit(kloss_func, f_slip_arr, fz_arr, p0=[np.max(fz_arr), 100.0], bounds=([0, 10], [np.max(fz_arr)*3, 500]))
        s_fit = np.linspace(10, 800, 300)
        fit_curve = kloss_func(s_fit, *popt)
        ax1.plot(s_fit, fit_curve, color='black', linestyle='--', lw=2.0, 
                 label=rf"Fit Dispersivo Kloss ($f_{{\mathrm{{crit}}}} = {popt[1]:.1f}$ Hz, $F_{{\mathrm{{max}}}} = {popt[0]:.2f}\,\mu$N)")
    except Exception as e:
        print(f"  [AVVISO] Fit di Kloss non convergente: {e}")
        
    ax1.set_title("Curva di Dispersione Elettromeccanica: Spinta di Lorentz $\langle F_z \rangle$ vs Scorrimento\n$f_{\mathrm{slip}} = |f - (p \cdot n / 60)|$ con $p = 3$", fontsize=10, fontweight='bold')
    ax1.set_xlabel("Frequenza di Scorrimento $f_{\mathrm{slip}}$ [Hz]", fontsize=9)
    ax1.set_ylabel(r"Spinta Assiale Media $\langle F_z \rangle$ [$\mu$N]", fontsize=9)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=8)
    
    # Subplot 2: Efficienza Elettrodinamica eta_F = <Fz> / P_J [uN / W] vs f_slip
    for idx, f in enumerate(freqs):
        f_pts = [p for p in slip_data if p['frequency_Hz'] == f]
        f_pts.sort(key=lambda x: x['f_slip_Hz'])
        s_vals = [p['f_slip_Hz'] for p in f_pts]
        eta_vals = [p['efficiency_uN_per_W'] for p in f_pts]
        ax2.plot(s_vals, eta_vals, marker='s', lw=1.6, color=colors[idx], label=f"f = {f} Hz")
        
    ax2.set_title(r"Figure of Merit Elettrodinamica: Efficienza di Spinta $\eta_F = \frac{\langle F_z \rangle}{P_J}$" + "\n(Spinta generata per unità di potenza Joule dissipata nel mantello)", fontsize=10, fontweight='bold')
    ax2.set_xlabel("Frequenza di Scorrimento $f_{\mathrm{slip}}$ [Hz]", fontsize=9)
    ax2.set_ylabel(r"Efficienza $\eta_F$ [$\mu$N / W]", fontsize=9)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='upper right', fontsize=8)
    
    plt.tight_layout()
    out_path = FIG_DIR / "fig_02_curva_dispersione_vs_slip.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Figura salvata: {out_path}")

def plot_fig_03_confronto_onde(nom_res, opt_res):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), dpi=300, gridspec_kw={'height_ratios': [2, 1]})
    
    # 1. Forma d'onda nel ciclo normalizzato t / T in [0, 1]
    t_norm = np.linspace(0, 1, len(nom_res["Fz_series_mN"]))
    fz_nom = np.array(nom_res["Fz_series_mN"]) * 1e3  # in uN
    fz_opt = np.array(opt_res["Fz_series_mN"]) * 1e3  # in uN
    
    ax1.plot(t_norm, fz_opt, color='red', lw=2.2, 
             label=rf"Punto Ottimale di Risonanza ({opt_res['frequency_Hz']} Hz, {opt_res['rpm']} RPM): $\langle F_z \rangle = {opt_res['mean_Fz_uN']:+.2f}\,\mu$N")
    ax1.plot(t_norm, fz_nom, color='navy', lw=1.8, linestyle='--', 
             label=rf"Punto Nominale Baseline (100 Hz, 1200 RPM): $\langle F_z \rangle = {nom_res['mean_Fz_uN']:+.2f}\,\mu$N")
    
    ax1.axhline(opt_res['mean_Fz_uN'], color='red', linestyle=':', alpha=0.8, lw=1.5)
    ax1.axhline(nom_res['mean_Fz_uN'], color='navy', linestyle=':', alpha=0.8, lw=1.5)
    
    ax1.set_title("Dinamica Temporale Comparativa della Spinta di Lorentz $F_z(t)$\nNormalizzata sul Periodo Elettrico $T = 1/f$", fontsize=11, fontweight='bold')
    ax1.set_ylabel(r"Forza Assiale Istantanea $F_z$ [$\mu$N]", fontsize=9)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=8.5)
    ax1.set_xlim(0, 1.0)
    
    # 2. Istogramma comparativo delle metriche chiave
    metrics = ["Spinta Media", "Picco Massimo", "Perdite Joule", "Efficienza"]
    v_nom = [nom_res['mean_Fz_uN'], nom_res['peak_Fz_uN'], nom_res['mean_Poule_W'], nom_res['efficiency_uN_per_W']]
    v_opt = [opt_res['mean_Fz_uN'], opt_res['peak_Fz_uN'], opt_res['mean_Poule_W'], opt_res['efficiency_uN_per_W']]
    
    # Normalizzati a nominale = 100%
    ratios_pct = [(v_opt[i] / (v_nom[i] + 1e-12)) * 100.0 for i in range(4)]
    
    bars = ax2.bar(metrics, ratios_pct, color=['crimson', 'orange', 'mediumseagreen', 'royalblue'], width=0.45, edgecolor='black')
    ax2.axhline(100.0, color='black', linestyle='--', alpha=0.7, label='Baseline Nominale (100%)')
    
    for bar in bars:
        h = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., h + 3.0, f"{h:.1f}%", ha='center', va='bottom', fontsize=8.5, fontweight='bold')
        
    ax2.set_title("Fattore di Incremento Percentuale Risonanza vs Baseline (%)", fontsize=10, fontweight='bold')
    ax2.set_ylabel("Valore Relativo [%]", fontsize=9)
    ax2.set_ylim(0, max(ratios_pct) * 1.25)
    ax2.grid(True, linestyle=':', alpha=0.6, axis='y')
    ax2.legend(loc='upper right', fontsize=8)
    
    plt.tight_layout()
    out_path = FIG_DIR / "fig_03_confronto_forma_onda_ottimo_vs_nominale.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Figura salvata: {out_path}")

def main():
    print("=" * 80)
    print("Post-processing Scientifico: Analisi di Risonanza e Curve di Dispersione 2D")
    print("=" * 80)
    
    if not RAW_JSON.is_file():
        print(f"[ERRORE] File raw non trovato: {RAW_JSON}")
        return 1
        
    with open(RAW_JSON, "r", encoding="utf-8") as f:
        raw_data = json.load(f)
        
    print(f"Caricati {len(raw_data)} punti operativi da {RAW_JSON}.")
    
    # Costruzione matrice 5x5
    fz_matrix_uN = np.zeros((len(FREQUENCIES), len(RPMS)))
    pj_matrix_W = np.zeros((len(FREQUENCIES), len(RPMS)))
    eta_matrix = np.zeros((len(FREQUENCIES), len(RPMS)))
    fslip_matrix = np.zeros((len(FREQUENCIES), len(RPMS)))
    
    dict_lookup = {}
    for pt in raw_data:
        f = pt['frequency_Hz']
        n = pt['rpm']
        dict_lookup[(f, n)] = pt
        if f in FREQUENCIES and n in RPMS:
            i = FREQUENCIES.index(f)
            j = RPMS.index(n)
            fz_matrix_uN[i, j] = pt['mean_Fz_uN']
            pj_matrix_W[i, j] = pt['mean_Poule_W']
            eta_matrix[i, j] = pt['efficiency_uN_per_W']
            fslip_matrix[i, j] = pt['f_slip_Hz']
            
    # Identificazione punto nominale
    nom_res = dict_lookup.get((100, 1200))
    if not nom_res:
        # Fallback al più vicino
        nom_res = raw_data[0]
        
    nom_point = (100, 1200, nom_res['mean_Fz_uN'])
    
    # Identificazione massimo assoluto di spinta assiale
    max_idx = np.unravel_index(np.argmax(fz_matrix_uN), fz_matrix_uN.shape)
    opt_f = FREQUENCIES[max_idx[0]]
    opt_rpm = RPMS[max_idx[1]]
    opt_res = dict_lookup[(opt_f, opt_rpm)]
    opt_point = (opt_f, opt_rpm, opt_res['mean_Fz_uN'])
    
    boost_pct = ((opt_res['mean_Fz_uN'] - nom_res['mean_Fz_uN']) / abs(nom_res['mean_Fz_uN'])) * 100.0
    
    print("\n" + "=" * 80)
    print("ANALISI DI RISONANZA ELETTROMECCANICA IDENTIFICATA:")
    print(f"  Punto Nominale:     f = {nom_point[0]} Hz, n = {nom_point[1]} RPM -> <Fz> = {nom_point[2]:+.2f} uN")
    print(f"  Punto Ottimo Assoluto: f = {opt_f} Hz, n = {opt_rpm} RPM -> <Fz> = {opt_res['mean_Fz_uN']:+.2f} uN (Picco Istantaneo: {opt_res['peak_Fz_uN']:+.2f} uN)")
    print(f"  Incremento Lift:    +{boost_pct:.1f}% di spinta media assiale")
    print(f"  Scorrimento Ottimo: f_slip = {opt_res['f_slip_Hz']:.1f} Hz")
    print(f"  Potenza Joule:      {opt_res['mean_Poule_W']:.3f} W (Efficienza: {opt_res['efficiency_uN_per_W']:.2f} uN/W)")
    print("=" * 80)
    
    # Generazione Figure a 300 DPI
    print("\n--- Generazione Figure a 300 DPI ---")
    plot_fig_01_superficie_2d(FREQUENCIES, RPMS, fz_matrix_uN, opt_point, nom_point)
    plot_fig_02_dispersione_slip(raw_data, FREQUENCIES)
    plot_fig_03_confronto_onde(nom_res, opt_res)
    
    # Esportazione JSON finale strutturato
    final_json = {
        "metadata": {
            "study": "sweep_paramedrico_2d_frequenza_vs_rpm_risonanza",
            "variant": "rotore_centrato_z0",
            "excitation": "Regime B co-rotante a 60 gradi",
            "mantle": "rete stirata anisotropa 30 gradi MATC cilindrico",
            "date": "2026-09-22",
            "author": "Alessandro Brescacin",
            "license": "CERN-OHL-S-2.0"
        },
        "resonance_optimal_operating_point": {
            "optimal_frequency_Hz": opt_f,
            "optimal_rpm": opt_rpm,
            "optimal_f_slip_Hz": opt_res['f_slip_Hz'],
            "mean_Fz_uN": opt_res['mean_Fz_uN'],
            "mean_Fz_mN": opt_res['mean_Fz_mN'],
            "peak_Fz_uN": opt_res['peak_Fz_uN'],
            "peak_Fz_mN": opt_res['peak_Fz_mN'],
            "mean_Joule_W": opt_res['mean_Poule_W'],
            "efficiency_uN_per_W": opt_res['efficiency_uN_per_W'],
            "nominal_benchmark_Fz_uN": nom_res['mean_Fz_uN'],
            "lift_boost_percentage": boost_pct
        },
        "full_scale_1e7_projection": {
            "scale_factor_current": 100.0,
            "scale_factor_Force": 10000.0,
            "projected_mean_Fz_N": opt_res['mean_Fz_mN'] * 1e-3 * 10000.0,
            "projected_peak_Fz_N": opt_res['peak_Fz_mN'] * 1e-3 * 10000.0,
            "projected_nominal_Fz_N": nom_res['mean_Fz_mN'] * 1e-3 * 10000.0
        },
        "grid_axes": {
            "frequencies_Hz": FREQUENCIES,
            "rpms": RPMS
        },
        "matrix_mean_Fz_uN": fz_matrix_uN.tolist(),
        "matrix_mean_Joule_W": pj_matrix_W.tolist(),
        "matrix_efficiency_uN_per_W": eta_matrix.tolist(),
        "matrix_f_slip_Hz": fslip_matrix.tolist(),
        "raw_points": raw_data
    }
    
    out_json_path = DATA_DIR / "sweep_risonanza_matrice.json"
    with open(out_json_path, "w", encoding="utf-8") as f:
        json.dump(final_json, f, indent=2)
        
    print(f"\nDataset finale salvato con successo: {out_json_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
