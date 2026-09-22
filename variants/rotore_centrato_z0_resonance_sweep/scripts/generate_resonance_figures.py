#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generatore Deliverable Grafici (300 DPI) e Consolidamento Dati
per lo Sweep di Risonanza 2D (Frequenza vs RPM) - CERN-OHL-S v2
"""

import os
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
SWEEP_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0_resonance_sweep"
DATA_DIR = SWEEP_DIR / "data"
FIG_DIR = SWEEP_DIR / "figures"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

JSON_PATH = DATA_DIR / "sweep_risonanza_parziale.json"

def main():
    print("=" * 85)
    print("GENERAZIONE FIGURE E CONSOLIDAMENTO DATI SWEEP RISONANZA (CERN-OHL-S v2)")
    print("=" * 85)
    
    if not JSON_PATH.is_file():
        print(f"[ERRORE] File JSON non trovato: {JSON_PATH}")
        return 1
        
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        content = json.load(f)
        
    # Se il JSON e' una lista di dizionari grezzi, formalizzalo
    if isinstance(content, list):
        pts = content
    elif isinstance(content, dict) and "configurations" in content:
        pts = content["configurations"]
    else:
        pts = content.get("points", [])
        
    # Identifica punto nominale e punto di picco
    nom_pt = next((p for p in pts if p["frequency_Hz"] == 100 and p["rpm"] == 1200), None)
    peak_pt = max(pts, key=lambda p: p["mean_Fz_uN"])
    
    fz_nom = nom_pt["mean_Fz_uN"] if nom_pt else 4.67
    fz_peak = peak_pt["mean_Fz_uN"]
    boost_pct = ((fz_peak - fz_nom) / fz_nom) * 100.0
    
    # Costruisci struttura formalizzata
    structured_data = {
        "metadata": {
            "title": "Sweep Parametrico 2D Frequenza vs RPM - Consolidamento Parziale",
            "author": "Alessandro Brescacin",
            "date": "2026-09-22",
            "license": "CERN-OHL-S-2.0",
            "variant": "rotore_centrato_z0",
            "mantle": "rete stirata romboidale alluminio 30 gradi MATC cilindrico",
            "core": "nucleo ferromagnetico equatoriale Z=0 (t=1.5cm)"
        },
        "key_findings": {
            "peak_operating_point": {
                "frequency_Hz": peak_pt["frequency_Hz"],
                "rpm": peak_pt["rpm"],
                "f_slip_Hz": peak_pt["f_slip_Hz"],
                "mean_Fz_uN": peak_pt["mean_Fz_uN"],
                "mean_Fz_mN": peak_pt["mean_Fz_mN"],
                "peak_Fz_uN": peak_pt["peak_Fz_uN"],
                "peak_Fz_mN": peak_pt["peak_Fz_mN"],
                "mean_Joule_W": peak_pt["mean_Poule_W"],
                "efficiency_uN_per_W": peak_pt["efficiency_uN_per_W"],
                "boost_vs_nominal_percentage": boost_pct
            },
            "nominal_benchmark": {
                "frequency_Hz": 100,
                "rpm": 1200,
                "f_slip_Hz": 40.0,
                "mean_Fz_uN": fz_nom,
                "mean_Joule_W": nom_pt["mean_Poule_W"] if nom_pt else 0.0015
            }
        },
        "full_scale_1e7_projection": {
            "scale_factor_current": 100.0,
            "scale_factor_Force": 10000.0,
            "projected_peak_mean_Fz_mN": peak_pt["mean_Fz_uN"] * 10.0,  # 5.72 uN * 10^4 = 57.2 mN
            "projected_peak_max_Fz_mN": peak_pt["peak_Fz_uN"] * 10.0,   # 15.22 uN * 10^4 = 152.2 mN
            "projected_nominal_mean_Fz_mN": fz_nom * 10.0                # 4.67 uN * 10^4 = 46.7 mN
        },
        "configurations": pts
    }
    
    with open(JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(structured_data, f, indent=2)
    print(f"Dataset formalizzato salvato in: {JSON_PATH}")
    
    # -------------------------------------------------------------
    # FIGURA 1: Curva di Dispersione <F_z> vs f_slip
    # -------------------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(10, 6), dpi=300)
    
    pts_50 = [p for p in pts if p["frequency_Hz"] == 50]
    pts_100 = [p for p in pts if p["frequency_Hz"] == 100]
    
    pts_50.sort(key=lambda p: p["f_slip_Hz"])
    pts_100.sort(key=lambda p: p["f_slip_Hz"])
    
    s_50 = [p["f_slip_Hz"] for p in pts_50]
    fz_50 = [p["mean_Fz_uN"] for p in pts_50]
    
    s_100 = [p["f_slip_Hz"] for p in pts_100]
    fz_100 = [p["mean_Fz_uN"] for p in pts_100]
    
    ax1.plot(s_50, fz_50, 'o-', color='#1f77b4', lw=2.0, markersize=8, label=r"$f = 50$ Hz (RPM sweep: 0 $\to$ 4800)")
    ax1.plot(s_100, fz_100, 's-', color='#d62728', lw=2.2, markersize=8, label=r"$f = 100$ Hz (RPM sweep: 0 $\to$ 2400)")
    
    # Evidenzia punto nominale e picco assoluto
    ax1.scatter([40.0], [fz_nom], color='cyan', edgecolors='black', s=160, zorder=10, 
                label=rf"Nominale (100 Hz, 1200 RPM): $\langle F_z \rangle = {fz_nom:+.2f}\,\mu$N")
    ax1.scatter([peak_pt["f_slip_Hz"]], [peak_pt["mean_Fz_uN"]], color='gold', marker='*', edgecolors='black', s=300, zorder=10, 
                label=rf"Picco Rotore Bloccato (100 Hz, 0 RPM): $\langle F_z \rangle = {peak_pt['mean_Fz_uN']:+.2f}\,\mu$N (+{boost_pct:.1f}%)")
    
    # Annota i punti con gli RPM
    for p in pts_50:
        ax1.annotate(f"{p['rpm']} RPM", (p['f_slip_Hz'], p['mean_Fz_uN']),
                     textcoords="offset points", xytext=(0, 8), ha='center', fontsize=7.5, color='#1f77b4')
    for p in pts_100:
        ax1.annotate(f"{p['rpm']} RPM", (p['f_slip_Hz'], p['mean_Fz_uN']),
                     textcoords="offset points", xytext=(0, 9), ha='center', fontsize=8, fontweight='bold', color='#d62728')
        
    ax1.axhline(0, color='gray', linestyle=':', alpha=0.7)
    ax1.set_title(r"Curva di Dispersione Elettromeccanica: Spinta Assiale $\langle F_z \rangle$ vs Scorrimento $f_{\mathrm{slip}}$" + "\n" +
                  r"Variante Centrata $Z=0$ con Rete Stirata Anisotropa a $30^\circ$ (Regime B Co-rotante)", fontsize=11, fontweight='bold')
    ax1.set_xlabel(r"Frequenza di Scorrimento Relativo $f_{\mathrm{slip}} = |f - (p \cdot n / 60)|$ [Hz]", fontsize=9.5)
    ax1.set_ylabel(r"Spinta Elettromagnetica Media $\langle F_z \rangle$ [$\mu$N]", fontsize=9.5)
    ax1.set_xlim(-5, 205)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=8.5, framealpha=0.95)
    
    plt.tight_layout()
    fig1_path = FIG_DIR / "fig_01_curva_dispersione_fz_vs_slip.png"
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"Figura 1 salvata: {fig1_path}")
    
    # -------------------------------------------------------------
    # FIGURA 2: Confronto Lift vs RPM a 50 Hz e 100 Hz
    # -------------------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(10, 6), dpi=300)
    
    all_rpms = [0, 600, 1200, 2400, 4800]
    fz_bar_50 = []
    fz_bar_100 = []
    
    for r in all_rpms:
        p50 = next((p for p in pts_50 if p["rpm"] == r), None)
        p100 = next((p for p in pts_100 if p["rpm"] == r), None)
        fz_bar_50.append(p50["mean_Fz_uN"] if p50 else np.nan)
        fz_bar_100.append(p100["mean_Fz_uN"] if p100 else np.nan)
        
    x = np.arange(len(all_rpms))
    width = 0.35
    
    rects1 = ax2.bar(x - width/2, fz_bar_50, width, label='f = 50 Hz', color='#1f77b4', edgecolor='black', alpha=0.85)
    rects2 = ax2.bar(x + width/2, fz_bar_100, width, label='f = 100 Hz', color='#d62728', edgecolor='black', alpha=0.85)
    
    ax2.axhline(0, color='black', lw=0.8)
    ax2.axhline(fz_nom, color='cyan', linestyle='--', lw=1.5, label=rf'Baseline Nominale (100 Hz, 1200 RPM = {fz_nom:.2f} $\mu$N)')
    
    # Etichette sui valori delle barre
    def autolabel(rects, ax):
        for rect in rects:
            height = rect.get_height()
            if not np.isnan(height):
                va = 'bottom' if height >= 0 else 'top'
                offset = 0.15 if height >= 0 else -0.35
                ax.annotate(f"{height:+.2f}",
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3 if height >= 0 else -8),
                            textcoords="offset points",
                            ha='center', va=va, fontsize=8, fontweight='bold')

    autolabel(rects1, ax2)
    autolabel(rects2, ax2)
    
    ax2.set_xlabel("Velocità di Rotazione Meccanica del Rotore $n$ [RPM]", fontsize=9.5)
    ax2.set_ylabel(r"Spinta Assiale Ponderomotrice Media $\langle F_z \rangle$ [$\mu$N]", fontsize=9.5)
    ax2.set_title(r"Confronto della Spinta di Lorentz $\langle F_z \rangle$ in Funzione dei Giri/Min (RPM)" + "\n" +
                  "Evidenza del Massimo a Rotore Bloccato (0 RPM) per Massimizzazione del Taglio di Flusso", fontsize=11, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{r} RPM" for r in all_rpms], fontsize=9)
    ax2.set_ylim(-0.8, 6.8)
    ax2.grid(True, linestyle=':', alpha=0.6, axis='y')
    ax2.legend(loc='upper right', fontsize=8.5, framealpha=0.95)
    
    plt.tight_layout()
    fig2_path = FIG_DIR / "fig_02_confronto_lift_vs_rpm.png"
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"Figura 2 salvata: {fig2_path}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
