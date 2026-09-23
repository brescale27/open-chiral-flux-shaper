#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generazione Deliverable Grafici ad Alta Risoluzione (300 DPI):
Variante con Mantello Chiuso a Barattolo.

Figure:
1. fig_01_bilancio_4settori_joule_e_lift.png:
   - Forme d'onda temporali F_z(t) dei 4 settori (+30° vs -30°)
   - Ripartizione delle perdite Joule P_J nei 4 settori
   - Decoupling di F_bias e F_chiral per ciascun settore
2. fig_02_confronto_chiuso_vs_aperto_decoupled.png:
   - Confronto diretto Tubo Aperto vs Barattolo Chiuso
   - Confronto spinta grezza vs bias di mesh vs lift chirale puro
   - Abbattimento termico delle perdite Joule (-98.9%)

Autore: Alessandro Brescacin
Licenza: CERN-OHL-S-2.0
"""

import sys
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_DIR = SCRIPT_DIR.parent
DATA_DIR = VARIANT_DIR / "data"
FIGURES_DIR = VARIANT_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

JSON_FILE = DATA_DIR / "risultati_mantello_chiuso_bias_chiral.json"

def main():
    print("=" * 80)
    print("GENERAZIONE GRAFICI COMPARATIVI 300 DPI: MANTELLO CHIUSO A BARATTOLO")
    print("=" * 80)
    
    if not JSON_FILE.is_file():
        print(f"[ERRORE] File risultati non trovato: {JSON_FILE}")
        sys.exit(1)
        
    with open(JSON_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    time_ms = np.array(data["time_series_ms"])
    plus = data["nominal_plus30deg"]
    minus = data["specular_minus30deg"]
    dec = data["decoupling_by_sector"]
    comp = data["comparison_closed_vs_open"]
    
    plt.style.use('default')
    
    # -------------------------------------------------------------------------
    # FIGURA 1: Bilancio Elettrodinamico nei 4 Settori
    # -------------------------------------------------------------------------
    fig1 = plt.figure(figsize=(16, 5.5), dpi=300)
    fig1.patch.set_facecolor('#ffffff')
    
    # Panel A: Forme d'onda F_z(t)
    ax1 = fig1.add_subplot(1, 3, 1)
    ax1.set_facecolor('#fafafa')
    ax1.plot(time_ms, plus["lateral"]["Fz_series_uN"], color='#1f77b4', lw=1.8, label=r'Parete Laterale (+30°)')
    ax1.plot(time_ms, plus["top_cap"]["Fz_series_uN"], color='#2ca02c', lw=1.8, linestyle='--', label=r'Coperchio Sup (+30°)')
    ax1.plot(time_ms, plus["bottom_cap"]["Fz_series_uN"], color='#d62728', lw=1.8, linestyle=':', label=r'Coperchio Inf (+30°)')
    ax1.plot(time_ms, plus["inner_core"]["Fz_series_uN"], color='#9467bd', lw=1.8, label=r'Rotore/Nucleo Int (+30°)')
    ax1.plot(time_ms, plus["total"]["Fz_series_uN"], color='black', lw=2.2, label=r'Totale Macchina (+30°)')
    ax1.axhline(0, color='gray', linestyle='-', lw=0.8, alpha=0.6)
    ax1.set_title(r"$\mathbf{(a)\ Forme\ d'Onda\ Temporali\ F_z(t)\ per\ Settore}$", fontsize=11, pad=8)
    ax1.set_xlabel("Tempo [ms]", fontsize=10)
    ax1.set_ylabel(r"Forza di Lorentz $F_z\ [\mu\mathrm{N}]$", fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(loc='lower left', fontsize=8.5)
    
    # Panel B: Ripartizione Perdite Joule P_J
    ax2 = fig1.add_subplot(1, 3, 2)
    ax2.set_facecolor('#fafafa')
    sector_keys = ["lateral", "top_cap", "bottom_cap", "inner_core"]
    sector_names = ["Parete Lat.", "Cop. Sup.", "Cop. Inf.", "Nucleo Int."]
    pj_vals_uW = [dec[k]["Pj_plus30_mW"] * 1000.0 for k in sector_keys]
    bar_colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd']
    
    x_pos = np.arange(len(sector_keys))
    bars = ax2.bar(x_pos, pj_vals_uW, color=bar_colors, width=0.55, edgecolor='black', lw=0.8, alpha=0.85)
    for rect, val in zip(bars, pj_vals_uW):
        h = rect.get_height()
        ax2.annotate(f'{val:.2f} µW\n({val/(sum(pj_vals_uW)+1e-12)*100:.1f}%)',
                     xy=(rect.get_x() + rect.get_width() / 2, h),
                     xytext=(0, 4), textcoords="offset points",
                     ha='center', va='bottom', fontsize=8.5, fontweight='bold')
                     
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(sector_names, fontsize=9.5)
    ax2.set_ylabel(r"Perdite Joule $[\mu\mathrm{W}]$", fontsize=10)
    ax2.set_title(r"$\mathbf{(b)\ Distribuzione\ Dissipazione\ Termica\ P_J}$", fontsize=11, pad=8)
    ax2.grid(True, axis='y', linestyle='--', alpha=0.5)
    ax2.set_ylim(0, max(pj_vals_uW) * 1.35)
    
    # Panel C: Decoupling Mesh Bias vs Lift Chirale
    ax3 = fig1.add_subplot(1, 3, 3)
    ax3.set_facecolor('#fafafa')
    sec_all = ["lateral", "top_cap", "bottom_cap", "can_total", "inner_core", "total"]
    sec_labels_all = ["Parete", "Cop. Sup", "Cop. Inf", "Barattolo", "Nucleo", "TOTALE"]
    bias_vals = [dec[k]["F_bias_uN"] for k in sec_all]
    chiral_vals = [dec[k]["F_chiral_uN"] for k in sec_all]
    
    x = np.arange(len(sec_all))
    width = 0.35
    rects1 = ax3.bar(x - width/2, bias_vals, width, label=r'Mesh Bias $F_{\mathrm{bias}}$', color='#ff7f0e', edgecolor='black', lw=0.8, alpha=0.85)
    rects2 = ax3.bar(x + width/2, chiral_vals, width, label=r'Lift Chirale $F_{\mathrm{chiral}}$', color='#2ca02c', edgecolor='black', lw=0.8, alpha=0.85)
    ax3.axhline(0, color='black', lw=1.0)
    ax3.set_xticks(x)
    ax3.set_xticklabels(sec_labels_all, fontsize=8.5, rotation=15)
    ax3.set_ylabel(r"Forza Assiale $[\mu\mathrm{N}]$", fontsize=10)
    ax3.set_title(r"$\mathbf{(c)\ Decoupling:\ Mesh\ Bias\ vs\ Lift\ Chirale}$", fontsize=11, pad=8)
    ax3.grid(True, axis='y', linestyle='--', alpha=0.5)
    ax3.legend(loc='lower left', fontsize=8.5)
    
    fig1.suptitle(r"$\mathbf{Variante\ Mantello\ Chiuso\ a\ Barattolo\ (Z=\pm H/2)\ -\ Bilancio\ Elettrodinamico\ nei\ 4\ Settori}$",
                  fontsize=12, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig1_path = FIGURES_DIR / "fig_01_bilancio_4settori_joule_e_lift.png"
    plt.savefig(fig1_path, dpi=300, facecolor=fig1.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"[FIG 01] Salvata in: {fig1_path}")
    
    # -------------------------------------------------------------------------
    # FIGURA 2: Confronto Diretto Tubo Aperto vs Barattolo Chiuso
    # -------------------------------------------------------------------------
    fig2 = plt.figure(figsize=(15, 5.5), dpi=300)
    fig2.patch.set_facecolor('#ffffff')
    
    # Panel A: Confronto Lift Grezzo vs Bias vs Lift Chirale
    ax21 = fig2.add_subplot(1, 3, 1)
    ax21.set_facecolor('#fafafa')
    categories = [r'Grezzo $\langle F_z\rangle$', r'Mesh Bias $F_{\mathrm{bias}}$', r'Lift Puro $F_{\mathrm{chiral}}$']
    open_forces = [comp["open_tube_baseline"]["mean_Fz_raw_uN"],
                   comp["open_tube_baseline"]["F_bias_uN"],
                   comp["open_tube_baseline"]["F_chiral_pure_uN"]]
    closed_forces = [comp["closed_can_variant"]["mean_Fz_raw_uN"],
                     comp["closed_can_variant"]["F_bias_uN"],
                     comp["closed_can_variant"]["F_chiral_pure_uN"]]
                     
    x = np.arange(len(categories))
    w = 0.35
    r1 = ax21.bar(x - w/2, open_forces, w, label='Tubo Aperto (1200 RPM)', color='#1f77b4', edgecolor='black', lw=0.8, alpha=0.85)
    r2 = ax21.bar(x + w/2, closed_forces, w, label='Barattolo Chiuso (1200 RPM)', color='#e377c2', edgecolor='black', lw=0.8, alpha=0.85)
    ax21.axhline(0, color='black', lw=1.0)
    
    for r, v in zip(r1, open_forces):
        va = 'bottom' if v >= 0 else 'top'
        ax21.annotate(f'{v:+.2f} µN', xy=(r.get_x() + r.get_width()/2, v), xytext=(0, 3 if v>=0 else -10),
                      textcoords='offset points', ha='center', va=va, fontsize=8.5, fontweight='bold')
    for r, v in zip(r2, closed_forces):
        va = 'bottom' if v >= 0 else 'top'
        ax21.annotate(f'{v:+.2f} µN', xy=(r.get_x() + r.get_width()/2, v), xytext=(0, 3 if v>=0 else -10),
                      textcoords='offset points', ha='center', va=va, fontsize=8.5, fontweight='bold')
                      
    ax21.set_xticks(x)
    ax21.set_xticklabels(categories, fontsize=9.5)
    ax21.set_ylabel(r"Forza Assiale $[\mu\mathrm{N}]$", fontsize=10)
    ax21.set_title(r"$\mathbf{(a)\ Confronto\ Componenti\ di\ Spinta}$", fontsize=11, pad=8)
    ax21.grid(True, axis='y', linestyle='--', alpha=0.5)
    ax21.legend(loc='lower left', fontsize=8.5)
    
    # Panel B: Confronto Perdite Joule (Abbattimento 98.9%)
    ax22 = fig2.add_subplot(1, 3, 2)
    ax22.set_facecolor('#fafafa')
    cases_pj = ['Tubo Aperto', 'Barattolo Chiuso']
    pj_vals_mW = [comp["open_tube_baseline"]["mean_Pj_mW"], comp["closed_can_variant"]["mean_Pj_total_mW"]]
    bar_pj = ax22.bar(cases_pj, pj_vals_mW, color=['#ff7f0e', '#2ca02c'], width=0.45, edgecolor='black', lw=0.8, alpha=0.85)
    for rect, val in zip(bar_pj, pj_vals_mW):
        h = rect.get_height()
        ax22.annotate(f'{val:.3f} mW', xy=(rect.get_x() + rect.get_width()/2, h),
                      xytext=(0, 4), textcoords='offset points', ha='center', va='bottom', fontsize=9.5, fontweight='bold')
    ax22.annotate(r'$\mathbf{-98.9\%}$' + '\nAbbattimento Termico',
                  xy=(1, pj_vals_mW[1]), xytext=(0.55, 0.75),
                  textcoords='data', ha='center', va='center',
                  fontsize=10.5, fontweight='bold', color='darkgreen',
                  arrowprops=dict(arrowstyle='->', lw=1.5, color='darkgreen'))
    ax22.set_ylabel(r"Perdite Joule Totali $[\mathrm{mW}]$", fontsize=10)
    ax22.set_title(r"$\mathbf{(b)\ Abbattimento\ Perdite\ Joule\ P_J}$", fontsize=11, pad=8)
    ax22.grid(True, axis='y', linestyle='--', alpha=0.5)
    ax22.set_ylim(0, max(pj_vals_mW) * 1.25)
    
    # Panel C: Impatto Geometrico dei Coperchi (Schermatura Magnetica)
    ax23 = fig2.add_subplot(1, 3, 3)
    ax23.set_facecolor('#fafafa')
    elements_can = ['Parete Laterale', 'Coperchio Sup.', 'Coperchio Inf.']
    f_elements = [dec['lateral']['F_bias_uN'], dec['top_cap']['F_bias_uN'], dec['bottom_cap']['F_bias_uN']]
    p_elements = [dec['lateral']['Pj_plus30_mW'], dec['top_cap']['Pj_plus30_mW'], dec['bottom_cap']['Pj_plus30_mW']]
    
    x = np.arange(len(elements_can))
    ax23.bar(x, p_elements, color='#9467bd', width=0.45, edgecolor='black', lw=0.8, alpha=0.85)
    for i, v in enumerate(p_elements):
        ax23.annotate(f'{v*1000:.1f} µW', xy=(i, v), xytext=(0, 4), textcoords='offset points',
                      ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax23.set_xticks(x)
    ax23.set_xticklabels(elements_can, fontsize=9.5)
    ax23.set_ylabel(r"Dissipazione $[\mathrm{mW}]$", fontsize=10)
    ax23.set_title(r"$\mathbf{(c)\ Dissipazione\ Settori\ del\ Barattolo}$", fontsize=11, pad=8)
    ax23.grid(True, axis='y', linestyle='--', alpha=0.5)
    ax23.set_ylim(0, max(p_elements)*1.35 if max(p_elements)>0 else 0.01)
    
    fig2.suptitle(r"$\mathbf{Open\ Chiral\ Flux\ Shaper\ -\ Impatto\ dei\ Coperchi:\ Tubo\ Aperto\ vs\ Barattolo\ Chiuso}$",
                  fontsize=12, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig2_path = FIGURES_DIR / "fig_02_confronto_chiuso_vs_aperto_decoupled.png"
    plt.savefig(fig2_path, dpi=300, facecolor=fig2.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"[FIG 02] Salvata in: {fig2_path}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
