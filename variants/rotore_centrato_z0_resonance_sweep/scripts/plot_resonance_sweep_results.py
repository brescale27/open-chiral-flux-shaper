#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASE 3: Generazione Figure ad Alta Risoluzione (300 DPI) per lo Sweep 2D
e la Validazione di Integrità Numerica (Mesh Bias & Maxwell Stress Tensor).

Figure generate in variants/rotore_centrato_z0_resonance_sweep/figures/:
1. fig_01_superficie_risonanza_lift_2d.png:
   Superficie 3D e mappa di contorno 2D di <F_z>(f, n) con evidenziazione del picco globale.
2. fig_02_curva_dispersione_vs_slip.png:
   Curva di dispersione <F_z> in funzione di f_slip = |f - 0.05*n| con fit di risonanza induttiva.
3. fig_03_validazione_bias_e_tensore_maxwell.png:
   Validazione rigorosa del bias di mesh e confronto tra Lorentz (J x B), MST e test chirale a +-30°.

Autore: Alessandro Brescacin
Licenza: CERN-OHL-S-2.0
"""

import os
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).resolve().parent
SWEEP_DIR = SCRIPT_DIR.parent
DATA_DIR = SWEEP_DIR / "data"
FIGURES_DIR = SWEEP_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def plot_fig03_validation_bias_mst():
    """Figura 3: Validazione Bias e Confronto Maxwell Stress Tensor vs Lorentz"""
    mst_json = DATA_DIR / "validazione_mst_chiral_bias.json"
    if not mst_json.is_file():
        print(f"[SKIP Fig 03] File {mst_json} non ancora disponibile.")
        return

    with open(mst_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    plus30 = data["nominal_plus30deg"]
    minus30 = data["specular_minus30deg"]
    bias_data = data["bias_decoupling"]["total_active_domain"]
    mst_val = data["maxwell_stress_tensor_validation"]

    time_ms = np.linspace(0.5, 10.0, len(plus30["fz_lorentz_active_uN"]))

    plt.style.use('default')
    fig = plt.figure(figsize=(16, 5.2), dpi=300)
    fig.patch.set_facecolor('#ffffff')

    # Panel A: Forme d'onda Lorentz vs MST
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.set_facecolor('#fafafa')
    ax1.plot(time_ms, plus30["fz_lorentz_active_uN"], 'b-', lw=2.2, label=r'Lorentz $\int (\vec{J}\times\vec{B})_z dV$')
    ax1.plot(time_ms, plus30["fz_mst_uN"], 'r--', lw=2.0, label=r'MST $\oint_{\partial\Omega} T_z dA$ ($R=8\mathrm{cm}$)')
    ax1.axhline(plus30["mean_Fz_active_uN"], color='navy', linestyle=':', lw=1.5,
                label=rf'$\langle F_{{z,\mathrm{{Lor}}}}\rangle = {plus30["mean_Fz_active_uN"]:+.2f}\,\mu\mathrm{{N}}$')
    ax1.axhline(plus30["mean_Fz_mst_uN"], color='darkred', linestyle=':', lw=1.5,
                label=rf'$\langle F_{{z,\mathrm{{MST}}}}\rangle = {plus30["mean_Fz_mst_uN"]:+.2f}\,\mu\mathrm{{N}}$')
    ax1.axhline(0, color='gray', linestyle='-', lw=0.8, alpha=0.6)
    ax1.set_title(r"$\mathbf{(a)\ Lorentz\ vs\ Maxwell\ Stress\ Tensor\ (MST)}$", fontsize=11, pad=8)
    ax1.set_xlabel("Tempo [ms]", fontsize=10)
    ax1.set_ylabel(r"Forza Assiale $F_z\ [\mu\mathrm{N}]$", fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(loc='lower right', fontsize=8.5)

    # Panel B: Inversione Chirale (+30° vs -30°)
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.set_facecolor('#fafafa')
    ax2.plot(time_ms, plus30["fz_lorentz_active_uN"], 'b-', lw=2.0, label=r'$\theta = +30^\circ$ (Nominale)')
    ax2.plot(time_ms, minus30["fz_lorentz_active_uN"], 'm-', lw=2.0, label=r'$\theta = -30^\circ$ (Speculare)')
    f_bias = bias_data["F_bias_uN"]
    f_chiral = bias_data["F_chiral_pure_uN"]
    ax2.axhline(f_bias, color='darkorange', linestyle='--', lw=1.8,
                label=rf'Mesh Bias $F_{{\mathrm{{bias}}}} = {f_bias:+.2f}\,\mu\mathrm{{N}}$')
    ax2.axhline(f_chiral, color='green', linestyle='-.', lw=1.8,
                label=rf'Lift Chirale Puro $F_{{\mathrm{{chiral}}}} = {f_chiral:+.2f}\,\mu\mathrm{{N}}$')
    ax2.axhline(0, color='gray', linestyle='-', lw=0.8, alpha=0.6)
    ax2.set_title(r"$\mathbf{(b)\ Inversione\ di\ Parit\grave{a}\ Chirale\ (\pm 30^\circ)}$", fontsize=11, pad=8)
    ax2.set_xlabel("Tempo [ms]", fontsize=10)
    ax2.set_ylabel(r"Forza Assiale $F_z\ [\mu\mathrm{N}]$", fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='lower right', fontsize=8.5)

    # Panel C: Decoupling Bar Chart
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.set_facecolor('#fafafa')
    bars = [
        (r'$\langle F_z(+30^\circ)\rangle$', plus30["mean_Fz_active_uN"], '#1f77b4'),
        (r'$\langle F_z(-30^\circ)\rangle$', minus30["mean_Fz_active_uN"], '#9467bd'),
        (r'Mesh Bias $F_{\mathrm{bias}}$', f_bias, '#ff7f0e'),
        (r'Lift Puro $F_{\mathrm{chiral}}$', f_chiral, '#2ca02c'),
        (r'MST $\langle F_{z}\rangle$', plus30["mean_Fz_mst_uN"], '#d62728')
    ]
    labels = [b[0] for b in bars]
    vals = [b[1] for b in bars]
    colors = [b[2] for b in bars]
    x_pos = np.arange(len(bars))
    rects = ax3.bar(x_pos, vals, color=colors, width=0.55, edgecolor='black', lw=0.8, alpha=0.85)
    ax3.axhline(0, color='black', lw=1.0)
    for rect, val in zip(rects, vals):
        h = rect.get_height()
        va = 'bottom' if h >= 0 else 'top'
        ax3.annotate(f'{val:+.2f} µN',
                     xy=(rect.get_x() + rect.get_width() / 2, h),
                     xytext=(0, 3 if h >= 0 else -10),
                     textcoords="offset points",
                     ha='center', va=va, fontsize=9, fontweight='bold')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(labels, fontsize=8.5, rotation=15)
    ax3.set_ylabel(r"Spinta Assiale $[\mu\mathrm{N}]$", fontsize=10)
    ax3.set_title(r"$\mathbf{(c)\ Decoupling\ Componenti\ e\ Bilancio\ Forze}$", fontsize=11, pad=8)
    ax3.grid(True, axis='y', linestyle='--', alpha=0.5)

    fig.suptitle(r"$\mathbf{Open\ Chiral\ Flux\ Shaper\ -\ Validazione\ Integrit\grave{a}\ Numerica\ e\ Tensore\ di\ Maxwell\ (MST)}$",
                 fontsize=12, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig3_path = FIGURES_DIR / "fig_03_validazione_bias_e_tensore_maxwell.png"
    plt.savefig(fig3_path, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"[FIG 03] Salvata in: {fig3_path}")

def plot_fig01_fig02_sweep_matrix():
    """Figure 1 e 2: Superficie di risonanza 2D e curva di dispersione vs slip"""
    matrix_json = DATA_DIR / "sweep_risonanza_matrice.json"
    if not matrix_json.is_file():
        # Fallback al parziale se matrice completa non ancora generata
        matrix_json = DATA_DIR / "sweep_risonanza_parziale.json"
        if not matrix_json.is_file():
            print(f"[SKIP Fig 01/02] Nessun dataset sweep trovato in {DATA_DIR}.")
            return

    with open(matrix_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    configs = data.get("configurations", data if isinstance(data, list) else [])
    if len(configs) == 0:
        print("[SKIP Fig 01/02] Nessuna configurazione nel dataset.")
        return

    # Estrai f, rpm, fz, pj, slip
    f_list = sorted(list(set(c["frequency_Hz"] for c in configs)))
    rpm_list = sorted(list(set(c["rpm"] for c in configs)))

    Fz_grid = np.full((len(rpm_list), len(f_list)), np.nan)
    Slip_grid = np.full((len(rpm_list), len(f_list)), np.nan)
    Pj_grid = np.full((len(rpm_list), len(f_list)), np.nan)

    slip_points = []
    fz_points = []
    f_colors = []

    for c in configs:
        f_val = c["frequency_Hz"]
        n_val = c["rpm"]
        fz = c.get("mean_Fz_uN", c.get("mean_Fz_raw_uN", 0.0))
        pj = c.get("mean_Poule_W", c.get("mean_Pj_W", 0.0))
        slip = c.get("f_slip_Hz", abs(f_val - 0.05 * n_val))

        i = rpm_list.index(n_val)
        j = f_list.index(f_val)
        Fz_grid[i, j] = fz
        Slip_grid[i, j] = slip
        Pj_grid[i, j] = pj

        slip_points.append(slip)
        fz_points.append(fz)
        f_colors.append(f_val)

    slip_points = np.array(slip_points)
    fz_points = np.array(fz_points)
    f_colors = np.array(f_colors)

    # -------------------------------------------------------------
    # FIGURA 1: Superficie di Risonanza 2D (f x RPM)
    # -------------------------------------------------------------
    plt.style.use('default')
    fig = plt.figure(figsize=(15, 6.2), dpi=300)
    fig.patch.set_facecolor('#ffffff')

    F_mesh, RPM_mesh = np.meshgrid(f_list, rpm_list)

    # Subplot 1: 3D Surface
    ax_3d = fig.add_subplot(1, 2, 1, projection='3d')
    ax_3d.set_facecolor('#ffffff')
    # Riempi eventuali nan per interpolazione griglia
    from scipy.interpolate import griddata
    valid_mask = ~np.isnan(Fz_grid)
    if np.sum(valid_mask) >= 6:
        pts_valid = np.column_stack([F_mesh[valid_mask], RPM_mesh[valid_mask]])
        vals_valid = Fz_grid[valid_mask]
        F_dense = np.linspace(min(f_list), max(f_list), 40)
        RPM_dense = np.linspace(min(rpm_list), max(rpm_list), 40)
        FD, RD = np.meshgrid(F_dense, RPM_dense)
        ZD = griddata(pts_valid, vals_valid, (FD, RD), method='cubic')
        # Fallback linear per bordi nan
        ZD_lin = griddata(pts_valid, vals_valid, (FD, RD), method='linear')
        ZD[np.isnan(ZD)] = ZD_lin[np.isnan(ZD)]

        surf = ax_3d.plot_surface(FD, RD, ZD, cmap=cm.viridis, edgecolor='none', alpha=0.88)
        ax_3d.scatter(pts_valid[:, 0], pts_valid[:, 1], vals_valid, color='red', s=35, zorder=5, label='Punti Simulati FEM')
        # Evidenzia picco
        idx_peak = np.nanargmax(vals_valid)
        ax_3d.scatter([pts_valid[idx_peak, 0]], [pts_valid[idx_peak, 1]], [vals_valid[idx_peak]],
                      color='gold', s=120, edgecolors='black', lw=1.5, zorder=10,
                      label=rf'Picco: ${vals_valid[idx_peak]:.2f}\,\mu\mathrm{{N}}$')

    ax_3d.set_title(r"$\mathbf{(a)\ Superficie\ 3D\ di\ Risonanza:\ }\langle F_z \rangle(f, \mathrm{RPM})$", fontsize=11, pad=12)
    ax_3d.set_xlabel("Frequenza $f$ [Hz]", fontsize=9.5)
    ax_3d.set_ylabel("Velocità $n$ [RPM]", fontsize=9.5)
    ax_3d.set_zlabel(r"$\langle F_z \rangle\ [\mu\mathrm{N}]$", fontsize=9.5)
    ax_3d.view_init(elev=28, azim=-125)
    ax_3d.legend(loc='upper left', fontsize=8.5)

    # Subplot 2: Contour Heatmap 2D
    ax_2d = fig.add_subplot(1, 2, 2)
    ax_2d.set_facecolor('#fafafa')
    if np.sum(valid_mask) >= 6:
        cnt = ax_2d.contourf(FD, RD, ZD, levels=25, cmap=cm.viridis)
        cbar = fig.colorbar(cnt, ax=ax_2d, pad=0.03)
        cbar.set_label(r"Spinta Assiale Media $\langle F_z \rangle\ [\mu\mathrm{N}]$", fontsize=9.5)
        # Linee di contorno iso-forza
        lines = ax_2d.contour(FD, RD, ZD, levels=10, colors='white', alpha=0.35, linewidths=0.8)
        ax_2d.clabel(lines, inline=True, fontsize=8, fmt='%.1f')

        ax_2d.scatter(pts_valid[:, 0], pts_valid[:, 1], c='black', s=25, zorder=5, label='Simulazioni')
        ax_2d.scatter([pts_valid[idx_peak, 0]], [pts_valid[idx_peak, 1]],
                      color='gold', s=130, edgecolors='black', lw=1.5, zorder=10,
                      label=rf'Massimo Assoluto ({pts_valid[idx_peak, 0]:.0f} Hz, {pts_valid[idx_peak, 1]:.0f} RPM)')

    ax_2d.set_title(r"$\mathbf{(b)\ Mappa\ di\ Contorno\ Isotropica\ di\ Lift}$", fontsize=11, pad=8)
    ax_2d.set_xlabel("Frequenza Elettrica $f$ [Hz]", fontsize=10)
    ax_2d.set_ylabel("Velocità Meccanica $n$ [RPM]", fontsize=10)
    ax_2d.grid(True, linestyle='--', alpha=0.4)
    ax_2d.legend(loc='upper right', fontsize=8.5)

    fig.suptitle(r"$\mathbf{Open\ Chiral\ Flux\ Shaper\ -\ Mappatura\ 2D\ della\ Spinta\ Assiale\ (Frequenza\ vs\ RPM)}$",
                 fontsize=12, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig1_path = FIGURES_DIR / "fig_01_superficie_risonanza_lift_2d.png"
    plt.savefig(fig1_path, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"[FIG 01] Salvata in: {fig1_path}")

    # -------------------------------------------------------------
    # FIGURA 2: Curva di Dispersione vs Slip Frequency f_slip
    # -------------------------------------------------------------
    fig2 = plt.figure(figsize=(10, 6.2), dpi=300)
    fig2.patch.set_facecolor('#ffffff')
    ax = fig2.add_subplot(1, 1, 1)
    ax.set_facecolor('#fafafa')

    # Fit di risonanza asintotica: F(s) = A * s / (1 + (s/s0)^2) o sigmoidale
    def induction_model(s, A, s0):
        return (A * s) / (1.0 + (s / s0)**2)

    sort_idx = np.argsort(slip_points)
    s_sorted = slip_points[sort_idx]
    f_sorted = fz_points[sort_idx]

    try:
        popt, _ = curve_fit(induction_model, s_sorted, f_sorted, p0=[0.1, 80.0], maxfev=5000)
        s_dense = np.linspace(0, max(slip_points) * 1.05, 200)
        f_fit = induction_model(s_dense, *popt)
        ax.plot(s_dense, f_fit, 'r-', lw=2.2, label=rf'Modello Risonanza Induttiva ($s_0 = {popt[1]:.1f}\,\mathrm{{Hz}}$)')
        peak_fit_slip = popt[1]
        peak_fit_f = induction_model(peak_fit_slip, *popt)
        ax.axvline(peak_fit_slip, color='darkred', linestyle='--', lw=1.2, alpha=0.7,
                   label=rf'Picco Teorico $f_{{\mathrm{{slip,opt}}}} = {peak_fit_slip:.1f}\,\mathrm{{Hz}}$')
    except Exception:
        pass

    # Scatter punti con colorbar frequenza
    sc = ax.scatter(slip_points, fz_points, c=f_colors, cmap=cm.plasma, s=65, edgecolors='black', lw=0.9, zorder=5)
    cbar = fig2.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Frequenza Elettrica $f$ [Hz]", fontsize=9.5)

    ax.axhline(0, color='gray', linestyle='-', lw=0.8, alpha=0.7)
    ax.set_title(r"$\mathbf{Curva\ di\ Dispersione\ Elettrodinamica:\ }\langle F_z \rangle\ \mathbf{in\ Funzione\ dello\ Scorrimento\ } f_{\mathrm{slip}}$",
                 fontsize=12, pad=10)
    ax.set_xlabel(r"Frequenza di Scorrimento Relativo $f_{\mathrm{slip}} = |f - 0.05\cdot n|\ [\mathrm{Hz}]$", fontsize=10)
    ax.set_ylabel(r"Spinta Assiale Ponderomotrice $\langle F_z \rangle\ [\mu\mathrm{N}]$", fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower right', fontsize=9.5)

    plt.tight_layout()
    fig2_path = FIGURES_DIR / "fig_02_curva_dispersione_vs_slip.png"
    plt.savefig(fig2_path, dpi=300, facecolor=fig2.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"[FIG 02] Salvata in: {fig2_path}")

def main():
    print("=" * 80)
    print("GENERAZIONE DELIVERABLE GRAFICI 300 DPI (FASE 3)")
    print("=" * 80)
    plot_fig03_validation_bias_mst()
    plot_fig01_fig02_sweep_matrix()

if __name__ == "__main__":
    main()
