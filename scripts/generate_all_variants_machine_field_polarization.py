#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================================================================
OPEN CHIRAL FLUX SHAPER - ALL VARIANTS MACHINE, FIELD & POLARIZATION GENERATOR
========================================================================================
Genera per CIASCUNA variante del framework la tavola diagnostica dedicata (300 DPI):
  variants/<variant_id>/figures/fig_machine_field_polarization.png

Ciascuna tavola comprende 3 pannelli coordinati:
  1. Pannello A: Geometria 3D della Macchina (Rotori, Bobine Statoriche, Mantello/Collimatore)
  2. Pannello B: Distribuzione del Campo Magnetico B (Linee di flusso e Mappa d'Induzione)
  3. Pannello C: Odografo di Polarizzazione Trasversa (B_theta vs B_phi, Stokes s3, AR, CP)

Inoltre genera la Tavola Sinottica Master (Figura 43, 300 DPI):
  figures/fig_43_all_variants_machine_field_polarization_matrix.png
che mette a confronto sinottico tutte le 8 varianti in una matrice complessiva.

Licenza: CERN-OHL-S-2.0 / Apache 2.0
Autore: Alessandro Brescacin per Open Chiral Flux Shaper
========================================================================================
"""

import os
import sys
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

ROOT_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT_DIR / "figures"
VAR_DIR = ROOT_DIR / "variants"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Definizione metrologica delle 8 varianti primarie
VARIANTS_CONFIG = [
    {
        "id": "rotore_centrato_mantello_chiuso",
        "name": "Single Rotor Baseline (Z-axis)",
        "r_coils": 55.0,
        "n_coils": 4,
        "type": "single_rotor",
        "mantle": "Cilindro Passivo Chiuso (mu_r=1)",
        "collimator": False,
        "b_max_gap": 5.99,
        "purity_cp": 15.9,
        "s3": 0.159,
        "ar_db": 21.94,
        "status_cp": "FAIL (Dipolo Planare Lineare)",
        "color": "#94a3b8"
    },
    {
        "id": "gabbia_sferica_doppio_rotore_90deg",
        "name": "Dual Continuous 90° NPNPNP",
        "r_coils": 55.0,
        "n_coils": 24,
        "type": "dual_continuous",
        "mantle": "Sferico Triplo Strato (+-30°)",
        "collimator": False,
        "b_max_gap": 11.85,
        "purity_cp": 94.0,
        "s3": 0.940,
        "ar_db": 3.09,
        "status_cp": "NEAR-CIRCULAR (Continuo 3D)",
        "color": "#818cf8"
    },
    {
        "id": "gabbia_sferica_fibonacci_24x24",
        "name": "Fibonacci 24x24 (Pisano mod 9)",
        "r_coils": 55.0,
        "n_coils": 24,
        "type": "fibonacci",
        "mantle": "Spirale Aurea Modulata",
        "collimator": False,
        "b_max_gap": 10.91,
        "purity_cp": 90.0,
        "s3": 0.900,
        "ar_db": 4.06,
        "status_cp": "ELLIPTIC (Guida Discreta)",
        "color": "#f59e0b"
    },
    {
        "id": "gabbia_sferica_triskelion_esagramma_24pulse",
        "name": "Triskelion 3-Lobe Hexagram",
        "r_coils": 55.0,
        "n_coils": 24,
        "type": "triskelion",
        "mantle": "Esagramma 3 Lobi (120°)",
        "collimator": False,
        "b_max_gap": 10.12,
        "purity_cp": 77.9,
        "s3": 0.779,
        "ar_db": 6.40,
        "status_cp": "TRIFOLIO (Modo m=3)",
        "color": "#ec4899"
    },
    {
        "id": "gabbia_sferica_doppio_gruppo_90deg_48coils",
        "name": "Dual Orthogonal 90° (48 Coils)",
        "r_coils": 55.0,
        "n_coils": 48,
        "type": "dual_48coils",
        "mantle": "Gabbia Anisotropa Sferica (+-30°)",
        "collimator": False,
        "b_max_gap": 14.37,
        "purity_cp": 95.5,
        "s3": 0.955,
        "ar_db": 2.67,
        "status_cp": "PASS (IEEE AR <= 3 dB)",
        "color": "#38bdf8"
    },
    {
        "id": "gabbia_sferica_chiral_wpt_actuator",
        "name": "Chiral WPT Benchtop Actuator",
        "r_coils": 55.0,
        "n_coils": 48,
        "type": "wpt_benchtop",
        "mantle": "Nucleo PEEK + Mantello Calibrato",
        "collimator": False,
        "b_max_gap": 12.64,
        "purity_cp": 95.5,
        "s3": 0.955,
        "ar_db": 2.67,
        "status_cp": "PASS (Prototipo 18.5 W)",
        "color": "#22c55e"
    },
    {
        "id": "gabbia_sferica_chiral_diode_asymmetric_pulse",
        "name": "Chiral Diode (+45°/+15°/-22.5°)",
        "r_coils": 55.0,
        "n_coils": 48,
        "type": "chiral_diode",
        "mantle": "Gradiente Asimmetrico Triplo Strato",
        "collimator": False,
        "b_max_gap": 16.77,
        "purity_cp": 99.98,
        "s3": 0.9998,
        "ar_db": 0.15,
        "status_cp": "PASS (Ultra-Pure CP Record)",
        "color": "#f97316"
    },
    {
        "id": "gabbia_sferica_inner_coils_copper_collimator",
        "name": "Inner Coils (28mm) + Cu Collimator",
        "r_coils": 28.0,
        "n_coils": 24,
        "type": "copper_collimator",
        "mantle": "Gabbia + Tubo Rame (L=200mm)",
        "collimator": True,
        "b_max_gap": 48.65,
        "purity_cp": 96.8,
        "s3": 0.968,
        "ar_db": 2.25,
        "status_cp": "PASS (Fascio Guidato 103x)",
        "color": "#f43f5e"
    }
]

def render_3d_machine(ax, v):
    """Disegna la geometria 3D schematica della macchina per la specifica variante."""
    ax.set_facecolor('#0f172a')
    
    # 1. Rotore(i)
    u = np.linspace(0, 2 * np.pi, 25)
    v_ang = np.linspace(0, np.pi, 25)
    r_rot = 0.022
    
    # Rotore principale asse Z
    xr = r_rot * np.outer(np.cos(u), np.sin(v_ang))
    yr = r_rot * np.outer(np.sin(u), np.sin(v_ang))
    zr = r_rot * np.outer(np.ones(np.size(u)), np.cos(v_ang))
    ax.plot_surface(xr, yr, zr, color='#64748b', alpha=0.6, edgecolor='none')
    
    # Se doppio rotore, rotore ausiliario asse X
    if "Dual" in v["name"] or "doppio_rotore" in v["id"]:
        ax.plot_surface(zr, yr, xr, color='#475569', alpha=0.5, edgecolor='none')
        ax.plot([-0.035, 0.035], [0, 0], [0, 0], color='#cbd5e1', lw=1.5, ls='--')
    
    # Albero motore Z
    ax.plot([0, 0], [0, 0], [-0.045, 0.045], color='#cbd5e1', lw=2.0)
    
    # 2. Bobine Statoriche
    rc = v["r_coils"] * 1e-3
    n_c = v["n_coils"]
    t_coil = np.linspace(0, 2 * np.pi, 20)
    r_coil_w = 0.007
    
    if v["type"] == "single_rotor":
        # 4 bobine cilindriche a 90° su asse equatoriale
        for i in range(4):
            phi = i * np.pi / 2.0
            cx = rc * np.cos(phi)
            cy = rc * np.sin(phi)
            cz = np.linspace(-0.025, 0.025, 20)
            for z_val in [-0.015, 0.0, 0.015]:
                ax.plot(cx + r_coil_w * np.cos(t_coil), cy + r_coil_w * np.sin(t_coil),
                        np.full_like(t_coil, z_val), color='#e2e8f0', lw=1.2)
    elif v["type"] == "triskelion":
        # 3 lobi a 120° con grappoli di bobine
        for l in range(3):
            phi_base = l * 2.0 * np.pi / 3.0
            for k in range(8):
                d_phi = (k - 3.5) * 0.08
                phi = phi_base + d_phi
                cx = rc * np.cos(phi)
                cy = rc * np.sin(phi)
                cz = 0.018 * np.sin(k * np.pi / 4.0)
                ax.plot(cx + r_coil_w * np.cos(t_coil) * 0.5,
                        cy + r_coil_w * np.sin(t_coil) * 0.5,
                        cz + np.linspace(-0.005, 0.005, 20), color=v["color"], lw=1.0)
    elif v["type"] == "copper_collimator":
        # Bobine compatte a r = 28 mm
        for i in range(24):
            phi = i * 2.0 * np.pi / 24.0
            theta = np.pi/2.0 + 0.3 * np.sin(2.0 * phi)
            cx = rc * np.sin(theta) * np.cos(phi)
            cy = rc * np.sin(theta) * np.sin(phi)
            cz = rc * np.cos(theta)
            ax.plot([cx], [cy], [cz], marker='o', markersize=4, color='#f43f5e')
            
        # Disegno del Tubo di Rame Coassiale (da z = 55 a 255 mm, r = 38-43 mm)
        z_tube = np.linspace(0.055, 0.255, 15)
        th_t = np.linspace(0, 2*np.pi, 25)
        ZT, THT = np.meshgrid(z_tube, th_t)
        XT = 0.040 * np.cos(THT)
        YT = 0.040 * np.sin(THT)
        ax.plot_surface(XT, YT, ZT, color='#ea580c', alpha=0.35, edgecolor='#c2410c', lw=0.4)
        ax.plot([0, 0], [0, 0], [0.055, 0.260], color='#38bdf8', lw=1.8, ls=':')
        ax.text(0, 0, 0.270, "Uscita Tubo (255mm)", color='#38bdf8', fontsize=7, ha='center')
    else:
        # Gabbia sferica standard (24 o 48 bobine a r = 55 mm)
        n_pts = min(n_c, 36)
        for i in range(n_pts):
            phi = i * 2.0 * np.pi / n_pts
            theta = np.pi/2.0 + 0.45 * np.cos(3.0 * phi)
            cx = rc * np.sin(theta) * np.cos(phi)
            cy = rc * np.sin(theta) * np.sin(phi)
            cz = rc * np.cos(theta)
            c_coil = v["color"] if i % 2 == 0 else '#38bdf8'
            ax.plot([cx], [cy], [cz], marker='o', markersize=3.5, color=c_coil)
            
    # 3. Mantello / Profilo Esterno
    if not v["collimator"]:
        r_m = 0.050
        xm = r_m * np.outer(np.cos(u), np.sin(v_ang))
        ym = r_m * np.outer(np.sin(u), np.sin(v_ang))
        zm = r_m * np.outer(np.ones(np.size(u)), np.cos(v_ang))
        ax.plot_wireframe(xm, ym, zm, color=v["color"], alpha=0.15, lw=0.5)

    ax.set_xlim(-0.065, 0.065)
    ax.set_ylim(-0.065, 0.065)
    z_lim = 0.280 if v["collimator"] else 0.065
    ax.set_zlim(-0.065, z_lim)
    ax.view_init(elev=22, azim=45)
    ax.axis('off')

def render_b_field_map(ax, v):
    """Disegna la mappa 2D d'intensità e le linee di flusso magnetico B."""
    ax.set_facecolor('#0f172a')
    
    x = np.linspace(-60, 60, 100)
    y = np.linspace(-60, 60, 100)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2) + 1e-6
    Phi = np.arctan2(Y, X)
    
    # Calcolo forma d'onda del campo specifica per la variante
    b_scale = v["b_max_gap"]
    
    if v["type"] == "single_rotor":
        # Dipolo puro orientato lungo X
        Bx = b_scale * (3.0 * X * Y / (R**4 + 200))
        By = b_scale * ((2.0 * Y**2 - X**2) / (R**4 + 200))
        B_mag = np.sqrt(Bx**2 + By**2)
        B_mag = np.clip(B_mag * 120.0, 0, b_scale)
    elif v["type"] == "triskelion":
        # Cuspide a 3 lobi
        m3_mod = 1.0 + 0.45 * np.cos(3.0 * Phi)
        B_mag = b_scale * np.exp(-((R - 35)**2) / 250.0) * m3_mod
        Bx = -B_mag * np.sin(Phi + np.pi/6)
        By =  B_mag * np.cos(Phi + np.pi/6)
    elif v["type"] == "fibonacci":
        # Modulazione ellittica aurea
        ellip = 1.0 + 0.3 * np.cos(2.0 * (Phi - np.radians(32)))
        ripple = 1.0 + 0.08 * np.cos(24.0 * Phi)
        B_mag = b_scale * np.exp(-((R - 45)**2) / 300.0) * ellip * ripple
        Bx = -B_mag * np.sin(Phi)
        By =  B_mag * np.cos(Phi)
    elif v["type"] == "chiral_diode":
        # Asimmetria direzionale marcata
        asymm = 1.0 + 0.4 * np.sin(Phi + np.pi/4)
        B_mag = b_scale * np.exp(-((R - 42)**2) / 220.0) * asymm
        Bx = -B_mag * (np.sin(Phi) - 0.2 * np.cos(Phi))
        By =  B_mag * (np.cos(Phi) + 0.2 * np.sin(Phi))
    elif v["type"] == "copper_collimator":
        # Campo concentrato vicino al rotore (R < 35 mm) e fascio guidato
        B_mag = b_scale * np.exp(-((R - 26)**2) / 120.0)
        Bx = -B_mag * np.sin(Phi)
        By =  B_mag * np.cos(Phi)
    else:
        # Modo rotante continuo quasi isotropo
        B_mag = b_scale * np.exp(-((R - 45)**2) / 280.0)
        Bx = -B_mag * np.sin(Phi)
        By =  B_mag * np.cos(Phi)
        
    cf = ax.contourf(X, Y, B_mag, levels=30, cmap='plasma', alpha=0.9)
    # Streamlines del campo vettoriale
    ax.streamplot(X, Y, Bx, By, color='#ffffff', density=0.85, linewidth=0.7, arrowsize=0.9)
    
    # Contorno della macchina
    r_circ = v["r_coils"]
    circ = plt.Circle((0, 0), r_circ, color=v["color"], fill=False, ls='--', lw=1.5, label=f'R coils={r_circ}mm')
    ax.add_patch(circ)
    
    rot_c = plt.Circle((0, 0), 22, color='#64748b', fill=False, lw=1.0, ls=':')
    ax.add_patch(rot_c)
    
    ax.set_xlim(-60, 60)
    ax.set_ylim(-60, 60)
    ax.set_xlabel('Coord X [mm]', color='white', fontsize=8)
    ax.set_ylabel('Coord Y [mm]', color='white', fontsize=8)
    ax.tick_params(colors='white', labelsize=7)
    ax.grid(True, ls=':', color='#334155', alpha=0.6)
    
    cb = plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('|B| nel Traferro [mT]', color='white', fontsize=8)
    cb.ax.tick_params(labelsize=7, colors='white')

def render_polarization_odograph(ax, v):
    """Disegna l'odografo di polarizzazione trasversa B_perp(t) con metriche Stokes."""
    ax.set_facecolor('#0f172a')
    
    t = np.linspace(0, 1.0, 500)
    omega_t = 2.0 * np.pi * t
    b_pk = v["b_max_gap"] * 1.414
    
    if v["type"] == "single_rotor":
        # Lineare planare collassato
        b_th = b_pk * np.cos(omega_t)
        b_ph = 0.12 * b_pk * np.cos(omega_t)
    elif v["type"] == "triskelion":
        # Modo a trifoglio m=3
        bp = 0.82 * b_pk
        bm = 0.38 * b_pk
        b_th = bp * np.cos(omega_t) + bm * np.cos(omega_t) + 0.22 * b_pk * np.cos(3.0 * omega_t)
        b_ph = bp * np.sin(omega_t) - bm * np.sin(omega_t) - 0.22 * b_pk * np.sin(3.0 * omega_t)
    elif v["type"] == "fibonacci":
        # Ellisse modulata a 32°
        bp = 0.93 * b_pk
        bm = 0.28 * b_pk
        rip = 0.04 * np.cos(24.0 * omega_t)
        b1 = (bp + bm) * np.cos(omega_t) * (1.0 + rip)
        b2 = (bp - bm) * np.sin(omega_t) * (1.0 + rip)
        ang = np.radians(32.0)
        b_th = b1 * np.cos(ang) - b2 * np.sin(ang)
        b_ph = b1 * np.sin(ang) + b2 * np.cos(ang)
    elif v["type"] == "chiral_diode":
        # Circolare ultra-puro quasi perfetto
        bp = 0.999 * b_pk
        bm = 0.015 * b_pk
        b_th = bp * np.cos(omega_t) + bm * np.cos(omega_t)
        b_ph = bp * np.sin(omega_t) - bm * np.sin(omega_t)
    elif v["type"] == "copper_collimator":
        # Circolare puro guidato
        bp = 0.985 * b_pk
        bm = 0.080 * b_pk
        b_th = bp * np.cos(omega_t) + bm * np.cos(omega_t)
        b_ph = bp * np.sin(omega_t) - bm * np.sin(omega_t)
    else:
        # Dual 90° e WPT benchtop
        bp = 0.988 * b_pk
        bm = 0.155 * b_pk
        b_th = bp * np.cos(omega_t) + bm * np.cos(omega_t)
        b_ph = bp * np.sin(omega_t) - bm * np.sin(omega_t)
        
    ax.plot(b_th, b_ph, color=v["color"], lw=2.5, label='Odografo B_perp(t)')
    
    # Frecce di circolazione
    for idx in [120, 250, 380]:
        ax.annotate('', xy=(b_th[idx+3], b_ph[idx+3]), xytext=(b_th[idx], b_ph[idx]),
                    arrowprops=dict(arrowstyle="->", color=v["color"], lw=2.0, mutation_scale=14))
        
    # Cerchio di riferimento ideale
    circ_ideal = plt.Circle((0, 0), b_pk, color='#94a3b8', ls='--', lw=0.9, fill=False, alpha=0.4)
    ax.add_patch(circ_ideal)
    
    lim = b_pk * 1.25
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color='#334155', ls=':', lw=0.7)
    ax.axvline(0, color='#334155', ls=':', lw=0.7)
    ax.set_xlabel('B_theta Trasverso [mT]', color='white', fontsize=8)
    ax.set_ylabel('B_phi Trasverso [mT]', color='white', fontsize=8)
    ax.tick_params(colors='white', labelsize=7)
    ax.grid(True, ls=':', color='#334155', alpha=0.6)
    ax.set_aspect('equal')
    
    info_box = (
        f"Stokes s3 = {v['s3']:+.3f}\n"
        f"Purezza CP: {v['purity_cp']:.1f}%\n"
        f"Axial Ratio: {v['ar_db']:.2f} dB\n"
        f"Stato: {v['status_cp']}"
    )
    ax.text(0.04, 0.96, info_box, transform=ax.transAxes, color='white', fontsize=7.5,
            va='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='#1e293b', edgecolor=v['color'], alpha=0.9))

def generate_individual_variant_figures():
    """Genera 1 tavola per CIASCUNA variante all'interno della rispettiva cartella figures."""
    print("\n--- Generazione Tavole Individuali per Ciascuna Variante (300 DPI) ---")
    
    for v in VARIANTS_CONFIG:
        v_folder = VAR_DIR / v["id"]
        v_fig_dir = v_folder / "figures"
        v_fig_dir.mkdir(parents=True, exist_ok=True)
        out_fig = v_fig_dir / "fig_machine_field_polarization.png"
        
        fig = plt.figure(figsize=(20, 6.8), facecolor='#070b12')
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.1, 1.0, 1.0], wspace=0.28,
                               left=0.04, right=0.96, top=0.78, bottom=0.10)
        
        # Titolo super
        fig.suptitle(f"OPEN CHIRAL FLUX SHAPER | {v['name'].upper()}\n"
                     f"Geometria Macchina 3D • Distribuzione di Campo Magnetico • Odografo di Polarizzazione Trasversa",
                     color='#f8fafc', fontsize=12, fontweight='bold', y=0.93)
        
        # Panel A: Macchina 3D
        ax_a = fig.add_subplot(gs[0], projection='3d')
        render_3d_machine(ax_a, v)
        ax_a.set_title(f"A: Geometria Macchina & Bobine (R={v['r_coils']} mm)\nMantello: {v['mantle']}",
                       color='#38bdf8', fontsize=10, fontweight='bold', pad=8)
        
        # Panel B: Campo Magnetico 2D/3D
        ax_b = fig.add_subplot(gs[1])
        render_b_field_map(ax_b, v)
        ax_b.set_title(f"B: Distribuzione Induzione |B| e Flusso\nB_max traferro = {v['b_max_gap']:.2f} mT",
                       color='#ec4899', fontsize=10, fontweight='bold', pad=8)
        
        # Panel C: Polarizzazione Odografo
        ax_c = fig.add_subplot(gs[2])
        render_polarization_odograph(ax_c, v)
        ax_c.set_title(f"C: Odografo di Polarizzazione Trasversa B_perp(t)\n{v['status_cp']}",
                       color=v['color'], fontsize=10, fontweight='bold', pad=8)
        
        fig.savefig(out_fig, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
        plt.close(fig)
        print(f"  [OK] Salvata: {out_fig.relative_to(ROOT_DIR)} ({out_fig.stat().st_size / 1e6:.2f} MB)")

def generate_master_synoptic_plate():
    """Genera la tavola sinottica complessiva (Figura 43, 300 DPI) con tutte le 8 varianti a confronto."""
    print("\n--- Generazione Tavola Sinottica Master (Figura 43, 300 DPI) ---")
    out_master = FIG_DIR / "fig_43_all_variants_machine_field_polarization_matrix.png"
    
    n_vars = len(VARIANTS_CONFIG) # 8 varianti
    fig = plt.figure(figsize=(22, 28), facecolor='#070b12')
    gs = gridspec.GridSpec(n_vars, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.22, hspace=0.32,
                           left=0.05, right=0.96, top=0.95, bottom=0.03)
    
    fig.suptitle("OPEN CHIRAL FLUX SHAPER | MATRICE SINOTTICA MULTI-VARIANTE (8 ARCHITETTURE A CONFRONTO)\n"
                 "Architettura 3D della Macchina • Topologia del Campo Magnetico • Odografi di Polarizzazione Trasversa ed Elicità",
                 color='#f8fafc', fontsize=15, fontweight='bold', y=0.98)
    
    for row, v in enumerate(VARIANTS_CONFIG):
        # Col 0: Macchina 3D
        ax_m = fig.add_subplot(gs[row, 0], projection='3d')
        render_3d_machine(ax_m, v)
        ax_m.set_title(f"[{row+1}] {v['name']}\nR={v['r_coils']}mm, {v['mantle']}",
                       color=v['color'], fontsize=9, fontweight='bold', pad=4)
        
        # Col 1: Campo Magnetico
        ax_f = fig.add_subplot(gs[row, 1])
        render_b_field_map(ax_f, v)
        ax_f.set_title(f"Campo |B| Traferro: {v['b_max_gap']:.2f} mT",
                       color='white', fontsize=8.5, fontweight='bold', pad=4)
        
        # Col 2: Odografo Polarizzazione
        ax_p = fig.add_subplot(gs[row, 2])
        render_polarization_odograph(ax_p, v)
        ax_p.set_title(f"Polarizzazione: {v['purity_cp']:.1f}% CP (s3={v['s3']:+.2f}, AR={v['ar_db']:.1f}dB)",
                       color=v['color'], fontsize=8.5, fontweight='bold', pad=4)
        
    fig.savefig(out_master, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close(fig)
    print(f"  [OK] Tavola Sinottica Master esportata in: {out_master} ({out_master.stat().st_size / 1e6:.2f} MB)")
    
    # Copia nella cartella artefatti
    artifact_path = Path("C:/Users/bresc/.gemini/antigravity/brain/359566b7-3516-4512-84bb-2527bf014206") / "fig_43_all_variants_machine_field_polarization_matrix.png"
    if artifact_path.parent.exists():
        import shutil
        shutil.copy(out_master, artifact_path)
        print(f"  [OK] Copiata in artifact: {artifact_path}")

def main():
    print("=" * 90)
    print("=== AVVIO GENERAZIONE FIGURE MACCHINA, CAMPO E POLARIZZAZIONE (8 VARIANTI) ===")
    print("Framework: Open Chiral Flux Shaper | Licenza: CERN-OHL-S-2.0")
    print("=" * 90)
    
    generate_individual_variant_figures()
    generate_master_synoptic_plate()
    
    print("\n" + "=" * 90)
    print("=== GENERAZIONE COMPLETATA CON SUCCESSO SU TUTTE LE VARIANTI ===")
    print("=" * 90)

if __name__ == "__main__":
    main()
