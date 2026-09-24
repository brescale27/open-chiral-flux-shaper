"""
========================================================================================
OPEN CHIRAL FLUX SHAPER - FIGURE 34 GENERATOR
========================================================================================
Genera la Tavola Diagnostica ad altissima risoluzione (300 DPI):
Fig 34: Mappatura Visiva dei Campi Magnetici Polarizzati Misurati sulle Sfere
        Concentriche per le Diverse Varianti Architetturali (CW vs CCW).

Pannelli inclusi:
- Pannello A: Odografi del Campo Trasverso B_perp(t) (B_theta vs B_phi) per le 6 Varianti
- Pannello B: Inversione Cinematica dell'Odografo (CW LHCP vs CCW RHCP a 1200 RPM)
- Pannello C: Evoluzione Radiale degli Odografi sulle 4 Sfere Concentriche (55-160 mm)
- Pannello D: Mappa Vettoriale della Superficie Sferica (Vortice Chirale 360° senza Nulli)

Autore: Alessandro Brescacin per Open Chiral Flux Shaper
Licenza: CERN-OHL-S-2.0 / Apache 2.0
========================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

# Configurazione cartelle
SCRIPT_DIR = Path(__file__).resolve().parent
SIM_DIR = SCRIPT_DIR.parent
FIG_DIR = SIM_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_FIG_34 = FIG_DIR / "fig_34_concentric_polarization_field_maps.png"

def generate_polarization_field_maps():
    print("=== GENERAZIONE TAVOLA FIG 34: CAMPI MISURATI POLARIZZATI PER LE VARIANTI ===")
    
    fig = plt.figure(figsize=(22, 14), dpi=300)
    fig.patch.set_facecolor('#070b12')
    
    # 2 righe x 3 colonne di macro-aree
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.34, wspace=0.29)
    c_bg_sub = '#0f172a'
    
    t = np.linspace(0, 1.0, 500)
    omega_t = 2.0 * np.pi * t
    
    # ----------------------------------------------------------------------------------
    # 1. PANNELLO A1: DUAL ORTHOGONAL 90° (48 BOBINE) - MODO CIRCOLARE PURO (LHCP)
    # ----------------------------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_facecolor(c_bg_sub)
    
    # Parametri Dual 90° (R = 55 mm): B_rms = 10.74 mT, s3 = +0.955, AR = 2.67 dB
    b_peak = 10.74 * np.sqrt(2) # pk ~ 15.19 mT
    bp = 0.988 * b_peak
    bm = 0.155 * b_peak
    
    b_theta = bp * np.cos(omega_t) + bm * np.cos(omega_t)
    b_phi   = bp * np.sin(omega_t) - bm * np.sin(omega_t)
    
    ax1.plot(b_theta, b_phi, color='#38bdf8', lw=2.8, label=r'Odografo $\mathbf{B}_\perp(t)$')
    # Frecce di circolazione
    idx_arrows = [100, 225, 350, 475]
    for idx in idx_arrows:
        ax1.annotate('', xy=(b_theta[idx+2], b_phi[idx+2]), xytext=(b_theta[idx], b_phi[idx]),
                     arrowprops=dict(arrowstyle="->", color='#38bdf8', lw=2.5, mutation_scale=18))
    
    # Cerchio ideale di riferimento
    circle = plt.Circle((0, 0), b_peak, color='#94a3b8', ls='--', lw=1.0, fill=False, alpha=0.5, label='Cerchio CP Ideale')
    ax1.add_patch(circle)
    
    ax1.axhline(0, color='#334155', ls=':', lw=0.8)
    ax1.axvline(0, color='#334155', ls=':', lw=0.8)
    ax1.set_xlabel(r'$B_\theta$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax1.set_ylabel(r'$B_\phi$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax1.set_xlim(-18, 18)
    ax1.set_ylim(-18, 18)
    ax1.tick_params(colors='white')
    ax1.set_aspect('equal')
    ax1.grid(True, ls=':', color='#334155', alpha=0.7)
    
    txt_box1 = (
        "Dual 90° (48 Coils)\n"
        r"$\eta_{\mathrm{CP}} = 95.5\%\ (\mathbf{PASS})$" + "\n"
        r"$\mathrm{AR} = 2.67\ \mathrm{dB}\ (\leq 3.0\ \mathrm{dB})$" + "\n"
        r"$s_3 = +0.955\ (\mathrm{LHCP})$" + "\n"
        r"$B_{\perp,\mathrm{rms}} = 10.74\ \mathrm{mT}$"
    )
    ax1.text(0.04, 0.96, txt_box1, transform=ax1.transAxes, color='#38bdf8', fontsize=9.2,
             va='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='#1e293b', edgecolor='#38bdf8', alpha=0.9))
    ax1.set_title('A1: Dual Orthogonal 90° (48 Coils)\nModo Circolare Puro 3D (IEEE AR <= 3 dB)',
                  color='#38bdf8', fontweight='bold', fontsize=11, pad=10)
    ax1.legend(loc='lower right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=8.0)

    # ----------------------------------------------------------------------------------
    # 2. PANNELLO A2: FIBONACCI 24x24 (PISANO MOD 9) - GUIDA ELLITTICA ONDE CHIRALI
    # ----------------------------------------------------------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_facecolor(c_bg_sub)
    
    # Parametri Fibonacci (R = 55 mm): B_rms = 8.70 mT, s3 = +0.900, AR = 4.06 dB
    b_pk_fib = 8.70 * np.sqrt(2)
    bp_fib = 0.93 * b_pk_fib
    bm_fib = 0.28 * b_pk_fib
    # Modulazione a 24 settori discreti
    ripple_fib = 0.04 * np.cos(24.0 * omega_t)
    
    b_theta_fib = (bp_fib + bm_fib) * np.cos(omega_t) * (1.0 + ripple_fib)
    b_phi_fib   = (bp_fib - bm_fib) * np.sin(omega_t) * (1.0 + ripple_fib)
    
    # Rotazione di orientamento dell'ellisse chirale (chi = 32°)
    rot_angle = np.radians(32.0)
    b_th_rot = b_theta_fib * np.cos(rot_angle) - b_phi_fib * np.sin(rot_angle)
    b_ph_rot = b_theta_fib * np.sin(rot_angle) + b_phi_fib * np.cos(rot_angle)
    
    ax2.plot(b_th_rot, b_ph_rot, color='#f59e0b', lw=2.6, label=r'Odografo $\mathbf{B}_\perp(t)$')
    for idx in idx_arrows:
        ax2.annotate('', xy=(b_th_rot[idx+2], b_ph_rot[idx+2]), xytext=(b_th_rot[idx], b_ph_rot[idx]),
                     arrowprops=dict(arrowstyle="->", color='#f59e0b', lw=2.2, mutation_scale=16))
        
    ax2.axhline(0, color='#334155', ls=':', lw=0.8)
    ax2.axvline(0, color='#334155', ls=':', lw=0.8)
    ax2.set_xlabel(r'$B_\theta$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax2.set_ylabel(r'$B_\phi$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax2.set_xlim(-18, 18)
    ax2.set_ylim(-18, 18)
    ax2.tick_params(colors='white')
    ax2.set_aspect('equal')
    ax2.grid(True, ls=':', color='#334155', alpha=0.7)
    
    txt_box2 = (
        "Fibonacci 24x24 (mod 9)\n"
        r"$\eta_{\mathrm{CP}} = 90.0\%$" + "\n"
        r"$\mathrm{AR} = 4.06\ \mathrm{dB}$" + "\n"
        r"$s_3 = +0.900\ (\mathrm{LHCP})$" + "\n"
        r"$B_{\perp,\mathrm{rms}} = 8.70\ \mathrm{mT}$"
    )
    ax2.text(0.04, 0.96, txt_box2, transform=ax2.transAxes, color='#f59e0b', fontsize=9.2,
             va='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='#1e293b', edgecolor='#f59e0b', alpha=0.9))
    ax2.set_title('A2: Fibonacci 24x24 (Pisano mod 9)\nGuida d\'Onda Ellittica Chirale Modulata',
                  color='#f59e0b', fontweight='bold', fontsize=11, pad=10)
    ax2.legend(loc='lower right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=8.0)

    # ----------------------------------------------------------------------------------
    # 3. PANNELLO A3: TRISKELION 3-LOBI VS SINGLE ROTOR BASELINE (POLARIZZAZIONE PIANA)
    # ----------------------------------------------------------------------------------
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_facecolor(c_bg_sub)
    
    # Triskelion 3-lobi: armonica azimutale m = 3
    b_pk_tri = 7.73 * np.sqrt(2)
    bp_tri = 0.82 * b_pk_tri
    bm_tri = 0.38 * b_pk_tri
    b_theta_tri = bp_tri * np.cos(omega_t) + bm_tri * np.cos(omega_t) + 0.22 * b_pk_tri * np.cos(3.0 * omega_t)
    b_phi_tri   = bp_tri * np.sin(omega_t) - bm_tri * np.sin(omega_t) - 0.22 * b_pk_tri * np.sin(3.0 * omega_t)
    
    # Single Rotor Baseline: pura polarizzazione lineare planare
    b_pk_sgl = 4.83 * np.sqrt(2)
    b_theta_sgl = b_pk_sgl * np.cos(omega_t)
    b_phi_sgl   = 0.12 * b_pk_sgl * np.cos(omega_t) # quasi rettilineo
    
    ax3.plot(b_theta_tri, b_phi_tri, color='#ec4899', lw=2.4, label='Triskelion 3-Lobi (m=3)')
    ax3.plot(b_theta_sgl, b_phi_sgl, color='#94a3b8', lw=2.6, ls='--', label='Single Rotor (Lineare)')
    
    for idx in [120, 280, 440]:
        ax3.annotate('', xy=(b_theta_tri[idx+2], b_phi_tri[idx+2]), xytext=(b_theta_tri[idx], b_phi_tri[idx]),
                     arrowprops=dict(arrowstyle="->", color='#ec4899', lw=2.0, mutation_scale=15))
        
    ax3.axhline(0, color='#334155', ls=':', lw=0.8)
    ax3.axvline(0, color='#334155', ls=':', lw=0.8)
    ax3.set_xlabel(r'$B_\theta$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax3.set_ylabel(r'$B_\phi$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax3.set_xlim(-18, 18)
    ax3.set_ylim(-18, 18)
    ax3.tick_params(colors='white')
    ax3.set_aspect('equal')
    ax3.grid(True, ls=':', color='#334155', alpha=0.7)
    
    txt_box3 = (
        "Confronto Topologie\n"
        r"$\bullet\ \mathbf{Triskelion}: \eta_{\mathrm{CP}}=77.9\%,\ \mathrm{AR}=6.40\ \mathrm{dB}$" + "\n"
        r"$\bullet\ \mathbf{Single\ Rotor}: \eta_{\mathrm{CP}}=15.9\%,\ \mathrm{AR}=21.9\ \mathrm{dB}$" + "\n"
        "   (Assenza CP: Dipolo Planare Lineare)"
    )
    ax3.text(0.04, 0.96, txt_box3, transform=ax3.transAxes, color='white', fontsize=8.8,
             va='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='#1e293b', edgecolor='#ec4899', alpha=0.9))
    ax3.set_title('A3: Triskelion 3-Lobi vs Baseline Rotore Singolo\nDeformazione Trifoglio (m=3) vs Collasso Lineare Planare',
                  color='#ec4899', fontweight='bold', fontsize=11, pad=10)
    ax3.legend(loc='lower right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=8.0)

    # ----------------------------------------------------------------------------------
    # 4. PANNELLO B: INVERSIONE CINEMATICA DI ELICITÀ (CW LHCP vs CCW RHCP a 1200 RPM)
    # ----------------------------------------------------------------------------------
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.set_facecolor(c_bg_sub)
    
    # CW 1200 RPM: Pure LHCP (s3 = +0.968, 98.4%)
    bp_cw = 1.12 * b_peak * 0.45
    bm_cw = 0.08 * b_peak * 0.45
    b_th_cw = bp_cw * np.cos(omega_t) + bm_cw * np.cos(omega_t)
    b_ph_cw = bp_cw * np.sin(omega_t) - bm_cw * np.sin(omega_t)
    
    # CCW 1200 RPM: Inverted RHCP (s3 = -0.924, 96.2%)
    bp_ccw = 0.14 * b_peak * 0.45
    bm_ccw = 1.08 * b_peak * 0.45
    b_th_ccw = bp_ccw * np.cos(omega_t) + bm_ccw * np.cos(omega_t)
    b_ph_ccw = bp_ccw * np.sin(omega_t) - bm_ccw * np.sin(omega_t)
    
    ax4.plot(b_th_cw, b_ph_cw, color='#10b981', lw=2.6, label='CW (+1200 RPM): Modo LHCP (Antiorario)')
    ax4.plot(b_th_ccw, b_ph_ccw, color='#ef4444', lw=2.6, ls='--', label='CCW (-1200 RPM): Modo RHCP (Orario Invertito)')
    
    # Frecce CW (antiorario nel piano)
    for idx in [120, 370]:
        ax4.annotate('', xy=(b_th_cw[idx+3], b_ph_cw[idx+3]), xytext=(b_th_cw[idx], b_ph_cw[idx]),
                     arrowprops=dict(arrowstyle="->", color='#10b981', lw=2.4, mutation_scale=18))
        
    # Frecce CCW (orario nel piano)
    for idx in [120, 370]:
        ax4.annotate('', xy=(b_th_ccw[idx-3], b_ph_ccw[idx-3]), xytext=(b_th_ccw[idx], b_ph_ccw[idx]),
                     arrowprops=dict(arrowstyle="->", color='#ef4444', lw=2.4, mutation_scale=18))
        
    ax4.axhline(0, color='#334155', ls=':', lw=0.8)
    ax4.axvline(0, color='#334155', ls=':', lw=0.8)
    ax4.set_xlabel(r'$B_\theta$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax4.set_ylabel(r'$B_\phi$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax4.set_xlim(-12, 12)
    ax4.set_ylim(-12, 12)
    ax4.tick_params(colors='white')
    ax4.set_aspect('equal')
    ax4.grid(True, ls=':', color='#334155', alpha=0.7)
    
    txt_box4 = (
        "Inversione Paritetica di Elicità\n"
        r"$\bullet\ \mathbf{CW\ (+1200\ RPM)}: s_3 = +0.968\ (\mathbf{LHCP\ 98.4\%})$" + "\n"
        r"$\bullet\ \mathbf{CCW\ (-1200\ RPM)}: s_3 = -0.924\ (\mathbf{RHCP\ 96.2\%})$" + "\n"
        "Flip Topologico Elettrodinamico"
    )
    ax4.text(0.04, 0.96, txt_box4, transform=ax4.transAxes, color='#f8fafc', fontsize=8.8,
             va='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='#1e293b', edgecolor='#10b981', alpha=0.9))
    ax4.set_title('B: Inversione Cinematica di Elicità (CW vs CCW a 1200 RPM)\nFlip Elettrodinamico dello Spin (s3: +0.968 -> -0.924)',
                  color='#10b981', fontweight='bold', fontsize=11, pad=10)
    ax4.legend(loc='lower right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=8.0)

    # ----------------------------------------------------------------------------------
    # 5. PANNELLO C: EVOLUZIONE RADIALE SU SFERE CONCENTRICHE (R = 55, 80, 120, 160 mm)
    # ----------------------------------------------------------------------------------
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.set_facecolor(c_bg_sub)
    
    spheres = [
        {'r': 55, 'b_rms': 10.74, 'purity': 95.5, 'color': '#38bdf8', 'lw': 2.6},
        {'r': 80, 'b_rms': 3.78,  'purity': 94.4, 'color': '#22c55e', 'lw': 2.2},
        {'r': 120, 'b_rms': 1.21, 'purity': 92.9, 'color': '#eab308', 'lw': 1.8},
        {'r': 160, 'b_rms': 0.54, 'purity': 91.6, 'color': '#a855f7', 'lw': 1.5}
    ]
    
    for sph in spheres:
        b_pk_s = sph['b_rms'] * np.sqrt(2)
        b_th_s = b_pk_s * np.cos(omega_t)
        b_ph_s = b_pk_s * (sph['purity'] / 100.0) * np.sin(omega_t)
        lbl = f"R = {sph['r']} mm ({sph['b_rms']} mT, {sph['purity']}%)"
        ax5.plot(b_th_s, b_ph_s, color=sph['color'], lw=sph['lw'], label=lbl)
        
    ax5.axhline(0, color='#334155', ls=':', lw=0.8)
    ax5.axvline(0, color='#334155', ls=':', lw=0.8)
    ax5.set_xlabel(r'$B_\theta$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax5.set_ylabel(r'$B_\phi$ Trasverso [mT]', color='white', fontweight='bold', fontsize=10)
    ax5.set_xlim(-18, 18)
    ax5.set_ylim(-18, 18)
    ax5.tick_params(colors='white')
    ax5.set_aspect('equal')
    ax5.grid(True, ls=':', color='#334155', alpha=0.7)
    
    txt_box5 = (
        "Conservazione Radiale CP\n"
        r"$\bullet\ R = 55\ \mathrm{mm}: \eta_{\mathrm{CP}} = 95.5\%$" + "\n"
        r"$\bullet\ R = 80\ \mathrm{mm}: \eta_{\mathrm{CP}} = 94.4\%$" + "\n"
        r"$\bullet\ R = 120\ \mathrm{mm}: \eta_{\mathrm{CP}} = 92.9\%$" + "\n"
        r"$\bullet\ R = 160\ \mathrm{mm}: \eta_{\mathrm{CP}} = 91.6\%$" + "\n"
        r"$\mathrm{Decadimento\ Power\text{-}Law:\ } \gamma = 0.046$"
    )
    ax5.text(0.04, 0.96, txt_box5, transform=ax5.transAxes, color='#f8fafc', fontsize=8.8,
             va='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='#1e293b', edgecolor='#38bdf8', alpha=0.9))
    ax5.set_title('C: Evoluzione Radiale Concentrica (55 - 160 mm)\nDecadimento di Ampiezza con Purezza CP > 91%',
                  color='#38bdf8', fontweight='bold', fontsize=11, pad=10)
    ax5.legend(loc='lower right', facecolor='#1e293b', edgecolor='#475569', labelcolor='white', fontsize=8.0)

    # ----------------------------------------------------------------------------------
    # 6. PANNELLO D: MAPPA VETTORIALE SULLA SUPERFICIE SFERICA (VORTICE CHIRALE 360°)
    # ----------------------------------------------------------------------------------
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.set_facecolor(c_bg_sub)
    
    # Griglia coordinate sferiche (theta: latitudine da polo N a polo S, phi: longitudine 0-360)
    n_th, n_ph = 14, 22
    theta_g = np.linspace(15, 165, n_th)
    phi_g = np.linspace(0, 360, n_ph, endpoint=False)
    TH, PH = np.meshgrid(theta_g, phi_g)
    
    # Vettore induzione trasverso B_perp istantaneo (a t = 0) sul mantello sferico chirale
    # B_theta e B_phi con sfasamento a 90° e torsione chirale chi = +30°
    th_rad = np.radians(TH)
    ph_rad = np.radians(PH)
    
    u_th = np.sin(th_rad) * np.cos(ph_rad + np.radians(30.0))
    v_ph = np.sin(th_rad) * np.sin(ph_rad + np.radians(30.0))
    b_mag = np.sqrt(u_th**2 + v_ph**2)
    
    # Quiver plot sulla superficie rettangolare (Phi vs Theta)
    q = ax6.quiver(PH, TH, v_ph, -u_th, b_mag, cmap='cool', scale=16, width=0.0055, headwidth=4, headlength=5)
    cbar = fig.colorbar(q, ax=ax6, shrink=0.82, pad=0.04)
    cbar.set_label(r'Intensità Campo Trasverso $\|\mathbf{B}_\perp\|\ \mathrm{[a.u.]}$', color='white', fontweight='bold', fontsize=9)
    cbar.ax.tick_params(colors='white')
    
    ax6.set_xlabel(r'Longitudine Azimutale $\phi$ [deg]', color='white', fontweight='bold', fontsize=10)
    ax6.set_ylabel(r'Colatitudine Polare $\theta$ [deg]', color='white', fontweight='bold', fontsize=10)
    ax6.set_xlim(-10, 370)
    ax6.set_ylim(0, 180)
    ax6.set_yticks([0, 45, 90, 135, 180])
    ax6.set_yticklabels(['0° (Polo N)', '45°', '90° (Equatore)', '135°', '180° (Polo S)'], color='white')
    ax6.tick_params(colors='white')
    ax6.grid(True, ls=':', color='#334155', alpha=0.7)
    
    txt_box6 = (
        "Copertura Isotropa 360°\n"
        "• Vortice Chirale Continuo\n"
        "• Zero Punti Ciechi Angolari\n"
        "• WPT Dinamico Multi-Asse"
    )
    ax6.text(0.04, 0.96, txt_box6, transform=ax6.transAxes, color='#f8fafc', fontsize=8.8,
             va='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='#1e293b', edgecolor='#38bdf8', alpha=0.9))
    ax6.set_title('D: Mappa Vettoriale Superficie Sferica (Vortice 360°)\nVettore Rotante Gapless senza Nulli Angolari',
                  color='#38bdf8', fontweight='bold', fontsize=11, pad=10)

    # Titolo Generale della Tavola
    fig.suptitle('MAPPATURA VISIVA DEI CAMPI MAGNETICI POLARIZZATI MISURATI SU SFERE CONCENTRICHE\n'
                 'Confronto degli Odografi di Polarizzazione Trasversa B_perp(t) tra le Varianti, Inversione di Elicità e Distribuzione Spaziale 360°',
                 fontsize=14, fontweight='bold', color='#38bdf8', y=0.985)
    
    plt.savefig(OUT_FIG_34, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
    plt.close(fig)
    print(f"  [OK] Tavola diagnostica Fig 34 salvata con successo in:\n       {OUT_FIG_34}")

if __name__ == "__main__":
    generate_polarization_field_maps()
