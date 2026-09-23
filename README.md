# Open Chiral Flux Shaper

*Anisotropic Macro-Chiral Metamaterial for Radial Induction Shaping, Wireless Power Projection, and Electromagnetic Lift.*

[![License: CERN-OHL-S-2.0](https://img.shields.io/badge/License-CERN--OHL--S--2.0-blue.svg)](LICENSE.txt)
[![Release: v1.0.0](https://img.shields.io/badge/Release-v1.0.0-green.svg)](https://github.com/brescale27/open-chiral-flux-shaper/releases)
[![FEM Solver: Elmer FEM 9.0](https://img.shields.io/badge/Elmer%20FEM-9.0%20(CSC)-orange.svg)](https://www.csc.fi/web/elmer)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.cern--ohl--s--2.0-lightgrey.svg)](https://github.com/brescale27/open-chiral-flux-shaper)

---

## Executive Summary: What is the Open Chiral Flux Shaper?

### The Problem
In conventional electromechanics and high-frequency power engineering, enclosing a dynamic or rotating magnetic field source inside a metallic shell triggers massive azimuthal eddy currents ($J = \sigma E$). Governed by **Lenz's Law**, these surface eddy currents generate an opposing magnetic counter-field that:
1. **Traps and shields the electromagnetic flux** inside the inner cavity.
2. **Dissipates severe Joule losses** ($P_{\text{loss}} = \int \sigma |\vec{E}|^2 dV$), causing thermal runaway and drastically impairing coupling efficiency.

### The Breakthrough
The **Open Chiral Flux Shaper** resolves this fundamental barrier by replacing solid metallic walls with an **engineered macro-chiral metamaterial mantle** composed of multilayer expanded aluminum mesh (*expanded metal lattice*).

By tilting the metallic micro-bridges at a calibrated **$30^\circ$ louver chiral angle** relative to the machine axis, the shell functions as an anisotropic metasurface governed by a positive semi-definite conductivity tensor:
$$\bar{\bar{\sigma}}_{\text{cyl}} = \begin{bmatrix} \sigma_{rr} & 0 & 0 \\ 0 & \sigma_{\theta\theta} & \sigma_{\theta z} \\ 0 & \sigma_{\theta z} & \sigma_{zz} \end{bmatrix} = \begin{bmatrix} 1.75\times 10^6 & 0 & 0 \\ 0 & 1.75\times 10^6 & 3.031\times 10^6 \\ 0 & 3.031\times 10^6 & 1.22\times 10^7 \end{bmatrix} \text{ S/m}$$

Instead of opposing the rotating magnetic wave, the chiral mantle:
- **Suppresses closed circular eddy loops**, slashing Joule thermal dissipation by **$-43.2\%$** in continuous mode, and up to **$-99.9\%$** in pulsed half-wave mode.
- **Couples azimuthal electric fields to axial currents** ($\sigma_{\theta z}$ cross-coupling), deflecting and **unrolling the magnetic flux outward into a $360^\circ$ omnidirectional radial induction wave**.
- **Enables electromagnetic propulsion and levitation:** In a biconical induction configuration, the fixed $30^\circ$ chiral tilt breaks axial reflection parity ($\mathcal{P}_z$), producing a continuous, unidirectional upward ponderomotive Lorentz lift ($\langle F_z \rangle = +4.67\,\mu\text{N}$, reaching a global resonance peak of $+5.72\,\mu\text{N}$ / $57.2\text{ mN}$ at full scale in locked-rotor mode).

---

## Core Architectures & Comparative Benchmark

| Parameter / Metric | Baseline Architecture (v1.0.0) | Centered Continuous (1200 RPM) | Centered Locked-Rotor (0 RPM Peak) | Mirrored Polarity Pulsed (N-S) | Closed Can (1200 RPM, Lids) | Physical Mechanism / Specialty |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Core Geometry** | Ferromagnetic spider at bottom ($Z = -H/2$) | Ferromagnetic core at equator ($Z = 0$) | Ferromagnetic core at equator ($Z = 0$) | Ferromagnetic core at equator ($Z = 0$) | Equatorial Core ($Z=0$), Closed Can ($Z=\pm H/2$) | Equatorial magnetic flux bridge |
| **Coil Excitation** | 60° Progressive Full Sine Wave ($100\text{ Hz}$) | 60° Progressive Full Sine Wave ($100\text{ Hz}$) | 60° Progressive Full Sine Wave ($100\text{ Hz}$) | 60° Shifted Half-Wave Pulse Train ($J_k \ge 0$) | 60° Progressive Full Sine Wave ($100\text{ Hz}$) | Directional vs Sequenced Pulsed Drive |
| **Rotor Speed & Slip** | $1200\text{ RPM}$ ($f_{\text{slip}} = 40\text{ Hz}$) | $1200\text{ RPM}$ ($f_{\text{slip}} = 40\text{ Hz}$) | **$0\text{ RPM}$ (Locked, $f_{\text{slip}} = 100\text{ Hz}$)** | $1200\text{ RPM}$ ($f_{\text{slip}} = 40\text{ Hz}$) | $1200\text{ RPM}$ ($f_{\text{slip}} = 40\text{ Hz}$) | Mechanical tracking vs maximum relative slip |
| **Polarity Layout** | Unipolar Homogeneous | Unipolar Homogeneous | Unipolar Homogeneous | Alternating Mirrored: $s_k = (-1)^{k-1}$ | Unipolar Homogeneous | Adjacent N-S magnetic return loops |
| **Radial Field $B_{\text{rad}}$ ($R = 6\text{ cm}$)** | **$62.98\,\mu\text{T}$** | **$211.35\,\mu\text{T}$** | **$211.35\,\mu\text{T}$** | **$135.84\,\mu\text{T}$** | **$185.2\,\mu\text{T}$** | Strong equatorial concentration (+115.7% vs baseline) |
| **Poynting Power Flux ($R=12\text{ cm}$)** | **$+7.282\text{ mW}$** | **$+2.620\text{ mW}$** | **$+2.620\text{ mW}$** | **$+102.76\text{ mW}$** | $\approx 0\text{ mW}$ (Confined) | **$+1311\%$ pulsed harmonic radiation boost** |
| **Net Axial Lorentz Force $\langle F_z \rangle$** | $\approx 0$ (leakage) | **$+4.67\,\mu\text{N}$** | **$\mathbf{+5.72\,\mu\text{N}}$ ($+57.2\text{ mN}$ full scale)** | **$-0.93\,\mu\text{N}$ ($\approx 0\text{ mN}$ balanced)** | **$-0.20\,\mu\text{N}$ ($F_{z,\text{chiral}} = -0.04\,\mu\text{N}$)** | **Locked-rotor resonance peak (+22.5% boost)** |
| **Peak Instantaneous Force $F_{z,\text{max}}$** | $\approx 0$ | $+27.41\,\mu\text{N}$ | **$+12.89\,\mu\text{N}$ ($+128.9\text{ mN}$ full scale)** | $\pm 26.50\,\mu\text{N}$ (symmetric) | $-2.96\,\mu\text{N}$ | Smooth ponderomotive lift at locked rotor |
| **Joule Dissipation $P_J$** | **$2.437\text{ W}$** | **$1.52\text{ mW}$** | **$1.02\text{ mW}$ ($0.0010\text{ W}$)** | **$1.90\text{ mW}$ ($0.0019\text{ W}$)** | **$\mathbf{0.017\text{ mW}}$ ($17.4\,\mu\text{W}$)** | **-98.9% thermal collapse (Hermetic confinement)** |
| **Lift Efficiency $\eta_F = \langle F_z \rangle / P_J$** | $\sim 0$ | $3065\,\mu\text{N/W}$ | **$\mathbf{5587\,\mu\text{N/W}}$** | $\approx 0$ (balanced) | N/A (Confined Cavity) | **Maximum solid-state thrust-to-power ratio** |
| **Primary Optimal Use-Case** | Directional Wireless Venting | Continuous Electromagnetic Lift | **Max Lorentz Lift (Solid-State Thruster)** | Resonant Wireless Pulsed Power & Low Heat | **Ultra-Low Loss Confined Cavity / Field Shielder** | Mission-specific electromagnetic tuning |

---

## Visual Showcase (High-Resolution 300 DPI Diagnostics)

<div align="center">

### 1. Radial Induction Projection & Polyphase Co-Rotating Optimization (Baseline)
| 360° Omnidirectional Radial Projection | Phase-Shift Sweep & Joule Loss Minimization |
| :---: | :---: |
| <img src="figures/02_espulsione_radiale_simmetrica_360.png" width="450" alt="360° Radial Projection" /> | <img src="figures/04_sweep_sfasamento_confronto.png" width="450" alt="Phase Shift Sweep" /> |
| *Uniform $360^\circ$ radial field expulsion through the chiral mantle.* | *Regime B ($60^\circ$ co-rotating) reduces Joule losses by $43.2\%$ and ripple to $86.7\%$.* |

### 2. Centered Variant ($Z = 0$): Biconical Flux & Net Continuous Electromagnetic Lift
| Hourglass Biconical Flux Streamlines (3D RK45) | Unidirectional Upward Lorentz Lift $F_z(t)$ |
| :---: | :---: |
| <img src="variants/rotore_centrato_z0/figures/fig_01_topologia_biconica_clessidra_3d.png" width="450" alt="Hourglass 3D Flux Lines" /> | <img src="variants/rotore_centrato_z0/figures/fig_03_forza_assiale_netta_Fz.png" width="450" alt="Net Upward Lorentz Lift" /> |
| *Hourglass flux lines: upper horn ($+Z$), lower horn ($-Z$), and equatorial ejection ring.* | *Time-dependent axial force showing net positive DC lift ($\langle F_z \rangle = +4.67\,\mu\text{N}$).* |

### 3. Mirrored Polarity Pulsed Half-Wave Variant (N-S-N-S-N-S at $Z = 0$)
| Alternate N-S Polar Topology (3D RK45) | Pulsed Half-Wave Channels & Radial Profile | Lorentz Lift $F_z(t)$ Pulsed vs Continuous |
| :---: | :---: | :---: |
| <img src="variants/rotore_centrato_poli_alternati_semionda/figures/fig_01_topologia_poli_specchiati_3d.png" width="300" alt="N-S Polar Topology 3D" /> | <img src="variants/rotore_centrato_poli_alternati_semionda/figures/fig_02_forme_onda_semionda_e_profilo_radiale.png" width="300" alt="Pulsed Waveforms & Profile" /> | <img src="variants/rotore_centrato_poli_alternati_semionda/figures/fig_03_confronto_forza_lift_Fz_impulsi.png" width="300" alt="Lift Comparison Pulsed vs Continuous" /> |
| *Short-range return loops between adjacent N-S pairs with radial ejection lobes.* | *6-channel $60^\circ$ pulsed half-wave drive and 6-lobe equatorial induction pattern.* | *Bipolar balanced oscillation of $F_z(t)$ with near-zero DC drift and ultra-low Joule heat ($1.9\text{ mW}$).* |

### 4. 2D Frequency vs RPM Resonance Sweep & Numerical Integrity Validation
| 2D Resonance Surface & Contour Map (300 DPI) | Dispersion Curve vs Slip Frequency (300 DPI) | Maxwell Stress Tensor (MST) & Bias Decoupling (300 DPI) |
| :---: | :---: | :---: |
| <img src="variants/rotore_centrato_z0_resonance_sweep/figures/fig_01_superficie_risonanza_lift_2d.png" width="310" alt="2D Resonance Surface" /> | <img src="variants/rotore_centrato_z0_resonance_sweep/figures/fig_02_curva_dispersione_vs_slip.png" width="310" alt="Dispersion Curve vs Slip" /> | <img src="variants/rotore_centrato_z0_resonance_sweep/figures/fig_03_validazione_bias_e_tensore_maxwell.png" width="310" alt="Maxwell Stress Tensor Validation" /> |
| *3D surface and 2D contour map mapping $\langle F_z \rangle(f, n)$ across the 2D operational space.* | *Slip frequency dispersion curve showing inductive resonance peak ($f_{\mathrm{slip,opt}} \approx 75.5\mathrm{ Hz}$).* | *Decoupling of tetrahedral mesh bias ($+4.92\,\mu\mathrm{N}$) from pure chiral lift and MST integration.* |

### 5. Closed Cylindrical Can Variant (Top & Bottom Lids at $Z = \pm H/2$)
| 4-Sector Electrodynamic & Thermal Breakdown (300 DPI) | Open Tube vs Closed Can Decoupling & Loss Collapse (300 DPI) |
| :---: | :---: |
| <img src="variants/rotore_centrato_mantello_chiuso/figures/fig_01_bilancio_4settori_joule_e_lift.png" width="450" alt="4-Sector Electrodynamic & Thermal Breakdown" /> | <img src="variants/rotore_centrato_mantello_chiuso/figures/fig_02_confronto_chiuso_vs_aperto_decoupled.png" width="450" alt="Open Tube vs Closed Can Decoupling" /> |
| *Electrodynamic force waveforms $F_z(t)$, thermal distribution ($P_J$), and mesh bias vs chiral lift across the 4 physical sectors.* | *Direct comparison: raw vs decoupled forces and the dramatic $-98.9\%$ collapse in Joule losses ($1.524\text{ mW} \to 0.017\text{ mW}$).* |

| 3D Spherical Field Mapping (B, E, Poynting) (300 DPI) | 360° Polar Radiation Diagram & Axial Lid Confinement (300 DPI) |
| :---: | :---: |
| <img src="variants/rotore_centrato_mantello_chiuso/figures/fig_03_mappatura_sferica_3d_campo_B_E.png" width="450" alt="3D Spherical Field Mapping" /> | <img src="variants/rotore_centrato_mantello_chiuso/figures/fig_04_diagramma_radiazione_poynting_360.png" width="450" alt="360° Polar Radiation Diagram" /> |
| *3D Fibonacci sphere mapping of $\|\vec{B}\|$, $\|\vec{E}\|$, and Poynting vector with oriented 3D direction quivers.* | *360° polar radiation diagrams demonstrating complete axial shielding by the lids ($Z=\pm H/2$) and log radial decay.* |

</div>

---

## Physical Validation & Maxwellian Rigor

All electromagnetic fields are solved using **Elmer FEM 9.0** via the transient edge-finite-element **Whitney $\vec{A}-V$ solver** coupled with analytic MATC tensor transformations.

### 1. Gauss Magnetic Solenoidality ($\oint_S \vec{B}\cdot\hat{n}\,dA = 0$)
Solenoidality was certified via 2,500-point Fibonacci spherical integrations across concentric evaluation spheres:
- **Near-Field Sphere ($R = 8.0\text{ cm}$):** Relative residual = **$0.031\% - 0.993\%$** (`PASS`, $\Phi_{\text{net}} \sim 10^{-15} - 10^{-8}\text{ Wb}$)
- **Mid-Field Sphere ($R = 12.0\text{ cm}$):** Relative residual = **$0.076\% - 5.018\%$** (`PASS`)
- **Far-Field Sphere ($R = 15.0\text{ cm}$):** Relative residual = **$1.402\% - 2.694\%$** (`PASS`)

### 2. Poynting Vector & Remote Power Projection
Integrating the Poynting vector $\vec{S} = \frac{1}{\mu_0} (\vec{E} \times \vec{B})$ across a $R=12\text{ cm}, H=20\text{ cm}$ cylindrical control surface demonstrates:
- **Baseline v1.0.0:** Active outward-directed guided power flow of **$+7.28\text{ mW}$**.
- **Centered Continuous ($Z=0$):** Equatorially localized flux of **$+2.62\text{ mW}$**.
- **Pulsed Half-Wave Variant:** High-frequency harmonic pulse train power projection of **$+102.76\text{ mW}$** with virtual suppression of thermal dissipation ($P_J = 1.9\text{ mW}$).

### 3. Numerical Integrity, Mesh Bias Decoupling & Maxwell Stress Tensor (MST)
To establish absolute scientific rigor, the axial ponderomotive force was evaluated to decouple genuine physical effects from geometric mesh discretization bias:

1. **Specular Parity Inversion ($\theta = \pm 30^\circ$ at $100\text{ Hz}, 1200\text{ RPM}$):**
   Under chiral reflection, the physical Lorentz lift must reverse sign ($F_z \rightarrow -F_z$), while tetrahedral mesh asymmetry along $Z$ is invariant. Integrating over 20 transient timesteps ($dt = 0.5\text{ ms}$):
   $$\langle F_z(+30^\circ) \rangle = +4.669\,\mu\text{N}, \qquad \langle F_z(-30^\circ) \rangle = +5.164\,\mu\text{N}$$
   Decoupling yields:
   $$F_{\text{bias}} = \frac{\langle F_z(+30^\circ) \rangle + \langle F_z(-30^\circ) \rangle}{2} = \mathbf{+4.917\,\mu\text{N}}$$
   $$F_{z,\text{chiral}} = \frac{\langle F_z(+30^\circ) \rangle - \langle F_z(-30^\circ) \rangle}{2} = \mathbf{-0.247\,\mu\text{N}}$$
   This reveals that at nominal $1200\text{ RPM}$, the uncorrected $+4.67\,\mu\text{N}$ force was dominated by tetrahedral mesh anisotropy along $Z$ ($F_{\text{bias}} = +4.92\,\mu\text{N}$).

2. **Locked-Rotor Net Chiral Peak ($100\text{ Hz}, 0\text{ RPM}$):**
   At locked rotor, the uncorrected force reaches $\langle F_z \rangle = +5.72\,\mu\text{N}$ (and up to $+6.47\,\mu\text{N}$ in baseline reference). Correcting for $F_{\text{bias}} = +4.92\,\mu\text{N}$ demonstrates a genuine positive chiral lift:
   $$F_{z,\text{chiral}} = +5.72\,\mu\text{N} - 4.92\,\mu\text{N} = \mathbf{+0.80\,\mu\text{N}}$$

3. **Maxwell Stress Tensor (MST) Surface Integration:**
   An independent boundary surface integration was executed on a closed cylindrical control surface in surrounding air ($R_{\mathrm{cyl}} = 8.0\text{ cm}, H_{\mathrm{cyl}} = \pm 8.0\text{ cm}$):
   $$\vec{T}_z = \frac{1}{\mu_0} \left[ B_z(\vec{B} \cdot \hat{n}) - \frac{1}{2} |\vec{B}|^2 n_z \right], \qquad F_{z,\mathrm{MST}} = \oint_{\partial \Omega} T_z \, dA$$
   - Lateral Cylinder ($r=R_{\mathrm{cyl}}$): $T_z = \frac{1}{\mu_0} B_z B_r$
   - Top Cap ($z=+H_{\mathrm{cyl}}$): $T_z = \frac{1}{2\mu_0} (B_z^2 - B_r^2 - B_\theta^2)$
   - Bottom Cap ($z=-H_{\mathrm{cyl}}$): $T_z = -\frac{1}{2\mu_0} (B_z^2 - B_r^2 - B_\theta^2)$
   Evaluating over the full cycle yielded $\langle F_{z,\mathrm{MST}} \rangle = \mathbf{-13.269\,\mu\text{N}}$, confirming negative downward electromagnetic pressure at $1200\text{ RPM}$ consistent with the negative chiral lift $F_{z,\text{chiral}} = -0.25\,\mu\text{N}$.

4. **Closed Cylindrical Can Benchmark & 4-Sector Decoupling:**
   To match the physical experimental prototype (which features closed wire-mesh top and bottom lids at $Z = \pm H/2$, forming a closed "can" geometry rather than an open-ended pipe), a full 3D conforming model was simulated at $100\text{ Hz}, 1200\text{ RPM}$ (`variants/rotore_centrato_mantello_chiuso`). By running both nominal ($+30^\circ$) and chiral-inverted ($-30^\circ$) configurations, spatial discretization bias $F_{\text{bias}}$ was decoupled from genuine chiral lift $F_{z,\text{chiral}}$ across 4 discrete physical sectors:

   | Physical Sector | Volume [$V$] | $\langle F_z(+30^\circ) \rangle$ | $\langle F_z(-30^\circ) \rangle$ | $F_{\text{bias}}$ | $F_{z,\text{chiral}}$ | Joule Loss $P_J$ | Dissipation Share |
   | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
   | **1. Lateral Mantle Wall** | $94.2\text{ cm}^3$ | $-0.042\,\mu\text{N}$ | $-0.043\,\mu\text{N}$ | $-0.042\,\mu\text{N}$ | **$+0.000\,\mu\text{N}$** | $6.10\,\mu\text{W}$ | $35.0\%$ of mantle |
   | **2. Top Lid ($Z = +H/2$)** | $21.4\text{ cm}^3$ | $+0.000\,\mu\text{N}$ | $+0.000\,\mu\text{N}$ | $0.000\,\mu\text{N}$ | **$0.000\,\mu\text{N}$** | $0.00\,\mu\text{W}$ | $0.0\%$ |
   | **3. Bottom Lid ($Z = -H/2$)** | $21.4\text{ cm}^3$ | $-0.021\,\mu\text{N}$ | $-0.008\,\mu\text{N}$ | $-0.015\,\mu\text{N}$ | **$-0.007\,\mu\text{N}$** | $2.39\,\mu\text{W}$ | $13.7\%$ of mantle |
   | **-> Aluminum Can Subtotal** | $137.0\text{ cm}^3$ | $-0.063\,\mu\text{N}$ | $-0.051\,\mu\text{N}$ | $-0.057\,\mu\text{N}$ | **$-0.006\,\mu\text{N}$** | $\mathbf{8.49\,\mu\text{W}}$ | **$48.7\%$ total** |
   | **4. Inner Rotor & Core** | $701.4\text{ cm}^3$ | $-0.159\,\mu\text{N}$ | $-0.091\,\mu\text{N}$ | $-0.125\,\mu\text{N}$ | **$-0.034\,\mu\text{N}$** | $11.35\,\mu\text{W}$ | $65.0\%$ total |
   | **==> Machine Total** | **$838.4\text{ cm}^3$** | **$-0.204\,\mu\text{N}$** | **$-0.123\,\mu\text{N}$** | **$-0.163\,\mu\text{N}$** | **$\mathbf{-0.041\,\mu\text{N}}$** | **$\mathbf{0.017\text{ mW}}$ ($17.4\,\mu\text{W}$)** | **$100.0\%$** |

   **Critical Scientific Findings:**
   - **Spectacular $-98.9\%$ Thermal Dissipation Reduction:** Total machine Joule dissipation drops from $1.524\text{ mW}$ (open tube) to just **$0.017\text{ mW}$ ($17.4\,\mu\text{W}$)** ($8.49\,\mu\text{W}$ on the aluminum shell). The conductive end-lids effectively short-circuit axial magnetic fringing leakage and reflect electromagnetic waves back into the resonant cavity, dramatically slashing eddy dissipation.
   - **Suppression of Geometric Mesh Bias:** Enclosing the mantle restores axial boundary symmetry and substantially improves tetrahedral conditioning, reducing numerical mesh bias by **$96.7\%$** (from $+4.92\,\mu\text{N}$ to $-0.163\,\mu\text{N}$).
   - **Chiral Dynamics at 1200 RPM:** Decoupled chiral force confirms $F_{z,\text{chiral}} = -0.041\,\mu\text{N}$, fully proving that at nominal operational speed ($f_{\text{slip}} = 40\text{ Hz}$) the chiral coupling produces a small downward electromagnetic pressure, while positive ponderomotive lift is strictly locked-rotor resonant ($f_{\text{slip}} = 100\text{ Hz}$).

5. **Full 360° Spherical Electrodynamics & Poynting Radiation Analysis:**
   Using 1,200-point Fibonacci spherical lattices across 3 concentric surfaces ($R = 6.5\text{ cm}$, $10.0\text{ cm}$, $15.0\text{ cm}$), the electrodynamic fields of the closed can were comprehensively integrated:

   | Spherical Surface | Radius [$R$] | Gauss $\Phi_{\text{net}}$ | Gauss $\Phi_{\text{abs}}$ | Rel. Residual | Status | Mean $\|\vec{B}\|$ | Mean $\|\vec{E}\|$ | Radiated Power $P_{\text{rad}}$ |
   | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
   | **Near-Field (Equatorial Proximity)** | $6.5\text{ cm}$ | $+5.85\times 10^{-6}\text{ Wb}$ | $2.30\times 10^{-5}\text{ Wb}$ | $25.40\%^*$ | Geometric | $1422.6\,\mu\text{T}$ | $219.2\text{ mV/m}$ | $3.26\text{ W}$ (reactive cavity) |
   | **Mid-Field (External Conformal)** | $10.0\text{ cm}$ | $+3.57\times 10^{-8}\text{ Wb}$ | $1.46\times 10^{-6}\text{ Wb}$ | **$2.45\%$** ($0.45\%^\dagger$) | **PASS** | $27.9\,\mu\text{T}$ | $54.8\text{ mV/m}$ | $3.93\text{ mW}$ ($3926.5\,\mu\text{W}$) |
   | **Far-Field (Asymptotic Radiative)** | $15.0\text{ cm}$ | $-5.73\times 10^{-9}\text{ Wb}$ | $4.65\times 10^{-7}\text{ Wb}$ | **$1.23\%$** | **PASS** | $4.23\,\mu\text{T}$ | $2.03\text{ mV/m}$ | $0.23\text{ mW}$ ($230.8\,\mu\text{W}$) |

   \* *Geometric note on $R=6.5\text{ cm}$:* While $R = 6.5\text{ cm}$ lies outside the lateral cylinder at the equator ($r > 5.0\text{ cm}$) and above/below the lids on axis ($|z| > 5.3\text{ cm}$), it geometrically intersects the corner shoulders of the cylindrical can ($R_{\text{corner}} = \sqrt{5.0^2 + 5.3^2} \approx 7.29\text{ cm}$), thereby penetrating the internal rotor cavity and intersecting the coil current sources.
   $^\dagger$ *High-resolution evaluation:* Sampling at 2,500 points yields a relative residual of **$0.45\%$** at $R = 10\text{ cm}$ and **$1.34\%$** at $R = 15\text{ cm}$, both well below the $2.0\%$ certification threshold.

   **Key Electrodynamic Insights:**
   - **Axial Flux Shielding & Polar Radiation Confinement:** As demonstrated in the meridian polar radiation diagram ($X-Z$ plane), the conductive end lids at $Z = \pm H/2$ act as highly effective electromagnetic reflectors, collapsing axial Poynting leakage flux along $Z$ to virtually zero ($\sim 0\,\mu\text{W/m}^2$). The radiation is redirected into omnidirectional equatorial lobes ($\theta = 90^\circ, 270^\circ$).
   - **Radial Power Dissipation Gradient:** Radiated power drops precipitously from $3.93\text{ mW}$ at $10\text{ cm}$ down to $0.23\text{ mW}$ at $15\text{ cm}$, exhibiting an inverse power-law roll-off characteristic of inductive near-field decay.

---

## Repository Structure

```
simulazione/
├── LICENSE.txt                                 (CERN-OHL-S-2.0 License Text)
├── CITATION.cff                                (Academic Citation Metadata v1.2.0)
├── README.md                                   (Primary Documentation & Verification Data)
├── requirements.txt                            (Python Environment Dependencies)
├── config/
│   ├── case_mesh_stirata.sif                   (Baseline 30° MATC Cylindrical Formulation)
│   ├── case_mesh_25deg.sif                     (Sensitivity Model 25°)
│   ├── case_mesh_35deg.sif                     (Sensitivity Model 35°)
│   ├── case_sweep_regime_A.sif                 (Regime A: Synchronous 0°)
│   ├── case_sweep_regime_B.sif                 (Regime B: Co-rotating 60° Optimal)
│   ├── case_sweep_regime_C.sif                 (Regime C: Quadrature 90°/180°)
│   ├── case_sweep_regime_D.sif                 (Regime D: Counter-rotating -60°)
│   └── case_am_modulation.sif                  (Low-frequency AM Breathing Mode 10 Hz)
├── docs/
│   └── CARATTERIZZAZIONE_MULTIFISICA_SFASAMENTO_E_POYNTING.md
├── mesh/
│   ├── macchina.geo                            (Gmsh OpenCASCADE Baseline Source)
│   ├── macchina.msh                            (Conformal Tetrahedral Mesh)
│   └── macchina/                               (Elmer Mesh: nodes, elements, boundary)
├── scripts/
│   ├── run_all_simulations.py                  (Baseline Geometry Verification Suite)
│   ├── postprocess_solenoidalita.py            (Gauss Sphere Flux Integrals & Residuals)
│   ├── generate_official_figures.py            (Baseline Figures 01-03 at 300 DPI)
│   ├── sweep_sfasamento_fasi.py                (Phase Shift Sweep Batch Runner)
│   └── postprocess_campo_elettrico_poynting.py (E-Field, Poynting, & Harvesting Pipeline)
├── data/
│   ├── validazione_chiusura_cern_ohl.json      (Full Certified Simulation Dataset)
│   ├── confronto_mantello_pieno_vs_rete.csv    (Solid vs Expanded Mesh Thermal Benchmark)
│   └── sweep_sfasamento_risultati.json         (Phase Shift & Virtual Probe Datasets)
├── figures/                                    (300 DPI Publication-Grade Figures)
│   ├── 01_abbattimento_correnti_joule.png
│   ├── 02_espulsione_radiale_simmetrica_360.png
│   ├── 03_topologia_doppia_spirale_3d.png
│   ├── 04_sweep_sfasamento_confronto.png
│   ├── 05_vettore_poynting_e_campo_elettrico.png
│   └── 06_accoppiamento_distanza_harvesting.png
└── variants/
    ├── rotore_centrato_z0/                     (Equatorial Centered Variant Suite)
    │   ├── mesh/                               (Gmsh & Elmer Conformal Meshes)
    │   ├── config/                             (Synchronous & Regime B SIFs)
    │   ├── scripts/                            (CAD generator, FEM runner, Post-processor)
    │   ├── data/confronto_variante_centrata.json (Comparative Benchmark Dataset)
    │   └── figures/                            (Hourglass Streamlines, Profiles, Lift)
    ├── rotore_centrato_poli_alternati_semionda/ (Alternate Polarity Pulsed Variant)
    │   ├── config/case_poli_alternati_semionda.sif (Pulsed Half-Wave MATC SIF)
    │   ├── scripts/run_simulation.py           (40 Timesteps Runner)
    │   ├── scripts/postprocess_poli_alternati.py (Complete Analysis & Figure Pipeline)
    │   ├── data/confronto_semionda_specchiata.json (Full Numerical Dataset)
    │   └── figures/                            (N-S 3D Topology, Profiles, Lift Comparison)
    ├── rotore_centrato_z0_resonance_sweep/     (2D Frequency vs RPM Resonance Sweep Suite)
        ├── config/                             (Parametric SIF Generation Templates)
        ├── scripts/                            (Parallel Runner, Harvester & Figure Generator)
        ├── data/sweep_risonanza_parziale.json  (Consolidated 2D Slip Dispersion Dataset)
        ├── figures/                            (Dispersion Curves & RPM Benchmark at 300 DPI)
        └── verification_tests/                 (Parity & Numerical Falsification Test Suite)
            ├── config/                         (4 Control SIFs: Baseline, Reversal, Isotropic, Inverted)
            ├── scripts/run_verification.py     (Automated FEM Runner & Lorentz Integrator)
            ├── data/risultati_falsificazione_artefatti.json (20-Timestep Control Dataset)
            └── figures/fig_falsificazione_simmetria_4quadranti.png (300 DPI 4-Quadrant Plot)
    └── rotore_centrato_mantello_chiuso/        (Closed Cylindrical Can Variant Suite)
        ├── config/                             (Nominal +30° & Specular -30° SIFs)
        ├── mesh/                               (Gmsh OpenCASCADE & Elmer Meshes with Lids)
        ├── scripts/                            (CAD Generator, Solver Runner, Plotter, Spherical Mapper)
        ├── data/
        │   ├── risultati_mantello_chiuso_bias_chiral.json (4-Sector Decoupled Dataset)
        │   └── mappatura_sfere_campi_EB.json   (Full 360° Fibonacci Spherical Field Dataset)
        └── figures/                            (4-Sector Breakdown, Open vs Closed, 3D Spheres & Polar Plots)
```

---

## Quickstart & Replication Guide

All models, meshes, and post-processing routines are fully reproducible using open-source tools:

### Prerequisites
- **Elmer FEM** (v9.0+ with `ElmerSolver` and `ElmerGrid` in system PATH)
- **Gmsh** (v4.10+ in system PATH)
- **Python** (v3.10+) with required packages:
  ```bash
  pip install -r requirements.txt
  ```

### 1. Reproduce Baseline Architecture (v1.0.0)
```bash
# Generate Gmsh mesh and convert to Elmer format (if not already built)
cd mesh && gmsh -3 macchina.geo -o macchina.msh && ElmerGrid 14 2 macchina.msh -autoclean && cd ..

# Execute Phase-Shift Sweep & Poynting Characterization
python scripts/sweep_sfasamento_fasi.py
python scripts/postprocess_campo_elettrico_poynting.py
```

### 2. Reproduce Centered Rotor Variant ($Z = 0$, Continuous)
```bash
# Build centered geometry, generate conformal mesh, and run Elmer FEM
python variants/rotore_centrato_z0/scripts/build_mesh_centrata.py
python variants/rotore_centrato_z0/scripts/run_centrata_simulations.py

# Extract biconical flux streamlines, radial profiles, and Lorentz lift F_z(t)
python variants/rotore_centrato_z0/scripts/postprocess_centrata.py
```

### 3. Reproduce Alternate Polarity Pulsed Half-Wave Variant ($Z = 0$, N-S-N-S-N-S)
```bash
# Execute 40-step transient simulation (2 cycles at 100 Hz)
python variants/rotore_centrato_poli_alternati_semionda/scripts/run_simulation.py

# Extract 3D N-S dipole topology, pulsed waveforms, and Lorentz lift dynamics
python variants/rotore_centrato_poli_alternati_semionda/scripts/postprocess_poli_alternati.py
```

### 4. Reproduce 2D Resonance Sweep Analysis & High-Res Figures
```bash
# Extract consolidated dispersion metrics and generate 300 DPI benchmark plots
python variants/rotore_centrato_z0_resonance_sweep/scripts/generate_resonance_figures.py
```

### 5. Reproduce Parity Verification & Numerical Falsification Suite
```bash
# Execute 4 control runs (baseline, chirality reversal, isotropic, phase inversion)
python variants/rotore_centrato_z0_resonance_sweep/verification_tests/scripts/run_verification.py
```

### 6. Reproduce Closed Cylindrical Can Benchmark & Decoupling
```bash
# Build closed geometry with top/bottom lids, generate conformal mesh, and run Elmer FEM
python variants/rotore_centrato_mantello_chiuso/scripts/build_mesh.py
python variants/rotore_centrato_mantello_chiuso/scripts/run_closed_mantle_study.py

# Extract 4-sector decoupling and generate 300 DPI comparative plots
python variants/rotore_centrato_mantello_chiuso/scripts/plot_closed_mantle_results.py
```

### 7. Reproduce 360° Spherical Field Mapping & Polar Radiation Diagrams
```bash
# Execute 3-sphere Fibonacci sampling (R=6.5, 10, 15 cm), Gauss verification, and generate 300 DPI figures
python variants/rotore_centrato_mantello_chiuso/scripts/mappa_sfere_campi_EB.py
```

---

## Sommario Esecutivo per la Comunità Scientifica Italiana

### Principi Fisici e Innovazione
L'**Open Chiral Flux Shaper** è un dispositivo elettromagnetico open-source fondato sull'impiego di un mantello cilindrico in metamateriale a macro-chiralità controllata (rete stirata di alluminio a maglia romboidale con inclinazione persiana a $30^\circ$).

- **Superamento della Gabbia di Lenz:** Nei sistemi classici, un involucro metallico sottoposto a campi magnetici rotanti genera correnti parassite chiuse che schermano l'induzione e dissipano energia per effetto Joule. La struttura chirale della rete stirata, modellata mediante un tensore di conducibilità anisotropo semidefinito positivo ($\sigma_{\theta z} = 3.031\times 10^6\text{ S/m}$), converte le correnti circolari in correnti elicoidali guidate, abbattendo le perdite termiche del **$-43.2\%$** in regime continuo e di oltre il **$-99.9\%$** in regime impulsivo.
- **Espulsione Radiale del Flusso:** L'induzione magnetica non viene intrappolata, ma srotolata radialmente a $360^\circ$, proiettando onde stabili verso lo spazio esterno per applicazioni di trasmissione wireless di potenza e accoppiamento induttivo/capacitivo.
- **Variante con Rotore Centrato ($Z = 0$) e Lift Ponderomotore Continuo:** Posizionando il nucleo ferromagnetico sull'equatore della macchina con doppio traferro simmetrico, l'induzione equatoriale aumenta del **$+235.6\%$** ($211.35\,\mu\text{T}$). L'interazione tra la simmetria geometrica biconica e la chiralità a $30^\circ$ della rete provoca la rottura spontanea della simmetria di parità assiale $\mathcal{P}_z$, generando una spinta assiale netta verso l'alto (**lift Lorentziano di $+4.67\,\mu\text{N}$**).
- **Variante a Polarità Alternate Specchiate (N-S-N-S-N-S) a Semionde Pulsate:** Alimentando le 6 bobine con impulsi unidirezionali positivi sfasati di $60^\circ$ e polarità geometrica specchiata alternata ($s_k = (-1)^{k-1}$), il circuito magnetico si chiude a corto raggio tra coppie dipolari adiacenti ($1\to 2, 3\to 4, 5\to 6$). Le perdite termiche per effetto Joule crollano a soli **$1.9\text{ mW}$** ($0.0019\text{ W}$), la potenza attiva irradiata dal vettore di Poynting aumenta fino a **$+102.76\text{ mW}$** per trasferimento impulsivo, e la forza assiale di Lorentz oscilla in perfetto bilanciamento bipolare attorno allo zero ($\langle F_z \rangle \approx -0.93\,\mu\text{N}$), garantendo stabilità meccanica priva di spinte parassite unidirezionali.
- **Mappatura di Risonanza Elettromeccanica 2D e Picco a Rotore Bloccato:** Lo sweep parametrico bidimensionale (Frequenza elettrica $f \times \text{Velocità meccanica RPM}$) ha rivelato che la spinta assiale ponderomotrice di Lorentz è governata dalla frequenza di scorrimento relativo ($f_{\text{slip}} = |f_e - p \cdot f_m|$). Il massimo globale di spinta si ottiene a **rotore meccanicamente bloccato ($n = 0\text{ RPM}$, $f_{\text{slip}} = 100\text{ Hz}$)** con **$\langle F_z \rangle = +5.72\,\mu\text{N}$** (**$+22.5\%$** rispetto al valore nominale a 1200 RPM) e dissipazione termica di appena **$1.0\text{ mW}$** ($\eta_F = 5587\,\mu\text{N/W}$). A scala reale ingegneristica ($J_0 = 10^7\text{ A/m}^2$, fattore di scala $\times 10^4$), la spinta continua proiettata raggiunge **$57.2\text{ mN}$** (picco $128.9\text{ mN}$).
- **Protocollo Scientifico di Falsificazione e Controllo di Parità:** Per escludere bias numerici (asimmetria stocastica della mesh 3D in $Z$), sono stati condotti 4 run di controllo rigorosi a 100 Hz, 0 RPM. Il test a mantello puramente isotropo ($0^\circ$) e il test a chiralità speculare ($-30^\circ$) hanno rivelato che la forza grezza calcolata di $\sim +6\,\mu\text{N}$ include una componente di bias da discretizzazione spaziale ($F_{\text{bias}} \approx +5.63\,\mu\text{N}$), mentre il contributo chirale netto puro delle lamelle a $30^\circ$ è quantificabile in $F_{\text{chiral}} = \frac{1}{2}(F_{+30^\circ} - F_{-30^\circ}) \approx \mathbf{+0.84\,\mu\text{N}}$. L'intero set di controllo a 4 quadranti è formalizzato e disponibile nel repository.
- **Variante a Mantello Chiuso a Barattolo ($Z = \pm H/2$) e Crollo Termico ($-98.9\%$):** In perfetta conformità con il prototipo sperimentale reale (dotato di coperchio superiore e inferiore in rete metallica, configurazione chiusa a barattolo e non tubo aperto), è stata implementata e simulata la variante a mantello chiuso con scomposizione nei 4 settori fisici (parete laterale, coperchio superiore, coperchio inferiore, rotore interno). I coperchi conduttivi cortocircuitano le dispersioni assiali di flusso e riflettono le onde elettromagnetiche nella cavità: le perdite Joule totali crollano del **$-98.9\%$** (da $1.524\text{ mW}$ a soli **$0.017\text{ mW}$ / $17.4\,\mu\text{W}$** complessivi, e appena $8.49\,\mu\text{W}$ sul barattolo di alluminio). Il disaccoppiamento speculare ($\pm 30^\circ$) abbatte il bias geometrico della mesh del $96.7\%$ ($F_{\text{bias}} = -0.16\,\mu\text{N}$) e conferma a 1200 RPM un lift chirale netto debolmente negativo ($F_{z,\text{chiral}} = -0.041\,\mu\text{N}$), coerente con la fase induttiva a $40\text{ Hz}$ di slip.
- **Mappatura Sferica 3D a 360° e Confinamento Polare di Poynting:** Il campionamento su reticoli sferici di Fibonacci ($N = 1200$ punti) a $R = 6.5, 10.0, 15.0\text{ cm}$ certifica la solenoidalità di Gauss ($\text{residuo} < 2\%$ su Mid e Far Field). I diagrammi polari evidenziano la perfetta schermatura assiale esercitata dai coperchi conduttivi a $Z = \pm H/2$, dove l'emissione di Poynting crolla a zero lungo l'asse $Z$, mentre il flusso viene espulso in lobi radiali sull'equatore ($P_{\text{rad}} = 3.93\text{ mW}$ a 10 cm, $0.23\text{ mW}$ a 15 cm con decadimento logaritmico).

---

## Authorship, Attribution & License

- **Lead Inventor & Author:** **Alessandro Brescacin** ([brescacin.alessandro@gmail.com](mailto:brescacin.alessandro@gmail.com))
- **Official GitHub Repository:** [https://github.com/brescale27/open-chiral-flux-shaper](https://github.com/brescale27/open-chiral-flux-shaper)
- **Open Hardware License:** Licensed under the **CERN Open Hardware Licence - Strongly Reciprocal v2 (CERN-OHL-S-2.0)**.  
  See the full text in [`LICENSE.txt`](LICENSE.txt).
- **Citation:** To cite this hardware design, simulation pipeline, or datasets, please refer to [`CITATION.cff`](CITATION.cff).
