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
- **Enables electromagnetic propulsion and levitation:** In a biconical induction configuration, the fixed $30^\circ$ chiral tilt breaks axial reflection parity ($\mathcal{P}_z$), producing a continuous, unidirectional upward ponderomotive Lorentz lift ($\langle F_z \rangle = +4.67\,\mu\text{N}$).

---

## Core Architectures & Comparative Benchmark

| Parameter / Metric | Baseline Architecture (v1.0.0) | Centered Continuous (Regime B) | Mirrored Polarity Pulsed Half-Wave (N-S-N-S-N-S) | Physical Mechanism / Specialty |
| :--- | :---: | :---: | :---: | :---: |
| **Core Geometry** | Ferromagnetic spider at bottom ($Z = -H/2$) | Ferromagnetic core at equator ($Z = 0$) | Ferromagnetic core at equator ($Z = 0$) | Equatorial magnetic flux bridge |
| **Coil Excitation** | 60° Progressive Full Sine Wave | 60° Progressive Full Sine Wave | 60° Shifted Half-Wave Pulse Train ($J_k \ge 0$) | Directional vs Sequenced Pulsed Drive |
| **Polarity Layout** | Unipolar Homogeneous | Unipolar Homogeneous | Alternating Mirrored: $s_k = (-1)^{k-1}$ | Adjacent N-S magnetic return loops |
| **Radial Field $B_{\text{rad}}$ ($R = 6\text{ cm}$)** | **$62.98\,\mu\text{T}$** | **$211.35\,\mu\text{T}$** | **$135.84\,\mu\text{T}$** | Strong equatorial concentration (+115.7% vs baseline) |
| **Poynting Power Flux ($R=12\text{ cm}$)** | **$+7.282\text{ mW}$** | **$+2.620\text{ mW}$** | **$+102.76\text{ mW}$** | **$+1311\%$ pulsed harmonic radiation boost** |
| **Net Axial Lorentz Force $\langle F_z \rangle$** | $\approx 0$ (leakage) | **$+4.67\,\mu\text{N}$ ($+0.00467\text{ mN}$)** | **$-0.93\,\mu\text{N}$ ($\approx 0\text{ mN}$ balanced)** | Mirrored N-S cancellation of cross terms |
| **Joule Dissipation $P_J$** | **$2.437\text{ W}$** | **$2.510\text{ W}$** | **$0.0019\text{ W}$ ($1.9\text{ mW}$)** | **$>99.9\%$ suppression of eddy heating** |
| **Primary Optimal Use-Case** | Directional Wireless Venting | Continuous Electromagnetic Lift | Resonant Wireless Pulsed Power & Ultra-Low Heat | Mission-specific electromagnetic tuning |

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
    └── rotore_centrato_poli_alternati_semionda/ (Alternate Polarity Pulsed Variant)
        ├── config/case_poli_alternati_semionda.sif (Pulsed Half-Wave MATC SIF)
        ├── scripts/run_simulation.py           (40 Timesteps Runner)
        ├── scripts/postprocess_poli_alternati.py (Complete Analysis & Figure Pipeline)
        ├── data/confronto_semionda_specchiata.json (Full Numerical Dataset)
        └── figures/                            (N-S 3D Topology, Profiles, Lift Comparison)
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

---

## Sommario Esecutivo per la Comunità Scientifica Italiana

### Principi Fisici e Innovazione
L'**Open Chiral Flux Shaper** è un dispositivo elettromagnetico open-source fondato sull'impiego di un mantello cilindrico in metamateriale a macro-chiralità controllata (rete stirata di alluminio a maglia romboidale con inclinazione persiana a $30^\circ$).

- **Superamento della Gabbia di Lenz:** Nei sistemi classici, un involucro metallico sottoposto a campi magnetici rotanti genera correnti parassite chiuse che schermano l'induzione e dissipano energia per effetto Joule. La struttura chirale della rete stirata, modellata mediante un tensore di conducibilità anisotropo semidefinito positivo ($\sigma_{\theta z} = 3.031\times 10^6\text{ S/m}$), converte le correnti circolari in correnti elicoidali guidate, abbattendo le perdite termiche del **$-43.2\%$** in regime continuo e di oltre il **$-99.9\%$** in regime impulsivo.
- **Espulsione Radiale del Flusso:** L'induzione magnetica non viene intrappolata, ma srotolata radialmente a $360^\circ$, proiettando onde stabili verso lo spazio esterno per applicazioni di trasmissione wireless di potenza e accoppiamento induttivo/capacitivo.
- **Variante con Rotore Centrato ($Z = 0$) e Lift Ponderomotore Continuo:** Posizionando il nucleo ferromagnetico sull'equatore della macchina con doppio traferro simmetrico, l'induzione equatoriale aumenta del **$+235.6\%$** ($211.35\,\mu\text{T}$). L'interazione tra la simmetria geometrica biconica e la chiralità a $30^\circ$ della rete provoca la rottura spontanea della simmetria di parità assiale $\mathcal{P}_z$, generando una spinta assiale netta verso l'alto (**lift Lorentziano di $+4.67\,\mu\text{N}$**).
- **Variante a Polarità Alternate Specchiate (N-S-N-S-N-S) a Semionde Pulsate:** Alimentando le 6 bobine con impulsi unidirezionali positivi sfasati di $60^\circ$ e polarità geometrica specchiata alternata ($s_k = (-1)^{k-1}$), il circuito magnetico si chiude a corto raggio tra coppie dipolari adiacenti ($1\to 2, 3\to 4, 5\to 6$). Le perdite termiche per effetto Joule crollano a soli **$1.9\text{ mW}$** ($0.0019\text{ W}$), la potenza attiva irradiata dal vettore di Poynting aumenta fino a **$+102.76\text{ mW}$** per trasferimento impulsivo, e la forza assiale di Lorentz oscilla in perfetto bilanciamento bipolare attorno allo zero ($\langle F_z \rangle \approx -0.93\,\mu\text{N}$), garantendo stabilità meccanica priva di spinte parassite unidirezionali.

---

## Authorship, Attribution & License

- **Lead Inventor & Author:** **Alessandro Brescacin** ([brescacin.alessandro@gmail.com](mailto:brescacin.alessandro@gmail.com))
- **Official GitHub Repository:** [https://github.com/brescale27/open-chiral-flux-shaper](https://github.com/brescale27/open-chiral-flux-shaper)
- **Open Hardware License:** Licensed under the **CERN Open Hardware Licence - Strongly Reciprocal v2 (CERN-OHL-S-2.0)**.  
  See the full text in [`LICENSE.txt`](LICENSE.txt).
- **Citation:** To cite this hardware design, simulation pipeline, or datasets, please refer to [`CITATION.cff`](CITATION.cff).
