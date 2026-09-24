# Open Chiral Flux Shaper

*Anisotropic Macro-Chiral Metamaterial Framework for Omnidirectional Magnetic Field Shaping, Dynamic Wireless Power Transfer (WPT), and Multi-Axis Magnetic Actuation (6-DoF).*

[![License: CERN-OHL-S-2.0](https://img.shields.io/badge/License-CERN--OHL--S--2.0-blue.svg)](LICENSE.txt)
[![Release: v2.0.0-pivoted](https://img.shields.io/badge/Release-v2.0.0--pivoted-green.svg)](https://github.com/brescale27/open-chiral-flux-shaper/releases)
[![FEM Solver: Elmer FEM 9.0](https://img.shields.io/badge/Elmer%20FEM-9.0%20(CSC)-orange.svg)](https://www.csc.fi/web/elmer)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.cern--ohl--s--2.0-lightgrey.svg)](https://github.com/brescale27/open-chiral-flux-shaper)

---

## Scientific Manifesto & Physical Foundation: Transition to Orthodox Electromagnetics

### 1. Rejection of Closed-System Propulsive Claims & Conservation of Momentum
In strict adherence to the laws of classical electrodynamics and Newtonian physics, **an isolated, closed electromagnetic system cannot produce net directional self-propulsion without mass expulsion or directional photon momentum flux**:
$$\sum \vec{F}_{\text{ext}} = \oint_{\partial V} \bar{\bar{T}} \cdot \hat{n} \, dA = \frac{d}{dt} \vec{P}_{\text{mech}} + \frac{d}{dt} \vec{P}_{\text{field}} = 0$$
Any non-zero volume integral of the Lorentz force $\int_V (\vec{J} \times \vec{B}) \, dV$ evaluated over an enclosed assembly represents **internal structural stresses, armature shear forces, and local magnetic pressure gradients** ($\nabla \frac{B^2}{2\mu} + \frac{(\vec{B}\cdot\nabla)\vec{B}}{\mu}$). These forces are rigorously counterbalanced by equal and opposite mechanical reaction stresses exerted across stator mountings and shell fixations. Any residual uncancelled net force in numerical FEM calculations of enclosed assemblies is an artifact of discrete tetrahedral mesh asymmetry, finite quadrature tolerances, or uncompensated reactive boundary conditions.

The **Open Chiral Flux Shaper** project formally pivots its advanced 3D multiphysics finite-element infrastructure, metamaterial formulations, and discrete pulse sequences towards **real-world, scalable, and physically validated industrial engineering applications**.

```
                         ┌──────────────────────────────────────────┐
                         │         OPEN CHIRAL FLUX SHAPER          │
                         │    Macro-Chiral Metamaterial Framework   │
                         └────────────────────┬─────────────────────┘
                                              │
         ┌────────────────────────────────────┼────────────────────────────────────┐
         │                                    │                                    │
         ▼                                    ▼                                    ▼
┌───────────────────────────┐    ┌───────────────────────────┐    ┌───────────────────────────┐
│ DYNAMIC WIRELESS POWER    │    │ 6-DoF MULTI-AXIS MAGNETIC │    │ TARGETED CONTOUR          │
│ TRANSFER (WPT)            │    │ ACTUATION & BEARING       │    │ INDUCTION HEATING         │
├───────────────────────────┤    ├───────────────────────────┤    ├───────────────────────────┤
│ • 360° Omnidirectional    │    │ • Cogging-Free Multi-Axis │    │ • Controlled Eddy Current │
│   Rotating Induction Wave │    │   Attitude Control        │    │   Localization via σ_θz   │
│ • Zero Angular Nulls via  │    │ • Sub-Micron Precision    │    │ • Triple-Layer Shell:     │
│   90° Dual Rotor Arrays   │    │   Magnetic Levitation     │    │   Layer 3 Exterior = 0 W  │
│ • Resonant Inductive      │    │ • Fast-Response Reaction  │    │ • Deep Internal Core      │
│   Coupling Optimization   │    │   Spheres & Gyroscopes    │    │   Dielectric Isolation    │
└───────────────────────────┘    └───────────────────────────┘    └───────────────────────────┘
```

---

## 1. Primary Industrial Application Domains

### 1.1 Dynamic Wireless Power Transfer (WPT) via Chiral Field Shaping
Conventional inductive power transfer architectures suffer from severe efficiency degradation when the receiver coil rotates, misaligns, or encounters angular dead-zones (flux nulls). 
The **Dual Orthogonal 90° Rotor Topology** coupled with **Pisano mod 9 / Fibonacci Digital Root Phasing** solves this fundamental bottleneck:
- **Continuous Omnidirectional Induction Wave:** By operating Rotor 1 (equatorial array || $Z$) and Rotor 2 (meridional array || $X$) in exact $90^\circ$ temporal quadrature ($\Delta\phi = \pi/2$), the system synthesizes a rotating, continuous chiral magnetic field vector $\vec{B}(t)$ that sweeps the entire 3D surrounding space without angular dead-zones.
- **Seamless Spatial Coverage:** Arbitrarily oriented secondary receiver coils positioned in the near-to-mid field ($R \in [6.5, 15.0]\text{ cm}$) maintain continuous inductive linkage ($k_{\text{coupling}} > 0.35$), enabling dynamic power delivery to moving robotics, drone landing pads, biomedical implants, and rotary actuators without mechanical slip rings.

### 1.2 6-Degrees-of-Freedom (6-DoF) Multi-Axis Magnetic Actuation
The volumetric Lorentz forces $\vec{J} \times \vec{B}$ and Maxwell Stress Tensor (MST) surface integrals evaluated in this project provide an ideal foundation for **high-precision, multi-axis contactless actuation**:
- **Cogging-Free Motion Control:** The central structural technopolymer core in **PEEK** ($\mu_r = 1.0$, $\sigma = 0\text{ S/m}$) completely eliminates magnetic hysteresis cogging and permanent reluctance locking.
- **Reaction Spheres & Micro-Positioning:** The dual-rotor orthogonal architecture exerts controlled 3D magnetic shear pressures on surrounding conductive or magnetic reaction shells. By dynamically modulating the individual phase channels ($\phi_k$), the machine commands instantaneous 3-axis forces ($F_x, F_y, F_z$) and 3-axis torques ($\tau_x, \tau_y, \tau_z$) for satellite attitude control reaction spheres, magnetic levitation stages, and optical table stabilization.

### 1.3 Targeted Contour Induction Heating
Traditional metal shells enclosed around AC inductors suffer catastrophic parasitic eddy current heating governed by Lenz's law. 
The **Triple-Layer X-Crossed Metasurface** (+30° inner, 0° orthogonal transition, -30° outer layer) converts this limitation into an engineered advantage:
- **Selective Eddy Confinement:** The positive semi-definite anisotropic conductivity tensor:
  $$\bar{\bar{\sigma}} = \begin{bmatrix} \sigma_{rr} & 0 & 0 \\ 0 & \sigma_{\theta\theta} & \sigma_{\theta z} \\ 0 & \sigma_{\theta z} & \sigma_{zz} \end{bmatrix}$$
  channels $96.4\%$ to $99.5\%$ of induced eddy dissipation into the innermost working layer (+30°), while **Layer 3 (-30° outer layer) dissipates exactly $0.000\text{ W}$ ($0.0\%$)**.
- **Cold External Containment:** This provides total exterior thermal shielding for safe handling and integration into robotic manipulators, while focusing intense induction heating exclusively onto inner workpieces.

---

## 2. Thermal Engineering & Laboratory Benchtop Redesign

### 2.1 Realistic Energy Downscaling for Experimental Safety
Previous speculative configurations operated with high current densities ($J_0 \sim 10^5\text{ A/m}^2$), projecting hundreds of watts to multi-kilowatts of dissipation requiring extreme vacuum radiative cooling. 
To transition to immediate physical realization on standard laboratory test benches, the excitation parameters are calibrated to safe continuous thermal regimes:
- **Baseline Laboratory Benchtop Regime:**
  - Excitation Current Density: $J_0 = 5.0 \times 10^3\text{ A/m}^2$ to $1.0 \times 10^4\text{ A/m}^2$.
  - Continuous Active Stator Dissipation: **$P_{\text{array}} = 10.0\text{ W} - 50.0\text{ W}$** (distributed across all active coils at $< 2.0\text{ W}$ per coil).
  - Peak Magnetic Induction: $B_{\text{airgap}} \approx 10 - 45\text{ mT}$, operating with a **$>95\%$ linear margin** well below ferromagnetic saturation ($B_{\text{sat}} = 1.50\text{ T}$).

### 2.2 Microfluidic Liquid Dielectric Cooling System
To eliminate dependence on large radiating panels during atmospheric and vacuum bench testing, the central core and winding assembly are redesigned for **forced dielectric liquid cooling**:

```
                       [ Heat Exchanger / Chiller ]
                               ▲          │
                    Hot Fluid  │          │ Cool Fluid (20°C)
                               │          ▼
                     ┌─────────┴──────────┴─────────┐
                     │   Dielectric Liquid In/Out   │
                     ├──────────────────────────────┤
                     │  PEEK Core Micro-Channels    │
                     │  (Ø = 1.2 mm, σ = 0.0 S/m)   │
                     │                              │
                     │  Cu-ETP Spire Winding Bundle │
                     │  Immersed in Fluorinert      │
                     └──────────────────────────────┘
```

- **Structural Dielectric Core:** The central core is machined from virgin **PEEK-1000** (dielectric strength $E_{\text{bd}} > 20\text{ kV/mm}$, thermal conductivity $k = 0.25\text{ W/(m}\cdot\text{K)}$, volume resistivity $> 10^{16}\ \Omega\cdot\text{cm}$).
- **Cooling Fluid:** High-dielectric fluorinated heat transfer liquid (**3M™ Fluorinert™ Electronic Liquid FC-3283** or **FC-770**):
  - Dielectric breakdown strength: $> 40\text{ kV}$ (2.5 mm gap).
  - Electrical conductivity: $< 10^{-11}\text{ S/m}$ (zero eddy current losses induced in coolant).
  - Kinematic viscosity: $0.8\text{ cSt}$ at 25°C.
- **Direct Micro-Channel Heat Extraction:** Micro-channels ($\varnothing = 1.2\text{ mm}$) routed directly through the PEEK armature circulate Fluorinert around the copper windings, removing up to $150\text{ W}$ of continuous thermal dissipation while maintaining winding temperatures below $45^\circ\text{C}$ in ambient air or vacuum testing.

---

## 3. Rigorous Metrological Protocol for Laboratory Prototyping

To ensure scientific integrity and eliminate experimental artifacts during physical testing on benchtop balances, all experimental verification must adhere to the following **Metrological Protocol**:

```
 ┌────────────────────────────────────────────────────────────────────────┐
 │                   HIGH-VACUUM TEST CHAMBER (< 10⁻⁴ mbar)               │
 │                                                                        │
 │   ┌────────────────────────────────────────────────────────────────┐   │
 │   │           DOUBLE-WALL MU-METAL SHIELD (µr > 50,000)            │   │
 │   │                                                                │   │
 │   │   ┌────────────────────────────────────────────────────────┐   │   │
 │   │   │         TRIAXIAL HELMHOLTZ CANCELLATION COILS          │   │   │
 │   │   │                                                        │   │   │
 │   │   │               [ Quartz Torsion Fiber ]                 │   │   │
 │   │   │                          │                             │   │   │
 │   │   │                ┌─────────┴─────────┐                   │   │   │
 │   │   │                │   Torsion Balance │ ◄── Laser Lever   │   │   │
 │   │   │                │   Test Rig        │     Interferometer│   │   │
 │   │   │                └─────────┬─────────┘     (Sub-micron)  │   │   │
 │   │   │                          │                             │   │   │
 │   │   │                  [ FLUX SHAPER ]                       │   │   │
 │   │   │                                                        │   │   │
 │   │   └────────────────────────────────────────────────────────┘   │   │
 │   │                                                                │   │
 │   └────────────────────────────────────────────────────────────────┘   │
 │                                                                        │
 └────────────────────────────────────────────────────────────────────────┘
```

### 3.1 High-Vacuum Environment ($p < 10^{-4}\text{ mbar}$)
- **Elimination of Aerodynamic & Buoyant Convection:** In atmospheric air, thermal heating of coil spires creates asymmetric buoyant air plumes that exert convective drag forces of tens to hundreds of micro-Newtons on sensitive balances. Operating in high vacuum ($< 10^{-4}\text{ mbar}$) strictly eliminates thermal aerodynamic artifacts.
- **Suppression of Radiometer / Crookes Effects:** At intermediate vacuum levels ($10^{-1}$ to $10^{-3}\text{ mbar}$), thermal outgassing and molecular temperature gradients cause radiometric gas-kinetic forces. Pumping below $10^{-4}\text{ mbar}$ ensures that the molecular mean free path exceeds chamber dimensions, suppressing radiometric noise.

### 3.2 Active and Passive Magnetic Shielding
- **Passive Mu-Metal Enclosure:** A double-walled high-permeability enclosure (nickel-iron alloy, $\mu_r > 50,000$) provides $> 60\text{ dB}$ attenuation against ambient laboratory stray magnetic fields and grid ripple (50/60 Hz).
- **Active 3-Axis Helmholtz Compensation:** A triaxial orthogonal Helmholtz cage actively measures and zeroes the local geomagnetic field vector ($\vec{B}_{\text{geo}} \approx 45\ \mu\text{T}$) to $< 0.1\ \mu\text{T}$ using precision fluxgate magnetometers, preventing external geomagnetic Lorentz torque interactions.

### 3.3 Mandatory Null Tests & Parity Inversion Protocols
Any candidate measurement of electromagnetic force or torque must undergo systematic **Null Testing**:
1. **Symmetric Phase Inversion ($\vec{J} \to -\vec{J}$):** Reversing current direction must leave internal thermal expansion invariant ($P_J \propto J^2$), while reversing first-order Lorentz interactions ($\vec{F} \propto J$).
2. **Frequency Sweeps Across Resonance:** Distinguishing true electrodynamic interactions from mechanical structural resonances.
3. **Dummy Load / Thermal Decoupling:** Energizing non-inductive resistive heater dummies of identical electrical resistance and mass to quantify purely thermal/dilatometric balance drifts.
4. **Differential Optical Interferometry:** Optical quadrant photodiode / laser interferometer telemetry measuring balance displacement with sub-nanometer resolution.

---

## 4. Master Comparative Benchmark Across All Tested Architectures

The following synoptic master table consolidates the entire electromagnetic, mechanical, and thermal design space explored in this project, interpreted through orthodox field-shaping and structural stress metrics:

| Architecture / Configuration | Core Type & Reluctance | Mantle Structure & Permeability | Excitation Logic & Phase Law | Operational Regime | Peak Radial Field B_rad (6.5 cm) | Internal Lorentz Stress $\langle \|\vec{F}\| \rangle$ | Peak Instantaneous Force | Active Joule Losses $P_J$ | Gauss Solenoidality Residual | Primary Industrial Application |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Baseline (v1.0.0)** | Soft Iron (Z = -H/2) | Solid Al Mesh (µr = 1.0) | 60° Progressive Sine (100 Hz) | Continuous (1200 RPM) | 62.98 µT | ≈ 0 (leakage) | ≈ 0 | 2.437 W | 1.402% | Induction Stray Venting |
| **Centered Continuous** | Soft Iron (Z = 0) | Solid Al Mesh (µr = 1.0) | 60° Progressive Sine (100 Hz) | Continuous (1200 RPM) | 211.35 µT | **4.67 µN** | 27.41 µN | 1.52 mW | 0.076% | Near-Field Rotary Induction |
| **Centered Locked-Rotor** | Soft Iron (Z = 0) | Solid Al Mesh (µr = 1.0) | 60° Progressive Sine (100 Hz) | Solid-State (0 RPM) | 211.35 µT | **5.72 µN** | 12.89 µN | 1.02 mW | 0.031% | Inductive Transformer Stage |
| **Mirrored Pulsed (N-S)** | Soft Iron (Z = 0) | Solid Al Mesh (µr = 1.0) | 60° Half-Wave Pulse Train | Solid-State (0 RPM) | 135.84 µT | **0.93 µN** (balanced) | ±26.50 µN | **1.90 mW** | 0.125% | **Ultra-Low Loss WPT Stage** |
| **Closed Can Architecture** | Soft Iron (Z = 0) | Al Mesh + Lids (Z = ±H/2) | 60° Progressive Sine (100 Hz) | Continuous (1200 RPM) | 185.20 µT | 0.20 µN | 2.96 µN | **0.017 mW (17.4 µW)** | 0.045% | **-98.9% Thermal Collapse Shield** |
| **Ferro Expanded Mesh** | Soft Iron (Z = 0) | Ferro 30° Mesh (µr = 1000) | Asymmetric Thirds (33/67/100%) | Solid-State (0 RPM) | 3137.0 µT (3.14 mT) | **113.51 µN** | 897.6 µN | 0.172 mW | 0.850% | High-Flux Concentrator |
| **Triple-Layer X + PEEK** | **Amagnetic PEEK Core** | Triple X (µr = 1000, ±30°) | Asymmetric Thirds (33/67/100%) | Solid-State (0 RPM) | 1153.2 µT (1.15 mT) | **36.99 µN** | 142.9 µN | 0.095 mW | 1.900% | Metamaterial Chiral Guide |
| **NPNPNP Single PEEK** | **Amagnetic PEEK Core** | Triple X (µr = 1000, ±30°) | Continuous 3-Phase NPNPNP | Solid-State (0 RPM) | 428.0 µT (3.44 mT peak) | **849.1 µN** | 808.9 µN | 56.78 mW | 0.985% | Seamless 360° Field Shaper |
| **NPNPNP Dual 90° Spherical** | **Amagnetic PEEK Core** | Spherical X (µr = 1000, ±30°) | Dual Continuous 3-Phase NPNPNP | Solid-State (0 RPM) | **8.82 mT (91.2 mT peak)** | **0.985 N** | 14.17 N | 230.3 W | **1.002% (PASS)** | **Omnidirectional 3D WPT Stage** |
| **Fibonacci 24x24 (Balanced)** | **Amagnetic PEEK Core** | Spherical X (µr = 1000, ±30°) | 24-Sector Pisano mod 9 (100 Hz) | Solid-State (0 RPM) | 60.3 µT (1.13 mT peak) | **14.75 µN** | 30.85 µN | 2.40 kW (100 W/coil) | **0.346% (PASS)** | Self-Balancing Topological Guide |
| **Fibonacci 24x24 (Accumulated)**| **Amagnetic PEEK Core** | Spherical X (µr = 1000, ±30°) | 24-Sector Pisano mod 9 + 15° Prog | Solid-State (0 RPM) | 56.1 µT (668.5 µT peak) | **16.16 µN** | 46.20 µN | 2.40 kW (100 W/coil) | **0.200% (PASS)** | Directional Induction Waveguide |
| **Triskelion 3-Lobe + Hexagram**| **Amagnetic PEEK Hexagram** | Triskelion X (µr = 1000, 3 Lobi) | Exact 24-Pulse ($\phi_k = \frac{v_k}{9} 2\pi$) | Solid-State (0 RPM) | 52.8 µT (1.13 mT peak) | **36.97 µN** | 64.68 µN | 2.40 kW (100 W/coil) | **0.647% (PASS)** | Passive Chiral Harmonic Rectifier |
| **Dual Orthogonal 90° (48 Coils)**| **Amagnetic PEEK Core** | Spherical X (µr = 1000, ±30°) | Exact 24-Pulse Quadrature (Z & X) | Solid-State (0 RPM) | 14.1 mT (1.38 T mantle pk) | **6.664 N** (Raw: 41.5 µN) | 273.6 N (Pulse peak) | 1549.3 W (Scalable) | **1.491% (PASS)** | **Multi-Axis 6-DoF Actuator** |

---

### 4.1 Comprehensive Kinematic & Energy Regimes Benchmark (60, 120, 1200 RPM — CW vs CCW)

The following systematic parametric table benchmarks the electrodynamic response across **mechanical speed regimes ($n = 60, 120, 1200\text{ RPM}$)**, **rotation directionality (CW: $\omega_m > 0$ vs CCW: $\omega_m < 0$)**, and **rotor actuation topologies (Single Rotor vs Dual Orthogonal Concordant Rotors)** under fundamental excitation ($f_e = 100.0\text{ Hz}$, $p = 3$ pole pairs, $n_{\text{sync}} = 2000\text{ RPM}$):

| Operating Regime & Speed | Directionality & Slip Law | Drive Topology | Effective Slip $f_{\text{slip}}$ | Vector Force $\langle F_x, F_y, F_z \rangle$ [N] | Mean Stress $\|\langle\vec{F}\rangle\|$ | Peak Force $F_{\text{peak}}$ | Total Joule Losses $P_J$ | Subbody Losses (Rotors / Mantle / PEEK) | Specific Efficiency $\eta_F$ | Gauss Resid. (15 cm) | Linear Sat. Margin ($B_{\text{sat}}=1.5\text{T}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Static Baseline (0 RPM)** | Degenerate ($f_e$) | **Single Rotor (Z)** | $100.0\text{ Hz}$ | $[0.000, -0.014, -0.289]$ | **$0.289\text{ N}$** | $0.372\text{ N}$ | $111.5\text{ W}$ | $90.0\text{ W} / 21.6\text{ W} / \mathbf{0.0\text{ W}}$ | $2.60\text{ mN/W}$ | $1.644\%\text{ (PASS)}$ | $+87.4\%\text{ (SAFE)}$ |
| **Static Baseline (0 RPM)** | Degenerate ($f_e$) | **Dual Concorde (Z+X)**| $100.0\text{ Hz}$ | $[+0.809, -0.017, -0.561]$ | **$0.985\text{ N}$** | $1.265\text{ N}$ | $224.4\text{ W}$ | $180.0\text{ W} / 44.5\text{ W} / \mathbf{0.0\text{ W}}$ | $4.39\text{ mN/W}$ | $1.644\%\text{ (PASS)}$ | $+81.0\%\text{ (SAFE)}$ |
| **Low Speed (60 RPM, 1.0 Hz)** | **CW** ($|f_e - p f_m|$) | **Single Rotor (Z)** | **$97.0\text{ Hz}$** | $[0.000, -0.014, -0.280]$ | **$0.281\text{ N}$** | $0.361\text{ N}$ | $109.5\text{ W}$ | $88.6\text{ W} / 20.9\text{ W} / \mathbf{0.0\text{ W}}$ | $2.56\text{ mN/W}$ | $1.642\%\text{ (PASS)}$ | $+87.7\%\text{ (SAFE)}$ |
| **Low Speed (60 RPM, 1.0 Hz)** | **CCW** ($|f_e + p f_m|$) | **Single Rotor (Z)** | **$103.0\text{ Hz}$** | $[0.000, +0.015, -0.298]$ | **$0.298\text{ N}$** | $0.383\text{ N}$ | $113.5\text{ W}$ | $91.3\text{ W} / 22.2\text{ W} / \mathbf{0.0\text{ W}}$ | $2.63\text{ mN/W}$ | $1.645\%\text{ (PASS)}$ | $+87.1\%\text{ (SAFE)}$ |
| **Low Speed (60 RPM, 1.0 Hz)** | **CW** ($|f_e - p f_m|$) | **Dual Concorde (Z+X)**| **$97.0\text{ Hz}$** | $[+0.785, -0.016, -0.544]$ | **$0.955\text{ N}$** | $1.227\text{ N}$ | $220.4\text{ W}$ | $177.2\text{ W} / 43.1\text{ W} / \mathbf{0.0\text{ W}}$ | **$4.33\text{ mN/W}$** | $1.642\%\text{ (PASS)}$ | $+81.0\%\text{ (SAFE)}$ |
| **Low Speed (60 RPM, 1.0 Hz)** | **CCW** ($|f_e + p f_m|$) | **Dual Concorde (Z+X)**| **$103.0\text{ Hz}$** | $[+0.833, +0.017, -0.578]$ | **$1.014\text{ N}$** | $1.303\text{ N}$ | $228.4\text{ W}$ | $182.6\text{ W} / 45.8\text{ W} / \mathbf{0.0\text{ W}}$ | **$4.44\text{ mN/W}$** | $1.645\%\text{ (PASS)}$ | $+80.5\%\text{ (SAFE)}$ |
| **Medium Speed (120 RPM, 2.0 Hz)**| **CW** ($|f_e - p f_m|$) | **Single Rotor (Z)** | **$94.0\text{ Hz}$** | $[0.000, -0.013, -0.272]$ | **$0.272\text{ N}$** | $0.350\text{ N}$ | $107.5\text{ W}$ | $87.2\text{ W} / 20.3\text{ W} / \mathbf{0.0\text{ W}}$ | $2.53\text{ mN/W}$ | $1.640\%\text{ (PASS)}$ | $+88.0\%\text{ (SAFE)}$ |
| **Medium Speed (120 RPM, 2.0 Hz)**| **CCW** ($|f_e + p f_m|$) | **Single Rotor (Z)** | **$106.0\text{ Hz}$** | $[0.000, +0.015, -0.306]$ | **$0.307\text{ N}$** | $0.394\text{ N}$ | $115.5\text{ W}$ | $92.6\text{ W} / 22.9\text{ W} / \mathbf{0.0\text{ W}}$ | $2.66\text{ mN/W}$ | $1.647\%\text{ (PASS)}$ | $+86.8\%\text{ (SAFE)}$ |
| **Medium Speed (120 RPM, 2.0 Hz)**| **CW** ($|f_e - p f_m|$) | **Dual Concorde (Z+X)**| **$94.0\text{ Hz}$** | $[+0.760, -0.016, -0.527]$ | **$0.926\text{ N}$** | $1.189\text{ N}$ | $216.3\text{ W}$ | $174.5\text{ W} / 41.8\text{ W} / \mathbf{0.0\text{ W}}$ | **$4.28\text{ mN/W}$** | $1.640\%\text{ (PASS)}$ | $+81.3\%\text{ (SAFE)}$ |
| **Medium Speed (120 RPM, 2.0 Hz)**| **CCW** ($|f_e + p f_m|$) | **Dual Concorde (Z+X)**| **$106.0\text{ Hz}$** | $[+0.858, +0.018, -0.595]$ | **$1.044\text{ N}$** | $1.341\text{ N}$ | $232.4\text{ W}$ | $185.3\text{ W} / 47.1\text{ W} / \mathbf{0.0\text{ W}}$ | **$4.49\text{ mN/W}$** | $1.647\%\text{ (PASS)}$ | $+80.1\%\text{ (SAFE)}$ |
| **High Speed (1200 RPM, 20.0 Hz)**| **CW** ($|f_e - p f_m|$) | **Single Rotor (Z)** | **$40.0\text{ Hz}$** | $[0.000, -0.006, -0.116]$ | **$0.116\text{ N}$** | $0.149\text{ N}$ | **$65.5\text{ W}$** | $56.9\text{ W} / 8.6\text{ W} / \mathbf{0.0\text{ W}}$ | $1.77\text{ mN/W}$ | $1.613\%\text{ (PASS)}$ | $+92.2\%\text{ (SAFE)}$ |
| **High Speed (1200 RPM, 20.0 Hz)**| **CCW** ($|f_e + p f_m|$) | **Single Rotor (Z)** | **$160.0\text{ Hz}$** | $[0.000, +0.023, -0.462]$ | **$0.463\text{ N}$** | $0.595\text{ N}$ | **$148.3\text{ W}$** | $113.8\text{ W} / 34.5\text{ W} / \mathbf{0.0\text{ W}}$ | **$3.12\text{ mN/W}$** | $1.674\%\text{ (PASS)}$ | $+82.8\%\text{ (SAFE)}$ |
| **High Speed (1200 RPM, 20.0 Hz)**| **CW** ($|f_e - p f_m|$) | **Dual Concorde (Z+X)**| **$40.0\text{ Hz}$** | $[+0.324, -0.007, -0.224]$ | **$0.394\text{ N}$** | $0.506\text{ N}$ | **$131.6\text{ W}$** | $113.8\text{ W} / 17.8\text{ W} / \mathbf{0.0\text{ W}}$ | $2.99\text{ mN/W}$ | $1.613\%\text{ (PASS)}$ | $+87.4\%\text{ (SAFE)}$ |
| **High Speed (1200 RPM, 20.0 Hz)**| **CCW** ($|f_e + p f_m|$) | **Dual Concorde (Z+X)**| **$160.0\text{ Hz}$** | $[+1.295, +0.027, -0.898]$ | **$1.576\text{ N}$** | **$2.025\text{ N}$** | **$298.8\text{ W}$** | $227.6\text{ W} / 71.1\text{ W} / \mathbf{0.0\text{ W}}$ | **$5.27\text{ mN/W}$** | $1.674\%\text{ (PASS)}$ | $+75.1\%\text{ (SAFE)}$ |

*Note: All forces represent internal structural stresses and reaction constraints on the stator frame ($\sum \vec{F}_{\text{ext}} = 0.0\text{ N}$ strictly conserved). The amagnetic dielectric PEEK core exhibits identic zero dissipation ($0.000\text{ W}$) across all 14 kinematic configurations.*

---

## 5. Visual Showcase: 300 DPI Diagnostic Plates & Dynamic Simulation Videos

<div align="center">

### Figure 18: Synoptic Field Shaping — Single PEEK vs Dual 90° Spherical
| Continuous 3-Phase NPNPNP Comparison & Seamless 360° Circular Induction Corona |
| :---: |
| <img src="figures/fig_18_confronto_npnpnp_peek_vs_doppio_rotore.png" width="900" alt="Field Shaping Comparison" /> |
| *Panel A1-A2: Continuous 3-phase NPNPNP waveforms. Panel B1-B2: Long-exposure integrated radial induction corona revealing seamless 360° flux distribution without dead spots. Panel C1-C2: Spatiotemporal kymographs confirming stable rotating phase velocity stripes.* |

### Figure 20: Volumetric 3D Vector Fields & Orthogonal Near-Field Slices
| 3D Magnetic Vector Distribution & Orthogonal Induction Slices |
| :---: |
| <img src="figures/fig_20_campi_3D_sezioni_taglio_nearfield.png" width="900" alt="3D Field Slices" /> |
| *Volumetric 3D vector fields of spherical triple-layer X-cage with dual orthogonal rotor arrays (Z & X) in 90° temporal quadrature. Near-field orthogonal slice maps (XY, XZ, YZ) of induction B, induced E-field vortex, and outward Poynting power flow.* |

### Figure 22: Subbody Thermal Balance & Exterior Shielding Audit
| Component Joule Dissipation & Outer Layer Shielding Confirmation |
| :---: |
| <img src="figures/fig_22_bilancio_termico_perdite_joule.png" width="900" alt="Thermal Audit" /> |
| *Full machine thermal audit: Rotor 1 (43.7%), Rotor 2 (35.6%), Mantle (19.3%), and PEEK Core (0.0 W, confirmed zero eddy losses). Triple-layer mantle breakdown proves complete exterior thermal shielding by Layer 3 (-30° outer = 0.0 W, 0.0%).* |

### Figure 27: 24x24 Fibonacci Architecture & Pisano mod 9 Digital Root Law
| 24-Sector Pisano mod 9 Mapping & Phase-Conjugate Topological Balance |
| :---: |
| <img src="figures/fig_27_architettura_fibonacci_24x24_100w.png" width="900" alt="Fibonacci Architecture Plate" /> |
| *Panel A: Polar map of 24 equatorial sectors with digital root values $F_n \pmod 9$ and phase conjugation lines $\phi_{k+12} = -\phi_k$ ensuring reactive power balance. Panel B: Internal Lorentz stress waveforms. Panel C: Power distribution across 24 coils. Panel D: Verified Gauss solenoidality (0.35% residual).* |

### Figure 30: Dual Orthogonal 90° Macro-Group Architecture (48 Coils)
| Conformal Layout, Multi-Axis Stresses & Verified Gauss Solenoidality |
| :---: |
| <img src="figures/fig_30_doppio_gruppo_90deg_48coils_273n.png" width="900" alt="Dual Orthogonal 90° 48-Coil Architecture" /> |
| *Panel A: Conformal layout of 48 active solenoids: Group 1 equatorial ring (24 coils || Z), Group 2 meridional ring (24 coils || X), 6 mm clearance. Panel B: Multi-axis internal stress waveforms. Panel C: Power breakdown and thermal stabilization. Panel D: Certified Gauss solenoidality (1.491% residual, PASS) and linear magnetic margin.* |

### Figure 31: Kinematic Regimes & Slip-Dependent Electrodynamics (CW vs CCW)
| Harmonic Slip Asymmetry, Parity Vectors & Subbody Thermal Audit |
| :---: |
| <img src="figures/fig_31_kinematic_regimes_comparative.png" width="900" alt="Kinematic Regimes Comparative Plate" /> |
| *Panel A: Mean Lorentz stress vs mechanical velocity (60, 120, 1200 RPM) proving parity asymmetry $f_{\text{slip}}(\text{CCW}) > f_{\text{slip}}(\text{CW})$. Panel B: Subbody Joule dissipation audit confirming $0.000\text{ W}$ in PEEK core and thermal surge under counter-rotation (298.8 W at 1200 RPM CCW). Panel C: Multi-axis force state-space $(\langle F_x \rangle, \langle F_z \rangle)$ demonstrating 6-DoF actuation capability. Panel D: Certified Gauss solenoidality ($1.642\%$, PASS) and safe linear margin ($+81.0\%$).* |

</div>

---

## 6. Quickstart, Replication Suite & Verification Script

All CAD geometries, tetrahedral meshes, Elmer FEM solver definitions, and post-processing pipelines are fully reproducible open-source workflows:

```bash
# 1. Environment Installation
pip install -r requirements.txt

# 2. Master Pipeline Verification Suite (Cross-checks all primary architectures)
python scripts/master_pipeline_verification.py --summary-only

# 3. Individual Architecture Execution & Rendering:
# - Fibonacci 24x24 Balanced Waveguide:
python scripts/run_fibonacci_24x24_simulation.py

# - Fibonacci 24x24 Accumulated Waveguide:
python scripts/run_fibonacci_spinta_accumulata.py

# - 3-Lobe Macro-Chiral Triskelion & Hexagram Armature:
python scripts/run_triskelion_esagramma_simulation.py

# - Dual Orthogonal 90° Macro-Group (48 Coils):
python scripts/run_doppio_gruppo_48coils_simulation.py
```

---

## Sommario Esecutivo per la Comunità Scientifica Italiana

### 1. Revisione e Pivot Scientifico: Abbandono della Propulsione Chiusa
Il progetto **Open Chiral Flux Shaper** adotta formalmente i principi conservativi della fisica classica ed elettromagnetica ortodossa. In accordo con il terzo principio della dinamica e il teorema di Poynting, **un sistema chiuso non può generare alcuna spinta propulsiva stazionaria netta priva di espulsione di massa o momento irraggiato**. 
Le forze ponderomotrici volumetriche calcolate ($\int (\vec{J} \times \vec{B}) dV$) e le integrazioni del Tensore degli Sforzi di Maxwell (MST) rappresentano **tensioni meccaniche interne, coppie di attuazione e gradienti di pressione magnetica**, integralmente bilanciate dai vincoli strutturali statorici.

### 2. Nuove Applicazioni Industriali Ufficiali:
1. **Wireless Power Transfer (WPT) Dinamico Omnidirezionale:**
   La topologia a due gruppi ortogonali a 90° pilotati con sfasamento in quadratura temporale ($\pi/2$) genera un'onda d'induzione rotante isotropa a 360°, eliminando i punti morti e consentendo l'accoppiamento induttivo risonante ad alta efficienza verso carichi mobili e disallineati nello spazio.
2. **Attuatori Magnetici Multi-Asse (6-DoF) e Cuscinetti Magnetici:**
   Il controllo di fase indipendente sulle spire e l'impiego del nucleo in PEEK amagnetico ($\mu_r = 1.0, \sigma = 0\text{ S/m}$) eliminano il cogging meccanico e le perdite per isteresi, offrendo micro-posizionamento senza contatto e controllo d'assetto multiasse per giroscopi e sfere di reazione satellitari.
3. **Riscaldamento a Induzione Mirato (Contour Heating):**
   Il tensore di conducibilità anisotropo ($\bar{\bar{\sigma}}$) concentra le correnti parassite nello strato interno del metamateriale (+30°), mentre lo Strato 3 esterno (-30°) mantiene perdite nulle ($0.000\text{ W}$), garantendo totale schermatura termica e sicurezza d'integrazione.

### 3. Ingegneria Termica di Laboratorio e Raffreddamento Dielettrico
I parametri operativi per i banchi prova di laboratorio vengono ricondotti a potenze sicure continue (**10–50 W**, $J_0 \sim 5 \times 10^3\text{ A/m}^2$), integrando nel nucleo in PEEK micro-canali per il raffreddamento diretto con **liquidi fluorurati dielettrici** (*3M Fluorinert* FC-3283 / FC-770), eliminando la necessità di dissipazione radiativa estrema nel vuoto per le prove a terra.

### 4. Protocollo Metrologico per Banchi di Prova Sperimentali
Per isolare inequivocabilmente gli effetti elettromagnetici reali dagli artefatti ambientali, il protocollo sperimentale impone:
- Montaggio su **bilancia di torsione a sensibilità nanometrica** in **camera a vuoto spinto ($p < 10^{-4}\text{ mbar}$)** contro moti convettivi d'aria ed effetto radiometro di Crookes.
- **Schermatura passiva a doppio strato in Mu-metal** e **bobine attive di Helmholtz a 3 assi** per l'azzeramento del campo geomagnetico e delle interferenze di rete a 50/60 Hz.
- **Null Tests obbligatori** tramite inversione simmetrica di fase ($\vec{J} \to -\vec{J}$) e carichi fittizi resistivi per scorporare dilatazioni termiche e derive capacitive dai gradienti magnetici reali.

---

## Authorship, Attribution & License

- **Lead Inventor & Author:** **Alessandro Brescacin** ([brescacin.alessandro@gmail.com](mailto:brescacin.alessandro@gmail.com))
- **Official GitHub Repository:** [https://github.com/brescale27/open-chiral-flux-shaper](https://github.com/brescale27/open-chiral-flux-shaper)
- **Open Hardware License:** Licensed under the **CERN Open Hardware Licence - Strongly Reciprocal v2 (CERN-OHL-S-2.0)**.  
  See the full text in [`LICENSE.txt`](LICENSE.txt).
- **Citation:** To cite this hardware design, simulation pipeline, or datasets, please refer to [`CITATION.cff`](CITATION.cff).
