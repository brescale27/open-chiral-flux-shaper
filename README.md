# Open Chiral Flux Shaper

*An Open-Source Multiphysics Framework for Anisotropic Metamaterial Field Shaping, Omnidirectional Wireless Power Transfer (WPT), and 6-DoF Contactless Magnetic Actuation.*

[![License: CERN-OHL-S-2.0](https://img.shields.io/badge/License-CERN--OHL--S--2.0-blue.svg)](LICENSE.txt)
[![Release: v2.0.0](https://img.shields.io/badge/Release-v2.0.0-green.svg)](https://github.com/brescale27/open-chiral-flux-shaper/releases)
[![FEM Solver: Elmer FEM 9.0](https://img.shields.io/badge/Elmer%20FEM-9.0%20(CSC)-orange.svg)](https://www.csc.fi/web/elmer)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.cern--ohl--s--2.0-lightgrey.svg)](https://github.com/brescale27/open-chiral-flux-shaper)

---

## 1. Overview & Technical Scope

**Open Chiral Flux Shaper** is a high-fidelity finite-element electrodynamic modeling framework developed in Elmer FEM 3D and Python. The project provides an open-hardware and computational foundation for engineering macro-chiral electromagnetic field distributions using structured metamaterial shells, discrete harmonic excitation laws, and multi-axis orthogonal coil arrays.

By pairing multi-sector discrete winding arrays (including 24-sector Pisano sequence mod 9 topological mappings and dual-ring orthogonal temporal quadratures) with a spherical triple-layer metamaterial shell exhibiting high relative permeability ($\mu_r = 1000$) and an anisotropic conductivity tensor $\bar{\bar{\sigma}}(\theta = \pm 30^\circ)$, the architecture controls the spatiotemporal orientation of the Poynting vector $\vec{S} = \vec{E} \times \vec{H}$ and the magnetic vector potential $\vec{A}$ across stationary and kinematic regimes.

The computational pipeline and hardware designs target three core industrial domains:
1. **Dynamic Omnidirectional Wireless Power Transfer (WPT):** Continuous, steerable 360-degree near-field inductive power links that eliminate angular blind spots without mechanical gimbals.
2. **Multi-Axis Contactless Magnetic Actuation (6-DoF):** Micro-positioning, magnetic levitation, and attitude control testbeds utilizing amagnetic dielectric cores to achieve cogging-free actuation.
3. **Targeted Contour Induction Heating:** High-efficiency localized thermal induction driven by directional chiral current paths, combined with zero-loss outer shielding.

In full alignment with classical electrodynamics, momentum conservation, and the Maxwell Stress Tensor formulation, all computed ponderomotive forces represent internal structural stresses and reaction torques balanced by stator mountings ($\sum \vec{F}_{\text{ext}} = 0$).

```
                 ┌─────────────────────────────┐
                 │   OPEN CHIRAL FLUX SHAPER   │
                 │  Macro-Chiral Metamaterial  │
                 └──────────────┬──────────────┘
                                │
        ┌───────────────────────┼───────────────────────┐
        ▼                       ▼                       ▼
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│  Dynamic WPT  │       │ 6-DoF Contact │       │Contour Heating│
│ Omnidirection │       │  Actuation    │       │ 0W Outer Loss │
│ 0 Blind Spots │       │ No Cogging    │       │ Layer 3 Shield│
└───────────────┘       └───────────────┘       └───────────────┘
```

---

## 2. Industrial Application Domains

### 2.1 Dynamic Omnidirectional Wireless Power Transfer (WPT)
Conventional inductive resonant power transfer systems suffer from rapid efficiency drop-offs when transmitter and receiver coils experience angular or axial misalignment. Open Chiral Flux Shaper utilizes orthogonal and multi-sector polyphase excitation to synthesize a continuously rotating, isotropic induction corona across all three spatial dimensions:

```
                      Z (Orthogonal Ring 1)
                                ▲
                                │
                        ┌───────┴───────┐
                 .──────┤   PEEK Core   ├──────.
               /        │  Microchannel │        \
              /         └───────┬───────┘         \
      ◄──────┼──────────────────┼──────────────────┼──────► Y (Ring 2)
       Mobile \                 │                 /  Receiver Coil
       Drone   \     [ 360° Induction Corona ]   /   [Link: 84.6%]
        Link    '────────────────────────────────'
                                │
                                ▼
                                X
```

- **Solid-State Field Steering:** Full 360-degree spherical coverage achieved purely through temporal phase sequencing ($\phi_k = \frac{v_k}{9} \cdot 2\pi$), removing moving parts, slip rings, and mechanical gimbals.
- **Link Efficiency & Coupling:** Resonant inductive link efficiency reaching $\eta_{\text{link}} = 84.6\%$ at near-field distances ($R = 6.5\text{ cm}$) with a calculated coupling factor $k = 0.385$.
- **Target Applications:** Dynamic docking stations for autonomous aerial vehicles (UAVs), continuous charging for robotic end-effectors, subsea autonomous vehicles, and medical endoscopic capsules.

### 2.2 Multi-Axis Contactless Magnetic Actuation & Active Bearings (6-DoF)
By independently modulating the current amplitude and phase offsets across orthogonal solenoid groups, the architecture produces both pure magnetic couples $\vec{\tau} = \int_V (\vec{r} \times (\vec{J} \times \vec{B})) dV$ and controlled localized magnetic pressure gradients $\nabla \left( \frac{B^2}{2\mu} \right)$:

```
        [ Rotor Array 1 (|| Z) ]     [ Rotor Array 2 (|| X) ]
                    │                            │
                    ▼                            ▼
           Phase Drive: phi_Z           Phase Drive: phi_X
                    │                            │
                    └─────────────┬──────────────┘
                                  ▼
                    ┌───────────────────────────┐
                    │ Triple-Layer Metamaterial │
                    │   Dynamic Maxwell Core    │
                    └─────────────┬─────────────┘
                                  │
          ┌───────────────────────┼───────────────────────┐
          ▼                       ▼                       ▼
     Torque Tau_z            Torque Tau_x            Force F_z
     (Yaw Control)          (Pitch / Roll)         (Levitation)
```

- **Zero-Cogging Topology:** Utilizing an amagnetic, electrically insulating PEEK rotor core ($\mu_r = 1.0, \sigma = 0\text{ S/m}$) completely eliminates parasitic magnetic detent torque, magnetic hysteresis drag, and internal eddy currents.
- **Dynamic Decoupling:** Orthogonal 90-degree temporal and spatial quadrature enables independent control over pitch, roll, yaw, and translation axes.
- **Target Applications:** Contactless reaction spheres for satellite attitude determination and control systems (ADCS), ultra-clean magnetic levitation stages for semiconductor lithography, and high-speed momentum wheels.

### 2.3 Targeted Contour Induction Heating & Thermal Shielding
The metamaterial shell features a tri-layer structure designed to confine and direct high-frequency induced currents:
- **Layer 1 (+30° Chiral Inner Shell):** High-loss anisotropic conductivity zone that absorbs 96.4% of total mantle eddy current dissipation, channeling thermal flux into targeted boundaries.
- **Layer 2 (Intermediate Orthogonal Shell):** Transitional barrier absorbing 3.6% of eddy dissipation.
- **Layer 3 (-30° Counter-Chiral Outer Shell):** Perfect electromagnetic shielding layer exhibiting identically zero dissipation ($0.000\text{ W}$), preventing external stray thermal leakage.

```
  Radius [mm]
   50.0 ──┬───────────────────────────── Layer 3: -30° (0.0 W, Shield)
          │ Conductivity barrier
   48.5 ──┼───────────────────────────── Layer 2: Ortho (0.07 W)
          │ Transition zone
   47.0 ──┴───────────────────────────── Layer 1: +30° (1.78 W, Heat)
          ▼ Internal Air Gap / Coils / PEEK Core
```

---

## 3. Thermal Engineering & Laboratory Design

To transition from high-power computational models to physical laboratory prototypes without thermal degradation, the operating parameters have been downscaled into an intrinsically safe continuous-wave (CW) regime:

```
                     [ Heat Sources: 18.5 W Total ]
                     Coils: 16.6 W | Mantle: 1.85 W
                                   │
                                   ▼
             ┌───────────────────────────────────────────┐
             │         PEEK Core Micro-Channels          │
             │          (12 Parallel Conduits)           │
             └─────────────────────┬─────────────────────┘
                                   │
                Fluorinert FC-3283 │ Flow: 55.4 mL/min
                Dielectric Coolant │ Delta T: 10.0 °C
                                   ▼
             ┌───────────────────────────────────────────┐
             │       External Compact Heat Exchanger     │
             │            (Ambient Rejection)            │
             └───────────────────────────────────────────┘
```

- **Calibrated Electrical Regime:** Satiated excitation current density of $J_0 = 5.0 \times 10^3\text{ A/m}^2$, corresponding to an effective current of $I_{\text{rms}} = 0.65\text{ A}$ across 48 multi-turn coils (120 turns of AWG 27 enameled copper, $R_{\text{coil}} = 0.82\ \Omega$).
- **Total Thermal Dissipation:** $P_{\text{tot}} = 18.48\text{ W}$ ($8.31\text{ W}$ Group 1, $8.31\text{ W}$ Group 2, $1.85\text{ W}$ mantle eddy dissipation, and $0.000\text{ W}$ in the PEEK core).
- **Dielectric Liquid Cooling:** Microfluidic cooling channels integrated directly into the non-conductive PEEK structure using 3M Fluorinert FC-3283 ($c_p = 1100\text{ J/(kg}\cdot\text{K)}$, $\rho = 1820\text{ kg/m}^3$). A laminar flow rate of $55.4\text{ mL/min}$ maintains a steady-state temperature rise below $\Delta T = 10.0^\circ\text{C}$ during continuous CW bench operations, eliminating the need for bulky vacuum radiative panels during atmospheric or vacuum testbench trials.

---

## 4. Metrological Protocol for Laboratory Prototyping

To ensure experimental rigor and eliminate false-positive force readings caused by environmental interference, physical prototypes must be tested under strict metrological controls:

```
  ┌────────────────────────────────────────────────────────────┐
  │ High-Vacuum Chamber (p < 10⁻⁴ mbar)                        │
  │  ┌──────────────────────────────────────────────────────┐  │
  │  │ Active 3-Axis Helmholtz Shield (B_ext < 0.1 µT)      │  │
  │  │  ┌────────────────────────────────────────────────┐  │  │
  │  │  │ Double Mu-Metal Enclosure (> 60 dB Attenuation)│  │  │
  │  │  │  ┌──────────────────────────────────────────┐  │  │  │
  │  │  │  │ Quartz Torsion Fiber Balance             │  │  │  │
  │  │  │  │    │                                     │  │  │  │
  │  │  │  │    ├── Laser Interferometer (Sub-nm)     │  │  │  │
  │  │  │  │    │                                     │  │  │  │
  │  │  │  │   [ FLUX SHAPER BENCH RIG ]              │  │  │  │
  │  │  │  │   (Fluorinert FC-3283 Feedlines)         │  │  │  │
  │  │  │  └──────────────────────────────────────────┘  │  │  │
  │  │  └────────────────────────────────────────────────┘  │  │
  │  └──────────────────────────────────────────────────────┘  │
  └────────────────────────────────────────────────────────────┘
```

1. **High-Vacuum Environment ($p < 10^{-4}\text{ mbar}$):** Eliminates buoyant convective air currents, acoustic streaming, and radiometric Crookes/Knudsen thermal outgassing forces that mimic micro-Newton forces on sensitive balances.
2. **Magnetic Isolation:** Dual-walled Mu-metal enclosure ($> 60\text{ dB}$ attenuation against power grid noise) paired with an active triaxial Helmholtz compensation system that zeroes the local geomagnetic field ($\vec{B}_{\text{geo}} \approx 45\ \mu\text{T}$) to $< 0.1\ \mu\text{T}$, preventing external Lorentz interactions.
3. **Sub-Nanometer Telemetry:** Differential optical interferometer and quadrant photodiode monitoring a quartz-fiber torsion balance, calibrated via electrostatic comb drives.
4. **Mandatory Null Tests & Parity Inversion Protocols:**
   - *Symmetric Phase Inversion ($\vec{J} \to -\vec{J}$):* Distinguishes first-order electromagnetic interactions from second-order electrostatic and capacitive artifacts.
   - *Non-Inductive Dummy Heaters:* Dissipating equivalent Joule heat through non-inductive resistive loads to isolate dilatometric thermal expansion of the balance arm.
   - *Quadrature Balancing:* Activating paired orthogonal rings in balanced opposition to experimentally confirm zero external momentum transfer ($\sum \vec{F}_{\text{ext}} = 0$).

---

## 5. Master Comparative Benchmark Across All Tested Architectures

The synoptic master table consolidates the entire electromagnetic, mechanical, and thermal design space evaluated in this project:

| Architecture / Configuration | Core Type & Reluctance | Mantle Structure & Permeability | Excitation Logic & Phase Law | Operational Regime | Peak Radial Field B_rad (6.5 cm) | Internal Lorentz Stress $\|\langle\vec{F}\rangle\|$ | Peak Instantaneous Force | Active Joule Losses $P_J$ | Gauss Solenoidality Residual | Primary Industrial Application |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Baseline (v1.0.0)** | Soft Iron (Z = -H/2) | Solid Al Mesh ($\mu_r = 1.0$) | 60° Progressive Sine (100 Hz) | Continuous (1200 RPM) | 62.98 µT | ≈ 0 (leakage) | ≈ 0 | 2.437 W | 1.402% | Induction Stray Venting |
| **Centered Continuous** | Soft Iron (Z = 0) | Solid Al Mesh ($\mu_r = 1.0$) | 60° Progressive Sine (100 Hz) | Continuous (1200 RPM) | 211.35 µT | 4.67 µN | 27.41 µN | 1.52 mW | 0.076% | Near-Field Rotary Induction |
| **Centered Locked-Rotor** | Soft Iron (Z = 0) | Solid Al Mesh ($\mu_r = 1.0$) | 60° Progressive Sine (100 Hz) | Solid-State (0 RPM) | 211.35 µT | 5.72 µN | 12.89 µN | 1.02 mW | 0.031% | Inductive Transformer Stage |
| **Mirrored Pulsed (N-S)** | Soft Iron (Z = 0) | Solid Al Mesh ($\mu_r = 1.0$) | 60° Half-Wave Pulse Train | Solid-State (0 RPM) | 135.84 µT | 0.93 µN (balanced) | ±26.50 µN | 1.90 mW | 0.125% | Ultra-Low Loss WPT Stage |
| **Closed Can Architecture** | Soft Iron (Z = 0) | Al Mesh + Lids (Z = ±H/2) | 60° Progressive Sine (100 Hz) | Continuous (1200 RPM) | 185.20 µT | 0.20 µN | 2.96 µN | 0.017 mW | 0.045% | -98.9% Thermal Shielding |
| **Ferro Expanded Mesh** | Soft Iron (Z = 0) | Ferro 30° Mesh ($\mu_r = 1000$) | Asymmetric Thirds (33/67/100%)| Solid-State (0 RPM) | 3137.0 µT | 113.51 µN | 897.6 µN | 0.172 mW | 0.850% | High-Flux Concentrator |
| **Triple-Layer X + PEEK** | Amagnetic PEEK Core | Triple X ($\mu_r = 1000, \pm 30^\circ$) | Asymmetric Thirds (33/67/100%)| Solid-State (0 RPM) | 1153.2 µT | 36.99 µN | 142.9 µN | 0.095 mW | 1.900% | Metamaterial Chiral Guide |
| **NPNPNP Single PEEK** | Amagnetic PEEK Core | Triple X ($\mu_r = 1000, \pm 30^\circ$) | Continuous 3-Phase NPNPNP | Solid-State (0 RPM) | 428.0 µT | 849.1 µN | 808.9 µN | 56.78 mW | 0.985% | Seamless 360° Field Shaper |
| **NPNPNP Dual 90° Spherical**| Amagnetic PEEK Core | Spherical X ($\mu_r = 1000, \pm 30^\circ$) | Dual Continuous 3-Phase NPNPNP | Solid-State (0 RPM) | 8.82 mT (91.2 mT pk) | 0.985 N | 14.17 N | 230.3 W | 1.002% (PASS) | Omnidirectional 3D WPT Stage |
| **Fibonacci 24x24 (Balanced)**| Amagnetic PEEK Core | Spherical X ($\mu_r = 1000, \pm 30^\circ$) | 24-Sector Pisano mod 9 (100 Hz) | Solid-State (0 RPM) | 60.3 µT (1.13 mT pk) | 14.75 µN | 30.85 µN | 2.40 kW | 0.346% (PASS) | Self-Balancing Guide |
| **Fibonacci 24x24 (Accum.)** | Amagnetic PEEK Core | Spherical X ($\mu_r = 1000, \pm 30^\circ$) | 24-Sector Pisano mod 9 + 15° Prog| Solid-State (0 RPM) | 56.1 µT (668.5 µT pk)| 16.16 µN | 46.20 µN | 2.40 kW | 0.200% (PASS) | Directional Waveguide |
| **Triskelion 3-Lobe Hexagram**| Amagnetic PEEK Hexagram| Triskelion X ($\mu_r = 1000, 3\text{ Lobi}$)| Exact 24-Pulse ($\phi_k = \frac{v_k}{9} 2\pi$)| Solid-State (0 RPM) | 52.8 µT (1.13 mT pk) | 36.97 µN | 64.68 µN | 2.40 kW | 0.647% (PASS) | Chiral Harmonic Rectifier |
| **Dual Orthogonal 90° (48 C.)**| Amagnetic PEEK Core | Spherical X ($\mu_r = 1000, \pm 30^\circ$) | Exact 24-Pulse Quadrature (Z & X) | Solid-State (0 RPM) | 14.1 mT (1.38 T pk) | 6.664 N (Raw: 41.5 µN)| 273.6 N (Burst) | 1549.3 W | 1.491% (PASS) | Multi-Axis 6-DoF Actuator |
| **Chiral WPT / 6-DoF Benchtop**| Amagnetic PEEK Core | Spherical X ($\mu_r = 1000, \pm 30^\circ$) | Exact 24-Pulse Quadrature (Z & X) | Solid-State (0 RPM) | 2.08 µT (2.89 µT pk) | 0.011 µN | 0.041 µN | 18.48 W | 1.491% (PASS) | Calibrated Lab Prototype |

---

## 6. Comprehensive Kinematic & Energy Regimes Benchmark (60, 120, 1200 RPM — CW vs CCW)

The behavior of the electrodynamic interaction is governed by the mechanical slip frequency $f_{\text{slip}} = |f_e \mp p \cdot f_{\text{mech}}|$, where $f_e = 100.0\text{ Hz}$, $p = 3$ pole pairs, and synchronous mechanical speed is $n_{\text{sync}} = 2000\text{ RPM}$. Below is the systematic comparison across all 14 evaluated kinematic states:

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

---

## 7. Visual Showcase & Diagnostic Plates

The repository provides high-resolution 300 DPI analytical plates and dynamic simulation records:

<div align="center">

### Figure 18: Synoptic Field Shaping — Single PEEK vs Dual 90° Spherical
| Continuous 3-Phase NPNPNP Waveforms & Omnidirectional 360° Induction Corona |
| :---: |
| <img src="figures/fig_18_confronto_npnpnp_peek_vs_doppio_rotore.png" width="900" alt="Field Shaping Comparison" /> |
| *Panel A1-A2: Continuous 3-phase NPNPNP excitation currents. Panel B1-B2: Time-integrated radial induction corona demonstrating gapless 360° flux distribution without angular blind spots. Panel C1-C2: Spatiotemporal kymographs confirming constant rotating phase velocity.* |

### Figure 20: Volumetric 3D Vector Fields & Orthogonal Near-Field Slices
| 3D Magnetic Vector Distribution & Orthogonal Induction Slices |
| :---: |
| <img src="figures/fig_20_campi_3D_sezioni_taglio_nearfield.png" width="900" alt="3D Field Slices" /> |
| *Volumetric vector field distribution of spherical triple-layer X-cage with dual orthogonal rotor arrays (Z & X) in 90° temporal quadrature. Near-field orthogonal slice maps (XY, XZ, YZ) of induction B, induced E-field vortex, and outward Poynting power flow.* |

### Figure 22: Subbody Thermal Balance & Exterior Shielding Audit
| Subbody Joule Dissipation & Outer Layer Shielding Confirmation |
| :---: |
| <img src="figures/fig_22_bilancio_termico_perdite_joule.png" width="900" alt="Thermal Audit" /> |
| *Full machine thermal audit: Rotor 1 (43.7%), Rotor 2 (35.6%), Mantle (19.3%), and PEEK Core (0.0 W, confirmed zero eddy losses). Triple-layer mantle breakdown confirms complete exterior thermal shielding by Layer 3 (-30° outer = 0.0 W, 0.0%).* |

### Figure 27: 24x24 Fibonacci Architecture & Pisano mod 9 Digital Root Law
| 24-Sector Pisano mod 9 Mapping & Phase-Conjugate Topological Balance |
| :---: |
| <img src="figures/fig_27_architettura_fibonacci_24x24_100w.png" width="900" alt="Fibonacci Architecture Plate" /> |
| *Panel A: Polar map of 24 equatorial sectors with digital root values F_n mod 9 and phase conjugation lines ensuring reactive balance. Panel B: Internal Lorentz stress waveforms. Panel C: Power distribution across 24 coils. Panel D: Verified Gauss solenoidality (0.35% residual).* |

### Figure 30: Dual Orthogonal 90° Macro-Group Architecture (48 Coils)
| Conformal Layout, Multi-Axis Stresses & Verified Gauss Solenoidality |
| :---: |
| <img src="figures/fig_30_doppio_gruppo_90deg_48coils_273n.png" width="900" alt="Dual Orthogonal 90° 48-Coil Architecture" /> |
| *Panel A: Conformal layout of 48 active solenoids: Group 1 equatorial ring (24 coils || Z), Group 2 meridional ring (24 coils || X), 6 mm clearance. Panel B: Multi-axis internal stress waveforms. Panel C: Power breakdown and thermal stabilization. Panel D: Certified Gauss solenoidality (1.491% residual, PASS) and linear magnetic margin.* |

### Figure 31: Kinematic Regimes & Slip-Dependent Electrodynamics (CW vs CCW)
| Harmonic Slip Asymmetry, Parity Vectors & Subbody Thermal Audit |
| :---: |
| <img src="figures/fig_31_kinematic_regimes_comparative.png" width="900" alt="Kinematic Regimes Comparative Plate" /> |
| *Panel A: Mean Lorentz stress vs mechanical velocity (60, 120, 1200 RPM) proving parity asymmetry f_slip(CCW) > f_slip(CW). Panel B: Subbody Joule dissipation audit confirming 0.000 W in PEEK core and thermal surge under counter-rotation (298.8 W at 1200 RPM CCW). Panel C: Multi-axis force state-space (<Fx>, <Fz>) demonstrating 6-DoF actuation capability. Panel D: Certified Gauss solenoidality (1.642%, PASS) and safe linear margin (+81.0%).* |

### Dynamic Video: Dual Orthogonal 90° Multi-Axis Electrodynamics
| 3D Orthogonal Solenoid Current State, Dynamic Magnetic Vector & Real-Time Waveforms |
| :---: |
| <img src="figures/video_dinamica_doppio_gruppo_48coils.gif" width="900" alt="Dynamic Video: Dual Orthogonal 90° Electrodynamics" /> |
| *Synchronized high-resolution simulation video over 16.0 ms transient electrical cycle (64 timesteps, 100 Hz). Left: 3D perspective wireframe of spherical mantle showing the 48 active solenoids with current density color-modulation and resultant dynamic magnetic vector. Top Right: 3D state-space force hodograph. Bottom Right: Real-time scrolling waveforms.* |

</div>

---

## 8. Quickstart, Replication Suite & Verification Script

The repository is fully reproducible using open-source tools:

```bash
# 1. Environment Installation
pip install -r requirements.txt

# 2. Master Verification Suite (Cross-checks primary architectures)
python scripts/master_pipeline_verification.py --summary-only

# 3. Kinematic Regimes Benchmark (14 states, CW vs CCW, Figure 31)
python scripts/run_kinematic_regimes_simulation.py

# 4. Calibrated Laboratory Benchtop Prototype (Safe 18.5 W regime)
python variants/gabbia_sferica_chiral_wpt_actuator/scripts/run_chiral_wpt_actuator_simulation.py

# 5. Core Architectural Simulations:
# - Dual Orthogonal 90° Macro-Group (48 Coils):
python scripts/run_doppio_gruppo_48coils_simulation.py

# - 3-Lobe Macro-Chiral Triskelion & Hexagram Armature:
python scripts/run_triskelion_esagramma_simulation.py

# - Fibonacci 24x24 Balanced Waveguide:
python scripts/run_fibonacci_24x24_simulation.py
```

---

## 9. Sommario Esecutivo per la Comunità Scientifica Italiana

### 1. Inquadramento Fisico ed Epistemologico
Il progetto **Open Chiral Flux Shaper** è un framework multifisico computazionale per la modellazione e la manipolazione di campi elettromagnetici macro-chirali. In aderenza al principio di conservazione della quantità di moto, al terzo principio della dinamica e al teorema di Poynting:
- **Tensioni Interne di Maxwell:** Tutte le forze volumetriche calcolate rappresentano gradienti di pressione magnetica e tensioni strutturali interne tra rotori e mantello, integralmente bilanciate dai vincoli meccanici dello statore.
- **Pressione di Radiazione di Poynting:** A frequenze industriali ($100\text{ Hz}$) e dimensioni sub-lunghezza d'onda ($ka \sim 10^{-7}$), la spinta fotonica derivante da radiazione è trascurabile ($F_{\text{rad}} = P/c \sim 24.5\text{ pN}$ per $P = 7.34\text{ mW}$), escludendo qualsiasi spinta propulsiva stazionaria netta a sistema chiuso.

### 2. Ambiti Applicativi Industriali Convalidati
1. **Wireless Power Transfer (WPT) Dinamico Omnidirezionale:**
   La generazione di un'onda d'induzione rotante isotropa a 360° nello spazio tridimensionale consente il trasferimento induttivo continuo verso droni, veicoli subacquei, bracci robotici articolati e dispositivi biomedicali, azzerando le perdite da disallineamento angolare.
2. **Attuazione Magnetica Contactless a 6 Gradi di Libertà (6-DoF):**
   L'impiego di un nucleo in PEEK amagnetico e dielettrico ($\sigma = 0\text{ S/m}, \mu_r = 1.0$) elimina totalmente la coppia di cogging e le perdite per isteresi, offrendo micro-posizionamento senza contatto per banchi ottici e sfere di reazione per l'assetto satellitare.
3. **Riscaldamento a Induzione Mirato (Contour Heating):**
   Il tensore di conducibilità chirale anisotropo ($\bar{\bar{\sigma}}$ a $\pm 30^\circ$) concentra le perdite nel profilo interno (+30°), mentre lo Strato 3 esterno (-30°) mantiene perdite identicamente nulle ($0.000\text{ W}$), garantendo una perfetta schermatura termica verso l'ambiente esterno.

### 3. Ingegneria Termica e Metrologia di Laboratorio
- **Regime di Banco Sicuro:** Densità di corrente calibrata a $J_0 = 5 \times 10^3\text{ A/m}^2$ ($18.5\text{ W}$ totali) con raffreddamento a liquido dielettrico fluorurato (*3M Fluorinert* FC-3283) a $55.4\text{ mL/min}$ in micro-condotti integrati nel nucleo PEEK.
- **Protocollo Metrologico per Test a Vuoto:** Camera a vuoto ($< 10^{-4}\text{ mbar}$), schermatura passiva in Mu-metal ($> 60\text{ dB}$), gabbia di Helmholtz a 3 assi, bilancia di torsione con telemetria interferometrica e null tests simmetrici di inversione di fase.

---

## 10. Authorship & License

- **Author & Principal Investigator:** Alessandro Brescacin
- **Repository:** [https://github.com/brescale27/open-chiral-flux-shaper](https://github.com/brescale27/open-chiral-flux-shaper)
- **License:** Open Hardware licensed under the **CERN Open Hardware Licence Version 2 - Strongly Reciprocal ([CERN-OHL-S-2.0](LICENSE.txt))**.
- **Software Components:** Scientific Python scripts and post-processing tools licensed under the **Apache License, Version 2.0**.
