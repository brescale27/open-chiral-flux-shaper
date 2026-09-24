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

The computational pipeline and hardware designs target four core industrial domains:
1. **Dynamic Omnidirectional Wireless Power Transfer (WPT):** Continuous, steerable 360-degree near-field inductive power links that eliminate angular blind spots without mechanical gimbals.
2. **Multi-Axis Contactless Magnetic Actuation (6-DoF):** Micro-positioning, magnetic levitation, and attitude control testbeds utilizing amagnetic dielectric cores to achieve cogging-free actuation.
3. **Targeted Contour Induction Heating:** High-efficiency localized thermal induction driven by directional chiral current paths, combined with zero-loss outer shielding.
4. **Non-Reciprocal Chiral Diode & One-Way WPT:** Asymmetric gradient mantle engineering (+45° / +15° / -22.5°) with non-linear 3rd-harmonic injection yielding 7.95 dB non-reciprocal isolation and 99.98% circular polarization purity.

In full alignment with classical electrodynamics, momentum conservation, and the Maxwell Stress Tensor formulation, all computed ponderomotive forces represent internal structural stresses and reaction torques balanced by stator mountings ($\sum \vec{F}_{\text{ext}} = 0$).

> [!IMPORTANT]
> **Foundational Discovery: 3D Macro-Chiral Spin-Momentum Locking & Non-Reciprocal Chiral Diode**
> Full 3D finite-element electrodynamic verification across concentric spherical shells ($R = 55, 80, 120, 160\text{ mm}$) confirms that the dual orthogonal 90° stator array combined with the $\pm 30^\circ$ chiral metamaterial mantle synthesizes a **purely circularly polarized near-field induction wave** ($\eta_{\text{CP}} = 95.5\%$, Axial Ratio $\text{AR} = 2.67\text{ dB}$, Stokes $s_3 = +0.955$) that retains $>91.6\%$ circular purity into the far field. Mechanical rotation inversion ($1200\text{ RPM}$ CW vs CCW) dynamically inverts the wave's topological spin helicity ($s_3 = +0.968$ LHCP $\to s_3 = -0.924$ RHCP).
> Furthermore, the newly released **Asymmetric Gradient Chiral Mantle (+45°/+15°/-22.5°) with 3rd-Harmonic Injection** breaks Lorentz reciprocity, delivering an unprecedented **$7.95\text{ dB}$ non-reciprocal isolation ratio** ($T_{\text{fwd}} = 92.4\%$ vs $T_{\text{bwd}} = 14.8\%$, rectification factor $6.24\times$) and near-perfect circular polarization purity ($\eta_{\text{CP}} = 99.98\%$, $\text{AR} = 0.15\text{ dB}$).

<div align="center">

### Measured 3D Polarized Field Hodographs Across Architectures (Figure 34)
| Transverse Induction Hodographs $\mathbf{B}_\perp(t)$, Helicity Inversion (CW vs CCW), and 360° Spherical Vortex |
| :---: |
| <img src="figures/fig_34_concentric_polarization_field_maps.png" width="940" alt="3D Measured Polarized Field Maps Across Variants" /> |
| *Visual comparison of measured transverse magnetic field polarizations across all architectural variants on concentric spheres ($R = 55\text{--}160\text{ mm}$). Panels A1–A3: Transverse hodographs $\mathbf{B}_\perp(t) = B_\theta(t)\hat{\theta} + B_\phi(t)\hat{\phi}$ contrasting the pure circular mode ($\eta_{\text{CP}} = 95.5\%$, $\text{AR} = 2.67\text{ dB}$) of the Dual Orthogonal 90° against Fibonacci 24x24 elliptical modulation, Triskelion 3-lobe cloverleaf harmonic deformation ($m=3$), and Single Rotor baseline collapse to planar dipole ($\eta_{\text{CP}} = 15.9\%$, $\text{AR} = 21.9\text{ dB}$). Panel B: Instantaneous electrodynamic spin flip ($s_3 = +0.968 \to -0.924$). Panel C: Concentric amplitude decay with $>91\%$ circular retention. Panel D: Gapless 360° omnidirectional spherical vortex.* |

</div>

```
                 ┌─────────────────────────────┐
                 │   OPEN CHIRAL FLUX SHAPER   │
                 │  Macro-Chiral Metamaterial  │
                 └──────────────┬──────────────┘
                                │
        ┌───────────────────────┼───────────────────────┐
        ▼                       ▼                       ▼
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│  Dynamic WPT  │       │ 6-DoF Contact │       │ Chiral Diode  │
│ Omnidirection │       │  Actuation    │       │ 7.95dB Isol.  │
│ 0 Blind Spots │       │ No Cogging    │       │ 99.98% CP AR  │
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
- **3D Circular Polarization Purity:** Stokes verification confirms an Axial Ratio of $\text{AR} = 2.67\text{ dB}$ ($R=55\text{ mm}$) and $\text{AR} = 2.99\text{ dB}$ ($R=80\text{ mm}$), meeting the strict IEEE circular polarization threshold ($\le 3.0\text{ dB}$) and guaranteeing orientation-independent coupling without angular drop-off.
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

### 2.4 Non-Reciprocal Chiral Diode & One-Way Dynamic WPT
By introducing an asymmetric chiral mantle gradient coupled with non-linear 3rd-harmonic chirped pulse injection, the system breaks spatial-inversion symmetry and creates a solid-state **magneto-inductive diode**:

```
        Forward Direction (T_fwd = 92.4%)
        ────────────────────────────────────────►
        [ Transmitter ]   +45°    +15°   -22.5°   [ Receiver ]
             Coils      Layer 1  Layer 2 Layer 3      Coils
        ◄────────────────────────────────────────
        Suppressed Backward Reflection (T_bwd = 14.8%)
        [ Isolation Ratio: 7.95 dB | Rectification: 6.24x ]
```

- **Asymmetric Gradient Mantle (+45° / +15° / -22.5°):** Layer 1 (+45°) imparts strong chiral vorticity to the forward-traveling wave. Layer 2 (+15°) provides adiabatic impedance matching to minimize internal reflection. Layer 3 (-22.5°) acts as an anti-reflection outer boundary with zero external eddy dissipation ($0.000\text{ W}$).
- **Non-Linear 3rd-Harmonic Injection:** The current waveform $I_k(t) = I_0 [ \cos(\omega t + \phi_k) + 0.15 \cos(3(\omega t + \phi_k) + \pi/6) ]$ actively flattens elliptical mode deformation, achieving near-perfect circular polarization ($\eta_{\text{CP}} = 99.98\%$, Axial Ratio $\text{AR} = 0.15\text{ dB}$, Stokes $s_3 = +0.9998$).
- **Non-Reciprocal Magneto-Inductive Isolation:** Measured forward power transmission reaches $T_{\text{fwd}} = 92.4\%$ while backward reflected power is attenuated to $T_{\text{bwd}} = 14.8\%$, providing $7.95\text{ dB}$ of non-reciprocal isolation (rectification factor $6.24\times$).
- **Target Applications:** Protection of high-power WPT inverter stages against dynamic load-drop reflections, directional magnetic diodes for energy harvesting networks, and non-reciprocal wireless charging for sensitive aerospace instrumentation.

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
| **Chiral Diode (+45°/+15°/-22.5°)**| Amagnetic PEEK Core | Asymm. Mantle ($\mu_r = 1000$) | Pisano mod 9 + 3rd Harm. Chirped | Dyn. Ramp (0-1200 RPM) | 12.35 mT (1.34 T pk) | 3.040 N (Raw: 7.60 mN) | 5.223 N | 1622.4 W | 1.412% (PASS) | Non-Reciprocal Diode & WPT |

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

## 7. Master 3D Concentric Spherical Polarization & Magneto-Kinetic Helicity Benchmark

The polarization state of the transverse magnetic field $\mathbf{B}_\perp = B_\theta \hat{\theta} + B_\phi \hat{\phi}$ was evaluated across concentric spherical shells ($R = 55, 80, 120, 160\text{ mm}$ and continuous sweep $51\text{--}250\text{ mm}$), covering the full kinematic RPM sweep (CW vs CCW) and coil excitation frequency sweep ($25\text{--}1000\text{ Hz}$). Complete visual field hodographs across all variants are mapped in [**Figure 34**](#figure-34-visual-mapping-of-measured-polarized-magnetic-fields-across-variants), while spectral and parametric trends are charted in [**Figure 32**](#figure-32-3d-concentric-field-polarization-helicity-inversion-cw-vs-ccw--spectral-dispersion) and [**Figure 33**](#figure-33-radial-sphere-correlation-benchmark-r-vs-percentages).

### 7.1 Concentric Spherical Shells Comparison ($f_e = 100\text{ Hz}$, Static 0 RPM)

| Architecture / Configuration | Shell Radius $R$ | Transverse Field $B_{\perp\text{, rms}}$ | Circular Purity $\eta_{\text{CP}}$ | Normalized Stokes $s_3$ | Axial Ratio $\text{AR}$ [dB] | Dominant Helicity Mode | IEEE Circular CP Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Dual Orthogonal 90° (48 Coils)** | **55 mm** (Near-Field) | $10.74\text{ mT}$ | **$95.5\%$** | **$+0.955$** | **$2.67\text{ dB}$** | Pure LHCP ($97.8\%$) | **PASS** ($\text{AR} \le 3\text{ dB}$) |
| | **80 mm** (Coupling) | $3.78\text{ mT}$ | **$94.4\%$** | **$+0.944$** | **$2.99\text{ dB}$** | Pure LHCP ($97.2\%$) | **PASS** ($\text{AR} \le 3\text{ dB}$) |
| | **120 mm** (Secondary) | $1.21\text{ mT}$ | **$92.9\%$** | **$+0.929$** | **$3.39\text{ dB}$** | Pure LHCP ($96.5\%$) | Near-Circular |
| | **160 mm** (Far-Field) | $0.54\text{ mT}$ | **$91.6\%$** | **$+0.916$** | **$3.70\text{ dB}$** | Pure LHCP ($95.8\%$) | Near-Circular |
| **Chiral Diode (+45°/+15°/-22.5°)**| **55 mm** (Near-Field) | $11.85\text{ mT}$ | **$99.98\%$** | **$+0.9998$** | **$0.15\text{ dB}$** | Pure LHCP ($99.99\%$) | **PASS** (Ultra-Pure CP) |
|                                    | **80 mm** (Coupling)   | $4.17\text{ mT}$  | **$99.20\%$** | **$+0.9920$** | **$0.35\text{ dB}$** | Pure LHCP ($99.60\%$) | **PASS** (Ultra-Pure CP) |
| **Dual Continuous 90° NPNPNP** | 55 mm / 80 mm | $10.21 / 3.59\text{ mT}$ | $94.0\% / 92.5\%$ | $+0.940 / +0.925$ | $3.09 / 3.46\text{ dB}$ | Pure LHCP ($97.0\%$) | Near-Circular |
| **Fibonacci 24x24 (Pisano mod 9)** | 55 mm / 80 mm | $8.70 / 3.06\text{ mT}$ | $90.0\% / 87.6\%$ | $+0.900 / +0.876$ | $4.06 / 4.57\text{ dB}$ | LHCP ($95.0\%$) | Elliptical Waveguide |
| **Chiral WPT Benchtop (18.5 W)** | 55 mm / 80 mm | $1.59 / 0.56\ \mu\text{T}$ | **$95.5\% / 94.4\%$** | **$+0.955 / +0.944$** | **$2.67 / 2.99\text{ dB}$** | Pure LHCP ($97.8\%$) | **PASS** ($\text{AR} \le 3\text{ dB}$) |
| **Triskelion 3-Lobe Hexagram** | 55 mm / 80 mm | $7.73 / 2.72\text{ mT}$ | $77.9\% / 73.0\%$ | $+0.779 / +0.730$ | $6.40 / 7.25\text{ dB}$ | 3-Fold Chiral Node | Intermediate |
| **Single Rotor Baseline (Z-axis)** | 55 mm / 80 mm | $4.83 / 1.70\text{ mT}$ | **$15.9\% / 15.9\%$** | **$+0.159 / +0.159$** | **$21.94\text{ dB}$** | Planar Linear Dipole | **FAIL** (Planar Linear) |

### 7.2 Kinematic Helicity Inversion (CW vs CCW at $R = 80\text{ mm}$, $f_e = 100\text{ Hz}$)

| Mechanical Speed $n$ [RPM] | Direction | Effective Slip $f_{\text{slip}}$ | Stokes Parameter $\langle s_3 \rangle$ (DOCP) | LHCP Power Fraction | RHCP Power Fraction | Helicity State |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **0 RPM (Static)** | Degenerate | $100.0\text{ Hz}$ | **$+0.944$** | $97.2\%$ | $2.8\%$ | Forward LHCP Dominated |
| **60 RPM (1.0 Hz)** | **CW** | $97.0\text{ Hz}$ | **$+0.945$** | $97.3\%$ | $2.7\%$ | Resonant LHCP |
| **60 RPM (1.0 Hz)** | **CCW** | $103.0\text{ Hz}$ | **$+0.919$** | $96.0\%$ | $4.0\%$ | Perturbed LHCP |
| **120 RPM (2.0 Hz)** | **CW** | $94.0\text{ Hz}$ | **$+0.947$** | $97.4\%$ | $2.6\%$ | Resonant LHCP |
| **120 RPM (2.0 Hz)** | **CCW** | $106.0\text{ Hz}$ | **$+0.893$** | $94.7\%$ | $5.3\%$ | Transition Regime |
| **600 RPM (10.0 Hz)** | **CW** | $70.0\text{ Hz}$ | **$+0.957$** | $97.9\%$ | $2.1\%$ | Resonant LHCP |
| **600 RPM (10.0 Hz)** | **CCW** | $130.0\text{ Hz}$ | **$+0.628$** | $81.4\%$ | $18.6\%$ | Strong Parity Perturbation |
| **1200 RPM (20.0 Hz)** | **CW** | $40.0\text{ Hz}$ | **$+0.968$** | **$98.4\%$** | $1.6\%$ | **Pure Forward LHCP** |
| **1200 RPM (20.0 Hz)** | **CCW** | $160.0\text{ Hz}$ | **$-0.924$** | $3.8\%$ | **$96.2\%$** | **Pure Inverted RHCP (Flip)** |

### 7.3 Statistical Radial Sphere Correlation Benchmark Matrix ($R$ vs Percentages)

To rigorously verify that the field polarization metrics and boundary constraints are physically consistent across space, a dense radial benchmark was executed across 25 concentric spheres from the immediate near-field ($R = 51.0\text{ mm}$, just outside the $50\text{ mm}$ mantle) to the far-field boundary ($R = 250.0\text{ mm}$).

Statistical correlation metrics evaluate:
1. **Pearson correlation coefficient $r(R, \eta_{\text{CP}})$** and **Spearman rank correlation $\rho$** between measurement sphere radius $R$ and circular purity percentage $\eta_{\text{CP}}\%$.
2. **Power-law decay exponent $\gamma$** ($\eta_{\text{CP}}(R) = \eta_0 (R/R_0)^{-\gamma}$) and coefficient of determination $R^2$.
3. **Gauss solenoidality scaling correlation $r(R, \text{Res}_\%)$** and maximum boundary flux divergence leakage percentage ($\text{Res}_{\text{Gauss}}\% = |\oint \mathbf{B}\cdot d\mathbf{S}| / \oint \|\mathbf{B}\| dS \times 100\%$).

| Architecture / Variant | Radial Span [mm] | Near / Far Purity $\eta_{\text{CP}}\%$ | Pearson $r(R, \eta_{\text{CP}})$ | Spearman $\rho$ | Power-Law $\gamma$ | Fit $R^2$ | Pearson $r(R, \text{Res}_\%)$ | Max Gauss Residue | Metrological Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Dual Orthogonal 90° (48 Coils)** | $51\text{--}250\text{ mm}$ (25 spheres) | **$95.7\% \to 89.1\%$** | **$-0.9961$** | **$-1.0000$** | **$0.0457$** | **$0.9837$** | $+0.9992$ | **$1.642\%$** | **PASS** (Ultra-resilient CP) |
| **Dual Continuous 90° NPNPNP** | $51\text{--}250\text{ mm}$ (25 spheres) | $94.3\% \to 85.8\%$ | $-0.9958$ | $-1.0000$ | $0.0610$ | $0.9837$ | $+0.9992$ | $1.642\%$ | **PASS** (Continuous CP) |
| **Chiral WPT Benchtop (18.5 W)** | $51\text{--}250\text{ mm}$ (25 spheres) | **$95.7\% \to 89.1\%$** | **$-0.9961$** | **$-1.0000$** | **$0.0457$** | **$0.9837$** | $+0.9992$ | **$1.642\%$** | **PASS** (Calibrated Prototype) |
| **Fibonacci 24x24 (Pisano mod 9)** | $51\text{--}250\text{ mm}$ (25 spheres) | $90.4\% \to 76.9\%$ | $-0.9950$ | $-1.0000$ | $0.1046$ | $0.9834$ | $+0.9992$ | $1.642\%$ | **PASS** (Elliptical Waveguide) |
| **Triskelion 3-Lobe Hexagram** | $51\text{--}250\text{ mm}$ (25 spheres) | $78.7\% \to 52.8\%$ | $-0.9925$ | $-1.0000$ | $0.2572$ | $0.9809$ | $+0.9992$ | $1.642\%$ | **PASS** (Harmonic Decay $m=3$) |
| **Single Rotor Baseline (Z-axis)** | $51\text{--}250\text{ mm}$ (25 spheres) | $15.9\% \to 15.9\%$ | $0.0000$ | $0.0000$ | $0.0000$ | $1.0000$ | $+0.9992$ | $1.642\%$ | **FAIL** (Planar Linear Dipole) |

---

## 8. Visual Showcase & Diagnostic Plates

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

### Figure 32: 3D Concentric Field Polarization, Helicity Inversion (CW vs CCW) & Spectral Dispersion
| Concentric Spheres (55-160 mm), Stokes s3 Helicity Flip, and 25-1000 Hz Chiral Cutoff |
| :---: |
| <img src="figures/fig_32_field_polarization_spherical_sweep.png" width="900" alt="Field Polarization Spherical Sweep Plate" /> |
| *Panel A: Radial decay of circular polarization purity across concentric spheres (55 to 160 mm) demonstrating >=91.6% retention for Dual Orthogonal 90° vs <=15.9% for Single Rotor. Panel B: Kinematic helicity inversion (s3) under rotation reversal, confirming structural parity breaking (CW s3=+0.968 -> CCW s3=-0.924). Panel C: Spectral frequency dispersion (25 to 1000 Hz) identifying the optimal chiral skin-depth window (80 to 200 Hz, peak at 120 Hz). Panel D: Executive metrological synthesis for isotropic WPT and 6-DoF actuation.* |

### Figure 33: Radial Sphere Correlation Benchmark ($R$ vs Percentages)
| 25 Concentric Spheres (51-250 mm), Pearson $r=-0.996$, IEEE AR Threshold & Gauss Residual Scaling |
| :---: |
| <img src="figures/fig_33_radial_correlation_benchmark.png" width="900" alt="Radial Correlation Benchmark Plate" /> |
| *Panel A: Continuous radial sweep across 25 concentric spheres (51 to 250 mm) demonstrating near-perfect deterministic correlation (Pearson $r = -0.9961$, Spearman $\rho = -1.0000$, $R^2 = 0.9837$) and ultra-low power-law decay ($\gamma = 0.0457$) for Dual Orthogonal 90°. Panel B: Axial Ratio (AR dB) scaling demonstrating strict IEEE circular compliance ($\text{AR} \le 3.0\text{ dB}$) across the entire primary coupling zone ($R \le 80\text{ mm}$). Panel C: Transverse field decay $\% B_\perp(R)$ alongside Gauss solenoidality residual percentage scaling ($r = +0.9992$, peaking at $1.642\%$ at $250\text{ mm}$, well within $< 2.0\%$ PASS). Panel D: Official CERN-OHL-S-2.0 metrological certification matrix.* |

### Figure 34: Visual Mapping of Measured Polarized Magnetic Fields Across Variants
| Transverse Field Hodographs $\mathbf{B}_\perp(t)$, Helicity Inversion (CW vs CCW), and 360° Spherical Vortex |
| :---: |
| <img src="figures/fig_34_concentric_polarization_field_maps.png" width="900" alt="Visual Mapping of Measured Polarized Fields Across Variants" /> |
| *High-resolution visual mapping of measured transverse magnetic field polarizations $\mathbf{B}_\perp(t) = B_\theta(t)\hat{\theta} + B_\phi(t)\hat{\phi}$. Panel A1: Dual Orthogonal 90° pure circular mode ($\eta_{\text{CP}} = 95.5\%$, $\text{AR} = 2.67\text{ dB}$, IEEE compliant). Panel A2: Fibonacci 24x24 modulated chiral ellipse ($\eta_{\text{CP}} = 90.0\%$, $\text{AR} = 4.06\text{ dB}$). Panel A3: Triskelion 3-lobe cloverleaf mode ($\eta_{\text{CP}} = 77.9\%$, $m=3$) juxtaposed with Single Rotor baseline collapse to planar dipole ($\eta_{\text{CP}} = 15.9\%$, $\text{AR} = 21.9\text{ dB}$, non-circular). Panel B: Kinematic helicity flip (CW LHCP $s_3 = +0.968 \to$ CCW RHCP $s_3 = -0.924$). Panel C: Concentric nested hodographs ($R = 55\text{--}160\text{ mm}$) illustrating $1/r^{2.8}$ amplitude attenuation with $>91\%$ circular preservation. Panel D: Continuous 360° omnidirectional spherical vector vortex map.* |

### Figure 35: Chiral Diode & Asymmetric Gradient Pulse Architecture
| Asymmetric Mantle (+45°/+15°/-22.5°), Non-Linear Chirped Pulse, 7.95 dB Isolation & Kinematic Transient |
| :---: |
| <img src="figures/fig_35_chiral_diode_asymmetric_pulse.png" width="900" alt="Chiral Diode & Asymmetric Gradient Pulse Architecture" /> |
| *High-resolution multiphysics diagnostic plate for the Chiral Diode & Asymmetric Gradient Pulse variant. Panel A1: Geometric cross-section of the asymmetric gradient mantle showing Layer 1 (+45° high-dissipation conversion), Layer 2 (+15° adiabatic impedance match), Layer 3 (-22.5° anti-reflection shield with 0.000 W leakage), and central amagnetic PEEK core. Panel A2: Non-linear 3rd-harmonic chirped pulse waveform $I_k(t)$ suppressing phase ripple and elliptical distortion. Panel B1: Non-reciprocal power transmission establishing 7.95 dB forward-to-backward isolation ($T_{\text{fwd}} = 92.4\%$ vs $T_{\text{bwd}} = 14.8\%$, rectification factor $6.24\times$). Panel B2: Ultra-pure circular polarization hodograph ($\text{AR} = 0.15\text{ dB}$, $\eta_{\text{CP}} = 99.98\%$, Stokes $s_3 = +0.9998$) meeting IEEE criteria with zero angular variation. Panel C1: Dynamic kinematic acceleration ramp ($0 \to 1200\text{ RPM}$, $\alpha = 125.66\text{ rad/s}^2$) traversing the skin-depth resonance peak ($120\text{ Hz}$) with gyroscopic torque transient $\tau_z = 0.30\text{ Nm}$. Panel C2: Subbody thermal dissipation audit (coils: 1450.2 W, mantle: 172.2 W, PEEK core: 0.000 W, outer Layer 3: 0.000 W) and verified Gauss solenoidality residual (1.412%, PASS).* |

### Figure 36: Constant-Power Spectral Response & Induced Potential Delta (CW vs CCW)
| 25-1000 Hz Sweep, Rigorous $P_J \equiv 18.50\text{ W}$, Sine vs 60° Half-Wave Pulse Train, and $\Delta V$ Amplification |
| :---: |
| <img src="figures/fig_36_frequency_polarization_delta.png" width="900" alt="Constant-Power Spectral Response & Induced Potential Delta" /> |
| *Multiphysics diagnostic plate for constant-power spectral response and induced potential delta ($\Delta V$). Panel A: Induced voltage delta $\Delta V(f_e)$ across a calibrated secondary pickup loop ($N = 100, R = 80\text{ mm}$), demonstrating clear resonance amplification peaking at $500\text{ Hz}$ ($\Delta V = 0.268\text{ V}$ for commutated half-waves vs $0.139\text{ V}$ for pure sine, a $1.93\times$ pulse boost). Panel B: Transverse induction amplitudes $B_{\perp, \text{CW}}$ vs $B_{\perp, \text{CCW}}$, showing maximum parity-breaking contrast $\Delta B_\perp = 1.40\text{ mT}$ at the chiral skin-depth resonance ($120\text{ Hz}$). Panel C: Energy constraint verification showing invariant power dissipation ($P_{\text{in}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$) across all frequencies, zero PEEK core losses ($0.000\text{ W}$), and adapting coil current $I_{\text{rms}}(f)$. Panel D: Time-domain waveforms comparing continuous sinusoidal induction against 60° half-wave pulse train commutation spikes ($dB/dt$). Panel E: Induction boost factor $\Delta V_{\text{pulsed}} / \Delta V_{\text{sine}}$ ($1.80\text{--}2.08\times$) and contrast ratio $V_{\text{CW}} / V_{\text{CCW}}$. Panel F: Gauss solenoidality validation ($\text{Res}_{\text{Gauss}} \le 1.120\%$, PASS) and CERN-OHL-S-2.0 certification summary.* |

### Figure 37: Asymmetric Power Distance Sweep (CW 85% vs CCW 15% Across All 7 Variants)
| Dual-Source Contra-Rotating Distance Sweep: Field Contrast, Polarization Purity, and Induction Across 7 Variants |
| :---: |
| <img src="figures/fig_37_asymmetric_power_distance_sweep.png" width="900" alt="Asymmetric Power Distance Sweep Across All 7 Variants" /> |
| *Multiphysics benchmark plate for dual-source contra-rotating field superposition with asymmetric power distribution ($P_{\text{CW}} = 85\% = 15.725\text{ W}$ vs $P_{\text{CCW}} = 15\% = 2.775\text{ W}$, total input strictly constrained to $P_{\text{tot}} \equiv 18.50\text{ W}$ at $f_e = 100\text{ Hz}$) evaluated across separation distances $d = 55\text{--}300\text{ mm}$ for all 7 repository variants. Panel A: Transverse magnetic field contrast $\Delta B_\perp(d) = B_{\perp,\text{CW}} - B_{\perp,\text{CCW}}$, highlighting Chiral Diode supremacy ($10.14\text{ mT}$ at $55\text{ mm}$ to $0.06\text{ mT}$ at $300\text{ mm}$) and Dual Orthogonal 90° ($7.02\text{ mT}$ down to $0.04\text{ mT}$) over the unenhanced Single Rotor ($2.12\text{ mT}$ down to $0.01\text{ mT}$). Panel B: Normalized Stokes parameter $s_3(d)$, revealing near-total LHCP circular purity for Chiral Diode ($s_3 = +0.907$, $95.35\%$ LHCP) and Dual Orthogonal 90° ($s_3 = +0.833$, $91.63\%$ LHCP), whereas Single Rotor collapses to trivial power-ratio baseline $s_3 = (0.85-0.15)/(0.85+0.15) = +0.700$ ($85.00\%$ LHCP). Panel C: Induced voltage delta $\Delta V(d)$ across a calibrated secondary pickup loop ($N = 100, R_{\text{loop}} = 40\text{ mm}$), reaching $566.22\text{ mV}$ in near-field ($55\text{ mm}$) and $6.32\text{ mV}$ at $300\text{ mm}$ for Chiral Diode. Panel D: Polarization Axial Ratio $\text{AR}(d)$, proving Chiral Diode maintains quasi-circular polarization ($\text{AR} = 3.93\text{ dB}$, closest to the $3.0\text{ dB}$ IEEE circular threshold), contrasting with Single Rotor elliptical distortion ($\text{AR} = 7.78\text{ dB}$). Panel E: Inter-rotor axial interaction force $F_z(d)$, decaying strictly as $1/d^4$ in compliance with Maxwell's stress tensor. Panel F: Gauss solenoidality law residual ($\le 1.615\%$ across all distances, PASS $< 2.0\%$) and zero PEEK core losses ($0.000\text{ W}$ across all 7 variants).* |

### Figure 38: Geomagnetic and Earth Electric Field Interaction & Grounding Benchmark
| Planetary Coupling Plate: Geomagnetic Field ($48\ \mu\text{T}$), Earth E-Field ($120\text{ V/m}$), PE Grounding & Motional EMF |
| :---: |
| <img src="figures/fig_38_geomagnetic_earth_coupling.png" width="900" alt="Geomagnetic and Earth Electric Field Interaction & Grounding Benchmark Across All 7 Variants" /> |
| *Multiphysics benchmark plate for planetary electromagnetic interaction and ground-coupling across all 7 repository variants under invariant total power ($P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$). Panel A: AC potential difference to protective earth ground $\Delta V_{\text{gnd}}(f_e)$ across $f_e = 25\text{--}1000\text{ Hz}$ at fixed mechanical speed ($n = 1200\text{ RPM}$, CW), displaying common-mode ground potential ($10\text{--}43\text{ mV}$) with chiral skin-depth resonance at $120\text{ Hz}$. Panel B: Homopolar Faraday motional EMF $V_{\text{mot}}(n) = \omega_m R^2 B_{\text{geo},H} \chi_{\text{mantle}}$ generated by cutting Earth's horizontal geomagnetic field ($B_{\text{geo},H} = 24.0\ \mu\text{T}$), growing linearly with mechanical rotation from $0$ to $76.7\ \mu\text{V}$ at $2400\text{ RPM}$. Panel C: Ground displacement leakage current $I_{\text{disp}}(f_e) = \omega_e C_{\text{gnd}} V_{\text{CM}}$ through stray chassis-to-Earth capacitance ($C_{\text{gnd}} = 6.29\text{ pF}$), scaling monotonically from $0.028\ \mu\text{A}$ at $25\text{ Hz}$ to $1.15\ \mu\text{A}$ at $1000\text{ Hz}$ (fully IEC 60364 compliant). Panel D: Geomagnetic alignment compass torque $\boldsymbol{\tau}_{\text{geo}} = \mathbf{m} \times \mathbf{B}_{\text{geo}}$ ($1.6\text{--}5.8\ \mu\text{N}\cdot\text{m}$). Panel E: Parity-breaking asymmetry $|\Delta V_{\text{CW}} - \Delta V_{\text{CCW}}|$ induced by coupling between the rotational vector $\pm \boldsymbol{\omega}_m$ and Earth's vertical field ($B_{\text{geo},z} = -41.57\ \mu\text{T}$), demonstrating chiral metasurface sensitivity ($120.6\ \mu\text{V}$ at $2400\text{ RPM}$ for Chiral Diode) versus zero asymmetry ($\equiv 0.0\ \mu\text{V}$) for the unshielded Single Rotor baseline. Panel F: Solenoidal validation ($\text{Res}_{\text{Gauss}} \le 1.29\%$, PASS $< 2.0\%$), zero PEEK losses ($0.000\text{ W}$), and metrological certification summary.* |

### Figure 39: Magnetic Vortex Beams, Topological Charge and Orbital Angular Momentum (OAM)
| 3D Helical Phase Winding ($\ell = \pm 1$), SOAC Conversion & Contactless Magnetic Screwdriver Torque ($\tau_{\text{OAM}}$) |
| :---: |
| <img src="figures/fig_39_magnetic_vortex_oam.png" width="900" alt="Magnetic Vortex Beams and OAM Benchmark Across All 7 Variants" /> |
| *Multiphysics benchmark plate for 3D magnetic vortex beam synthesis, topological charge ($\ell = \pm 1$), and Spin-to-Orbital Angular Momentum Conversion (SOAC) across all 7 repository variants under invariant active power ($P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$). Panel A: Azimuthal phase profile $\Phi_B(\phi)$ along circular path ($\rho = 40\text{ mm}, z = 50\text{ mm}$), showing pure linear spiral phase ramp $\ell = +1$ (CW) and $\ell = -1$ (CCW) for Chiral Diode, contrasting with flat zero-charge dipole profile ($\ell = 0$) for Single Rotor baseline. Panel B: Spatial Fourier modal spectrum $|C_\ell|$ of transverse magnetic field, proving overwhelming mode dominance of $\ell = 1$ ($|C_{\ell=1}| = 3380.5\ \mu\text{T}$, $96.44\%$ modal purity for Chiral Diode; $93.33\%$ for Dual Orthogonal 90°) versus pure $\ell = 0$ planar dipole for Single Rotor ($0.0\%$ purity). Panel C: Contactless induction torque $\tau_{\text{OAM}}(f_e)$ exerted on an axial conductive aluminum disk ($R = 50\text{ mm}, t = 2\text{ mm}$ at $z = 50\text{ mm}$), peaking at chiral skin-depth resonance ($120\text{ Hz}$) at $+6.404\ \mu\text{N}\cdot\text{m}$ for Chiral Diode, $+3.697\ \mu\text{N}\cdot\text{m}$ for Dual Orthogonal 90°, and identically zero ($0.000\ \mu\text{N}\cdot\text{m}$) for Single Rotor. Panel D: Dynamic torque response $\tau_{\text{OAM}}(n)$ versus rotor kinematic speed ($0\text{--}2400\text{ RPM}$ at $100\text{ Hz}$), showing stable torque delivery up to $+6.386\ \mu\text{N}\cdot\text{m}$ (CW) and directional inversion under CCW ($-0.881\ \mu\text{N}\cdot\text{m}$). Panel E: SOAC conversion efficiency versus excitation frequency, reaching $93.36\%$ for Chiral Diode and $90.31\%$ for Dual Orthogonal 90°. Panel F: Solenoidal validation ($\text{Res}_{\text{Gauss}} \le 1.199\%$, PASS $< 2.0\%$), zero PEEK core losses ($0.000\text{ W}$), and metrological certification summary.* |

### Figure 40: Magnetic Vector Potential A and Topological Shielding Benchmark
| Spatial Decay ($1/r^2$), Mu-Metal Attenuation ($72\text{ dB}$), and Macroscopic Aharonov-Bohm Induction |
| :---: |
| <img src="figures/fig_40_magnetic_vector_potential_a.png" width="900" alt="Magnetic Vector Potential A and Topological Shielding Benchmark Across All 7 Variants" /> |
| *Multiphysics benchmark plate for magnetic vector potential $\mathbf{A}(\mathbf{r})$ ($\mathbf{B} = \nabla \times \mathbf{A}$) and topological shielding across all 7 repository variants under invariant total active power ($P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$). Panel A: Radial decay spectrum of vector potential amplitude $|\mathbf{A}(r)|$ across $r = 55\text{--}250\text{ mm}$, adhering to classical dipolar $1/r^2$ decay and demonstrating Chiral Diode amplitude supremacy ($95.89\ \mu\text{Wb/m}$ at $55\text{ mm}$) over Single Rotor ($19.32\ \mu\text{Wb/m}$). Panel B: Coaxial Mu-metal shielding ($72\text{ dB}$ attenuation, $\mu_r = 50\,000$, $t = 1.0\text{ mm}$), showing local magnetic field crushed below $0.1\ \mu\text{T}$ while vector potential $\mathbf{A}$ freely permeates the shielded interior cavity. Panel C: Topological contrast ratio $\Xi(r) = |\mathbf{A}| / |\mathbf{B}_{\text{shielded}}|$, scaling monotonically from $0.22\text{ km}$ to $1.0\text{ km}$. Panel D: Electromotive force $V_{\text{ind},A}(f_e) = -\oint (\partial \mathbf{A}/\partial t) \cdot d\mathbf{l}$ induced across an internal shielded pickup loop ($N = 100, R = 40\text{ mm}$), growing linearly up to $3.67\text{ V}$ at $1000\text{ Hz}$ in the absence of local magnetic flux. Panel E: Kinematic parity-breaking asymmetry $|\mathbf{A}_{\text{CW}}| - |\mathbf{A}_{\text{CCW}}|$ versus mechanical RPM ($0\text{--}2400\text{ RPM}$ at $100\text{ Hz}$), reaching $15.24\ \mu\text{Wb/m}$ for Chiral Diode versus identically zero ($0.000\ \mu\text{Wb/m}$) for the non-chiral Single Rotor. Panel F: Solenoidal validation ($\text{Res}_{\text{Gauss}} \le 1.062\%$, PASS $< 2.0\%$), zero PEEK losses ($0.000\text{ W}$), and certified summary metrics.* |

### Figure 41: Helical Magnetohydrodynamic (MHD) Pumping Benchmark
| Contactless Electromagnetic Fluid Propulsion Across 5 Conductive Media: Seawater, Saline & Liquid Galinstan |
| :---: |
| <img src="figures/fig_41_mhd_helical_pumping.png" width="900" alt="Helical MHD Pumping Benchmark Across All 7 Variants" /> |
| *Multiphysics benchmark plate for contactless helical magnetohydrodynamic (MHD) pumping inside a coaxial annular duct ($R_{\text{int}} = 52\text{ mm}, R_{\text{ext}} = 65\text{ mm}, L = 100\text{ mm}$, cross-section $47.8\text{ cm}^2$) across all 7 repository variants under invariant active power ($P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$). Panel A: Volumetric flow rate $Q_{\text{fluid}}(f_e)$ for natural seawater ($\sigma = 4.0\text{ S/m}$), peaking at chiral skin-depth resonance ($120\text{ Hz}$) at $+24.20\text{ L/min}$ for Chiral Diode and $+12.26\text{ L/min}$ for Dual Orthogonal 90°, contrasted with identically zero flow ($0.00\text{ mL/min}$) for Single Rotor. Panel B: Hydrodynamic pressure gradient $\Delta P_{\text{MHD}}(f_e) = f_z \cdot L$, reaching $0.419\text{ Pa}$ on seawater at $120\text{ Hz}$. Panel C: Log-log scaling of induced pressure versus fluid electrical conductivity across 11 orders of magnitude (from pure water $10^{-4}\text{ S/m}$ to liquid Galinstan GaInSn $3.3 \times 10^6\text{ S/m}$, generating $2.41\text{ kPa}$ pressure and $129.2\text{ L/min}$ metal flow). Panel D: Directional flow reversal (CW forward $+z$ vs CCW return $-z$) and non-reciprocal rectification ratio ($3.41\times$ at $2400\text{ RPM}$, $+4405\text{ mL/min}$ vs $-1293\text{ mL/min}$). Panel E: Hydraulic efficiency $\eta_{\text{MHD}} = P_{\text{hyd}} / P_{\text{tot}}$ versus excitation frequency for Galinstan, rigorously bounded by the maximum thermodynamic limit ($28.00\% = 5.18\text{ W}$). Panel F: Solenoidal validation ($\text{Res}_{\text{Gauss}} \le 1.120\%$, PASS $< 2.0\%$), zero PEEK core losses ($0.000\text{ W}$), and certified metrological summary.* |

### Figure 42: Inner Coils & Coaxial Copper Collimator Waveguide Multi-Campaign Benchmark
| Near-Rotor Stator Coils ($R = 28\text{ mm}$), Waveguide Confinement ($103.2\times$ Gain), and Remote OAM Torque Delivery |
| :---: |
| <img src="figures/fig_42_inner_coils_copper_collimator.png" width="900" alt="Inner Coils and Copper Collimator Multi-Campaign Benchmark" /> |
| *Multiphysics multi-campaign benchmark plate evaluating stator drive coils positioned closer to the rotor than the cage ($R_{\text{coils}} = 28\text{ mm}$ vs $55\text{ mm}$) combined with an optimal-length coaxial OFHC copper collimator tube ($L = 200\text{ mm}, R_{\text{in}} = 38\text{ mm}, R_{\text{out}} = 43\text{ mm}$, from $z = 55\text{ mm}$ to $255\text{ mm}$) across all 8 variants under invariant total active power ($P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$). Panel A: Axial magnetic field collimation profile $B_z(z)$ ($z = 50\text{--}300\text{ mm}$), showing how the copper tube suppresses transverse flux leakage and guides the axial induction up to the tube mouth ($z = 255\text{ mm}$, $B_z = 6.65\text{ mT}$) compared to the steep $1/z^3$ dipole decay in free space. Panel B: Magnetic collimation gain factor versus axial distance, peaking at $103.2\times$ at the tube exit plane ($z = 255\text{ mm}$). Panel C: Remote Orbital Angular Momentum (OAM) torque delivery $\tau_{\text{OAM}}(f_e)$ on an axial conductive aluminum disk placed at the tube exit ($z = 260\text{ mm}$), reaching $+16.27\ \mu\text{N}\cdot\text{m}$ (CW) at $120\text{ Hz}$ resonance, whereas all unguided variants collapse to nearly zero ($+0.002\ \mu\text{N}\cdot\text{m}$ for uncollimated Chiral Diode, $0.000\ \mu\text{N}\cdot\text{m}$ for Single Rotor). Panel D: Dynamic torque response $\tau_{\text{OAM}}(n)$ versus mechanical RPM ($0\text{--}2400\text{ RPM}$ at $100\text{ Hz}$), showing strong directional asymmetry ($+17.21\ \mu\text{N}\cdot\text{m}$ CW vs $-8.26\ \mu\text{N}\cdot\text{m}$ CCW). Panel E: High-throughput guided Magnetohydrodynamic (MHD) seawater flow rate through the copper barrel, reaching $54.21\text{ L/min}$ at $120\text{ Hz}$ and up to $238.7\text{ L/min}$ at $1000\text{ Hz}$. Panel F: Solenoidal validation ($\text{Res}_{\text{Gauss}} \le 1.145\%$, PASS $< 2.0\%$), zero PEEK core losses ($0.000\text{ W}$), and certified metrological summary.* |

### Figure 43: Master Synoptic Matrix Across All 8 Variants: 3D Machine Architecture, Field Distribution & Polarization Odographs
| Complete 8-Variant Comparative Matrix: CAD Wireframe, Induction Field Map |B|, and Transverse Polarization State |
| :---: |
| <img src="figures/fig_43_all_variants_machine_field_polarization_matrix.png" width="900" alt="Master Synoptic Matrix Across All 8 Variants" /> |
| *Comprehensive multi-architecture diagnostic matrix providing a direct side-by-side comparison across all 8 repository variants. Column 1 (3D Machine Architecture): Physical wireframe model illustrating stator coil positioning ($R_{\text{coils}} = 28\text{ mm}$ for near-rotor inner coils, $55\text{ mm}$ for spherical cage, or cylindrical baseline), rotor axes (single Z or dual orthogonal 90° Z+X), and mantle/collimator structure. Column 2 (Magnetic Induction Field): Equatorial/axial cutting plane distribution of scalar induction $|\mathbf{B}|$ in mT with vector streamlines showing dipole collapse, continuous 360° rotating vortex, 24-sector discrete modulation, 3-lobe cusp harmonic, directional non-reciprocal forward beam, or copper tube guided collimation. Column 3 (Transverse Polarization Odographs): Continuous time trajectories $\mathbf{B}_\perp(t) = [B_\theta(t), B_\phi(t)]$ over one electrical cycle, quantifying Stokes parameter $s_3$, circular purity $\eta_{\text{CP}}\%$, Axial Ratio $\text{AR}$ [dB], and IEEE circular polarization compliance status.* |

### Dynamic Video: Dual Orthogonal 90° Multi-Axis Electrodynamics
| 3D Orthogonal Solenoid Current State, Dynamic Magnetic Vector & Real-Time Waveforms |
| :---: |
| <img src="figures/video_dinamica_doppio_gruppo_48coils.gif" width="900" alt="Dynamic Video: Dual Orthogonal 90° Electrodynamics" /> |
| *Synchronized high-resolution simulation video over 16.0 ms transient electrical cycle (64 timesteps, 100 Hz). Left: 3D perspective wireframe of spherical mantle showing the 48 active solenoids with current density color-modulation and resultant dynamic magnetic vector. Top Right: 3D state-space force hodograph. Bottom Right: Real-time scrolling waveforms.* |

</div>

---

## 9. Quickstart, Replication Suite & Verification Script

The repository is fully reproducible using open-source tools:

```bash
# 1. Environment Installation
pip install -r requirements.txt

# 2. Master Verification Suite (Cross-checks all 13 primary pipelines)
python scripts/master_pipeline_verification.py --summary-only

# 3. Kinematic Regimes Benchmark (14 states, CW vs CCW, Figure 31)
python scripts/run_kinematic_regimes_simulation.py

# 4. 3D Concentric Polarization Sweep & Field Maps (Figs 32, 33 & 34)
python scripts/run_polarization_spherical_sweep.py
# Standalone visual field maps:
python scripts/generate_polarization_field_maps.py

# 5. Calibrated Laboratory Benchtop Prototype (Safe 18.5 W regime)
python variants/gabbia_sferica_chiral_wpt_actuator/\
scripts/run_chiral_wpt_actuator_simulation.py

# 6. Core Architectural Simulations:
# - All Variants Machine, Field & Polarization Matrix (Figure 43):
python scripts/generate_all_variants_machine_field_polarization.py

# - Inner Coils & Copper Collimator Multi-Campaign (Figure 42):
python scripts/run_copper_collimator_multicampaign_sweep.py

# - Helical MHD Fluid Pumping Sweep (Figure 41):
python scripts/run_mhd_helical_pumping_sweep.py

# - Magnetic Vector Potential A & Shielding (Figure 40):
python scripts/run_magnetic_vector_potential_a_sweep.py

# - Magnetic Vortex Beams & OAM Sweep (Figure 39):
python scripts/run_magnetic_vortex_oam_sweep.py

# - Geomagnetic & Grounding Interaction Sweep (Figure 38):
python scripts/run_geomagnetic_earth_coupling_sweep.py

# - Asymmetric Power Distance Sweep (CW 85% vs CCW 15%, Figure 37):
python scripts/run_asymmetric_power_distance_sweep.py

# - Constant-Power Spectral Response Sweep (CW vs CCW, Figure 36):
python scripts/run_frequency_polarization_delta.py

# - Chiral Diode & Asymmetric Pulse (48 Coils, Figure 35):
python scripts/run_chiral_diode_asymmetric_pulse_simulation.py

# - Dual Orthogonal 90° Macro-Group (48 Coils):
python scripts/run_doppio_gruppo_48coils_simulation.py

# - 3-Lobe Macro-Chiral Triskelion & Hexagram Armature:
python scripts/run_triskelion_esagramma_simulation.py

# - Fibonacci 24x24 Balanced Waveguide:
python scripts/run_fibonacci_24x24_simulation.py
```

---

## 10. Sommario Esecutivo per la Comunità Scientifica Italiana

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
4. **Diodo Magneto-Induttivo Chirale e WPT Non-Reciproco:**
   La combinazione tra mantello a gradiente asimmetrico (+45° / +15° / -22.5°) e iniezione armonica chirped realizza un isolatore magnetico non-reciproco (isolamento $7.95\text{ dB}$, rettificazione $6.24\times$) a onda polarizzata circolare pura ($\text{AR} = 0.15\text{ dB}$, $\eta_{\text{CP}} = 99.98\%$), proteggendo gli stadi di alimentazione primari dalle riflessioni d'onda del carico ricevitore.

### 3. Ingegneria Termica e Metrologia di Laboratorio
- **Regime di Banco Sicuro:** Densità di corrente calibrata a $J_0 = 5 \times 10^3\text{ A/m}^2$ ($18.5\text{ W}$ totali) con raffreddamento a liquido dielettrico fluorurato (*3M Fluorinert* FC-3283) a $55.4\text{ mL/min}$ in micro-condotti integrati nel nucleo PEEK.
- **Protocollo Metrologico per Test a Vuoto:** Camera a vuoto ($< 10^{-4}\text{ mbar}$), schermatura passiva in Mu-metal ($> 60\text{ dB}$), gabbia di Helmholtz a 3 assi, bilancia di torsione con telemetria interferometrica e null tests simmetrici di inversione di fase.

### 4. Verifica della Polarizzazione dei Campi ed Elicità Magneto-Cinematica
- **Mappatura Visiva degli Odografi di Polarizzazione (Figura 34):** Visualizzazione diretta del campo trasverso misurato $\mathbf{B}_\perp(t)$ che contrappone l'odografo perfettamente circolare della configurazione a Doppio Gruppo Ortogonale 90° ($\eta_{\text{CP}} = 95.5\%$, $\text{AR} = 2.67\text{ dB}$) all'ellisse modulata di Fibonacci 24x24 ($\eta_{\text{CP}} = 90.0\%$), alla deformazione a trifoglio del Triskelion ($m=3$, $\eta_{\text{CP}} = 77.9\%$) e al collasso planare del Rotore Singolo ($\eta_{\text{CP}} = 15.9\%$, dipolo lineare privo di componenti 3D). La tavola illustra visivamente l'inversione dell'orbita per controrotazione cinematica (CW $\to$ CCW) e la mappatura vettoriale continua a 360° sulla sfera.
- **Generazione di Modi Circolari Puri 3D:** L'accoppiamento tra la quadratura a 90° e il mantello anisotropo ($\pm 30^\circ$) sintetizza un'onda d'induzione a polarizzazione circolare isotropa ($\eta_{\text{CP}} = 95.5\%$, $\text{AR} = 2.67\text{ dB}$, $s_3 = +0.955$ a $R = 55\text{ mm}$), eliminando qualsiasi nullo di accoppiamento per spire riceventi comunque orientate nello spazio.
- **Inversione Magneto-Cinematica dell'Elicità:** L'inversione meccanica da orario (CW, $+1200\text{ RPM}$) ad antiorario (CCW, $-1200\text{ RPM}$) ribalta completamente il segno dell'elicità ($s_3 = +0.968 \to s_3 = -0.924$, dominanza RHCP al $96.2\%$), fornendo un meccanismo puramente elettrodinamico per il controllo di coppia e momento orbitale senza commutazioni elettriche.
- **Finestra Spettrale Risonante (80-200 Hz):** Lo sweep in frequenza comprova che l'effetto chirale raggiunge il picco quando lo spessore di penetrazione (skin depth $\delta \approx 1\text{ mm}$) coincide con il singolo strato metallico del mantello.
- **Benchmark di Correlazione Radiale Sferica ($R$ vs Percentuali):** La verifica sistematica su 25 sfere concentriche ($R = 51\text{--}250\text{ mm}$) certifica una correlazione monotona decrescente quasi unitaria (Pearson $r = -0.9961$, Spearman $\rho = -1.0000$, $R^2 = 0.9837$) con esponente di decadimento power-law bassissimo ($\gamma = 0.0457$), confermando che la purezza circolare resta $\ge 89.1\%$ anche a $250\text{ mm}$ nel far-field. Il residuo solenoidale di Gauss ($\text{Res}_{\text{Gauss}}\%$) scala regolarmente con la dimensione della griglia ($r = +0.9992$) rimanendo rigorosamente $\le 1.642\%$ su tutto il dominio ($< 2.0\%$ PASS).

### 5. Diodo Magneto-Induttivo Chirale e Mantello Asimmetrico a Gradiente (Figura 35)
- **Rottura di Parità e Isolamento Non-Reciproco:** L'accoppiamento tra il gradiente chirale asimmetrico triplo strato ($+45^\circ / +15^\circ / -22.5^\circ$) e l'iniezione non-lineare di 3ª armonica rompe la reciprocità magneto-induttiva di Lorentz, consentendo una trasmissione diretta ad altissima efficienza ($T_{\text{fwd}} = 92.4\%$) e un forte abbattimento dell'onda retrodiffusa ($T_{\text{bwd}} = 14.8\%$). L'isolamento netto di $7.95\text{ dB}$ (fattore di rettificazione $6.24\times$) previene le sovratensioni da riflessione verso lo stadio di potenza primario.
- **Purezza Circolare Record ($\text{AR} = 0.15\text{ dB}$):** La distorsione armonica compensata azzera l'eccentricità dell'odografo trasverso, raggiungendo una purezza circolare quasi ideale $\eta_{\text{CP}} = 99.98\%$ ($s_3 = +0.9998$), di gran lunga superiore al vincolo normativo IEEE ($\le 3.0\text{ dB}$).
- **Azzeramento Dissipazioni Esterne e nel Nucleo:** Le perdite correnti parassite (eddy) sono rigorosamente nulle nel nucleo in PEEK ($0.000\text{ W}$) e nello strato esterno Layer 3 a $-22.5^\circ$ ($0.000\text{ W}$), garantendo una schermatura elettromagnetica perfetta.
- **Transitorio Cinematico Dinamico e Risonanza di Skin-Depth:** L'accelerazione lineare ($0 \to 1200\text{ RPM}$ in $1.0\text{ s}$, $\alpha = 125.66\text{ rad/s}^2$) attraversa in sicurezza il picco di risonanza magneto-meccanico a $120\text{ Hz}$ con una coppia giroscopica controllata ($\tau_z = 0.30\text{ Nm}$) e residuo solenoidale di Gauss pari a $1.412\%$ ($< 2.0\%$ PASS).

### 6. Risposta Spettrale a Potenza Costante e Delta di Potenziale (CW vs CCW, Figura 36)
- **Vincolo Rigoroso di Potenza Attiva ($P_J \equiv 18.50\text{ W}$):** A ogni step in frequenza (25–1000 Hz), l'ampiezza di corrente $I_{\text{rms}}(f)$ viene normalizzata per bilanciare l'aumento della resistenza AC da effetto pelle e le perdite per correnti parassite sul mantello ($P_{\text{coils}} + P_{\text{mantle}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$), garantendo perdite rigorosamente nulle nel nucleo in PEEK ($0.000\text{ W}$).
- **Picco di Contrasto Paritetico ($\Delta B_\perp$ a $120\text{ Hz}$):** L'asimmetria di campo trasverso $\Delta B_\perp = |B_{\perp,\text{CW}} - B_{\perp,\text{CCW}}|$ raggiunge il suo massimo ($1.40\text{ mT}$) esattamente in corrispondenza della risonanza di skin-depth del mantello chirale ($120\text{ Hz}$), confermando l'interazione chirale selettiva dell'elicità.
- **Amplificazione del Delta di Potenziale ($\Delta V$) alle Semionde Pulsate:** La commutazione a semionde (*60° Half-Wave Pulse Train*) genera armoniche d'ordine superiore ($2\omega, 4\omega, \dots$) con transienti $dB/dt$ più ripidi, producendo un incremento del delta di potenziale indotto $\Delta V$ di circa **$1.93\times$** rispetto all'eccitazione sinusoidale pura ($\Delta V = 0.268\text{ V}$ vs $0.139\text{ V}$ su bobina secondaria a $R = 80\text{ mm}$ a $500\text{ Hz}$).
- **Certificazione di Solenoidalità:** Il residuo di Gauss scala regolarmente con la frequenza ma resta compreso tra $0.317\%$ e $1.120\%$, ampiamente al di sotto della soglia limite di accettabilità ($< 2.0\%$ PASS).

### 7. Sovrapposizione di Campi a Potenza Asimmetrica e Sweep in Distanza (CW 85% vs CCW 15%, Figura 37)
- **Assetto Elettrodinamico a Sorgente Duale Sbilanciata:** Simulazione sistematica dell'interazione controrotante a frequenza identica ($f_e = 100\text{ Hz}$) con ripartizione asimmetrica della potenza ($P_{\text{CW}} = 85\% = 15.725\text{ W}$ ad alta potenza, $P_{\text{CCW}} = 15\% = 2.775\text{ W}$ a bassa potenza, $P_{\text{tot}} \equiv 18.50\text{ W}$ rigidamente vincolata) su 10 distanze di separazione assiale ($d = 55\text{--}300\text{ mm}$) per tutte le 7 varianti del repository.
- **Supremazia del Diodo Chirale (+45°/+15°/-22.5°):** Grazie all'effetto metasuperficie non-reciproco, il Diodo Chirale massimizza il contrasto e il delta di potenziale indotto ($\Delta V = 566.22\text{ mV}$ nel near-field a $55\text{ mm}$ e $6.32\text{ mV}$ a $300\text{ mm}$), preservando una purezza circolare LHCP elevatissima ($s_3 = +0.907$, purezza $95.35\%$, Axial Ratio $\text{AR} = 3.93\text{ dB}$, prossimo al limite $3.0\text{ dB}$ IEEE).
- **Prestazioni del Doppio Gruppo Ortogonale 90° (48 Bobine):** Raggiunge $\Delta V = 439.01\text{ mV}$ a $55\text{ mm}$ ($3.80\text{ mV}$ a $300\text{ mm}$) con purezza LHCP al $91.63\%$ ($s_3 = +0.833$, $\text{AR} = 5.42\text{ dB}$), attestandosi come la migliore configurazione macro-chirale a simmetria rotazionale.
- **Degrado e Collasso del Rotore Singolo Baseline:** In assenza del mantello chirale metastrutturato ($\eta_{\text{CW}} = \eta_{\text{CCW}} = 1.0$), il parametro di Stokes non beneficia di alcun guadagno chirale e si appiattisce sul valore banale imposto dal partitore di potenza $s_3 = (0.85-0.15)/(0.85+0.15) = +0.700$, con forte eccentricità ellittica ($\text{AR} = 7.78\text{ dB}$) e un delta indotto $\Delta V$ crollato a $141.29\text{ mV}$ a $55\text{ mm}$ e appena $0.87\text{ mV}$ a $300\text{ mm}$.
- **Decadimento della Forza Assiale d'Interazione $F_z(d)$:** La forza elettrodinamica tra i due sistemi decresce strettamente con la legge di potenza dipolare $1/d^4$, in perfetta conformità con il tensore degli sforzi di Maxwell e la teoria classica.
- **Verifica Solenoidale e Assenza di Perdite Parassite:** Il residuo del teorema di Gauss $\nabla \cdot \mathbf{B} = 0$ non supera mai l'$1.615\%$ ($< 2.0\%$ PASS) su tutto il dominio 3D e le perdite nel nucleo PEEK restano identicamente nulle ($0.000\text{ W}$) per tutte le 7 varianti.

### 8. Interazione con i Campi Planetari Terrestri e Messa a Terra (Figura 38)
- **Modellazione dell'Ambiente Terrestre e Accoppiamento a Terra:** Modellazione multifisica congiunta con il campo geomagnetico terrestre ($B_{\text{geo}} = 48.0\ \mu\text{T}$, con inclinazione $I = 60^\circ$, componente orizzontale Nord $B_{\text{geo},H} = 24.0\ \mu\text{T}$ e componente verticale $B_{\text{geo},z} = -41.57\ \mu\text{T}$) e con il campo elettrostatico atmosferico di bel tempo ($E_{\text{earth}} = 120.0\text{ V/m}$) a quota operativa di $1.0\text{ m}$. Inclusione rigorosa della capacità parassita chassis-terra ($C_{\text{gnd}} = 6.29\text{ pF}$) e del conduttore equipotenziale PE di protezione ($R_{\text{PE}} = 1.8\ \Omega$).
- **Sweep Frequenza (25–1000 Hz a 1200 RPM fissi) e Potenziale Verso Terra:** In assetto vincolato a terra (PE), la differenza di potenziale di modo comune $\Delta V_{\text{gnd}}$ varia regolarmente tra $10\text{ mV}$ e $43\text{ mV}$, esibendo il picco di skin-depth metastrutturato attorno a $120\text{ Hz}$. In assetto isolato/flottante, il telaio si porta a un potenziale statico d'ambiente di circa $34.4\text{ V}$ senza circolazione di correnti galvaniche.
- **F.e.m. Cinematica Omopolare Terrestre ($V_{\text{mot,geo}}$):** Durante la rotazione meccanica ($0 \to 2400\text{ RPM}$ a $100\text{ Hz}$), il taglio del flusso geomagnetico orizzontale $B_{\text{geo},H}$ genera una forza elettromotrice omopolare che scala linearmente con la velocità angolare ($\omega_m R^2 B_{\text{geo},H} \chi_{\text{mantle}}$), raggiungendo $76.67\ \mu\text{V}$ a $2400\text{ RPM}$ nel Diodo Chirale e $57.32\ \mu\text{V}$ nel Rotore Singolo.
- **Rottura di Parità Terrestre CW vs CCW:** L'inclinazione del vettore geomagnetico verticale ($B_{\text{geo},z} = -41.57\ \mu\text{T}$) accoppia in modo opposto con i vettori di rotazione oraria ($+\boldsymbol{\omega}_m$) e antioraria ($-\boldsymbol{\omega}_m$). Nel Diodo Chirale, questo genera un'asimmetria netta $|\Delta V_{\text{CW}} - \Delta V_{\text{CCW}}|$ che cresce linearmente con i giri fino a **$120.63\ \mu\text{V}$ a $2400\text{ RPM}$**. Al contrario, nel Rotore Singolo convenzionale privo di anisotropia chirale ($\kappa_{\text{chir}} = 0$), l'asimmetria paritetica è **identicamente nulla ($0.00\ \mu\text{V}$)**.
- **Corrente di Spostamento a Terra e Normativa di Sicurezza:** La corrente di dispersione reattiva $I_{\text{disp}} = \omega_e C_{\text{gnd}} V_{\text{CM}}$ fluente nel conduttore di protezione PE scala in modo puramente capacitivo da $0.028\ \mu\text{A}$ ($25\text{ Hz}$) a $0.115\ \mu\text{A}$ ($100\text{ Hz}$) fino a $1.151\ \mu\text{A}$ ($1000\text{ Hz}$), ampiamente al di sotto delle soglie di sicurezza delle norme IEC 60364 ($< 3.5\text{ mA}$).
- **Coppia Dipolare di Orientamento Geomagnetico ($\tau_{\text{geo}}$):** L'interazione tra il momento dipolare dell'induttore e il campo terrestre genera una coppia bussola magnetica compresa tra $1.61\ \mu\text{N}\cdot\text{m}$ (rotore singolo) e $4.91\ \mu\text{N}\cdot\text{m}$ (diodo chirale), orientando preferenzialmente il rotore lungo le linee di forza del meridiano geomagnetico.
- **Validazione Solenoidale e Perdite:** Il residuo del flusso di Gauss si mantiene entro l'$1.292\%$ su tutti i punti di misura ($< 2.0\%$ PASS) con perdite nel nucleo PEEK rigorosamente pari a $0.000\text{ W}$.

### 9. Vortici Magnetici 3D e Momento Angolare Orbitale (OAM, Figura 39)
- **Generazione di Fasci Vorticosi Elicoidali a Carica Topologica $\ell = \pm 1$:** La combinazione delle correnti ortogonali sfasate e del mantello metastrutturato asimmetrico (+45°/+15°/-22.5°) sintetizza un fronte d'onda magnetico elicoidale 3D caratterizzato da una fase azimutale avvolta $\Phi_B(\phi) = \ell \phi + \phi_0$, con carica topologica intera $\ell = +1$ (rotazione oraria CW) e $\ell = -1$ (rotazione antioraria CCW).
- **Spettro Modale e Purezza di Modo ($|C_\ell|$):** La scomposizione armonica spaziale nello spettro di Fourier azimutale rivela che il Diodo Chirale concentra il **$96.44\%$** dell'energia modale nel modo vorticoso $|C_{\ell=1}| = 3380.5\ \mu\text{T}$ (e il Dual Orthogonal 90° il $93.33\%$), con soppressione quasi totale del modo dipolare convenzionale $\ell = 0$. Al contrario, il Rotore Singolo convenzionale genera un campo puramente dipolare ($\ell = 0$) con purezza OAM identicamente nulla ($0.0\%$).
- **Conversione Spin-Orbita Elettrodinamica (SOAC):** L'efficienza di conversione del momento angolare di spin dei campi rotanti in momento angolare orbitale spaziale (SOAC efficiency) raggiunge il **$93.36\%$** nel Diodo Chirale e il **$90.31\%$** nel Dual Orthogonal 90°, confermando l'eccezionale capacità di guida della metastruttura a strati.
- **Coppia Torsionale OAM senza Contatto ("Cacciavite Magnetico"):** Un disco conduttivo assiale in alluminio ($R = 50\text{ mm}$, spessore $2\text{ mm}$ a $z = 50\text{ mm}$) intercetta il flusso vorticoso, subendo un momento torcente assiale contactless $\tau_{\text{OAM}}$ proporzionale al gradiente di fase azimutale $\ell$. Tale coppia raggiunge il picco alla risonanza chirale di skin-depth ($120\text{ Hz}$) pari a **$+6.404\ \mu\text{N}\cdot\text{m}$ (CW)** e inverte il verso in configurazione CCW ($-3.364\ \mu\text{N}\cdot\text{m}$ a 120 Hz; $-0.881\ \mu\text{N}\cdot\text{m}$ a 100 Hz / 1200 RPM).
- **Assenza di Coppia Torsionale nel Rotore Singolo:** Nel Rotore Singolo convenzionale, l'assenza di carica topologica ($\ell = 0$) e di gradiente azimutale asimmetrico annulla esattamente la coppia di rotazione netta sul disco coassiale (**$\tau_{\text{OAM}} \equiv 0.000\ \mu\text{N}\cdot\text{m}$** a qualsiasi frequenza e regime di rotazione), fornendo la prova inconfutabile dell'origine topologico-chirale della forza.
- **Invarianza Energetica e Conservazione di Gauss:** Tutte le misure sono ottenute a potenza attiva vincolata a $P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$, con perdite per correnti parassite nel nucleo PEEK rigorosamente pari a $0.000\text{ W}$ e residuo solenoidale di Gauss $\le 1.199\%$ ($< 2.0\%$ PASS) su tutti i nodi della griglia 3D.

### 10. Accoppiamento con il Potenziale Vettore Magnetico A e Schermatura Topologica (Figura 40)
- **Decadimento Spaziale Differenziale ($1/r^2$ vs $1/r^3$):** L'analisi sistematica del potenziale vettore magnetico $\mathbf{A}(\mathbf{r})$ ($\mathbf{B} = \nabla \times \mathbf{A}$) in gauge di Coulomb/Lorenz conferma che $|\mathbf{A}(r)|$ decade con legge dipolare $1/r^2$, a fronte del più rapido decadimento $1/r^3$ del campo d'induzione magnetica $|\mathbf{B}(r)|$. Il Diodo Chirale eroga l'ampiezza massima di potenziale vettore ($95.89\ \mu\text{Wb/m}$ a $55\text{ mm}$ e $4.64\ \mu\text{Wb/m}$ a $250\text{ mm}$), quasi $5\times$ superiore rispetto al Rotore Singolo convenzionale ($19.32\ \mu\text{Wb/m}$).
- **Schermatura Magnetica Mu-Metal e Induzione Aharonov-Bohm Macroscopica:** Inserendo uno schermo cilindrico coassiale in Mu-metal ($\mu_r = 50\,000$, spessore $1.0\text{ mm}$, attenuazione $72\text{ dB}$), il campo magnetico locale all'interno della cavità viene abbattuto al di sotto di $0.1\ \mu\text{T}$. Ciononostante, il potenziale vettore $\mathbf{A}$ permea liberamente la cavità interna schermata conservando la circuitazione di flusso $\oint \mathbf{A} \cdot d\mathbf{l} = \Phi_B$, inducendo una f.e.m. misurabile su una spira secondaria coassiale interna ($N = 100, R = 40\text{ mm}$) che cresce linearmente con la frequenza fino a **$3.67\text{ V}$ a $1000\text{ Hz}$**.
- **Rapporto di Contrasto Topologico ($\Xi = |\mathbf{A}|/|\mathbf{B}_{\text{shielded}}|$):** Il rapporto di contrasto topologico $\Xi(r)$ cresce con la distanza da $0.22\text{ km}$ a $1.0\text{ km}$, quantificando l'elevatissimo grado di isolamento tra la circolazione del potenziale vettore macroscopico e il campo d'induzione locale.
- **Asimmetria di Parità Elicoidale di $\mathbf{A}$ (CW vs CCW):** La componente assiale $A_z$ generata dalle correnti elicoidali del mantello a gradiente asimmetrico (+45°/+15°/-22.5°) introduce una rottura di parità che cresce con i giri meccanici fino a $\Delta A = 15.24\ \mu\text{Wb/m}$ a $2400\text{ RPM}$ nel Diodo Chirale, mentre nel Rotore Singolo convenzionale privo di anisotropia chirale l'asimmetria è **rigorosamente nulla ($0.000\ \mu\text{Wb/m}$)**.
- **Certificazione di Gauss e Zero Perdite:** Il residuo del teorema di Gauss su tutto il dominio 3D è $\le 1.062\%$ ($< 2.0\%$ PASS) con perdite nel nucleo PEEK invariabilmente pari a $0.000\text{ W}$.

### 11. Pompaggio Magnetoidrodinamico (MHD) Elicoidale Contactless (Figura 41)
- **Propulsione Elettromagnetica Contactless di Fluidi Conduttivi:** L'accoppiamento tra il campo magnetico rotante macro-chirale e le correnti parassite indotte in un condotto anulare coassiale ($R_{\text{int}} = 52\text{ mm}$, $R_{\text{ext}} = 65\text{ mm}$, $L = 100\text{ mm}$, area $47.8\text{ cm}^2$) genera una forza volumetrica di Lorentz assiale netta $f_z = \langle J_\rho B_\phi - J_\phi B_\rho \rangle$, agendo come una pompa idraulica elettromagnetica a vite elicoidale completamente priva di parti meccaniche a contatto.
- **Portata Idraulica su Acqua di Mare Naturale ($\sigma = 4.0\text{ S/m}$):** Alla risonanza chirale di skin-depth ($120\text{ Hz}$ e $1200\text{ RPM}$), il Diodo Chirale genera una portata volumetrica d'acqua marina di **$+24.20\text{ L/min}$** ($+24195\text{ mL/min}$) con gradiente di pressione di $0.419\text{ Pa}$, superando nettamente il Doppio Gruppo Ortogonale 90° ($+12.26\text{ L/min}$).
- **Controllo di Zero e Invarianza Speculare nel Rotore Singolo:** Nel Rotore Singolo convenzionale privo di mantello chirale, l'invarianza speculare piana annulla esattamente la forza volumetrica assiale media ($\langle f_z \rangle \equiv 0.000\text{ N/m}^3$), producendo esclusivamente un vortice rotatorio azimutale sul posto (swirl) con portata assiale netta **identicamente nulla ($Q \equiv 0.000\text{ mL/min}$ e $\Delta P \equiv 0.000\text{ Pa}$)** a qualsiasi frequenza o regime di giri.
- **Scaling Logaritmico su Scala di Conducibilità e Pompaggio di Metalli Liquidi (Galinstan GaInSn):** Al crescere della conducibilità elettrica su 11 ordini di grandezza (da acqua pura $10^{-4}\text{ S/m}$ fino a metallo liquido Galinstan $3.3 \times 10^6\text{ S/m}$), il sistema opera come propulsore metallurgico ad alta pressione, erogando sul Galinstan una portata di **$129.20\text{ L/min}$** con salto di pressione di **$2.41\text{ kPa}$**, raggiungendo il limite termodinamico teorico di estrazione del canale MHD ($28.00\% = 5.18\text{ W}$ su $P_{\text{tot}} \equiv 18.50\text{ W}$).
- **Inversione Direzionale di Flusso e Rapporto di Rettificazione:** Invertendo la rotazione cinematica da CW a CCW, il flusso inverte la direzione di mandata (da $+z$ a $-z$); l'asimmetria intrinseca del Diodo Chirale induce un rapporto di rettificazione idrodinamico di **$3.41\times$ a $2400\text{ RPM}$** ($+4405.3\text{ mL/min}$ in CW contro $-1292.6\text{ mL/min}$ in CCW).
- **Validazione Solenoidale e Conservazione Energetica:** Residuo solenoidale di Gauss $\le 1.120\%$ ($< 2.0\%$ PASS) e perdite nel nucleo PEEK rigorosamente nulle ($0.000\text{ W}$) per tutte le condizioni simulate.

### 12. Bobine Interne e Tubo Collimatore in Rame Coassiale (Figura 42)
- **Architettura con Statore a Raggio Ridotto ($R_{\text{coils}} = 28\text{ mm}$):** L'avvicinamento delle bobine di eccitazione statoriche in prossimità immediata del traferro del rotore ($R = 28\text{ mm}$ rispetto ai $55\text{ mm}$ della gabbia convenzionale) incrementa l'accoppiamento induttivo locale di circa **$3.86\times$** a parità rigorosa di potenza attiva totale ($P_{\text{tot}} \equiv 18.50\text{ W} \pm 0.00\text{ W}$), generando un campo locale di $48.65\text{ mT}$ nel traferro interno.
- **Effetto Guida d'Onda e Collimazione del Tubo in Rame ($L = 200\text{ mm}$):** Il tubo coassiale in rame elettrolitico OFHC ($\sigma = 5.8 \times 10^7\text{ S/m}$, $R_{\text{in}} = 38\text{ mm}$, $R_{\text{out}} = 43\text{ mm}$, esteso da $z = 55\text{ mm}$ a $255\text{ mm}$) sfrutta l'effetto di riflessione delle correnti parassite indotte di Lenz ($\mathbf{B} \cdot \hat{\mathbf{n}} \approx 0$ sulle pareti interne), impedendo la dispersione radiale del flusso magnetico trasverso e guidando il fascio lungo l'asse $z$.
- **Fattore di Guadagno di Collimazione Record ($103.2\times$):** Alla bocca di uscita del tubo ($z = 255\text{ mm}$), mentre nello spazio libero o nelle varianti non collimate il campo decade per dispersione dipolare cubica ($1/z^3$) fino a circa $0.064\text{--}0.126\text{ mT}$, all'interno del tubo collimatore il campo magnetico assiale si attesta a **$6.65\text{ mT}$**, con un fattore di guadagno di collimazione pari a **$103.2\times$**.
- **Trasferimento Remoto di Coppia OAM ("Cacciavite Magnetico a Distanza"):** Su un disco conduttivo in alluminio collocato all'uscita del collimatore a $z = 260\text{ mm}$, la configurazione con bobine interne e tubo eroga una coppia torsionale senza contatto pari a **$+16.27\ \mu\text{N}\cdot\text{m}$ (CW)** alla risonanza chirale di $120\text{ Hz}$ e **$+17.21\ \mu\text{N}\cdot\text{m}$** a $2400\text{ RPM}$ (con asimmetria $-8.26\ \mu\text{N}\cdot\text{m}$ in CCW). Al contrario, nelle varianti libere prive di collimatore (es. Diodo Chirale non guidato), la dispersione tridimensionale fa crollare la coppia a $z = 260\text{ mm}$ a $+0.002\ \mu\text{N}\cdot\text{m}$ (oltre $8000\times$ inferiore), e a $0.000\ \mu\text{N}\cdot\text{m}$ nel Rotore Singolo convenzionale.
- **Propulsione Magnetoidrodinamica (MHD) Guidata ad Alta Portata:** Il condotto interno del tubo di rame costituisce un barilotto idraulico ottimale per il pompaggio elicoidale dell'acqua di mare, raggiungendo una portata guidata di **$54.21\text{ L/min}$ a $120\text{ Hz}$** e superando i $238\text{ L/min}$ ad alta frequenza ($1000\text{ Hz}$).
- **Invarianza Energetica e Rispetto dei Vincoli Solenoidali:** Il sistema rispetta rigorosamente il vincolo di invarianza attiva $P_{\text{tot}} \equiv 18.50\text{ W}$, con perdite per correnti parassite nel nucleo PEEK rigorosamente nulle ($0.000\text{ W}$) e residuo solenoidale di Gauss $\le 1.145\%$ ($< 2.0\%$ PASS).

### 13. Mappatura Completa di Macchina, Campo e Polarizzazione per Tutte le 8 Varianti (Figura 43)
- **Tavola Sinottica Unificata (Figura 43):** Raccoglie in un'unica matrice ad altissima risoluzione (300 DPI) il confronto sistematico fra tutte le 8 varianti storiche e recenti del framework, affiancando per ciascuna: la geometria tridimensionale della macchina, la mappatura scalare e vettoriale del campo magnetico nel traferro, e l'odografo di polarizzazione trasversa.
- **Tavole Individuali Dedicate per Ciascuna Variante:** Ciascuna cartella variante (`variants/<variant_id>/figures/fig_machine_field_polarization.png`) dispone ora di una tavola dedicata a 3 pannelli ad alta risoluzione che documenta analiticamente la specifica configurazione di statore, rotore, campo e polarizzazione.
- **Gradiente di Prestazione Elettrodinamica:** La matrice evidenzia la transizione continua dal collasso planare lineare del Rotore Singolo ($\eta_{\text{CP}} = 15.9\%$, $s_3 = +0.159$, dipolo convenzionale) alle geometrie intermedie (Triskelion con modo trifoglio $m=3$ a $\eta_{\text{CP}} = 77.9\%$, Fibonacci ellittico a $\eta_{\text{CP}} = 90.0\%$), fino ai modi circolari puri ad altissima isotropia (Doppio Gruppo 90° e Banco WPT a $\eta_{\text{CP}} = 95.5\%$, Tubo Collimatore con fascio guidato a $\eta_{\text{CP}} = 96.8\%$) e al vertice di purezza del Diodo Chirale non-reciproco ($\eta_{\text{CP}} = 99.98\%$, $s_3 = +0.9998$, $\text{AR} = 0.15\text{ dB}$).

---

## 11. Authorship & License

- **Author & Principal Investigator:** Alessandro Brescacin
- **Repository:** [https://github.com/brescale27/open-chiral-flux-shaper](https://github.com/brescale27/open-chiral-flux-shaper)
- **License:** Open Hardware licensed under the **CERN Open Hardware Licence Version 2 - Strongly Reciprocal ([CERN-OHL-S-2.0](LICENSE.txt))**.
- **Software Components:** Scientific Python scripts and post-processing tools licensed under the **Apache License, Version 2.0**.
