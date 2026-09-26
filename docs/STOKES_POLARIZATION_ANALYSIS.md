# Stokes Polarization Analysis Tool

Independent post-processing utility for validating circular polarization purity
claims from Elmer FEM MagnetoDynamics VTU outputs.

## Purpose

Several breakthrough concepts in the Open Chiral Flux Shaper repository make
quantitative claims about circular polarization purity ($\eta_{\text{CP}}$), axial ratio
($\text{AR}$), and Stokes parameter $s_3$ on spherical shells of various radii around
the device:

| Claim | Radius | Reported $\eta_{\text{CP}}$ | Reported $\text{AR}$ | Reference |
|-------|--------|------------------|-------------|-----------|
| Dual Orthogonal 90° stator | $R = 55\text{ mm}$ | 95.5% | 2.67 dB | README §2.2 #1, Fig. 34 |
| Chiral Mode Shaper (Asym. Mantle) | $R = 55\text{ mm}$ | 99.98% (Compensated Mode Target) / 70.28% (Raw FEM Nodal) | 0.15 dB / 1.83 dB | README §2.2 #2, Fig. 35 |
| Apex Magnetic Cusp | apex | 97.2% | 15.43 dB | README §2.2 #5, Fig. 47 |
| Concentric shells decay | $R = 51\text{--}250\text{ mm}$ | $\ge 89.1\%$ | — | README §2.2 #1, Fig. 32 |

## Mathematical Formulation

For a transient simulation with $N_t$ timesteps at interval $\Delta t$, the time series
$B_\theta(t_n, \mathbf{r})$ and $B_\phi(t_n, \mathbf{r})$ on a spherical shell of radius $R$ are
extracted via discrete Fourier transform:

$$
\tilde{B}_\theta(\mathbf{r}) = \sum_n B_\theta(t_n, \mathbf{r}) e^{-i \omega_1 t_n}, \quad
\tilde{B}_\phi(\mathbf{r}) = \sum_n B_\phi(t_n, \mathbf{r}) e^{-i \omega_1 t_n}
$$

where $\omega_1 = \frac{2\pi}{N_t \Delta t}$ is the fundamental excitation frequency bin.

The four Stokes parameters (transverse-plane convention used in IEEE antenna
engineering) are:

$$
s_0 = |\tilde{B}_\theta|^2 + |\tilde{B}_\phi|^2
$$

$$
s_1 = |\tilde{B}_\theta|^2 - |\tilde{B}_\phi|^2
$$

$$
s_2 = 2 \operatorname{Re}(\tilde{B}_\theta \tilde{B}_\phi^*)
$$

$$
s_3 = 2 \operatorname{Im}(\tilde{B}_\theta \tilde{B}_\phi^*)
$$

Derived metrics:

$$
\eta_{\text{CP}} = \frac{s_0 + |s_3|}{2 s_0} \times 100\%
$$

$$
\text{AR} = 10 \log_{10}\left( \frac{s_0 + \sqrt{s_1^2 + s_2^2}}{s_0 - \sqrt{s_1^2 + s_2^2}} \right) \quad [\text{dB}]
$$

The helicity is determined by the sign of $s_3$: positive ($s_3 > 0$) = LHCP, negative ($s_3 < 0$) = RHCP.

## Usage

### Quick start (single shell)

```bash
python scripts/postprocess_stokes_polarization.py \
    --vtu-glob "variants/<variant_name>/work_dirs/run_*/results/macchina_out_t*.vtu" \
    --radius 0.055 \
    --frequency 100 \
    --dt 5e-4 \
    --output stokes_R55.json
```

### Multi-radius sweep

```bash
for R in 0.051 0.055 0.080 0.120 0.160 0.200 0.250; do
    python scripts/postprocess_stokes_polarization.py \
        --vtu-glob "variants/<variant>/results/macchina_out_t*.vtu" \
        --radius $R \
        --frequency 100 \
        --dt 5e-4 \
        --output stokes_R${R}.json \
        --quiet
done
```

### Batch multi-variant validation suite

```bash
python scripts/validate_polarization_claims.py --repo-root .
```

## Output Format

JSON output structure:

```json
{
  "metadata": {
    "tool": "postprocess_stokes_polarization.py",
    "radius_m": 0.055,
    "frequency_hz": 100.0
  },
  "stokes_parameters": {
    "s0": 1.234e-05,
    "s1": -1.234e-06,
    "s2": 2.345e-06,
    "s3": 5.001e-06
  },
  "metrics": {
    "circular_purity_pct": 70.28,
    "axial_ratio_dB": 1.83,
    "helicity": "LHCP (left-hand circular polarization)",
    "ieee_cp_pass": true
  },
  "per_point_stats": {
    "eta_cp_mean_pct": 77.01,
    "eta_cp_std_pct": 12.10,
    "eta_cp_min_pct": 52.44,
    "eta_cp_max_pct": 97.65
  }
}
```

## License

Apache License 2.0 (software) — compatible with CERN-OHL-S-2.0 (hardware).
