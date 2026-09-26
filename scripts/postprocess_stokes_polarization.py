#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stokes Polarization Analyzer for Elmer FEM MagnetoDynamics VTU outputs.

Computes the 4-component Stokes vector (s0, s1, s2, s3) and the standard
circular-purity metrics used in IEEE antenna characterization:

    s0 = |B_theta|^2 + |B_phi|^2           (total transverse power)
    s1 = |B_theta|^2 - |B_phi|^2           (linear horizontal/vertical)
    s2 = 2 Re(B_theta * conj(B_phi))      (linear +/-45 deg)
    s3 = 2 Im(B_theta * conj(B_phi))      (circular polarization)

    eta_CP = (s0 + |s3|) / (2 s0) * 100%   (circular polarization purity)
    AR_dB  = 10 log10( (s0 + sqrt(s1^2+s2^2)) / (s0 - sqrt(s1^2+s2^2)) )

For a transient Elmer simulation with N_t timesteps covering one electrical
period T = 1/f, the complex Fourier components at frequency f are extracted
via FFT, and the Stokes parameters are evaluated from those harmonics.

Usage:
    python postprocess_stokes_polarization.py \
        --vtu-glob "/path/to/run/results/macchina_out_t*.vtu" \
        --radius 0.055 \
        --frequency 100 \
        --dt 5e-4 \
        --output stokes_R55.json

License: Apache 2.0 / CERN-OHL-S-2.0
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import meshio


def load_vtu_bfield(vtu_path):
    """Load VTU file and return (points, B_field_at_points)."""
    m = meshio.read(str(vtu_path))
    pts = m.points
    bkey = None
    for cand in ('magnetic flux density', 'magnetic flux density e', 'B'):
        if cand in m.point_data:
            bkey = cand
            break
    if bkey is None:
        raise KeyError(f"No 'magnetic flux density' field in {vtu_path}. "
                       f"Available: {list(m.point_data.keys())}")
    B = m.point_data[bkey]
    return pts, B


def extract_shell_points(pts, R_target, tol=0.005):
    """Find indices of points lying on a spherical shell of radius R_target +/- tol."""
    r = np.linalg.norm(pts, axis=1)
    mask = np.abs(r - R_target) <= tol
    return np.where(mask)[0]


def to_spherical_components(B, pts):
    """Convert B from Cartesian (Bx, By, Bz) to spherical (B_r, B_theta, B_phi)."""
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r = np.linalg.norm(pts, axis=1)
    rho = np.hypot(x, y)
    rho_safe = np.where(rho < 1e-12, 1e-12, rho)
    r_safe = np.where(r < 1e-12, 1e-12, r)
    B_r     = (B[:, 0]*x + B[:, 1]*y + B[:, 2]*z) / r_safe
    B_theta = (B[:, 0]*x*z + B[:, 1]*y*z - B[:, 2]*rho*rho) / (r_safe * rho_safe)
    B_phi   = (-B[:, 0]*y + B[:, 1]*x) / rho_safe
    return B_r, B_theta, B_phi


def main():
    p = argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vtu-glob", required=True,
                   help="Glob pattern for VTU files")
    p.add_argument("--radius", type=float, required=True,
                   help="Radius of spherical shell for polarization analysis [m]")
    p.add_argument("--frequency", type=float, default=100.0,
                   help="Excitation frequency [Hz]")
    p.add_argument("--dt", type=float, default=5e-4,
                   help="Timestep [s]")
    p.add_argument("--tol", type=float, default=0.005,
                   help="Tolerance for shell point selection [m] (default 5 mm)")
    p.add_argument("--output", default=None,
                   help="Output JSON file path")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress informational output")
    args = p.parse_args()

    import glob
    vtu_files = sorted(glob.glob(args.vtu_glob))
    vtu_files = [Path(f) for f in vtu_files]
    if not vtu_files:
        vtu_files = sorted(Path("/").glob(args.vtu_glob.lstrip("/")))
    if not vtu_files:
        print(f"ERROR: No VTU files match {args.vtu_glob}", file=sys.stderr)
        sys.exit(1)
    if not args.quiet:
        print(f"Found {len(vtu_files)} VTU files")
        for f in vtu_files[:5]:
            print(f"  {f}")
        if len(vtu_files) > 5:
            print(f"  ... ({len(vtu_files) - 5} more)")

    R = args.radius
    tol = args.tol

    pts0, _ = load_vtu_bfield(vtu_files[0])
    shell_idx = extract_shell_points(pts0, R, tol)
    n_shell = len(shell_idx)
    if not args.quiet:
        print(f"\nShell at R = {R*1000:.1f} mm (+/-{tol*1000:.1f} mm): {n_shell} points")
    if n_shell < 10:
        print(f"WARNING: Only {n_shell} points on shell. Consider increasing --tol.", file=sys.stderr)

    n_t = len(vtu_files)
    B_theta_t = np.zeros((n_t, n_shell))
    B_phi_t = np.zeros((n_t, n_shell))
    B_r_t = np.zeros((n_t, n_shell))

    for i, vf in enumerate(vtu_files):
        pts, B = load_vtu_bfield(vf)
        if pts.shape != pts0.shape:
            raise ValueError(f"Mesh changed between timesteps: {pts.shape} vs {pts0.shape}")
        Br, Bt, Bp = to_spherical_components(B[shell_idx], pts[shell_idx])
        B_r_t[i] = Br
        B_theta_t[i] = Bt
        B_phi_t[i] = Bp
        if not args.quiet and (i + 1) % 5 == 0:
            print(f"  Loaded {i+1}/{n_t}")

    n_per_period = int(round(1.0 / (args.frequency * args.dt)))
    if not args.quiet:
        print(f"\nSampling: dt={args.dt*1e3:.2f} ms, f={args.frequency} Hz")
        print(f"Period = {1/args.frequency*1e3:.2f} ms = {n_per_period} steps")
        print(f"Have {n_t} samples ({n_t/n_per_period*100:.1f}% of one period)")

    fft_theta = np.fft.rfft(B_theta_t, axis=0)
    fft_phi = np.fft.rfft(B_phi_t, axis=0)

    bin_idx = 1 if n_t >= n_per_period else max(1, n_t // 4)
    if not args.quiet:
        print(f"Using FFT bin {bin_idx} for fundamental")

    B_theta_complex = fft_theta[bin_idx]
    B_phi_complex = fft_phi[bin_idx]

    s0 = np.abs(B_theta_complex)**2 + np.abs(B_phi_complex)**2
    s1 = np.abs(B_theta_complex)**2 - np.abs(B_phi_complex)**2
    s2 = 2.0 * np.real(B_theta_complex * np.conj(B_phi_complex))
    s3 = 2.0 * np.imag(B_theta_complex * np.conj(B_phi_complex))

    s0_mean = float(np.mean(s0))
    s1_mean = float(np.mean(s1))
    s2_mean = float(np.mean(s2))
    s3_mean = float(np.mean(s3))

    eta_CP = (s0_mean + abs(s3_mean)) / (2.0 * s0_mean + 1e-20) * 100.0
    linear_pol = np.sqrt(s1_mean**2 + s2_mean**2)
    if s0_mean > linear_pol:
        AR_dB = 10.0 * np.log10((s0_mean + linear_pol) / (s0_mean - linear_pol + 1e-20))
    else:
        AR_dB = float('inf')

    if s3_mean > 0:
        helicity = "LHCP (left-hand circular polarization)"
    elif s3_mean < 0:
        helicity = "RHCP (right-hand circular polarization)"
    else:
        helicity = "Linear (zero circular component)"

    ieee_pass = AR_dB <= 3.0

    print(f"\n{'='*70}")
    print(f"STOKES POLARIZATION ANALYSIS - Shell at R = {R*1000:.2f} mm")
    print(f"{'='*70}")
    print(f"  s0 (total transverse power)  = {s0_mean:.6e}")
    print(f"  s1 (linear H/V asymmetry)    = {s1_mean:+.6e}")
    print(f"  s2 (linear +/-45 asymmetry) = {s2_mean:+.6e}")
    print(f"  s3 (circular polarization)   = {s3_mean:+.6e}")
    print(f"  |s3|/s0                      = {abs(s3_mean)/s0_mean:.6f}")
    print(f"")
    print(f"  eta_CP (circular purity)     = {eta_CP:.4f} %")
    print(f"  AR  (axial ratio)            = {AR_dB:.4f} dB")
    print(f"  Helicity                     = {helicity}")
    print(f"  IEEE CP PASS (AR <= 3 dB)?   = {'YES' if ieee_pass else 'NO'}")
    print(f"{'='*70}")

    s3_per_point = 2.0 * np.imag(B_theta_complex * np.conj(B_phi_complex))
    s0_per_point = np.abs(B_theta_complex)**2 + np.abs(B_phi_complex)**2
    eta_per_point = (s0_per_point + np.abs(s3_per_point)) / (2.0 * s0_per_point + 1e-20) * 100.0
    eta_mean = float(np.mean(eta_per_point))
    eta_std = float(np.std(eta_per_point))
    eta_min = float(np.min(eta_per_point))
    eta_max = float(np.max(eta_per_point))
    print(f"\n  Per-point eta_CP stats:")
    print(f"    mean = {eta_mean:.4f} %")
    print(f"    std  = {eta_std:.4f} %")
    print(f"    min  = {eta_min:.4f} %")
    print(f"    max  = {eta_max:.4f} %")

    out = {
        "metadata": {
            "tool": "postprocess_stokes_polarization.py",
            "vtu_glob": str(args.vtu_glob),
            "n_vtu_files": n_t,
            "radius_m": R,
            "tol_m": tol,
            "n_shell_points": n_shell,
            "frequency_hz": args.frequency,
            "dt_s": args.dt,
            "samples_per_period": n_per_period,
            "fft_bin_used": bin_idx,
        },
        "stokes_parameters": {
            "s0": s0_mean, "s1": s1_mean, "s2": s2_mean, "s3": s3_mean,
        },
        "metrics": {
            "circular_purity_pct": eta_CP,
            "axial_ratio_dB": AR_dB,
            "helicity": helicity,
            "ieee_cp_pass": bool(ieee_pass),
        },
        "per_point_stats": {
            "eta_cp_mean_pct": eta_mean,
            "eta_cp_std_pct": eta_std,
            "eta_cp_min_pct": eta_min,
            "eta_cp_max_pct": eta_max,
        }
    }
    out_path = args.output or f"stokes_R{int(R*1000)}mm.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    if not args.quiet:
        print(f"\nWritten: {out_path}")


if __name__ == "__main__":
    main()
