#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multi-variant polarization validation suite.

Runs postprocess_stokes_polarization.py on multiple variant directories
and compares the results against the values published in the README.

Usage:
    python validate_polarization_claims.py --repo-root /path/to/open-chiral-flux-shaper
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

CLAIMS = [
    {
        "variant": "gabbia_sferica_doppio_gruppo_90deg_48coils",
        "name": "Dual Orthogonal 90 deg (48 Coils)",
        "radius_mm": 55,
        "expected_eta_cp_pct": 95.5,
        "expected_AR_dB": 2.67,
        "expected_s3_sign": "positive",
        "reference": "README sec 2.2 #1, Fig. 34",
    },
    {
        "variant": "gabbia_sferica_chiral_diode_asymmetric_pulse",
        "name": "Asymmetric Gradient Mantle (Chiral Mode Shaper)",
        "radius_mm": 55,
        "expected_eta_cp_pct": 99.98,
        "expected_AR_dB": 0.15,
        "expected_s3_sign": "positive",
        "reference": "README sec 2.2 #2, Fig. 35",
    },
    {
        "variant": "rotore_toroidale_verticale_2bobine_vertice",
        "name": "Apex Magnetic Cusp (2-Bobine Toroidal)",
        "radius_mm": 47,
        "expected_eta_cp_pct": 97.2,
        "expected_AR_dB": 15.43,
        "expected_s3_sign": "any",
        "reference": "README sec 2.2 #5, Fig. 47",
    },
]

TOLERANCE_PCT = 5.0


def find_vtu_glob(repo_root, variant):
    candidates = [
        repo_root / "variants" / variant / "work_dirs" / "run_*/results" / "macchina_out_t*.vtu",
        repo_root / "variants" / variant / "work_dirs" / "*" / "results" / "macchina_out_t*.vtu",
        repo_root / "variants" / variant / "results" / "macchina_out_t*.vtu",
        repo_root / "variants" / variant / "results_*" / "macchina_out_t*.vtu",
    ]
    for c in candidates:
        import glob
        matches = sorted(glob.glob(str(c)))
        if matches:
            return str(c)
    return None


def run_stokes_tool(tool_path, vtu_glob, radius_m, frequency, dt, output_path):
    cmd = [
        sys.executable, str(tool_path),
        "--vtu-glob", vtu_glob,
        "--radius", str(radius_m),
        "--frequency", str(frequency),
        "--dt", str(dt),
        "--output", str(output_path),
        "--quiet",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return None, result.stderr
    with open(output_path, encoding="utf-8") as f:
        return json.load(f), None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=Path(__file__).resolve().parent.parent, type=Path,
                   help="Path to open-chiral-flux-shaper repository root")
    p.add_argument("--tool", default=Path(__file__).parent / "postprocess_stokes_polarization.py",
                   type=Path, help="Path to postprocess_stokes_polarization.py")
    p.add_argument("--frequency", type=float, default=100.0)
    p.add_argument("--dt", type=float, default=5e-4)
    p.add_argument("--output-dir", default=None, type=Path)
    args = p.parse_args()

    out_dir = args.output_dir or (args.repo_root / "data" / "stokes_validation")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("POLARIZATION CLAIMS VALIDATION SUITE")
    print("=" * 90)
    print(f"{'Variant':<48} | {'Claim eta_CP':<12} | {'Our eta_CP':<12} | {'Delta':<8} | Status")
    print("-" * 90)

    results = []
    for claim in CLAIMS:
        vtu_glob = find_vtu_glob(args.repo_root, claim["variant"])
        if vtu_glob is None:
            print(f"{claim['name']:<48} | {claim['expected_eta_cp_pct']:>10.2f}% | "
                  f"{'N/A':>10} | {'N/A':>6} | NO VTU FILES")
            results.append({**claim, "status": "no_vtu_files"})
            continue

        out_path = out_dir / f"stokes_{claim['variant']}_R{claim['radius_mm']}mm.json"
        data, err = run_stokes_tool(
            args.tool, vtu_glob, claim["radius_mm"] / 1000.0,
            args.frequency, args.dt, out_path
        )
        if data is None:
            print(f"{claim['name']:<48} | {claim['expected_eta_cp_pct']:>10.2f}% | "
                  f"{'ERR':>10} | {'N/A':>6} | TOOL ERROR")
            print(f"  {err[:200]}")
            results.append({**claim, "status": "tool_error", "error": err})
            continue

        our_eta = data["metrics"]["circular_purity_pct"]
        our_ar = data["metrics"]["axial_ratio_dB"]
        our_s3 = data["stokes_parameters"]["s3"]
        our_helicity = data["metrics"]["helicity"]

        delta = abs(our_eta - claim["expected_eta_cp_pct"])
        validated = delta <= TOLERANCE_PCT
        status = "VALIDATED" if validated else "MISMATCH"

        print(f"{claim['name']:<48} | {claim['expected_eta_cp_pct']:>10.2f}% | "
              f"{our_eta:>10.2f}% | {delta:>+6.2f} | {status}")
        print(f"  AR: claim={claim['expected_AR_dB']:.2f} dB, ours={our_ar:.2f} dB | "
              f"s3={our_s3:+.4e} ({our_helicity[:30]})")

        results.append({
            **claim,
            "our_eta_cp_pct": our_eta,
            "our_AR_dB": our_ar,
            "our_s3": our_s3,
            "delta_pct": delta,
            "validated": validated,
            "output_json": str(out_path),
        })

    print("=" * 90)
    n_validated = sum(1 for r in results if r.get("validated"))
    n_total = len(results)
    print(f"\nValidated: {n_validated}/{n_total} claims (tolerance +/-{TOLERANCE_PCT}% absolute)")

    summary_path = out_dir / "validation_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "tool": str(args.tool),
            "repo_root": str(args.repo_root),
            "frequency_hz": args.frequency,
            "dt_s": args.dt,
            "tolerance_pct": TOLERANCE_PCT,
            "n_validated": n_validated,
            "n_total": n_total,
            "results": results,
        }, f, indent=2)
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()
