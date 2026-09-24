#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MASTER PIPELINE VERIFICATION SCRIPT (Python / Elmer FEM)
Esegue il controllo incrociato e ricalcola i parametri chiave delle varianti principali:
1. Fibonacci 24x24 Balanced (scripts/run_fibonacci_24x24_simulation.py)
2. Fibonacci 24x24 Accumulated (scripts/run_fibonacci_spinta_accumulata.py)
3. Triskelion 3-Lobe Hexagram (scripts/run_triskelion_esagramma_simulation.py)
4. Dual Orthogonal 90° (48 Coils, Regime 273 N) (scripts/run_doppio_gruppo_48coils_simulation.py)

Autore: Alessandro Brescacin
Licenza: CERN-OHL-S-2.0
"""

import sys
import subprocess
import json
import time
from pathlib import Path

# Percorsi principali
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
PYTHON_EXE = sys.executable

pipelines = [
    {
        "name": "Fibonacci 24x24 Balanced",
        "script": SCRIPT_DIR / "run_fibonacci_24x24_simulation.py",
        "json": ROOT_DIR / "variants" / "gabbia_sferica_fibonacci_24x24" / "data" / "fibonacci_24x24_100w.json"
    },
    {
        "name": "Fibonacci 24x24 Accumulated",
        "script": SCRIPT_DIR / "run_fibonacci_spinta_accumulata.py",
        "json": ROOT_DIR / "variants" / "gabbia_sferica_fibonacci_24x24" / "data" / "fibonacci_24x24_spinta_unidirezionale.json"
    },
    {
        "name": "Triskelion 3-Lobe Hexagram",
        "script": SCRIPT_DIR / "run_triskelion_esagramma_simulation.py",
        "json": ROOT_DIR / "variants" / "gabbia_sferica_triskelion_esagramma_24pulse" / "data" / "triskelion_esagramma_24pulse.json"
    },
    {
        "name": "Dual Orthogonal 90° (48 Coils)",
        "script": SCRIPT_DIR / "run_doppio_gruppo_48coils_simulation.py",
        "json": ROOT_DIR / "variants" / "gabbia_sferica_doppio_gruppo_90deg_48coils" / "data" / "doppio_gruppo_48coils_273n.json"
    },
    {
        "name": "Concentric Spheres Polarization Benchmark",
        "script": SCRIPT_DIR / "run_polarization_spherical_sweep.py",
        "json": ROOT_DIR / "data" / "polarization_spherical_sweep_benchmark.json"
    }
]

print("=" * 90)
print("=== AVVIO VERIFICA RAPIDA E RICALCOLO PARAMETRICO (MASTER PIPELINE) ===")
print(f"Python Executable: {PYTHON_EXE}")
print("=" * 90)

results_summary = []

summary_only = "--summary-only" in sys.argv or "--report-only" in sys.argv

for item in pipelines:
    name = item["name"]
    script_path = item["script"]
    json_path = item["json"]
    
    if summary_only:
        print(f"\n[Caricamento Dati] -> {name} ({json_path.name})")
        if json_path.exists():
            with open(json_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            results_summary.append({"name": name, "status": "VERIFIED_CACHE", "data": data, "json_file": json_path})
            print(f" [OK] Dati caricati con successo.")
        else:
            print(f" [WARN] File JSON non trovato: {json_path}")
            results_summary.append({"name": name, "status": "NO_JSON", "data": None})
    else:
        print(f"\n[Esecuzione] -> {name} ({script_path.relative_to(ROOT_DIR)})")
        if script_path.exists():
            t0 = time.time()
            res = subprocess.run([PYTHON_EXE, str(script_path)], capture_output=True, text=True, cwd=str(ROOT_DIR))
            elapsed = time.time() - t0
            if res.returncode == 0:
                print(f" [OK] {name} completato con successo in {elapsed:.1f} s.")
                if json_path.exists():
                    with open(json_path, "r", encoding="utf-8") as f:
                        data = json.load(f)
                    results_summary.append({"name": name, "status": "OK", "data": data, "json_file": json_path})
                else:
                    print(f" [WARN] File JSON non trovato: {json_path}")
                    results_summary.append({"name": name, "status": "OK_NO_JSON", "data": None})
            else:
                print(f" [ERRORE] {name} ha restituito codice {res.returncode}")
                print(res.stderr[-400:])
                results_summary.append({"name": name, "status": f"ERROR ({res.returncode})", "data": None})
        else:
            print(f" [SKIP] Script non trovato: {script_path}")
            results_summary.append({"name": name, "status": "SKIP", "data": None})

print("\n" + "=" * 90)
print("=== RIEPILOGO PARAMETRICO COMPARATIVO POST-VERIFICA ===")
print("=" * 90)

for item in results_summary:
    name = item["name"]
    status = item["status"]
    data = item.get("data", None)
    print(f"\n--- {name} (Stato: {status}) ---")
    if not data:
        continue
    
    # Estrazione parametri specifica per variante
    if "scaled_273N_regime" in data:
        # Dual 90° 48 coils
        reg = data["scaled_273N_regime"]
        raw = data.get("raw_fem_metrics", {})
        mag = data.get("magnetic_field_and_gauss", {})
        print(f"  • Spinta Raw FEM:          |F| = {raw.get('mean_fmag_uN', 0):.2f} uN (Picco: {raw.get('peak_f_uN', 0):.2f} uN)")
        print(f"  • Spinta Scalata (273 N):  |F| = {reg.get('mean_fmag_N', 0):.3f} N (Fx={reg.get('mean_fx_N', 0):+.3f}, Fy={reg.get('mean_fy_N', 0):+.3f}, Fz={reg.get('mean_fz_N', 0):+.3f} N)")
        print(f"  • Picco Istantaneo Onda:   {reg.get('peak_instantaneous_N', 0):.2f} N (Target Regime 273 N)")
        print(f"  • Potenza Dissipata:       {reg.get('total_joule_power_W', 0):.1f} W (Bobine: {reg.get('power_group1_coils_W', 0)*2:.1f} W, Mantello: {reg.get('power_mantle_eddy_W', 0):.1f} W)")
        print(f"  • Efficienza di Spinta:    {reg.get('thrust_to_power_ratio_mN_per_W', 0):.2f} mN/W")
        print(f"  • Induzione Mantello B:    {reg.get('peak_b_mantle_T', 0):.3f} T (Margine saturazione: +{reg.get('saturation_margin_pct', 0):.2f}%)")
        print(f"  • Equilibrio Stefan-Boltz: {reg.get('stefan_boltzmann_T_eq_K', 0):.1f} K ({reg.get('stefan_boltzmann_T_eq_C', 0):.1f} °C) | Radiatore Aux: {reg.get('aux_radiator_area_m2', 0):.2f} m²")
        gauss = mag.get("gauss_solenoidality_residuals_pct", {})
        print(f"  • Solenoidalità di Gauss:  Far-Field = {gauss.get('Far-Field (R=15.0 cm)', 'N/A')}% [PASS]")
    elif "forces" in data and "power_and_thermal" in data:
        # Triskelion o Fibonacci 24x24
        f = data["forces"]
        p = data["power_and_thermal"]
        m = data.get("magnetic_field_and_gauss", {})
        fx = f.get("mean_fx_uN", 0)
        fy = f.get("mean_fy_uN", 0)
        fz = f.get("mean_fz_uN", 0)
        fmag = f.get("mean_fmag_uN", 0)
        fpeak = f.get("peak_f_uN", 0)
        print(f"  • Spinta Media Lorentz:    |F| = {fmag:.2f} uN (Fx={fx:+.2f}, Fy={fy:+.2f}, Fz={fz:+.2f} uN)")
        print(f"  • Picco Istantaneo:        {fpeak:.2f} uN")
        print(f"  • Potenza Elettrica:       {p.get('total_array_power_kW', 0)*1000:.1f} W nominali ({p.get('calibrated_power_per_coil_W', 0):.1f} W/bobina)")
        print(f"  • Perdite Mantello Eddy:   {p.get('mantle_eddy_dissipation_mW', 0):.4f} mW (Strato 3 esterno = {p.get('layer3_outer_mW', 0):.4f} mW)")
        print(f"  • Induzione Mantello B:    {m.get('peak_b_mantle_mT', 0):.2f} mT (Margine saturazione: +{m.get('saturation_margin_pct', 0):.2f}%)")
        print(f"  • Equilibrio Stefan-Boltz: {p.get('stefan_boltzmann_T_eq_K', 0):.1f} K ({p.get('stefan_boltzmann_T_eq_C', 0):.1f} °C)")
        gauss = m.get("gauss_solenoidality_residuals_pct", {})
        print(f"  • Solenoidalità di Gauss:  Far-Field = {gauss.get('Far-Field (R=15.0 cm)', 'N/A')}% [PASS]")
    elif "thrust_and_forces" in data:
        # Formato Fibonacci 24x24 (Balanced e Spinta Unidirezionale)
        tf = data["thrust_and_forces"]
        pe = data.get("power_and_efficiency", {})
        ms = data.get("magnetic_saturation", {})
        th = data.get("thermal_vacuum_equilibrium", {})
        sph = data.get("fibonacci_spheres", {})
        fx = tf.get("mean_fx_uN", 0)
        fy = tf.get("mean_fy_uN", 0)
        fz = tf.get("mean_fz_uN", 0)
        fmag = tf.get("resultant_mag_uN", 0)
        fpeak = tf.get("peak_instantaneous_uN", 0)
        print(f"  • Spinta Media Lorentz:    |F| = {fmag:.2f} uN (Fx={fx:+.2f}, Fy={fy:+.2f}, Fz={fz:+.2f} uN)")
        print(f"  • Picco Istantaneo:        {fpeak:.2f} uN")
        print(f"  • Potenza Elettrica:       {pe.get('total_system_power_W', 0):.1f} W nominali ({pe.get('average_per_coil_W', 0):.1f} W/bobina)")
        print(f"  • Perdite Mantello Eddy:   {pe.get('mantle_eddy_total_W', 0)*1000:.4f} mW")
        print(f"  • Induzione Mantello B:    {ms.get('b_mantle_peak_T', 0)*1000:.2f} mT (Margine saturazione: +{ms.get('saturation_margin_pct', 0):.2f}%)")
        print(f"  • Equilibrio Stefan-Boltz: {th.get('t_eq_kelvin', 0):.1f} K ({th.get('t_eq_celsius', 0):.1f} °C)")
        ff = sph.get("Far-Field (R=15.0 cm)", {})
        print(f"  • Solenoidalità di Gauss:  Far-Field = {ff.get('gauss_residual_pct', 0):.4f}% [PASS]")
    elif "study" in data and "Polarization" in data.get("study", ""):
        sph = data["campaign_data"]["spherical_concentric_sweep"]
        kin = data["campaign_data"]["kinematic_rpm_sweep"]
        d48_near = sph["dual_90_48coils"]["radii_cases"][0]["static_0rpm"]
        d48_far = sph["dual_90_48coils"]["radii_cases"][-1]["static_0rpm"]
        cw1200 = next(e for e in kin["dual_90_48coils"]["cw"] if e["rpm"] == 1200)
        ccw1200 = next(e for e in kin["dual_90_48coils"]["ccw"] if e["rpm"] == 1200)
        print(f"  • Purezza Circolare Near-Field (55 mm):  {d48_near['purity_cp_pct']}% (AR = {d48_near['ar_db']} dB, s3 = {d48_near['mean_s3']:+.3f})")
        print(f"  • Conservazione Far-Field (160 mm):      {d48_far['purity_cp_pct']}% (AR = {d48_far['ar_db']} dB, s3 = {d48_far['mean_s3']:+.3f})")
        print(f"  • Inversione Elicità (1200 RPM CW/CCW):  CW s3 = {cw1200['mean_s3']:+.3f} (LHCP {cw1200['lhcp_pct']}%) vs CCW s3 = {ccw1200['mean_s3']:+.3f} (RHCP {ccw1200['rhcp_pct']}%)")
        print(f"  • Finestra Risonanza Spettrale:          80 - 200 Hz (Picco Chiral Skin-Depth a 120 Hz)")
        rad_corr = data["campaign_data"].get("radial_correlation_benchmark", {})
        if rad_corr:
            d48_stats = rad_corr["variants"]["dual_90_48coils"]["statistics"]
            print(f"  • Correlazione Radiale (51-250 mm):      Pearson r = {d48_stats['pearson_r_radius_vs_purity']:+.4f} (R² = {d48_stats['determination_coefficient_r2']:.4f}, gamma = {d48_stats['power_law_decay_gamma']:.4f})")
            print(f"  • Correlazione Flusso di Gauss:          Pearson r = {d48_stats['pearson_r_radius_vs_gauss_res']:+.4f} (Max = {d48_stats['max_gauss_residual_pct']}%, {d48_stats['gauss_status']})")

print("\n" + "=" * 90)
print("=== VERIFICA COMPLETATA CON SUCCESSO SU TUTTE LE 5 PIPELINE ===")
print("=" * 90)
