#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASE 1: Sanity Check contro il Mesh Bias, Inversione Chirale (-30°)
e Validazione con Tensore degli Sforzi di Maxwell (MST) su Superficie Chiusa.

Protocollo:
1. Esecuzione simulazione speculare theta = -30° al punto nominale (100 Hz, 1200 RPM, 20 timesteps, dt=0.5 ms).
2. Calcolo F_bias = (<F_z(+30°)> + <F_z(-30°)>) / 2 e spinta netta corretta F_chiral = (<F_z(+30°)> - <F_z(-30°)>) / 2.
3. Integrazione superficiale del Tensore di Maxwell (MST) su cilindro chiuso in aria (R=8 cm, H=+-8 cm):
   T_z = (1/mu0) * [ B_z(B · n) - 0.5 * |B|^2 * n_z ]
   F_{z,MST} = ∮ T_z dA
4. Confronto percentuale tra Lorentz volumetrico e MST.

Autore: Alessandro Brescacin
Licenza: CERN-OHL-S-2.0
"""

import os
import sys
import json
import time
import shutil
import subprocess
from pathlib import Path
import numpy as np
import meshio
from scipy.interpolate import LinearNDInterpolator

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).resolve().parent
SWEEP_DIR = SCRIPT_DIR.parent
ROOT_DIR = SWEEP_DIR.parent.parent
CONFIG_DIR = SWEEP_DIR / "config"
DATA_DIR = SWEEP_DIR / "data"
FIGURES_DIR = SWEEP_DIR / "figures"
WORK_DIR = SWEEP_DIR / "work_dirs" / "run_f100_rpm1200_minus30deg"
MESH_DB_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0" / "mesh"

CONFIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
WORK_DIR.mkdir(parents=True, exist_ok=True)

MU0 = 4.0 * np.pi * 1e-7
TIMESTEPS = 20
DT = 0.0005
F_HZ = 100
RPM = 1200

def find_elmersolver():
    cmd = shutil.which("ElmerSolver")
    if cmd:
        return cmd
    candidates = [
        r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe",
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerSolver.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerSolver.exe"
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return "ElmerSolver"

def generate_sif_minus30(mesh_db_posix, res_dir_posix):
    w_e = 2.0 * np.pi * F_HZ
    w_m = 2.0 * np.pi * (RPM / 60.0)
    
    sif = f"""! Elmer FEM: Specular Control Run theta = -30 deg (100 Hz, 1200 RPM)
! Physical Parity & Mesh Bias Falsification Check
$ function rotating_current_regime_B(tx) {{\\
  x = tx(0); y = tx(1); z = tx(2); t = tx(3);\\
  res = 0.0;\\
  if (abs(z) > 0.05) {{ res = 0.0; }} else {{\\
    rsq = x*x + y*y;\\
    if (rsq < 0.0007) {{ res = 0.0; }} else {{ if (rsq > 0.0019) {{ res = 0.0; }} else {{\\
      w = {w_m:.14f};\\
      rc = 0.035; rwsq = 4.9e-05; j0 = 1.0e+05;\\
      th0 = w * t;\\
      c0 = cos(th0); s0 = sin(th0);\\
      s60 = 0.8660254037844386;\\
      x1 = rc * c0; y1 = rc * s0;\\
      x2 = rc * (0.5*c0 - s60*s0); y2 = rc * (0.5*s0 + s60*c0);\\
      x3 = rc * (-0.5*c0 - s60*s0); y3 = rc * (-0.5*s0 + s60*c0);\\
      se = sin({w_e:.14f} * t); ce = cos({w_e:.14f} * t);\\
      s1 = se; s2 = 0.5*se - s60*ce; s3 = -0.5*se - s60*ce;\\
      d1 = (x-x1)*(x-x1) + (y-y1)*(y-y1);\\
      d2 = (x-x2)*(x-x2) + (y-y2)*(y-y2);\\
      d3 = (x-x3)*(x-x3) + (y-y3)*(y-y3);\\
      d4 = (x+x1)*(x+x1) + (y+y1)*(y+y1);\\
      d5 = (x+x2)*(x+x2) + (y+y2)*(y+y2);\\
      d6 = (x+x3)*(x+x3) + (y+y3)*(y+y3);\\
      if (d1 < rwsq) {{ res = j0 * 0.5 * (s1 + abs(s1)); }} else {{\\
      if (d2 < rwsq) {{ res = j0 * 0.5 * (s2 + abs(s2)); }} else {{\\
      if (d3 < rwsq) {{ res = j0 * 0.5 * (s3 + abs(s3)); }} else {{\\
      if (d4 < rwsq) {{ res = j0 * 0.5 * (abs(s1) - s1); }} else {{\\
      if (d5 < rwsq) {{ res = j0 * 0.5 * (abs(s2) - s2); }} else {{\\
      if (d6 < rwsq) {{ res = j0 * 0.5 * (abs(s3) - s3); }} else {{ res = 0.0; }}; }}; }}; }}; }}; }};\\
    }}; }};\\
  }};\\
  _rotating_current_regime_B = res;\\
}}

! Inversione Chirale: s_tz negativo (-3.031089e6 S/m)
$ function sigma_cyl(tx) {{\\
  x = tx(0); y = tx(1);\\
  r = sqrt(x*x + y*y);\\
  if (r < 1.0e-5) {{ r = 1.0e-5; }} else {{ r = r; }};\\
  c = x / r; s = y / r;\\
  s_rr = 1.75e6; s_tt = 1.75e6; s_zz = 1.22e7; s_tz = -3.031089e6;\\
  sxz = -s_tz * s;\\
  syz =  s_tz * c;\\
  a = zeros(3,3);\\
  a(0,0) = s_rr; a(1,1) = s_tt; a(2,2) = s_zz;\\
  a(0,2) = sxz;  a(2,0) = sxz;\\
  a(1,2) = syz;  a(2,1) = syz;\\
  _sigma_cyl = a;\\
}}

Header
  CHECK KEYWORDS Warn
  Mesh DB "{mesh_db_posix}" "macchina_centrata"
  Include Path ""
  Results Directory "{res_dir_posix}"
End

Simulation
  Max Output Level = 3
  Coordinate System = Cartesian 3D
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Transient
  Steady State Max Iterations = 1
  Output Intervals(1) = 1
  Timestep Intervals(1) = {TIMESTEPS}
  Timestep Sizes(1) = {DT:.8e}
End

Constants
  Permittivity of Vacuum = 8.8542e-12
  Permeability of Vacuum = 1.2566370614e-6
End

Equation 1
  Name = "Coupled Electromagnetics"
  Active Solvers(3) = 1 2 3
End

Solver 1
  Equation = "MGDynamics"
  Variable = "AV"
  Procedure = "MagnetoDynamics" "WhitneyAVSolver"
  Fix Input Current Density = Logical False
  Use Piola Transform = Logical False
  Linear System Solver = "Direct"
  Linear System Direct Method = Umfpack
End

Solver 2
  Equation = "MGDynamicsCalc"
  Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"
  Potential Variable = String "AV"
  Calculate Elemental Fields = Logical True
  Calculate Nodal Fields = Logical True
  Calculate Magnetic Flux Density = Logical True
  Calculate Magnetic Vector Potential = Logical True
  Calculate Current Density = Logical True
  Calculate Electric Field = Logical True
  Calculate Joule Heating = Logical True
  Calculate JxB = Logical True
  Linear System Solver = "Direct"
  Linear System Direct Method = Umfpack
End

Solver 3
  Equation = "ResultOutput"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = "macchina_out"
  Output Format = "vtu"
  Vtu Format = Logical True
  Save Geometry Ids = Logical True
  Discontinuous Bodies = Logical True
End

Material 1
  Name = "Air"
  Relative Permeability = 1.0
  Electric Conductivity = 0.0
End

Material 2
  Name = "ExpandedMeshMantle"
  Relative Permeability = 1.0
  Electric Conductivity(3,3) = Variable Coordinate 1, Coordinate 2
    Real MATC "sigma_cyl(tx)"
End

Material 3
  Name = "FerromagneticCore"
  Relative Permeability = 1000.0
  Electric Conductivity = 0.0
End

Body Force 1
  Current Density 3 = Variable Coordinate 1, Coordinate 2, Coordinate 3, time
    Real MATC "rotating_current_regime_B(tx)"
End

Body 1
  Target Bodies(1) = 1
  Name = "AirExterior"
  Equation = 1
  Material = 1
End

Body 2
  Target Bodies(1) = 2
  Name = "Aluminum"
  Equation = 1
  Material = 2
End

Body 3
  Target Bodies(1) = 3
  Name = "RotorAir"
  Equation = 1
  Material = 1
  Body Force = 1
End

Body 4
  Target Bodies(1) = 4
  Name = "FerromagneticCore"
  Equation = 1
  Material = 3
  Body Force = 1
End

Boundary Condition 1
  Target Boundaries(1) = 1
  Name = "FarFieldInfinity"
  AV {{e}} = Real 0.0
  AV = Real 0.0
End
"""
    return sif

def run_simulation(elmer_bin):
    res_dir = WORK_DIR / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    existing_vtus = sorted(list(res_dir.glob("macchina_out_t*.vtu")))
    
    if len(existing_vtus) >= TIMESTEPS and "--force" not in sys.argv:
        print(f"  [CACHE] Trovati gia {len(existing_vtus)} file VTU completi in {res_dir}. Skip solver.")
        return True

    sif_content = generate_sif_minus30(MESH_DB_DIR.as_posix(), res_dir.as_posix())
    (CONFIG_DIR / "case_100Hz_1200RPM_minus30deg.sif").write_text(sif_content, encoding="utf-8")
    (WORK_DIR / "case.sif").write_text(sif_content, encoding="utf-8")
    (WORK_DIR / "ELMERSOLVER_STARTINFO").write_text("case.sif\n1\n", encoding="utf-8")

    print(f"  [SOLVER] Avvio ElmerSolver per {TIMESTEPS} timestep (100 Hz, 1200 RPM, theta=-30°)...")
    t0 = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    proc = subprocess.run([elmer_bin, "case.sif"], cwd=str(WORK_DIR), env=env, capture_output=True, text=True)
    duration = time.time() - t0
    
    if proc.returncode != 0:
        print(f"  [ERRORE] ElmerSolver fallito con codice {proc.returncode}")
        print("STDERR tail:\n", proc.stderr[-1000:])
        return False
        
    vtus = list(res_dir.glob("macchina_out_t*.vtu"))
    print(f"  [SOLVER] Simulazione completata in {duration:.1f}s ({len(vtus)} VTU generati).")
    return True

def precalculate_mesh_and_mst_grid():
    ref_vtu = ROOT_DIR / "variants" / "rotore_centrato_z0" / "results_regime_B" / "macchina_out_t0001.vtu"
    m_ref = meshio.read(str(ref_vtu))
    pts = m_ref.points
    cells = m_ref.cells_dict['tetra']
    v0, v1, v2, v3 = pts[cells[:, 0]], pts[cells[:, 1]], pts[cells[:, 2]], pts[cells[:, 3]]
    elem_vols = np.abs(np.einsum('ij,ij->i', v1 - v0, np.cross(v2 - v0, v3 - v0))) / 6.0
    elem_com = (v0 + v1 + v2 + v3) / 4.0
    r_com = np.hypot(elem_com[:, 0], elem_com[:, 1])
    z_com = np.abs(elem_com[:, 2])

    active_mask = (r_com <= 0.051) & (z_com <= 0.051)
    mantle_mask = (r_com >= 0.046) & (r_com <= 0.051) & (z_com <= 0.051)

    # Griglia superficiale per Tensore di Maxwell MST (Cilindro chiuso R=8 cm, H=+-8 cm)
    Rc = 0.08
    Hc = 0.08
    nth = 64
    nz = 32
    nr = 24

    th = np.linspace(0, 2*np.pi, nth, endpoint=False)
    z = np.linspace(-Hc, Hc, nz)
    TH, Z = np.meshgrid(th, z)
    cos_th = np.cos(TH).ravel()
    sin_th = np.sin(TH).ravel()
    dA_lat = (2.0 * np.pi * Rc / nth) * (2.0 * Hc / (nz - 1))
    pts_lat = np.column_stack([(Rc * np.cos(TH)).ravel(), (Rc * np.sin(TH)).ravel(), Z.ravel()])

    r_arr = np.linspace(0.002, Rc, nr)
    R_top, TH_top = np.meshgrid(r_arr, th)
    xt = (R_top * np.cos(TH_top)).ravel()
    yt = (R_top * np.sin(TH_top)).ravel()
    pts_top = np.column_stack([xt, yt, np.full_like(xt, Hc)])
    pts_bot = np.column_stack([xt, yt, np.full_like(xt, -Hc)])
    dr = r_arr[1] - r_arr[0]
    dth = 2.0 * np.pi / nth
    dA_top = (R_top * dr * dth).ravel()

    mst_grid = {
        "Rc": Rc, "Hc": Hc,
        "pts_lat": pts_lat, "dA_lat": dA_lat, "cos_th": cos_th, "sin_th": sin_th,
        "pts_top": pts_top, "pts_bot": pts_bot, "dA_top": dA_top
    }

    return pts, cells, elem_vols, active_mask, mantle_mask, mst_grid

def evaluate_series(vtus, pts, cells, elem_vols, active_mask, mantle_mask, mst_grid):
    fz_lorentz_active = []
    fz_lorentz_mantle = []
    pj_mantle = []
    fz_mst_list = []

    pts_lat = mst_grid["pts_lat"]
    dA_lat = mst_grid["dA_lat"]
    cos_th = mst_grid["cos_th"]
    sin_th = mst_grid["sin_th"]
    pts_top = mst_grid["pts_top"]
    pts_bot = mst_grid["pts_bot"]
    dA_top = mst_grid["dA_top"]

    for idx, vf in enumerate(vtus):
        m = meshio.read(str(vf))
        jxb_z = m.point_data['jxb'][:, 2]
        pj = m.point_data['joule heating'].ravel()
        B = m.point_data['magnetic flux density']

        jxb_elem = np.mean(jxb_z[cells], axis=1)
        pj_elem = np.mean(pj[cells], axis=1)

        f_l_act = float(np.sum(elem_vols[active_mask] * jxb_elem[active_mask]))
        f_l_man = float(np.sum(elem_vols[mantle_mask] * jxb_elem[mantle_mask]))
        p_j_man = float(np.sum(elem_vols[mantle_mask] * pj_elem[mantle_mask]))

        fz_lorentz_active.append(f_l_act)
        fz_lorentz_mantle.append(f_l_man)
        pj_mantle.append(p_j_man)

        # Interpolazione per Tensore di Maxwell MST
        interp_Bx = LinearNDInterpolator(pts, B[:, 0], fill_value=0.0)
        interp_By = LinearNDInterpolator(pts, B[:, 1], fill_value=0.0)
        interp_Bz = LinearNDInterpolator(pts, B[:, 2], fill_value=0.0)

        # Laterale
        Bxl, Byl, Bzl = interp_Bx(pts_lat), interp_By(pts_lat), interp_Bz(pts_lat)
        Brl = Bxl * cos_th + Byl * sin_th
        fz_lat = np.sum((1.0 / MU0) * Bzl * Brl) * dA_lat

        # Top
        Bxt, Byt, Bzt = interp_Bx(pts_top), interp_By(pts_top), interp_Bz(pts_top)
        fz_top = np.sum((0.5 / MU0) * (Bzt**2 - Bxt**2 - Byt**2) * dA_top)

        # Bottom
        Bxb, Byb, Bzb = interp_Bx(pts_bot), interp_By(pts_bot), interp_Bz(pts_bot)
        fz_bot = np.sum(-(0.5 / MU0) * (Bzb**2 - Bxb**2 - Byb**2) * dA_top)

        fz_mst_list.append(float(fz_lat + fz_top + fz_bot))

    return {
        "fz_lorentz_active_uN": [x * 1e6 for x in fz_lorentz_active],
        "fz_lorentz_mantle_uN": [x * 1e6 for x in fz_lorentz_mantle],
        "pj_mantle_mW": [x * 1000.0 for x in pj_mantle],
        "fz_mst_uN": [x * 1e6 for x in fz_mst_list],
        "mean_Fz_active_uN": float(np.mean(fz_lorentz_active) * 1e6),
        "mean_Fz_mantle_uN": float(np.mean(fz_lorentz_mantle) * 1e6),
        "mean_Pj_mW": float(np.mean(pj_mantle) * 1000.0),
        "mean_Fz_mst_uN": float(np.mean(fz_mst_list) * 1e6)
    }

def main():
    print("=" * 85)
    print("FASE 1: SANITY CHECK MESH BIAS E TENSORE DEGLI SFORZI DI MAXWELL (MST)")
    print("Confronto Nominale (+30°) vs Speculare (-30°) a 100 Hz, 1200 RPM")
    print("=" * 85)

    elmer_bin = find_elmersolver()
    print(f"Solutore: {elmer_bin}")

    pts, cells, elem_vols, active_mask, mantle_mask, mst_grid = precalculate_mesh_and_mst_grid()
    print(f"Mesh caricata con successo: {len(elem_vols)} elementi.")

    # 1. Esecuzione simulazione speculare theta = -30°
    ok = run_simulation(elmer_bin)
    if not ok:
        print("[ERRORE] Impossibile completare la simulazione theta = -30°.")
        sys.exit(1)

    # 2. Caricamento VTU
    vtus_plus30 = sorted(list((ROOT_DIR / "variants" / "rotore_centrato_z0" / "results_regime_B").glob("*.vtu")))[:TIMESTEPS]
    vtus_minus30 = sorted(list((WORK_DIR / "results").glob("macchina_out_t*.vtu")))[:TIMESTEPS]

    print("\n[ANALISI] Elaborazione serie temporale per Baseline (+30°)...")
    res_plus30 = evaluate_series(vtus_plus30, pts, cells, elem_vols, active_mask, mantle_mask, mst_grid)

    print("\n[ANALISI] Elaborazione serie temporale per Speculare (-30°)...")
    res_minus30 = evaluate_series(vtus_minus30, pts, cells, elem_vols, active_mask, mantle_mask, mst_grid)

    # 3. Calcolo F_bias e F_chiral
    F_plus = res_plus30["mean_Fz_active_uN"]
    F_minus = res_minus30["mean_Fz_active_uN"]

    F_bias = (F_plus + F_minus) / 2.0
    F_chiral = (F_plus - F_minus) / 2.0

    # Mantello soltanto
    F_man_plus = res_plus30["mean_Fz_mantle_uN"]
    F_man_minus = res_minus30["mean_Fz_mantle_uN"]
    F_man_bias = (F_man_plus + F_man_minus) / 2.0
    F_man_chiral = (F_man_plus - F_man_minus) / 2.0

    # MST vs Lorentz
    mst_plus = res_plus30["mean_Fz_mst_uN"]
    diff_mst_lorentz = abs(mst_plus - F_plus)
    err_mst_pct = diff_mst_lorentz / (abs(F_plus) + 1e-12) * 100.0

    print("\n" + "=" * 85)
    print("SINTESI QUANTITATIVA: DECOUPLING BIAS DI MESH E LIFT CHIRALE REALE")
    print("=" * 85)
    print(f"Spinta Nominale (+30°):         <F_z(+30°)> = {F_plus:+8.3f} µN")
    print(f"Spinta Speculare (-30°):        <F_z(-30°)> = {F_minus:+8.3f} µN")
    print(f"---------------------------------------------------------------------")
    print(f"Componente Bias Asimmetrico:    F_bias      = {F_bias:+8.3f} µN")
    print(f"Spinta Chirale Pura Corretta:   F_chiral    = {F_chiral:+8.3f} µN")
    print(f"---------------------------------------------------------------------")
    print(f"Mantello (+30°):                <F_z,man>   = {F_man_plus:+8.3f} µN")
    print(f"Mantello (-30°):                <F_z,man>   = {F_man_minus:+8.3f} µN")
    print(f"Mantello Bias:                  F_man_bias  = {F_man_bias:+8.3f} µN")
    print(f"Mantello Chiral Pure:           F_man_chir  = {F_man_chiral:+8.3f} µN")
    print(f"---------------------------------------------------------------------")
    print(f"Tensore di Maxwell (MST +30°):  F_z,MST     = {mst_plus:+8.3f} µN")
    print(f"Differenza MST vs Lorentz:      Delta F_z   = {diff_mst_lorentz:8.3f} µN")
    print("=" * 85)

    # Salvataggio JSON
    output_data = {
        "metadata": {
            "title": "Fase 1: Sanity Check Mesh Bias e Validazione Tensore di Maxwell (MST)",
            "author": "Alessandro Brescacin",
            "license": "CERN-OHL-S-2.0",
            "date": "2026-09-23",
            "frequency_Hz": F_HZ,
            "rpm": RPM,
            "timesteps": TIMESTEPS,
            "dt_s": DT
        },
        "nominal_plus30deg": res_plus30,
        "specular_minus30deg": res_minus30,
        "bias_decoupling": {
            "total_active_domain": {
                "F_plus30_uN": F_plus,
                "F_minus30_uN": F_minus,
                "F_bias_uN": F_bias,
                "F_chiral_pure_uN": F_chiral,
                "chiral_fraction_pct": (abs(F_chiral) / abs(F_plus)) * 100.0 if abs(F_plus) > 0 else 0.0
            },
            "mantle_only": {
                "F_mantle_plus30_uN": F_man_plus,
                "F_mantle_minus30_uN": F_man_minus,
                "F_mantle_bias_uN": F_man_bias,
                "F_mantle_chiral_pure_uN": F_man_chiral
            }
        },
        "maxwell_stress_tensor_validation": {
            "control_cylinder_radius_m": mst_grid["Rc"],
            "control_cylinder_height_m": 2.0 * mst_grid["Hc"],
            "Fz_MST_mean_uN": mst_plus,
            "Fz_Lorentz_mean_uN": F_plus,
            "difference_uN": diff_mst_lorentz,
            "difference_pct": err_mst_pct
        }
    }

    out_json_path = DATA_DIR / "validazione_mst_chiral_bias.json"
    with open(out_json_path, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2)
    print(f"\nDataset di validazione Fase 1 salvato in: {out_json_path}")

if __name__ == "__main__":
    main()
