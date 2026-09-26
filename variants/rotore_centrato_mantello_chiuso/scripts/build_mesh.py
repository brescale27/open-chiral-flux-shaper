#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generazione e Conversione Mesh 3D per Variante con Mantello Chiuso a Barattolo.
Genera la mesh conforme 3D tramite API Gmsh OpenCASCADE ed esegue la conversione ElmerGrid.
"""

import sys
import subprocess
from pathlib import Path
import gmsh

SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_DIR = SCRIPT_DIR.parent
MESH_DIR = VARIANT_DIR / "mesh"
MSH_FILE = MESH_DIR / "macchina_mantello_chiuso.msh"
GEO_FILE = MESH_DIR / "macchina_mantello_chiuso.geo"
OUT_MESH_DIR = MESH_DIR / "macchina_mantello_chiuso"
def find_elmergrid():
    cmd = os.environ.get("ELMERGRID_BIN") or shutil.which("ElmerGrid")
    if cmd:
        return cmd
    candidates = [
        r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe",
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerGrid.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerGrid.exe",
        "/usr/local/bin/ElmerGrid",
        "/usr/bin/ElmerGrid",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return "ElmerGrid"

ELMERGRID_BIN = find_elmergrid()

def build():
    print("=" * 80)
    print("BUILD MESH: MANTELLO CHIUSO A BARATTOLO (Z=+H/2 E Z=-H/2)")
    print("=" * 80)
    
    MESH_DIR.mkdir(parents=True, exist_ok=True)
    
    gmsh.initialize()
    gmsh.model.add("macchina_mantello_chiuso")
    occ = gmsh.model.occ

    r_ext = 0.05
    r_int = 0.047
    h_alu = 0.1
    r_far = 0.2
    t_core = 0.015
    t_cap = 0.003

    v_rot_up = occ.addCylinder(0, 0, t_core/2, 0, 0, h_alu/2 - t_core/2, r_int)
    v_rot_down = occ.addCylinder(0, 0, -h_alu/2, 0, 0, h_alu/2 - t_core/2, r_int)
    v_core = occ.addCylinder(0, 0, -t_core/2, 0, 0, t_core, r_int)

    v_can_out = occ.addCylinder(0, 0, -h_alu/2 - t_cap, 0, 0, h_alu + 2*t_cap, r_ext)
    v_can_in = occ.addCylinder(0, 0, -h_alu/2, 0, 0, h_alu, r_int)
    v_mantle = occ.cut([(3, v_can_out)], [(3, v_can_in)])[0][0][1]

    v_sph = occ.addSphere(0, 0, 0, r_far)

    occ.fragment([(3, v_sph)], [(3, v_rot_up), (3, v_rot_down), (3, v_core), (3, v_mantle)])
    occ.synchronize()

    # Gruppi fisici conformi
    gmsh.model.addPhysicalGroup(3, [5], 1, 'AirExterior')
    gmsh.model.addPhysicalGroup(3, [4], 2, 'Aluminum')
    gmsh.model.addPhysicalGroup(3, [1, 2], 3, 'RotorAir')
    gmsh.model.addPhysicalGroup(3, [3], 4, 'FerromagneticCore')

    surfs = gmsh.model.getEntities(2)
    far_surf = None
    for dim, tag in surfs:
        bb = occ.getBoundingBox(dim, tag)
        r = max(abs(bb[0]), abs(bb[1]), abs(bb[3]), abs(bb[4]), abs(bb[2]), abs(bb[5]))
        if abs(r - r_far) < 1e-4:
            far_surf = tag
            break
    gmsh.model.addPhysicalGroup(2, [far_surf], 1, 'FarField')

    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.002)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', 0.025)
    gmsh.option.setNumber('Mesh.Algorithm', 6)
    gmsh.option.setNumber('Mesh.Algorithm3D', 1)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)

    print("Generazione mesh tetraedrica 3D...")
    gmsh.model.mesh.generate(3)
    gmsh.write(str(MSH_FILE))
    gmsh.finalize()
    print(f"Mesh Gmsh 2.2 generata: {MSH_FILE} ({MSH_FILE.stat().st_size / 1e6:.2f} MB)")

    print("\nConversione con ElmerGrid...")
    p = subprocess.run([ELMERGRID_BIN, '14', '2', str(MSH_FILE), '-autoclean'], cwd=str(MESH_DIR), capture_output=True, text=True)
    if p.returncode != 0:
        print("[ERRORE] ElmerGrid:", p.stderr)
        sys.exit(1)
        
    print(f"ElmerGrid completato con successo in: {OUT_MESH_DIR}")
    return 0

if __name__ == "__main__":
    sys.exit(build())
