#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generatore della mesh 3D conforme per la variante a rotore centrato a Z=0:
- Nucleo ferromagnetico a Z=0 (spessore 15 mm)
- Doppio traferro in aria (+H/2 e -H/2)
- Mantello in alluminio (rete stirata omogeneizzata, r_int=47mm, r_ext=50mm)
- Sfera far-field (R=200mm)
- Conversione ElmerGrid in mesh/macchina_centrata/
"""

import gmsh
import math
import os
import shutil
import subprocess
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
VARIANTS_DIR = ROOT_DIR / "variants" / "rotore_centrato_z0"
MESH_DIR = VARIANTS_DIR / "mesh"

def find_elmergrid():
    cmd = shutil.which("ElmerGrid")
    if cmd:
        return cmd
    candidates = [
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerGrid.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerGrid.exe"
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return "ElmerGrid"

def generate_mesh():
    print(f"\n{'='*75}")
    print("Costruzione Geometria 3D OpenCASCADE e Mesh Conforme con Gmsh")
    print(f"{'='*75}")
    
    gmsh.initialize()
    gmsh.clear()
    occ = gmsh.model.occ

    # Parametri di discretizzazione
    gmsh.option.setNumber("Mesh.Algorithm", 5) # Delaunay 2D (conforme su superfici periodiche)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1) # Delaunay 3D (Tetgen)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.005)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.035)

    r_ext = 0.05
    r_int = 0.047
    h_alu = 0.1
    r_far = 0.2
    t_core = 0.015


    # 1. Volumi rotorici interni: traferro superiore, inferiore e nucleo centrale
    v_rot_up = occ.addCylinder(0, 0, t_core/2, 0, 0, h_alu/2 - t_core/2, r_int)
    v_rot_down = occ.addCylinder(0, 0, -h_alu/2, 0, 0, h_alu/2 - t_core/2, r_int)
    v_core = occ.addCylinder(0, 0, -t_core/2, 0, 0, t_core, r_int)

    # 2. Mantello cilindrico in alluminio
    v_m1 = occ.addCylinder(0, 0, -h_alu/2, 0, 0, h_alu, r_ext)
    v_m2 = occ.addCylinder(0, 0, -h_alu/2, 0, 0, h_alu, r_int)
    v_mantle, _ = occ.cut([(3, v_m1)], [(3, v_m2)])

    # 3. Sfera Far-Field
    v_sph = occ.addSphere(0, 0, 0, r_far)

    occ.synchronize()

    # Frammentazione conforme
    out, out_map = occ.fragment([(3, v_sph)], [(3, v_rot_up), (3, v_rot_down), (3, v_core), v_mantle[0]])
    occ.synchronize()

    # Mappatura rigorosa dei corpi tramite out_map
    tag_air_ext = [t for dim, t in out_map[0] if t not in [out_map[1][0][1], out_map[2][0][1], out_map[3][0][1], out_map[4][0][1]]]
    tag_rot_up = [out_map[1][0][1]]
    tag_rot_down = [out_map[2][0][1]]
    tag_core = [out_map[3][0][1]]
    tag_mantle = [out_map[4][0][1]]

    # Campo di raffinamento Box: mesh fine nel rotore/mantello (7 mm) e transizione verso il farfield
    gmsh.model.mesh.field.add("Box", 1)
    gmsh.model.mesh.field.setNumber(1, "VIn", 0.007)
    gmsh.model.mesh.field.setNumber(1, "VOut", 0.035)
    gmsh.model.mesh.field.setNumber(1, "XMin", -0.055)
    gmsh.model.mesh.field.setNumber(1, "XMax", 0.055)
    gmsh.model.mesh.field.setNumber(1, "YMin", -0.055)
    gmsh.model.mesh.field.setNumber(1, "YMax", 0.055)
    gmsh.model.mesh.field.setNumber(1, "ZMin", -0.055)
    gmsh.model.mesh.field.setNumber(1, "ZMax", 0.055)
    gmsh.model.mesh.field.setNumber(1, "Thickness", 0.03)
    gmsh.model.mesh.field.setAsBackgroundMesh(1)


    print(f"AirExterior volume tag: {tag_air_ext}")
    print(f"Aluminum Mantle volume tag: {tag_mantle}")
    print(f"RotorAir volume tags: {tag_rot_up + tag_rot_down}")
    print(f"FerromagneticCore volume tag: {tag_core}")

    gmsh.model.addPhysicalGroup(3, tag_air_ext, 1)
    gmsh.model.setPhysicalName(3, 1, "AirExterior")

    gmsh.model.addPhysicalGroup(3, tag_mantle, 2)
    gmsh.model.setPhysicalName(3, 2, "Aluminum")

    gmsh.model.addPhysicalGroup(3, tag_rot_up + tag_rot_down, 3)
    gmsh.model.setPhysicalName(3, 3, "RotorAir")

    gmsh.model.addPhysicalGroup(3, tag_core, 4)
    gmsh.model.setPhysicalName(3, 4, "FerromagneticCore")

    # Identificazione superficie sferica FarField
    farfield_surfs = []
    for dim, tag in gmsh.model.getEntities(2):
        bb = occ.getBoundingBox(dim, tag)
        rmax = max(math.hypot(bb[0], bb[1]), math.hypot(bb[3], bb[4]))
        zmax = max(abs(bb[2]), abs(bb[5]))
        if abs(rmax - r_far) < 1e-4 or abs(zmax - r_far) < 1e-4:
            farfield_surfs.append(tag)

    print(f"FarField boundary surfaces: {farfield_surfs}")
    gmsh.model.addPhysicalGroup(2, farfield_surfs, 1)
    gmsh.model.setPhysicalName(2, 1, "FarField")

    # Generazione mesh 3D
    gmsh.model.mesh.generate(3)

    msh_file = MESH_DIR / "macchina_centrata.msh"
    gmsh.write(str(msh_file))
    gmsh.finalize()
    print(f"[OK] Mesh Gmsh salvata: {msh_file} ({os.path.getsize(msh_file)} bytes)")

    # Conversione tramite ElmerGrid
    elmergrid_bin = find_elmergrid()
    cmd = [elmergrid_bin, "14", "2", "macchina_centrata.msh", "-autoclean"]
    print(f"\nEsecuzione ElmerGrid: {' '.join(cmd)}")
    res = subprocess.run(cmd, cwd=MESH_DIR, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[ERRORE] ElmerGrid fallito:\n{res.stderr}\n{res.stdout}")
        return False


    print("[OK] ElmerGrid completato con successo!")
    elmer_out_dir = MESH_DIR / "macchina_centrata"
    names_file = elmer_out_dir / "mesh.names"
    if names_file.is_file():

        print(f"Contenuto di {names_file.name}:")
        with open(names_file, "r") as f:
            print(f.read())
    return True

if __name__ == "__main__":
    generate_mesh()
