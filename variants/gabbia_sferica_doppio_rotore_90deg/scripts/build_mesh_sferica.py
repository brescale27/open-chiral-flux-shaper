#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generazione e Conversione Mesh Conforme 3D per Variante:
Gabbia Sferica Metamateriale con Doppio Rotore Ortogonale a 90°.

Geometria CAD OpenCASCADE in Gmsh:
- Guscio sferico metamateriale (R_ext = 50 mm, R_int = 47 mm, t = 3 mm)
- Nucleo amagnetico sferico centrale in PEEK (R_core = 12 mm, mu_r = 1.0, sigma = 0.0 S/m)
- Doppio Rotore Ortogonale a 90°:
  * Rotore 1 (Equatoriale Z=0): 6 solenoidi paralleli all'asse Z (R_pos = 35 mm)
  * Rotore 2 (Trasversale X=0): 6 solenoidi paralleli all'asse X (R_pos = 35 mm, sfalsati di 30°)
- Sfera esterna di Far-Field (R_far = 200 mm)
- Meshing 3D Netgen conforme e conversione ElmerGrid.

Autore: Alessandro Brescacin
Licenza: CERN-OHL-S-2.0
"""

import sys
import subprocess
from pathlib import Path
import numpy as np
import gmsh

SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_DIR = SCRIPT_DIR.parent
MESH_DIR = VARIANT_DIR / "mesh"
MSH_FILE = MESH_DIR / "macchina_sferica_ortogonale.msh"
OUT_MESH_DIR = MESH_DIR / "macchina_sferica_ortogonale"
ELMERGRID_BIN = r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe"


def build():
    print("=" * 80)
    print("BUILD MESH 3D: GABBIA SFERICA CON DOPPIO ROTORE ORTOGONALE A 90°")
    print("=" * 80)

    MESH_DIR.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.model.add("macchina_sferica_ortogonale")
    occ = gmsh.model.occ

    r_far = 0.200    # 200 mm raggio sfera esterna
    r_ext = 0.050    # 50 mm raggio esterno mantello
    r_int = 0.047    # 47 mm raggio interno mantello (spessore 3 mm)
    rc = 0.035       # 35 mm raggio posizionamento bobine
    r_wire = 0.004   # 4 mm raggio cilindro bobina
    h_wire = 0.014   # 14 mm altezza bobina
    r_core = 0.012   # 12 mm raggio nucleo centrale PEEK

    # Sfere concentriche
    s_far = occ.addSphere(0, 0, 0, r_far)
    s_ext = occ.addSphere(0, 0, 0, r_ext)
    s_int = occ.addSphere(0, 0, 0, r_int)
    s_core = occ.addSphere(0, 0, 0, r_core)

    # Rotore 1: 6 solenoidi verticali paralleli a Z nel piano Z=0
    coils_1 = []
    for deg in [0, 60, 120, 180, 240, 300]:
        rad = np.radians(deg)
        x = rc * np.cos(rad)
        y = rc * np.sin(rad)
        cyl = occ.addCylinder(x, y, -h_wire / 2, 0, 0, h_wire, r_wire)
        coils_1.append(cyl)

    # Rotore 2: 6 solenoidi trasversali paralleli a X nel piano X=0 (Y-Z), sfalsati di 30°
    coils_2 = []
    for deg in [30, 90, 150, 210, 270, 330]:
        rad = np.radians(deg)
        y = rc * np.cos(rad)
        z = rc * np.sin(rad)
        cyl = occ.addCylinder(-h_wire / 2, y, z, h_wire, 0, 0, r_wire)
        coils_2.append(cyl)

    # Frammentazione conforme OpenCASCADE
    tool_vols = [(3, s_core)] + [(3, c) for c in coils_1] + [(3, c) for c in coils_2] + [(3, s_int), (3, s_ext)]
    occ.fragment([(3, s_far)], tool_vols)
    occ.synchronize()

    # Classificazione robusta dei solidi
    entities = gmsh.model.getEntities(3)
    classification = {
        'mantle': [],
        'coils1': [],
        'coils2': [],
        'core': [],
        'cavity_air': [],
        'air_ext': []
    }

    for dim, tag in entities:
        com = occ.getCenterOfMass(dim, tag)
        bb = occ.getBoundingBox(dim, tag)
        rmax = max(abs(bb[0]), abs(bb[1]), abs(bb[2]), abs(bb[3]), abs(bb[4]), abs(bb[5]))
        rcom = np.sqrt(com[0]**2 + com[1]**2 + com[2]**2)

        if rmax > 0.15:
            classification['air_ext'].append(tag)
        elif rmax > 0.048 and abs(rcom) < 0.005:
            classification['mantle'].append(tag)
        elif rmax <= 0.013 and abs(rcom) < 0.005:
            classification['core'].append(tag)
        elif abs(com[2]) < 0.005 and 0.030 < rcom < 0.040:
            classification['coils1'].append(tag)
        elif abs(com[0]) < 0.005 and 0.030 < rcom < 0.040:
            classification['coils2'].append(tag)
        else:
            classification['cavity_air'].append(tag)

    print("Classificazione corpi fisici 3D:")
    for k, v in classification.items():
        print(f"  - {k:12s}: {len(v)} solidi (tag {v})")

    # Assegnazione Gruppi Fisici 3D (Body 1 to 6)
    gmsh.model.addPhysicalGroup(3, classification['mantle'], 1, 'MantelloSferico')
    gmsh.model.addPhysicalGroup(3, classification['coils1'], 2, 'Rotore1Coils')
    gmsh.model.addPhysicalGroup(3, classification['coils2'], 3, 'Rotore2Coils')
    gmsh.model.addPhysicalGroup(3, classification['core'], 4, 'RotorePEEKCore')
    gmsh.model.addPhysicalGroup(3, classification['cavity_air'], 5, 'CavitaAriaInterna')
    gmsh.model.addPhysicalGroup(3, classification['air_ext'], 6, 'AriaEsterna')

    # Identificazione superficie 2D di FarField (R = 0.200 m)
    surfs = gmsh.model.getEntities(2)
    far_surf = None
    for dim, tag in surfs:
        bb = occ.getBoundingBox(dim, tag)
        r = max(abs(bb[0]), abs(bb[1]), abs(bb[3]), abs(bb[4]), abs(bb[2]), abs(bb[5]))
        if abs(r - r_far) < 1e-4:
            far_surf = tag
            break

    if far_surf is not None:
        gmsh.model.addPhysicalGroup(2, [far_surf], 101, 'FarField')
        print(f"  - FarField 2D: superficie tag {far_surf} -> Physical Tag 101")
    else:
        print("[ATTENZIONE] Superficie FarField non trovata univocamente, seleziono la più esterna.")
        max_r_surf = max(surfs, key=lambda st: max(abs(occ.getBoundingBox(st[0], st[1]))))
        gmsh.model.addPhysicalGroup(2, [max_r_surf[1]], 101, 'FarField')

    # Parametri meshing
    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.003)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', 0.030)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.option.setNumber('Mesh.Algorithm', 1)    # MeshAdapt 2D
    gmsh.option.setNumber('Mesh.Algorithm3D', 1)  # TetGen 3D

    print("\nGenerazione mesh tetraedrica 3D conforme (TetGen + MeshAdapt)...")
    gmsh.model.mesh.generate(3)
    gmsh.write(str(MSH_FILE))
    gmsh.finalize()
    print(f"Mesh Gmsh 2.2 generata con successo: {MSH_FILE} ({MSH_FILE.stat().st_size / 1e6:.2f} MB)")

    # Conversione con ElmerGrid
    print("\nConversione mesh con ElmerGrid...")
    cmd = [ELMERGRID_BIN, '14', '2', str(MSH_FILE), '-autoclean']
    p = subprocess.run(cmd, cwd=str(MESH_DIR), capture_output=True, text=True)
    if p.returncode != 0:
        print("[ERRORE] ElmerGrid fallito:", p.stderr)
        sys.exit(1)

    print(f"ElmerGrid completato con successo in: {OUT_MESH_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(build())
