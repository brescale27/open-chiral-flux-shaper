#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generazione e Conversione Mesh Conforme 3D per Variante:
Gabbia Sferica Ibrida con Architettura a 24 Gruppi (Mappa di Fibonacci & Radice Numerica).

Geometria CAD OpenCASCADE in Gmsh:
- Guscio sferico metamateriale (R_ext = 50 mm, R_int = 47 mm, t = 3 mm)
- Nucleo amagnetico sferico centrale in PEEK (R_core = 12 mm, mu_r = 1.0, sigma = 0.0 S/m)
- 24 Gruppi di Solenoidi distribuiti a Delta_theta = 15° (R_c = 35 mm, r_wire = 3 mm, h_wire = 14 mm)
- Sfera esterna di Far-Field (R_far = 200 mm)
- Meshing 3D Netgen/MeshAdapt conforme e conversione ElmerGrid.

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
MSH_FILE = MESH_DIR / "macchina_fibonacci_24x24.msh"
OUT_MESH_DIR = MESH_DIR / "macchina_fibonacci_24x24"
ELMERGRID_BIN = r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe"


def build():
    print("=" * 80)
    print("BUILD MESH 3D: GABBIA SFERICA CON ARRAY A 24 GRUPPI DI FIBONACCI")
    print("=" * 80)

    MESH_DIR.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.model.add("macchina_fibonacci_24x24")
    occ = gmsh.model.occ

    r_far = 0.200    # 200 mm raggio sfera esterna
    r_ext = 0.050    # 50 mm raggio esterno mantello
    r_int = 0.047    # 47 mm raggio interno mantello (spessore 3 mm)
    rc = 0.035       # 35 mm raggio cerchio bobine
    r_wire = 0.003   # 3 mm raggio cilindro bobina
    h_wire = 0.014   # 14 mm altezza bobina
    r_core = 0.012   # 12 mm raggio nucleo centrale PEEK

    # Sfere concentriche
    s_far = occ.addSphere(0, 0, 0, r_far)
    s_ext = occ.addSphere(0, 0, 0, r_ext)
    s_int = occ.addSphere(0, 0, 0, r_int)
    s_core = occ.addSphere(0, 0, 0, r_core)

    # 24 solenoidi distribuiti a 15°
    coils = []
    for k in range(24):
        ang = k * (2 * np.pi / 24)
        x = rc * np.cos(ang)
        y = rc * np.sin(ang)
        cyl = occ.addCylinder(x, y, -h_wire / 2, 0, 0, h_wire, r_wire)
        coils.append((3, cyl))

    # Frammentazione conforme OpenCASCADE
    tool_vols = [(3, s_core)] + coils + [(3, s_int), (3, s_ext)]
    occ.fragment([(3, s_far)], tool_vols)
    occ.synchronize()

    # Classificazione robusta dei solidi
    entities = gmsh.model.getEntities(3)
    classification = {
        'mantle': [],
        'coils': [],
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
            classification['coils'].append(tag)
        else:
            classification['cavity_air'].append(tag)

    print("Classificazione corpi fisici 3D:")
    for k, v in classification.items():
        print(f"  - {k:12s}: {len(v)} solidi (tag {v})")

    if len(classification['coils']) != 24:
        print(f"[ATTENZIONE] Trovati {len(classification['coils'])} solidi bobine anziché 24!")
        sys.exit(1)

    # Assegnazione Gruppi Fisici 3D (Body 1 to 5)
    gmsh.model.addPhysicalGroup(3, classification['mantle'], 1, 'MantelloSferico')
    gmsh.model.addPhysicalGroup(3, classification['coils'], 2, 'CoilsFibonacci24')
    gmsh.model.addPhysicalGroup(3, classification['core'], 3, 'RotorePEEKCore')
    gmsh.model.addPhysicalGroup(3, classification['cavity_air'], 4, 'CavitaAriaInterna')
    gmsh.model.addPhysicalGroup(3, classification['air_ext'], 5, 'AriaEsterna')

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
        max_r_surf = max(surfs, key=lambda st: max(abs(occ.getBoundingBox(st[0], st[1]))))
        gmsh.model.addPhysicalGroup(2, [max_r_surf[1]], 101, 'FarField')

    # Parametri meshing
    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.0035)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', 0.028)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.option.setNumber('Mesh.Algorithm', 1)    # MeshAdapt 2D
    gmsh.option.setNumber('Mesh.Algorithm3D', 1)  # TetGen 3D

    print("\nGenerazione mesh tetraedrica 3D conforme (TetGen + MeshAdapt)...")
    gmsh.model.mesh.generate(3)

    nodes = gmsh.model.mesh.getNodes()
    elements = gmsh.model.mesh.getElements()
    n_nodes = len(nodes[0])
    n_elem = sum(len(el) for el in elements[1])
    print(f"Mesh generata con successo: {n_nodes} nodi, {n_elem} elementi.")

    gmsh.write(str(MSH_FILE))
    print(f"Mesh salvata in: {MSH_FILE}")
    gmsh.finalize()

    # Conversione ElmerGrid
    print(f"\nConversione in formato ElmerGrid: {OUT_MESH_DIR}")
    cmd = [ELMERGRID_BIN, "14", "2", str(MSH_FILE), "-autoclean"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[ERRORE] ElmerGrid fallito (codice {res.returncode}):\n{res.stderr}\n{res.stdout}")
        sys.exit(1)
    print("Conversione ElmerGrid completata con successo.")
    print("=" * 80)


if __name__ == "__main__":
    build()
