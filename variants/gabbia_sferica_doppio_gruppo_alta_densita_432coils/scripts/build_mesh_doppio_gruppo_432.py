#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generazione e Conversione Mesh Conforme 3D per Variante:
Gabbia Sferica con Doppio Macro-Gruppo Ortogonale ad Alta Densità (432 Solenoidi Totali, 216 Bobine/Gruppo)
- Gruppo 1 (Equatoriale Z=0): 216 micro-solenoidi paralleli all'asse Z (R_c1 = 37 mm, Delta_theta = 1.667°)
- Gruppo 2 (Meridiano X=0): 216 micro-solenoidi paralleli all'asse X (R_c2 = 31 mm, Delta_psi = 1.667°)
- Clearance radiale di 6.0 mm tra i due anelli
- Guscio sferico metamateriale a triplo strato X (R_ext = 50 mm, R_int = 47 mm, spessore = 3 mm, mu_r = 1000.0)
- Nucleo amagnetico sferico centrale in PEEK (R_core = 12 mm, mu_r = 1.0, sigma = 0.0 S/m)
- Sfera esterna di Far-Field (R_far = 200 mm)
- Meshing 3D TetGen conforme e conversione ElmerGrid.

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
MSH_FILE = MESH_DIR / "macchina_doppio_gruppo_432.msh"
OUT_MESH_DIR = MESH_DIR / "macchina_doppio_gruppo_432"
ELMERGRID_BIN = r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe"


def build():
    print("=" * 80)
    print("BUILD MESH 3D: DOPPIO MACRO-GRUPPO AD ALTA DENSITÀ (432 SOLENOIDI, 216/GRUPPO)")
    print("=" * 80)

    MESH_DIR.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.model.add("macchina_doppio_gruppo_432")
    occ = gmsh.model.occ

    r_far = 0.200    # 200 mm sfera esterna Far-Field
    r_ext = 0.050    # 50 mm raggio esterno mantello
    r_int = 0.047    # 47 mm raggio interno mantello (spessore 3 mm)
    rc1 = 0.037      # 37 mm raggio cerchio Gruppo 1 (Z)
    rc2 = 0.031      # 31 mm raggio cerchio Gruppo 2 (X) -> clearance 6.0 mm
    r_wire = 0.00015 # 0.15 mm raggio filo solenoide elementare (diametro 0.30 mm)
    h_wire = 0.010   # 10 mm altezza solenoide
    r_core = 0.012   # 12 mm raggio nucleo centrale PEEK

    # Sfere concentriche
    s_far = occ.addSphere(0, 0, 0, r_far)
    s_ext = occ.addSphere(0, 0, 0, r_ext)
    s_int = occ.addSphere(0, 0, 0, r_int)
    s_core = occ.addSphere(0, 0, 0, r_core)

    # Gruppo 1: 216 solenoidi nel piano equatoriale XY (paralleli a Z, Delta_theta = 1.6667°)
    coils_1 = []
    for k in range(216):
        ang = k * (2.0 * np.pi / 216.0)
        x = rc1 * np.cos(ang)
        y = rc1 * np.sin(ang)
        cyl = occ.addCylinder(x, y, -h_wire / 2.0, 0, 0, h_wire, r_wire)
        coils_1.append((3, cyl))

    # Gruppo 2: 216 solenoidi nel piano meridiano YZ (paralleli a X, Delta_psi = 1.6667°)
    coils_2 = []
    for k in range(216):
        ang = k * (2.0 * np.pi / 216.0)
        y = rc2 * np.cos(ang)
        z = rc2 * np.sin(ang)
        cyl = occ.addCylinder(-h_wire / 2.0, y, z, h_wire, 0, 0, r_wire)
        coils_2.append((3, cyl))

    print(f"Creati {len(coils_1)} solenoidi nel Gruppo 1 (Z) e {len(coils_2)} nel Gruppo 2 (X). Totale: 432")

    # Frammentazione conforme OpenCASCADE
    print("Frammentazione conforme OpenCASCADE in corso...")
    tools = [(3, s_core)] + coils_1 + coils_2 + [(3, s_int), (3, s_ext)]
    occ.fragment([(3, s_far)], tools)
    occ.synchronize()

    # Classificazione robusta dei corpi fisici 3D
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
        elif abs(com[2]) < 0.005 and 0.035 < rcom < 0.039:
            classification['coils1'].append(tag)
        elif abs(com[0]) < 0.005 and 0.029 < rcom < 0.033:
            classification['coils2'].append(tag)
        else:
            classification['cavity_air'].append(tag)

    print("Classificazione corpi fisici 3D:")
    for k, v in classification.items():
        print(f"  - {k:12s}: {len(v)} solidi")

    if len(classification['coils1']) != 216:
        print(f"[ERRORE] Trovate {len(classification['coils1'])} bobine nel Gruppo 1 anziché 216!")
        sys.exit(1)
    if len(classification['coils2']) != 216:
        print(f"[ERRORE] Trovate {len(classification['coils2'])} bobine nel Gruppo 2 anziché 216!")
        sys.exit(1)
    if len(classification['mantle']) != 1:
        print(f"[ERRORE] Trovato {len(classification['mantle'])} solido mantello anziché 1!")
        sys.exit(1)

    # Assegnazione Gruppi Fisici 3D (Body 1 to 6)
    gmsh.model.addPhysicalGroup(3, classification['mantle'], 1, 'MantelloSferico')
    gmsh.model.addPhysicalGroup(3, classification['coils1'], 2, 'CoilsGruppo1_Z')
    gmsh.model.addPhysicalGroup(3, classification['coils2'], 3, 'CoilsGruppo2_X')
    gmsh.model.addPhysicalGroup(3, classification['core'], 4, 'RotorePEEKCore')
    gmsh.model.addPhysicalGroup(3, classification['cavity_air'], 5, 'CavitaAriaInterna')
    gmsh.model.addPhysicalGroup(3, classification['air_ext'], 6, 'AriaEsterna')

    # Identificazione superficie 2D di FarField (R = 0.200 m)
    surfs = gmsh.model.getEntities(2)
    far_surfs = [tag for dim, tag in surfs if abs(max(abs(x) for x in occ.getBoundingBox(dim, tag)) - r_far) < 1e-4]

    if far_surfs:
        gmsh.model.addPhysicalGroup(2, far_surfs, 101, 'FarField')
        print(f"  - FarField 2D: superfici tag {far_surfs} -> Physical Tag 101")
    else:
        max_r_surf = max(surfs, key=lambda st: max(abs(occ.getBoundingBox(st[0], st[1]))))
        gmsh.model.addPhysicalGroup(2, [max_r_surf[1]], 101, 'FarField')

    # Parametri meshing conformi
    gmsh.option.setNumber('Mesh.AngleToleranceFacetOverlap', 1e-4)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.0004)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', 0.020)
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
