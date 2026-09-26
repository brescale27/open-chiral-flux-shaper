#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generazione e Conversione Mesh Conforme 3D per Variante:
Gabbia Sferica Ibrida con Mantello a 3 Lobi Macro-Chirali (Triskelion),
Armatura Centrale a Esagramma e Array di 24 Gruppi di Solenoidi.

Geometria CAD OpenCASCADE in Gmsh:
- Guscio sferico metamateriale a 3 lobi macro-chirali a triskelion (R_ext = 50 mm, R_int = 47 mm, tilt +30°)
- Nucleo amagnetico centrale a stella / esagramma in PEEK (R_star = 14 mm, h = 16 mm, mu_r = 1.0, sigma = 0.0 S/m)
- 24 Gruppi di Solenoidi distribuiti a Delta_theta = 15° (R_c = 35 mm, r_wire = 3 mm, h_wire = 14 mm)
- Sfera esterna di Far-Field (R_far = 200 mm)
- Meshing 3D TetGen/MeshAdapt conforme e conversione ElmerGrid.

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
MSH_FILE = MESH_DIR / "macchina_triskelion_24pulse.msh"
OUT_MESH_DIR = MESH_DIR / "macchina_triskelion_24pulse"
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
    print("BUILD MESH 3D: GABBIA A 3 LOBI (TRISKELION) CON ARMATURA A ESAGRAMMA E 24 COILS")
    print("=" * 80)

    MESH_DIR.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.model.add("macchina_triskelion_24pulse")
    occ = gmsh.model.occ

    r_far = 0.200    # 200 mm sfera esterna Far-Field
    r_ext = 0.050    # 50 mm raggio base esterno mantello
    r_int = 0.047    # 47 mm raggio interno mantello
    rc = 0.035       # 35 mm raggio cerchio bobine
    r_wire = 0.003   # 3 mm raggio cilindro bobina
    h_wire = 0.014   # 14 mm altezza bobina
    r_star = 0.014   # 14 mm raggio esterno esagramma
    h_star = 0.016   # 16 mm altezza nucleo esagramma

    # 1. Nucleo centrale a Esagramma (due prismi triangolari incrociati a 60°)
    p1 = [occ.addPoint(r_star * np.cos(a), r_star * np.sin(a), -h_star / 2) for a in [0, 2 * np.pi / 3, 4 * np.pi / 3]]
    l1 = [occ.addLine(p1[0], p1[1]), occ.addLine(p1[1], p1[2]), occ.addLine(p1[2], p1[0])]
    s1 = occ.addPlaneSurface([occ.addWire(l1)])
    v1 = [x for x in occ.extrude([(2, s1)], 0, 0, h_star) if x[0] == 3]

    p2 = [occ.addPoint(r_star * np.cos(a), r_star * np.sin(a), -h_star / 2) for a in [np.pi / 3, np.pi, 5 * np.pi / 3]]
    l2 = [occ.addLine(p2[0], p2[1]), occ.addLine(p2[1], p2[2]), occ.addLine(p2[2], p2[0])]
    s2 = occ.addPlaneSurface([occ.addWire(l2)])
    v2 = [x for x in occ.extrude([(2, s2)], 0, 0, h_star) if x[0] == 3]
    core_fused, _ = occ.fuse(v1, v2)

    # 2. 24 Solenoidi equidistanziati a Delta_theta = 15°
    coils = []
    for k in range(24):
        ang = k * (2 * np.pi / 24)
        x = rc * np.cos(ang)
        y = rc * np.sin(ang)
        cyl = occ.addCylinder(x, y, -h_wire / 2, 0, 0, h_wire, r_wire)
        coils.append((3, cyl))

    # 3. Mantello a 3 Lobi Macro-Chirali (Triskelion)
    s_ext = occ.addSphere(0, 0, 0, r_ext)
    s_int = occ.addSphere(0, 0, 0, r_int)

    lobes = []
    tilt = np.radians(30.0)  # Inclinazione macro-chirale di 30°
    l_len = 0.040
    r_lobe = 0.007
    for k in range(3):
        ang = k * (2 * np.pi / 3)
        dx = -np.sin(ang) * np.sin(tilt) * l_len
        dy =  np.cos(ang) * np.sin(tilt) * l_len
        dz =  np.cos(tilt) * l_len
        x0 = 0.048 * np.cos(ang) - dx / 2
        y0 = 0.048 * np.sin(ang) - dy / 2
        z0 = -dz / 2
        cyl = occ.addCylinder(x0, y0, z0, dx, dy, dz, r_lobe)
        lobes.append((3, cyl))

    mantle_fused, _ = occ.fuse([(3, s_ext)], lobes)

    # 4. Sfera Far-Field esterna
    s_far = occ.addSphere(0, 0, 0, r_far)

    # 5. Frammentazione conforme OpenCASCADE
    tools = core_fused + coils + [(3, s_int)] + mantle_fused
    occ.fragment([(3, s_far)], tools)
    occ.synchronize()

    # 6. Classificazione dei corpi fisici 3D
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
        elif rmax > 0.048 and abs(rcom) < 0.008:
            classification['mantle'].append(tag)
        elif rmax <= 0.016 and abs(rcom) < 0.005:
            classification['core'].append(tag)
        elif abs(com[2]) < 0.005 and 0.030 < rcom < 0.040:
            classification['coils'].append(tag)
        else:
            classification['cavity_air'].append(tag)

    print("Classificazione corpi fisici 3D:")
    for k, v in classification.items():
        print(f"  - {k:12s}: {len(v)} solidi (tag {v})")

    if len(classification['coils']) != 24:
        print(f"[ERRORE] Trovati {len(classification['coils'])} solidi bobine anziché 24!")
        sys.exit(1)
    if len(classification['mantle']) != 1:
        print(f"[ERRORE] Trovati {len(classification['mantle'])} solidi mantello anziché 1!")
        sys.exit(1)
    if len(classification['core']) != 1:
        print(f"[ERRORE] Trovati {len(classification['core'])} solidi core anziché 1!")
        sys.exit(1)

    # Assegnazione Gruppi Fisici 3D (Body 1 to 5)
    gmsh.model.addPhysicalGroup(3, classification['mantle'], 1, 'MantelloTriskelion')
    gmsh.model.addPhysicalGroup(3, classification['coils'], 2, 'Coils24')
    gmsh.model.addPhysicalGroup(3, classification['core'], 3, 'CoreEsagramma')
    gmsh.model.addPhysicalGroup(3, classification['cavity_air'], 4, 'CavitaAriaInterna')
    gmsh.model.addPhysicalGroup(3, classification['air_ext'], 5, 'AriaEsterna')

    # Identificazione superficie 2D di FarField (R = 0.200 m)
    surfs = gmsh.model.getEntities(2)
    far_surfs = []
    for dim, tag in surfs:
        bb = occ.getBoundingBox(dim, tag)
        r = max(abs(bb[0]), abs(bb[1]), abs(bb[3]), abs(bb[4]), abs(bb[2]), abs(bb[5]))
        if abs(r - r_far) < 1e-4:
            far_surfs.append(tag)

    if far_surfs:
        gmsh.model.addPhysicalGroup(2, far_surfs, 101, 'FarField')
        print(f"  - FarField 2D: superfici tag {far_surfs} -> Physical Tag 101")
    else:
        max_r_surf = max(surfs, key=lambda st: max(abs(occ.getBoundingBox(st[0], st[1]))))
        gmsh.model.addPhysicalGroup(2, [max_r_surf[1]], 101, 'FarField')

    # Parametri meshing
    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.004)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', 0.030)
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
