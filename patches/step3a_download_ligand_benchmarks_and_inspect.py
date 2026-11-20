#!/usr/bin/env python3
"""
Phase 2 - Step 3A (Read-Only):
  1. Downloads RCSB PDB 4WNP.cif (Human ULK1 Kinase + Small-Molecule Inhibitor) and
     7DIY.cif (SARS-CoV-2 Nsp14-Nsp10 + Zn + Cofactor) into benchmarks/.
  2. Runs unmodified Step-2C unippi_cli.py on 4WNP and 7DIY to record BEFORE baselines.
  3. Tests our universal 3D covalent-graph planar ring detector & ligand SASA burial
     geometry on 5ZC9 (RCG600, ANP500), 4WNP (ULK1 inhibitor), and 7DIY (cofactor/metals)
     so Step 3B parameters are calibrated on real crystal geometry.
"""

import math
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, ".")
import unippi_cli as u

bench = Path("benchmarks")
bench.mkdir(exist_ok=True)

for pdb_id, desc in [
    ("4WNP", "Human ULK1 Kinase Domain + Small-Molecule Inhibitor (MARVAS ULK1/2 Reference)"),
    ("7DIY", "SARS-CoV-2 Nsp14-Nsp10 Complex + Zn2+ + Cofactor (MARVAS Nsp14-Nsp10 Reference)"),
]:
    dst = bench / f"{pdb_id}.cif"
    if not dst.exists():
        print(f"[DOWNLOAD] {pdb_id}.cif — {desc}")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.cif", dst)
        print(f"  Saved {dst} ({dst.stat().st_size:,} bytes)")
    else:
        print(f"[CACHED]   {dst} ({dst.stat().st_size:,} bytes)")

def detect_planar_rings(res):
    """Universal 5-/6-membered planar ring detector from 3D heavy-atom coordinates."""
    atoms = [a for a in res.atoms if a.element in ("C", "N", "O", "S")]
    n = len(atoms)
    adj = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            d = math.dist(atoms[i].coord, atoms[j].coord)
            if 1.15 <= d <= 1.68:
                adj[i].append(j)
                adj[j].append(i)
    rings = set()
    for start in range(n):
        stack = [(start, [start])]
        while stack:
            curr, path = stack.pop()
            if len(path) > 6:
                continue
            for nxt in adj[curr]:
                if nxt == start and len(path) in (5, 6):
                    rings.add(tuple(sorted(path)))
                elif nxt > start and nxt not in path and len(path) < 6:
                    stack.append((nxt, path + [nxt]))
    planar = []
    for r_idx in sorted(rings):
        pts = [atoms[k].coord for k in r_idx]
        m = len(pts)
        cx = sum(p[0] for p in pts) / m
        cy = sum(p[1] for p in pts) / m
        cz = sum(p[2] for p in pts) / m
        # Plane normal via cross product of 3 well-spaced ring vertices
        p1, p2, p3 = pts[0], pts[m // 3], pts[(2 * m) // 3]
        ux, uy, uz = p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]
        vx, vy, vz = p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        nl = math.sqrt(nx * nx + ny * ny + nz * nz)
        if nl < 1e-6:
            continue
        nx, ny, nz = nx / nl, ny / nl, nz / nl
        max_dev = max(abs((p[0] - cx) * nx + (p[1] - cy) * ny + (p[2] - cz) * nz) for p in pts)
        if max_dev <= 0.15:
            planar.append((len(r_idx), (cx, cy, cz), (nx, ny, nz), max_dev, [atoms[k].name for k in r_idx]))
    return planar

print("\n" + "=" * 92)
print("UNIVERSAL LIGAND & CHAIN GEOMETRY AUDIT (5ZC9, 4WNP, 7DIY)")
print("=" * 92)
for fname in ("5ZC9.cif", "4WNP.cif", "7DIY.cif"):
    atoms, meta = u.parse_structure_file(str(bench / fname))
    ordered, chains = u.build_residues(atoms)
    print(f"\n[{fname}] method={meta.get('method')} | Chains ({len(chains)}): " +
          ", ".join(f"{cid}({len(rlist)} {rlist[0].entity_type})" for cid, rlist in chains.items()))
    sasa_full = u.compute_residue_sasa(ordered)
    for cid, rlist in chains.items():
        if rlist[0].entity_type == "LIGAND":
            lig = rlist[0]
            sasa_free = u.compute_residue_sasa([lig])[lig.key]
            sasa_bound = sasa_full.get(lig.key, 0.0)
            f_bur_tot = max(0.0, (sasa_free - sasa_bound) / max(1.0, sasa_free))
            rings = detect_planar_rings(lig)
            halogens = [a.name for a in lig.atoms if a.element in ("F", "CL", "BR", "I")]
            # Check burial by each macromolecular chain separately
            by_chain = []
            for mcid, mrlist in chains.items():
                if mrlist[0].entity_type in ("PROTEIN", "NUCLEIC"):
                    if any(math.dist(lig.centroid, r.centroid) <= 14.0 for r in mrlist):
                        s_pair = u.compute_residue_sasa(mrlist + [lig])[lig.key]
                        fb = max(0.0, (sasa_free - s_pair) / max(1.0, sasa_free))
                        if fb >= 0.05:
                            by_chain.append(f"{mcid}:{fb*100:.1f}%")
            print(f"  LIGAND {cid} ({lig.res_name}, {len(lig.atoms)} heavy atoms): "
                  f"SASA_free={sasa_free:.1f} A^2 | Total_Burial={f_bur_tot*100:.1f}% ({', '.join(by_chain)}) | "
                  f"PlanarRings={len(rings)} ({[r[0] for r in rings]}) | Halogens={halogens}")

for title, inp, outp in [
    ("BEFORE BENCHMARK — 4WNP.cif (Human ULK1 Kinase + Small-Molecule Inhibitor)", "benchmarks/4WNP.cif", "/tmp/v2_before_4wnp"),
    ("BEFORE BENCHMARK — 7DIY.cif (SARS-CoV-2 Nsp14-Nsp10 Complex)", "benchmarks/7DIY.cif", "/tmp/v2_before_7diy"),
]:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", inp, "-o", outp], check=True)
