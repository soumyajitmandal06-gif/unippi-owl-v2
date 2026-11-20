#!/usr/bin/env python3
"""Read-only diagnostic: inspects 5ZC9/1URN ligands & chains, v1.8 JSON metrics, and unippi_cli.py anchors."""

import json
from collections import defaultdict
from pathlib import Path
import sys
sys.path.insert(0, ".")
import unippi_cli as u

for name in ("5ZC9.cif", "1URN.cif"):
    p = Path("benchmarks") / name
    atoms, meta = u.parse_structure_file(p)
    res_list, chains = u.build_residues(atoms)
    print("=" * 85)
    print(f"STRUCTURE LAYOUT: {name} (method={meta.get('method')})")
    for cid, rlist in chains.items():
        types = defaultdict(list)
        for r in rlist:
            types[r.entity_type].append(f"{r.res_name}{r.res_seq}")
        summary = " | ".join(f"{t}: {len(v)} ({', '.join(v[:4])}{'...' if len(v)>4 else ''})" for t, v in types.items())
        print(f"  Chain {cid}: {summary}")

for tag, jpath in [("5ZC9", "/tmp/v18_before_5zc9/unippi_full_report.json"), ("1URN", "/tmp/v18_before_1urn/unippi_full_report.json")]:
    data = json.loads(Path(jpath).read_text(encoding="utf-8"))
    print("\n" + "=" * 85)
    print(f"v1.8 JSON DETAILS ({tag}):")
    for r in data[:4]:
        print(f"  {r['chain_pair']} ({r['modality']}): UIS={r['UIS']} | dG_uni={r['dG_uni_kcal']} | dG_PRODIGY={r['dG_prodigy_kcal']} | dG_phys={r['dG_phys_kcal']} | "
              f"SB={r['salt_bridges']} HB={r['h_bonds']} HP={r['hydrophobic_contacts']} pi={r['pi_stacking']} clash={r['steric_clashes']} pro={r['helical_prolines']} | "
              f"wBSA={r['wBSA_A2']} mean_w={r['mean_w_intf']} pDockQ_coord={r['pDockQ_coord']}")

lines = Path("unippi_cli.py").read_text(encoding="utf-8").splitlines()
print("\n" + "=" * 85)
print("unippi_cli.py: build_residues() & atom-loop & UIS lines")
for title, a, b in [("BUILD_RESIDUES", 390, 445), ("ATOM_LOOP", 770, 845), ("UIS_RETURN", 1150, 1185)]:
    print(f"--- {title} (L{a}..L{b}) ---")
    for idx in range(a, min(b, len(lines)) + 1):
        print(f"L{idx:4d}: {lines[idx-1]}")
