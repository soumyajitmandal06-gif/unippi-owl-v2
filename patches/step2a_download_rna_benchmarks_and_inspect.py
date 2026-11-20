#!/usr/bin/env python3
"""
Phase 2 - Step 2A (Read-Only Inspection + Real Experimental RNA BEFORE Benchmark):
  1. Downloads RCSB PDB 5ZC9.cif (eIF4A1 + Polypurine RNA + Rocaglamide A) and
     1URN.cif (U1A RRM + Hairpin RNA, ITC ground truth dG ~ -13.8 kcal/mol) into benchmarks/.
  2. Runs unmodified v1.8 unippi_cli.py on both to record the exact BEFORE baseline.
  3. Prints exact anchor lines in unippi_cli.py for a zero-error Step 2B patch.
"""

import subprocess
import sys
import urllib.request
from pathlib import Path

bench_dir = Path("benchmarks")
bench_dir.mkdir(exist_ok=True)

for pdb_id, desc in [
    ("5ZC9", "Human eIF4A1 + Polypurine RNA + Rocaglamide A (MARVAS Rocaglate/Aglain Reference)"),
    ("1URN", "Human U1A RRM + Hairpin RNA (Experimental ITC Ground Truth: Kd ~ 0.03-0.1 nM, dG ~ -14.0 kcal/mol)"),
]:
    dst = bench_dir / f"{pdb_id}.cif"
    if not dst.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        print(f"[DOWNLOAD] {pdb_id}.cif — {desc}")
        urllib.request.urlretrieve(url, dst)
        print(f"  Saved {dst} ({dst.stat().st_size:,} bytes)")
    else:
        print(f"[CACHED]   {dst} ({dst.stat().st_size:,} bytes)")

print("\n" + "="*95)
print("BEFORE BENCHMARK (v1.8 Engine on 5ZC9.cif — eIF4A1 + RNA + Rocaglamide A)")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "-i", str(bench_dir / "5ZC9.cif"), "-o", "/tmp/v18_before_5zc9"], check=True)

print("\n" + "="*95)
print("BEFORE BENCHMARK (v1.8 Engine on 1URN.cif — U1A + Hairpin RNA)")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "-i", str(bench_dir / "1URN.cif"), "-o", "/tmp/v18_before_1urn"], check=True)

print("\n" + "="*95)
print("EXACT ANCHOR LINES IN unippi_cli.py FOR STEP 2B PATCH")
print("="*95)
lines = Path("unippi_cli.py").read_text(encoding="utf-8").splitlines()
for i, l in enumerate(lines, 1):
    if any(k in l for k in (
        "NUCLEOTIDES =", "AROMATIC_RING", "def get_aromatic_centroid",
        "def analyze_chain_pair", "if d <= 4.0 and", "if d <= 3.5 and",
        "pair_pi =", "dg_phys =", "modality =", "return {"
    )):
        print(f"L{i:4d}: {l}")
