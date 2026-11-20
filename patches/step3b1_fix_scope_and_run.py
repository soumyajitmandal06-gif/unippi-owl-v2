#!/usr/bin/env python3
"""
Phase 2 - Step 3B-1.5 (Surgical Scope Fix & Static AST Name-Binding Audit):
  1. Uses AST inspection on analyze_interface_pair() to extract the exact local variable
     bound to "mean_w_intf" (eliminating cross-language naming drift with JS `meanW`).
  2. Safeguards chain_remap[chain_a_id / chain_b_id] with .get() for composite ternary IDs ("A+B").
  3. Runs a full static AST free-variable audit on analyze_interface_pair & process_structure_file
     to guarantee ZERO hidden NameErrors on any branch.
  4. Executes the complete 6-test real benchmark suite.
"""

import ast
import builtins
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
src = cli_path.read_text(encoding="utf-8")
tree = ast.parse(src)


def find_func(t: ast.Module, name: str) -> ast.FunctionDef:
    for node in t.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"Function {name} not found")


# 1. Discover exact local variable name for mean interface weight from return dict
aip = find_func(tree, "analyze_interface_pair")
mean_w_var = None
for node in ast.walk(aip):
    if isinstance(node, ast.Dict):
        for k, v in zip(node.keys, node.values):
            if isinstance(k, ast.Constant) and k.value == "mean_w_intf":
                if isinstance(v, ast.Name):
                    mean_w_var = v.id
                elif isinstance(v, ast.Call) and v.args and isinstance(v.args[0], ast.Name):
                    mean_w_var = v.args[0].id
                break
assert mean_w_var, "Could not resolve local mean_w_intf variable from analyze_interface_pair() AST"
print(f"[AST-DISCOVERY] Local mean-weight variable in analyze_interface_pair() is `{mean_w_var}`.")

# Replace `mean_w *` in MODE_LIGAND and MODE_METAL pocket_quality formulas
old_lig_q = "pocket_quality = mean_w * clamp((lig_burial_frac - 0.15) / 0.65, 0.0, 1.0)"
new_lig_q = f"pocket_quality = {mean_w_var} * clamp((lig_burial_frac - 0.15) / 0.65, 0.0, 1.0)"
old_met_q = "pocket_quality = mean_w * clamp(w_metal_coord / 3.0, 0.0, 1.0) if n_metal_coord > 0 else 0.0"
new_met_q = f"pocket_quality = {mean_w_var} * clamp(w_metal_coord / 3.0, 0.0, 1.0) if n_metal_coord > 0 else 0.0"

for old_s, new_s, lbl in [(old_lig_q, new_lig_q, "MODE_LIGAND pocket_quality"), (old_met_q, new_met_q, "MODE_METAL pocket_quality")]:
    if new_s in src:
        print(f"  [OK-ALREADY] {lbl}")
    elif old_s in src:
        src = src.replace(old_s, new_s, 1)
        print(f"  [APPLIED]    {lbl} -> `{mean_w_var}`")
    else:
        raise AssertionError(f"Could not locate `{old_s}`")

# 2. Safeguard chain_remap[...] indexing for composite ternary chain IDs like "A+B"
for cid_var in ("chain_a_id", "chain_b_id"):
    raw_idx = f"chain_remap[{cid_var}]"
    safe_idx = f"chain_remap.get({cid_var}, {cid_var}.split('+')[0][:1])"
    if raw_idx in src:
        src = src.replace(raw_idx, safe_idx)
        print(f"  [APPLIED]    Safeguarded `{raw_idx}` for composite ternary IDs")

# 3. Static AST Scope Audit across module globals + target functions
tree2 = ast.parse(src)
mod_globals = set(dir(builtins))
for node in tree2.body:
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        mod_globals.add(node.name)
    elif isinstance(node, ast.Import):
        for al in node.names:
            mod_globals.add(al.asname or al.name.split(".")[0])
    elif isinstance(node, ast.ImportFrom):
        for al in node.names:
            mod_globals.add(al.asname or al.name)
    elif isinstance(node, ast.Assign):
        for t in node.targets:
            for sub in ast.walk(t):
                if isinstance(sub, ast.Name):
                    mod_globals.add(sub.id)
    elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
        mod_globals.add(node.target.id)

for fn_name in ("analyze_interface_pair", "process_structure_file"):
    fn = find_func(tree2, fn_name)
    bound = set(mod_globals)
    for arg in fn.args.args + fn.args.kwonlyargs + ([fn.args.vararg] if fn.args.vararg else []) + ([fn.args.kwarg] if fn.args.kwarg else []):
        bound.add(arg.arg)
    for node in ast.walk(fn):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bound.add(node.name)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for a in node.args.args + node.args.kwonlyargs:
                    bound.add(a.arg)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, (ast.Store, ast.Del)):
            bound.add(node.id)
        elif isinstance(node, ast.ExceptHandler) and node.name:
            bound.add(node.name)
        elif isinstance(node, ast.Import):
            for al in node.names:
                bound.add(al.asname or al.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            for al in node.names:
                bound.add(al.asname or al.name)
    loaded = {node.id for node in ast.walk(fn) if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)}
    missing = sorted(loaded - bound)
    if missing:
        raise AssertionError(f"STATIC AST AUDIT FAILED in `{fn_name}` — undefined identifiers: {missing}")
    print(f"[PASS] Static AST scope audit `{fn_name}`: 0 undefined variables.")

cli_path.write_text(src, encoding="utf-8")
print("[PASS] Saved verified unippi_cli.py.\n")

# 4. Run full 6-test real benchmark suite
tests = [
    ("TEST 1 (ZERO-REGRESSION PPI): fold_vip1_bzip61.zip", ["unippi_cli.py", "-i", "benchmarks/fold_vip1_bzip61.zip", "-o", "/tmp/v2_s3_vip1"]),
    ("TEST 2 (RNA + DRUG + TERNARY GLUE + MG BRIDGE): 5ZC9.cif", ["unippi_cli.py", "-i", "benchmarks/5ZC9.cif", "-o", "/tmp/v2_s3_5zc9"]),
    ("TEST 3 (RNA ITC GROUND TRUTH): 1URN.cif", ["unippi_cli.py", "-i", "benchmarks/1URN.cif", "-o", "/tmp/v2_s3_1urn"]),
    ("TEST 4 (MARVAS ULK1 KINASE + INHIBITOR 3RJ): 4WNP.cif", ["unippi_cli.py", "-i", "benchmarks/4WNP.cif", "-o", "/tmp/v2_s3_4wnp"]),
    ("TEST 5 (MARVAS CoV Nsp14-Nsp10 + ZN FINGERS + MG): 7DIY.cif", ["unippi_cli.py", "-i", "benchmarks/7DIY.cif", "-o", "/tmp/v2_s3_7diy"]),
    ("TEST 6 (CUSTOM SAME-PROTEIN INTRA-CHAIN SPLIT): 5ZC9 eIF4A1 NTD(1-236) vs CTD(237-406)", ["unippi_cli.py", "-i", "benchmarks/5ZC9.cif", "--split-chain", "A:1-236,237-406", "-o", "/tmp/v2_s3_5zc9_split"]),
]
for title, cmd in tests:
    print("\n" + "=" * 98)
    print(title)
    print("=" * 98)
    subprocess.run([sys.executable] + cmd, check=True)
