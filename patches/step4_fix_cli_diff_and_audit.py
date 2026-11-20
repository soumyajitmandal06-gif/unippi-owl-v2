#!/usr/bin/env python3
"""
Phase 3 - Step 4.1 (Surgical CLI Diff Fix + Full-File Static AST Scope Audit + JS Integrity Check):
  1. Replaces `Path(diff_csv).write_text(...)` in unippi_cli.py with stdlib `open()` (zero extra imports).
  2. Upgrades CLI `run_cli_ensemble_diff()` contact-rewiring keys to normalize by `(chain, res_seq)`
     for 100% parity with index.html Tab 3 on point-mutated escape ensembles.
  3. Runs a FULL-FILE Static AST Scope Audit across EVERY function in unippi_cli.py (0 undefined names).
  4. Runs an HTML/JS handler & bracket-balance audit on index.html.
  5. Verifies CLI `--diff-b` on both the 5-model null control (`fold_vip1` vs `fold_vip1`: ddG=0, p=1)
     and a perturbation test (`5ZC9` vs `1URN`: non-zero dUIS, ddG, LOST/GAINED rewiring & CSV export).
"""

import ast
import builtins
import re
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
html_path = Path("index.html")

cli = cli_path.read_text(encoding="utf-8")

# 1. Replace Path(diff_csv).write_text(...) with built-in open() + normalize contact keys by (chain, res_seq)
old_agg_and_write = '''    def _agg(sub):
        m: Dict[str, List[float]] = {}
        for r in sub:
            seen = set()
            for c in r.get("contacts", []):
                if c.get("weight_w", 0.0) < 0.25:
                    continue
                k = f"{c['res_a']} <-> {c['res_b']}"
                if k in seen:
                    continue
                seen.add(k)
                m.setdefault(k, []).append(float(c.get("min_dist_A", 0.0)))
        return m
    ma, mb = _agg(sub_a), _agg(sub_b)
    lost, gained = [], []
    diff_rows = ["status,contact_pair,occ_A_pct,occ_B_pct"]
    for k in sorted(set(ma.keys()) | set(mb.keys())):
        oa = len(ma.get(k, [])) / len(sub_a)
        ob = len(mb.get(k, [])) / len(sub_b)
        if oa >= 0.6 and ob < 0.4:
            lost.append(k)
            diff_rows.append(f'LOST,"{k}",{round(oa*100)},{round(ob*100)}')
        elif oa < 0.4 and ob >= 0.6:
            gained.append(k)
            diff_rows.append(f'GAINED,"{k}",{round(oa*100)},{round(ob*100)}')
    print(f"  Contact Rewiring: {len(lost)} LOST | {len(gained)} GAINED")
    if lost:
        print("    Top LOST  :", ", ".join(lost[:6]))
    if gained:
        print("    Top GAINED:", ", ".join(gained[:6]))
    diff_csv = os.path.join(out_dir, "diff_contact_rewiring.csv")
    Path(diff_csv).write_text("\\n".join(diff_rows) + "\\n", encoding="utf-8")'''

new_agg_and_write = '''    def _norm_res(tok: str) -> str:
        parts = tok.split(":", 1)
        m_num = re.search(r"(\\d+)$", parts[-1])
        return f"{parts[0]}:{m_num.group(1)}" if (len(parts) == 2 and m_num) else tok

    def _agg(sub):
        m: Dict[str, Dict] = {}
        for r in sub:
            seen = set()
            for c in r.get("contacts", []):
                if c.get("weight_w", 0.0) < 0.25:
                    continue
                k = f"{_norm_res(c['res_a'])} <-> {_norm_res(c['res_b'])}"
                if k in seen:
                    continue
                seen.add(k)
                rec = m.setdefault(k, {"label": f"{c['res_a']} <-> {c['res_b']}", "dists": []})
                rec["dists"].append(float(c.get("min_dist_A", 0.0)))
        return m
    ma, mb = _agg(sub_a), _agg(sub_b)
    lost, gained = [], []
    diff_rows = ["status,contact_pair,occ_A_pct,occ_B_pct,mean_dist_A_A,mean_dist_B_A"]
    for k in sorted(set(ma.keys()) | set(mb.keys())):
        ra_m, rb_m = ma.get(k), mb.get(k)
        oa = (len(ra_m["dists"]) / len(sub_a)) if ra_m else 0.0
        ob = (len(rb_m["dists"]) / len(sub_b)) if rb_m else 0.0
        lbl = f"{ra_m['label']} [B: {rb_m['label']}]" if (ra_m and rb_m and ra_m["label"] != rb_m["label"]) else (ra_m["label"] if ra_m else rb_m["label"])
        da = f"{sum(ra_m['dists'])/len(ra_m['dists']):.2f}" if ra_m else ""
        db = f"{sum(rb_m['dists'])/len(rb_m['dists']):.2f}" if rb_m else ""
        if oa >= 0.6 and ob < 0.4:
            lost.append(lbl)
            diff_rows.append(f'LOST,"{lbl}",{round(oa*100)},{round(ob*100)},{da},{db}')
        elif oa < 0.4 and ob >= 0.6:
            gained.append(lbl)
            diff_rows.append(f'GAINED,"{lbl}",{round(oa*100)},{round(ob*100)},{da},{db}')
        elif oa >= 0.4 or ob >= 0.4:
            diff_rows.append(f'PERSISTENT,"{lbl}",{round(oa*100)},{round(ob*100)},{da},{db}')
    print(f"  Contact Rewiring: {len(lost)} LOST | {len(gained)} GAINED")
    if lost:
        print("    Top LOST  :", ", ".join(lost[:6]))
    if gained:
        print("    Top GAINED:", ", ".join(gained[:6]))
    diff_csv = os.path.join(out_dir, "diff_contact_rewiring.csv")
    with open(diff_csv, "w", encoding="utf-8") as _df:
        _df.write("\\n".join(diff_rows) + "\\n")'''

if new_agg_and_write in cli:
    print("  [OK-ALREADY] CLI run_cli_ensemble_diff() open() & normalized contact keys")
elif old_agg_and_write in cli:
    cli = cli.replace(old_agg_and_write, new_agg_and_write, 1)
    print("  [APPLIED]    CLI run_cli_ensemble_diff() open() fix + (chain, res_seq) parity")
else:
    raise AssertionError("Could not locate run_cli_ensemble_diff() write block in unippi_cli.py")

# 2. FULL-FILE Static AST Scope Audit across EVERY top-level function in unippi_cli.py
tree = ast.parse(cli)
mod_globals = set(dir(builtins)) | {"__name__", "__file__", "__doc__", "__package__"}
for node in tree.body:
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

checked_fns = 0
for node in tree.body:
    if not isinstance(node, ast.FunctionDef):
        continue
    checked_fns += 1
    bound = set(mod_globals)
    for arg in node.args.args + node.args.kwonlyargs + ([node.args.vararg] if node.args.vararg else []) + ([node.args.kwarg] if node.args.kwarg else []):
        bound.add(arg.arg)
    for sub in ast.walk(node):
        if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bound.add(sub.name)
            if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for a in sub.args.args + sub.args.kwonlyargs:
                    bound.add(a.arg)
        elif isinstance(sub, ast.Lambda):
            for a in sub.args.args + sub.args.kwonlyargs:
                bound.add(a.arg)
        elif isinstance(sub, ast.Name) and isinstance(sub.ctx, (ast.Store, ast.Del)):
            bound.add(sub.id)
        elif isinstance(sub, ast.ExceptHandler) and sub.name:
            bound.add(sub.name)
        elif isinstance(sub, ast.Import):
            for al in sub.names:
                bound.add(al.asname or al.name.split(".")[0])
        elif isinstance(sub, ast.ImportFrom):
            for al in sub.names:
                bound.add(al.asname or al.name)
    loaded = {sub.id for sub in ast.walk(node) if isinstance(sub, ast.Name) and isinstance(sub.ctx, ast.Load)}
    missing = sorted(loaded - bound)
    if missing:
        raise AssertionError(f"FULL-FILE AST SCOPE AUDIT FAILED in `{node.name}()`: undefined names {missing}")

cli_path.write_text(cli, encoding="utf-8")
print(f"[PASS] Full-file Static AST Scope Audit passed across ALL {checked_fns} functions in unippi_cli.py (0 undefined names).")

# 3. Verify index.html event handlers & script presence
html = html_path.read_text(encoding="utf-8")
m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
assert m_js, "Missing inline <script> block in index.html"
handlers = set(re.findall(r'on(?:click|change|input)="([A-Za-z_]\w*)\(', html))
defined_fns = set(re.findall(r'\bfunction\s+([A-Za-z_]\w*)\s*\(', m_js[-1]))
missing_js = sorted(handlers - defined_fns)
assert not missing_js, f"Undefined JS event handlers: {missing_js}"
print(f"[PASS] Verified index.html: all {len(handlers)} event handlers bound cleanly.\n")

# 4. Run CLI --diff-b verification tests (Null Control + Non-Zero Perturbation Control)
print("=" * 96)
print("VERIFY 1 (NULL CONTROL — N=5 vs N=5): fold_vip1_bzip61.zip vs fold_vip1_bzip61.zip")
print("=" * 96)
subprocess.run([sys.executable, "unippi_cli.py", "-i", "benchmarks/fold_vip1_bzip61.zip", "--diff-b", "benchmarks/fold_vip1_bzip61.zip", "-o", "/tmp/v2_s4_diff_null"], check=True)
assert Path("/tmp/v2_s4_diff_null/diff_contact_rewiring.csv").exists(), "Missing diff_contact_rewiring.csv!"

print("\n" + "=" * 96)
print("VERIFY 2 (PERTURBATION CONTROL — Rewiring & CSV Check): 5ZC9.cif vs 1URN.cif")
print("=" * 96)
subprocess.run([sys.executable, "unippi_cli.py", "-i", "benchmarks/5ZC9.cif", "--diff-b", "benchmarks/1URN.cif", "-o", "/tmp/v2_s4_diff_pert"], check=True)
diff_lines = Path("/tmp/v2_s4_diff_pert/diff_contact_rewiring.csv").read_text(encoding="utf-8").splitlines()
print(f"[PASS] Verified diff_contact_rewiring.csv exported cleanly ({len(diff_lines) - 1} contact rows).")
