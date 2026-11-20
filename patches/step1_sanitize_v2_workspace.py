#!/usr/bin/env python3
"""
Phase 2 - Step 1: Safe Workspace Sanitization & v1.8 Baseline Freeze
  1. Verifies CWD is strictly 'UNI-PPI-OWL-v2-MARVAS' (never touches v1.8).
  2. Verifies index.html and unippi_cli.py exist before deleting anything.
  3. Detaches copied .git folder & purges stale v1.8 test folders/patch scripts.
  4. Organizes clean structure: baseline_v18_locked/, benchmarks/, patches/.
  5. Runs syntax + v1.8 baseline sanity check on unippi_cli.py.
"""

import ast
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

cwd = Path.cwd()
assert cwd.name == "UNI-PPI-OWL-v2-MARVAS", f"SAFETY ABORT: Run inside UNI-PPI-OWL-v2-MARVAS, not {cwd}"
assert (cwd / "index.html").exists(), "Missing index.html in UNI-PPI-OWL-v2-MARVAS!"
assert (cwd / "unippi_cli.py").exists(), "Missing unippi_cli.py in UNI-PPI-OWL-v2-MARVAS!"

# 1. Preserve essential files & regression zip
KEEP_ROOT = {
    "index.html",
    "unippi_cli.py",
    ".nojekyll",
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "requirements.txt",
    "baseline_v18_locked",
    "benchmarks",
    "patches",
    "step1_sanitize_v2_workspace.py",
}

(cwd / "baseline_v18_locked").mkdir(exist_ok=True)
(cwd / "benchmarks").mkdir(exist_ok=True)
(cwd / "patches").mkdir(exist_ok=True)

# Snapshot untouched v1.8 baseline for future diffs
shutil.copy2(cwd / "index.html", cwd / "baseline_v18_locked" / "index_v18.html")
shutil.copy2(cwd / "unippi_cli.py", cwd / "baseline_v18_locked" / "unippi_cli_v18.py")

# Rescue fold_vip1_bzip61.zip into benchmarks/ as our PPI zero-regression guard
reg_zip_dst = cwd / "benchmarks" / "fold_vip1_bzip61.zip"
if not reg_zip_dst.exists():
    for src in (cwd / "fold_vip1_bzip61.zip", cwd.parent / "UNI-PPI OWL" / "fold_vip1_bzip61.zip"):
        if src.exists():
            shutil.copy2(src, reg_zip_dst)
            break

def force_remove(path: Path):
    if path.is_dir() and not path.is_symlink():
        for root, dirs, files in os.walk(path, topdown=False):
            for f in files:
                fp = Path(root) / f
                try:
                    os.chmod(fp, stat.S_IWRITE | stat.S_IREAD)
                except OSError:
                    pass
            for d in dirs:
                dp = Path(root) / d
                try:
                    os.chmod(dp, stat.S_IRWXU)
                except OSError:
                    pass
        shutil.rmtree(path)
    else:
        try:
            os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
        except OSError:
            pass
        path.unlink()

removed_items = []
for item in sorted(cwd.iterdir()):
    if item.name in KEEP_ROOT:
        continue
    force_remove(item)
    removed_items.append(item.name)

# Move this setup script into patches/ so root stays spotless
shutil.move(str(cwd / "step1_sanitize_v2_workspace.py"), str(cwd / "patches" / "step1_sanitize_v2_workspace.py"))

# 2. Verify locked baseline integrity
cli_src = (cwd / "unippi_cli.py").read_text(encoding="utf-8")
html_src = (cwd / "index.html").read_text(encoding="utf-8")
ast.parse(cli_src)
assert "if (hasOrdCore && platConf(r) < 0.35) continue;" in html_src, "index.html is missing Step 11 parity lock!"
assert 'conf_contact_log = [c for c in residue_contact_log if c["weight_w"] >= 0.25]' in cli_src, "unippi_cli.py is missing Step 12 parity lock!"

print(f"[PASS] Sanitized UNI-PPI-OWL-v2-MARVAS (removed {len(removed_items)} stale items including .git if present).")
print("[PASS] Verified locked v1.8 parity baseline in both index.html and unippi_cli.py.")
print("\n[CLEAN WORKSPACE INVENTORY]")
for p in sorted(cwd.iterdir()):
    tag = "/" if p.is_dir() else f" ({p.stat().st_size:,} bytes)"
    print(f"  - {p.name}{tag}")
if reg_zip_dst.exists():
    print(f"  - benchmarks/{reg_zip_dst.name} ({reg_zip_dst.stat().st_size:,} bytes) [PPI Zero-Regression Guard]")

# 3. Quick CLI baseline check on fold_vip1_bzip61.zip (writes to /tmp so workspace stays clean)
if reg_zip_dst.exists():
    print("\n[BASELINE REGRESSION CHECK: v1.8 on fold_vip1_bzip61.zip]")
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(reg_zip_dst), "-o", "/tmp/unippi_v2_step1_check"], check=True)
