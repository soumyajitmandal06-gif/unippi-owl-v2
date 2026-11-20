#!/usr/bin/env python3
"""
Phase 4 - Step 7 (Final Pre-Release Master Polish — UniPPI-OWL v2.1-MARVAS):
  Closes the final 5 audit items before live GitHub Pages deployment:
    1. Fixes aromatic-exocyclic false positive (ANP N6) in get_ligand_cationic_nitrogens()
       by excluding aromatic ring carbons from the amidinium/guanidinium rule and enforcing
       sp3 C-N bond length (>=1.38 A) for aliphatic amines (both CLI & Browser).
    2. Surfaces Boltz-1 / AF3 summary-only ipTM in Leaderboard, Track 1 & CSV/TSV exports
       even when PAE matrices are absent (MODE_B_PRED_COORD + hasJsonIptm).
    3. Clears stale PDB residue maps & file handles when the user manually edits Tab 2
       sequence boxes (`onManualSeqEdit`), preventing silent mutation/merge mismatches.
    4. Adds halogen-bond (`xb`) credit to Hotspot scoring (CLI & Browser) and adds
       biochemical WT-role guidance for ALA / GLY / RNA / Ligand rows in the Hotspot table.
    5. Fixes 4+/5-component notice (`entities.length >= 3`), displays active pH & NaCl in
       Step 4 math derivation, removes dead `toggleSlotZBox()`, unifies `v2.1-MARVAS`
       branding, and creates a clean `.gitignore` for GitHub publication.
"""

import ast
import builtins
import json
import re
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
html_path = Path("index.html")


def ensure_replaced(text: str, old: str, new: str, label: str, count: int = 1) -> str:
    if new and new in text:
        print(f"  [OK-ALREADY] {label}")
        return text
    if old in text:
        print(f"  [APPLIED]    {label}")
        return text.replace(old, new, count)
    raise AssertionError(f"Anchor mismatch for [{label}]:\n{old[:200]!r}")


def ensure_removed(text: str, old: str, label: str) -> str:
    if old not in text:
        print(f"  [OK-ALREADY] {label}")
        return text
    print(f"  [APPLIED]    {label}")
    return text.replace(old, "", 1)


# ==============================================================================
# PART A: unippi_cli.py (TIGHTEN CATIONIC LIGAND DETECTOR + XB HOTSPOT CREDIT)
# ==============================================================================
print("=== PART A: POLISHING unippi_cli.py ===")
cli = cli_path.read_text(encoding="utf-8")

# A1. Require sp3 C-N >= 1.38 A for aliphatic amines & exclude aromatic ring carbons from guanidinium/amidinium rule
cli = ensure_replaced(
    cli,
    '        if not all(ot.element == "C" for ot in nbrs):',
    '        if not all(ot.element == "C" and math.dist(na.coord, ot.coord) >= 1.38 for ot in nbrs):',
    "CLI: Enforce sp3 C-N bond length (>=1.38 A) for aliphatic cationic amines",
)
cli = ensure_replaced(
    cli,
    '''    for ca in c_atoms:
        bn = [na for na in n_atoms if 1.25 <= math.dist(ca.coord, na.coord) <= 1.42]''',
    '''    for ca in c_atoms:
        if ca.name in arom_c_names:
            continue
        bn = [na for na in n_atoms if 1.25 <= math.dist(ca.coord, na.coord) <= 1.42]''',
    "CLI: Exclude aromatic ring carbons (e.g. adenine C6) from amidinium/guanidinium rule",
)

# A2. Credit cation-pi & halogen bonds in res_energy_contrib
cli = ensure_replaced(
    cli,
    '''            for r_key in (ra.key, rb.key):
                rec = res_energy_contrib.setdefault(r_key, {"sb": 0.0, "hb": 0.0, "hp": 0.0, "pi": 0.0})
                rec["sb"] += w_ij * eff_sb
                rec["hb"] += w_ij * eff_hb
                rec["hp"] += w_ij * eff_hp_weighted
                rec["pi"] += w_ij * pair_pi''',
    '''            for r_key in (ra.key, rb.key):
                rec = res_energy_contrib.setdefault(r_key, {"sb": 0.0, "hb": 0.0, "hp": 0.0, "pi": 0.0, "xb": 0.0})
                rec["sb"] += w_ij * eff_sb
                rec["hb"] += w_ij * eff_hb
                rec["hp"] += w_ij * eff_hp_weighted
                rec["pi"] += w_ij * (pair_pi + pair_catpi)
                rec["xb"] += w_ij * pair_xb''',
    "CLI: Track cation-pi & halogen-bond credit in res_energy_contrib",
)

# Full-file AST Scope Audit on unippi_cli.py
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

fn_cnt = 0
for node in tree.body:
    if not isinstance(node, ast.FunctionDef):
        continue
    fn_cnt += 1
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
        raise AssertionError(f"AST SCOPE AUDIT FAILED in `{node.name}()`: {missing}")

cli_path.write_text(cli, encoding="utf-8")
print(f"[PASS] unippi_cli.py polished & verified across ALL {fn_cnt} functions (0 undefined names).\n")


# ==============================================================================
# PART B: index.html (SYNC ALL 5 FIXES + VERSION UNIFICATION TO v2.1-MARVAS)
# ==============================================================================
print("=== PART B: POLISHING index.html ===")
html = html_path.read_text(encoding="utf-8")

# B1. Tighten getLigandCationicNitrogensJS() to 100% parity with CLI
html = ensure_replaced(
    html,
    '        if (nbrs.length < 1 || nbrs.length > 4 || !nbrs.every(ot => ot.element === "C")) continue;',
    '        if (nbrs.length < 1 || nbrs.length > 4 || !nbrs.every(ot => ot.element === "C" && Math.hypot(na.x-ot.x, na.y-ot.y, na.z-ot.z) >= 1.38)) continue;',
    "JS: Enforce sp3 C-N bond length (>=1.38 A) for aliphatic cationic amines",
)
html = ensure_replaced(
    html,
    '''      for (const ca of cAtoms) {
        const bn = nAtoms.filter(na => { const d = Math.hypot(ca.x-na.x, ca.y-na.y, ca.z-na.z); return d >= 1.25 && d <= 1.42; });''',
    '''      for (const ca of cAtoms) {
        if (aromCNames.has(ca.name)) continue;
        const bn = nAtoms.filter(na => { const d = Math.hypot(ca.x-na.x, ca.y-na.y, ca.z-na.z); return d >= 1.25 && d <= 1.42; });''',
    "JS: Exclude aromatic ring carbons (e.g. adenine C6) from amidinium/guanidinium rule",
)

# B2. Track xb in resContrib + entityType in mutCatalog + upgrade Hotspot table (ALA/GLY/RNA/XB aware)
html = ensure_replaced(
    html,
    '''              if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0 });
              const rc = resContrib.get(rk);
              rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*(pPi + pCatPi);''',
    '''              if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0, xb:0 });
              const rc = resContrib.get(rk);
              rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*(pPi + pCatPi); rc.xb += wij*pXb;''',
    "JS: Track halogen-bond (xb) credit inside resContrib",
)
html = ensure_replaced(
    html,
    '          const cInfo = resContrib.get(r.key) || { sb:0, hb:0, hp:0, pi:0 };',
    '          const cInfo = resContrib.get(r.key) || { sb:0, hb:0, hp:0, pi:0, xb:0 };',
    "JS: Default xb:0 on mutCatalog cInfo",
)
html = ensure_replaced(
    html,
    '          mutCatalog.push({ key: r.key, label: `${r.chainId}:${cname}${r.resSeq}`, chainId: r.chainId, resSeq: r.resSeq, wtName: cname, cInfo, fBuried, inHel, coord: r.cbCoord });',
    '          mutCatalog.push({ key: r.key, label: `${r.chainId}:${cname}${r.resSeq}`, chainId: r.chainId, resSeq: r.resSeq, wtName: cname, entityType: r.entityType, cInfo, fBuried, inHel, coord: r.cbCoord });',
    "JS: Pass entityType into mutCatalog for biochemical hotspot annotation",
)

html = ensure_replaced(
    html,
    '                <th>Est. ΔΔG<sub>trunc</sub> (kcal/mol)</th>',
    '                <th>Est. ΔΔG<sub>pert</sub> (kcal/mol)</th>',
    "HTML: Clarify Hotspot perturbation header (covers Protein, RNA & Ligand)",
)

old_hs_fn = '''    function computeHotspotList(r) {
      if (!r || !r.mutCatalog) return [];
      const rows = [];
      for (const m of r.mutCatalog) {
        const ci = m.cInfo || { sb:0, hb:0, hp:0, pi:0 };
        const ddgEst = 0.85 * ci.sb + 0.55 * ci.hb + 0.85 * ci.pi + 0.14 * ci.hp + 1.15 * clamp(m.fBuried, 0, 1);
        if (ddgEst < 0.35 && m.fBuried < 0.10) continue;
        const tier = ddgEst >= 2.0 ? "CRITICAL HOTSPOT" : (ddgEst >= 1.0 ? "STRONG HOTSPOT" : "SUPPORTING");
        rows.push({
          label: m.label, chainId: m.chainId, resSeq: m.resSeq, coord: m.coord,
          burPct: Number((m.fBuried * 100.0).toFixed(1)),
          sb: Number(ci.sb.toFixed(1)), hb: Number(ci.hb.toFixed(1)),
          hp: Number(ci.hp.toFixed(1)), pi: Number(ci.pi.toFixed(1)),
          ddgEst: Number(ddgEst.toFixed(2)), tier
        });
      }
      rows.sort((a, b) => b.ddgEst - a.ddgEst);
      return rows;
    }'''

new_hs_fn = '''    function computeHotspotList(r) {
      if (!r || !r.mutCatalog) return [];
      const rows = [];
      for (const m of r.mutCatalog) {
        const ci = m.cInfo || { sb:0, hb:0, hp:0, pi:0, xb:0 };
        const ddgEst = 0.85 * ci.sb + 0.55 * ci.hb + 0.85 * ci.pi + 0.55 * (ci.xb || 0) + 0.14 * ci.hp + 1.15 * clamp(m.fBuried, 0, 1);
        if (ddgEst < 0.35 && m.fBuried < 0.10) continue;
        const tier = ddgEst >= 2.0 ? "CRITICAL HOTSPOT" : (ddgEst >= 1.0 ? "STRONG HOTSPOT" : "SUPPORTING");
        const mutHint = m.entityType === "NUCLEIC" ? " [nt sub]" : (m.entityType === "LIGAND" ? " [ligand]" : (m.wtName === "ALA" ? " [WT=Ala→D/G]" : (m.wtName === "GLY" ? " [WT=Gly→D/W]" : "")));
        rows.push({
          label: `${m.label}${mutHint}`, chainId: m.chainId, resSeq: m.resSeq, coord: m.coord,
          burPct: Number((m.fBuried * 100.0).toFixed(1)),
          sb: Number(ci.sb.toFixed(1)), hb: Number(ci.hb.toFixed(1)),
          hp: Number(ci.hp.toFixed(1)), pi: Number((ci.pi + (ci.xb || 0)).toFixed(1)),
          ddgEst: Number(ddgEst.toFixed(2)), tier
        });
      }
      rows.sort((a, b) => b.ddgEst - a.ddgEst);
      return rows;
    }'''
html = ensure_replaced(html, old_hs_fn, new_hs_fn, "JS: Add XB credit & ALA/GLY/RNA/Ligand annotations to computeHotspotList()")

# B3. Surface Boltz-1 / summary-only JSON ipTM in Leaderboard, Track 1 & CSV/TSV
html = ensure_replaced(
    html,
    '        ipTM: Number(iptm.toFixed(3)), LIS: Number(lis.toFixed(3)), pDockQ_PAE: Number(pdPae.toFixed(3)), pDockQ_coord: Number(pdCoord.toFixed(3)),',
    '        ipTM: Number(iptm.toFixed(3)), hasJsonIptm: confData.iptm !== null, LIS: Number(lis.toFixed(3)), pDockQ_PAE: Number(pdPae.toFixed(3)), pDockQ_coord: Number(pdCoord.toFixed(3)),',
    "JS: Record hasJsonIptm flag for Boltz-1 & AF3 summary-only JSONs",
)
html = ensure_replaced(
    html,
    '          <td>${r.mode === "MODE_A_PRED_FULL" ? r.ipTM.toFixed(2) : "—"}</td>',
    '          <td>${(r.mode === "MODE_A_PRED_FULL" || r.hasJsonIptm) ? r.ipTM.toFixed(2) : "—"}</td>',
    "JS: Display Boltz-1/summary JSON ipTM in Leaderboard",
)
html = ensure_replaced(
    html,
    '              ? `ipTM = N/A (No JSON) · Mean Intf pLDDT = ${r.mean_intf_pLDDT.toFixed(1)} · Disordered (pLDDT &lt; 70) = ${r.pct_disordered ?? 0}%`',
    '              ? `ipTM = ${r.hasJsonIptm ? `${r.ipTM.toFixed(3)} (Summary JSON)` : "N/A (No JSON)"} · Mean Intf pLDDT = ${r.mean_intf_pLDDT.toFixed(1)} · Disordered (pLDDT &lt; 70) = ${r.pct_disordered ?? 0}%`',
    "JS: Display Boltz-1/summary JSON ipTM in Track 1",
)
html = ensure_replaced(
    html,
    '          if ((c === "ipTM" || c === "LIS") && r.mode !== "MODE_A_PRED_FULL") return "NA";',
    '          if (c === "LIS" && r.mode !== "MODE_A_PRED_FULL") return "NA";\n'
    '          if (c === "ipTM" && r.mode !== "MODE_A_PRED_FULL" && !r.hasJsonIptm) return "NA";',
    "JS: Preserve Boltz-1/summary JSON ipTM in CSV/TSV exports",
)

# B4. Display active pH & Ionic Strength in Step 4 math box
html = ensure_replaced(
    html,
    '''        const tempK = parseFloat(document.getElementById("tempInput")?.value) || 298.15;
        const lambdaPae = parseFloat(document.getElementById("lambdaInput")?.value) || 8.0;
        const cutoffA = parseFloat(document.getElementById("cutoffInput")?.value) || 5.5;''',
    '''        const tempK = parseFloat(document.getElementById("tempInput")?.value) || 298.15;
        const lambdaPae = parseFloat(document.getElementById("lambdaInput")?.value) || 8.0;
        const cutoffA = parseFloat(document.getElementById("cutoffInput")?.value) || 5.5;
        const phDisp = clamp(parseFloat(document.getElementById("phInput")?.value || "7.2") || 7.2, 4.0, 9.5);
        const saltDisp = clamp(parseFloat(document.getElementById("ionicInput")?.value || "150") || 150, 10, 1000);''',
    "JS: Read active pH & NaCl for Step 4 derivation display",
)
html = ensure_replaced(
    html,
    '<strong>Step 4 (Thermodynamic Free Energy Equations &amp; UIS Score at T = ${tempK} K):</strong><br>',
    '<strong>Step 4 (Thermodynamic Free Energy Equations &amp; UIS Score at T = ${tempK} K, pH = ${phDisp.toFixed(1)}, I = ${saltDisp} mM NaCl):</strong><br>',
    "JS: Show active pH & NaCl inside Step 4 math box",
)

# B5. Clear stale PDB handles & residue numbers when user manually edits Tab 2 sequences
html = ensure_replaced(
    html,
    '''oninput="refreshSlotBadge('X')"''',
    '''oninput="onManualSeqEdit('X')"''',
    "HTML: Wire onManualSeqEdit('X') on Slot X textarea",
)
html = ensure_replaced(
    html,
    '''oninput="refreshSlotBadge('Y')"''',
    '''oninput="onManualSeqEdit('Y')"''',
    "HTML: Wire onManualSeqEdit('Y') on Slot Y textarea",
)
old_toggle_z = '''    function toggleSlotZBox() {
      const en = document.getElementById("slotZEnable")?.checked;
      const box = document.getElementById("slotZFields");
      if (box) box.style.display = en ? "grid" : "none";
    }'''
new_manual_edit = '''    function onManualSeqEdit(slot) {
      STUDIO_SLOT_RESNUMS[slot] = null;
      STUDIO_SLOT_FILES[slot] = null;
      STUDIO_SLOT_PARSED[slot] = null;
      const pick = document.getElementById(`slot${slot}ChainPick`);
      if (pick) pick.style.display = "none";
      if (slot === "X") {
        const iw = document.getElementById("slotXIntraWrap");
        if (iw) iw.style.display = "none";
      }
      const mergeBtn = document.getElementById("btnMergeXYStructs");
      if (mergeBtn) mergeBtn.style.display = (STUDIO_SLOT_FILES.X && STUDIO_SLOT_FILES.Y) ? "inline-flex" : "none";
      refreshSlotBadge(slot);
    }'''
html = ensure_replaced(html, old_toggle_z, new_manual_edit, "JS: Replace dead toggleSlotZBox() with onManualSeqEdit() stale-PDB guard")

# B6. Fix 4+/5-component notice in Tab 2 (`entities.length >= 3`)
html = ensure_replaced(
    html,
    '''      if (entities.length === 3) {
        notes.push(`Ternary 3-component job configured (${entities.map(e=>`${e.slot}:${e.type.toUpperCase()}`).join(" + ")}) — Tab 1 will automatically trigger Ternary Molecular-Glue & Metal-Bridge auditing on the resulting 3D structure.`);
      }''',
    '''      if (entities.length >= 3) {
        notes.push(`Multi-component (${entities.length}-entity) complex configured (${entities.map(e=>`${e.slot}:${e.type.toUpperCase()}`).join(" + ")}) — Tab 1 will automatically trigger Ternary Molecular-Glue & Metal-Bridge auditing on the resulting 3D structure.`);
      }''',
    "JS: Show multi-component ternary glue/bridge notice for 3, 4, 5+ entities",
)

# B7. Unify version branding to v2.1-MARVAS across index.html
for old_v, new_v, lbl in [
    ('const MODEL_VERSION = "UniPPI-v2.0-MARVAS";', 'const MODEL_VERSION = "UniPPI-v2.1-MARVAS";', "JS: Bump MODEL_VERSION to UniPPI-v2.1-MARVAS"),
    ('UniPPI-OWL v2.0 (MARVAS)', 'UniPPI-OWL v2.1 (MARVAS)', "HTML: Update page title & metadata to v2.1"),
    ('<span>v2.0-MARVAS</span>', '<span>v2.1-MARVAS</span>', "HTML: Update nav badge to v2.1-MARVAS"),
    ('>v2.0-MARVAS</span></h1>', '>v2.1-MARVAS</span></h1>', "HTML: Update hero version to v2.1-MARVAS"),
    ('UniPPI-v2.0-MARVAS', 'UniPPI-v2.1-MARVAS', "HTML: Update model strings to UniPPI-v2.1-MARVAS"),
    ('UniPPI-OWL v2.0-MARVAS', 'UniPPI-OWL v2.1-MARVAS', "HTML: Update footer & scope strings to v2.1-MARVAS"),
]:
    if old_v in html:
        html = html.replace(old_v, new_v)
        print(f"  [APPLIED]    {lbl}")

html_path.write_text(html, encoding="utf-8")

# Create clean .gitignore for GitHub repo publication
gitignore_txt = """__pycache__/
*.pyc
*.pyo
.DS_Store
Thumbs.db
unippi_results/
*_out/
*.tmp
"""
Path(".gitignore").write_text(gitignore_txt, encoding="utf-8")
print("  [APPLIED]    Repo hygiene: wrote clean .gitignore")

# Verify ALL inline HTML JS event handlers are bound cleanly
m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
assert m_js, "Missing inline <script> block"
handlers = set(re.findall(r'on(?:click|change|input)="([A-Za-z_]\w*)\(', html))
defined_fns = set(re.findall(r'\bfunction\s+([A-Za-z_]\w*)\s*\(', m_js[-1]))
missing_js = sorted(handlers - defined_fns)
assert not missing_js, f"Undefined JS event handlers: {missing_js}"
print(f"[PASS] index.html polished cleanly: all {len(handlers)} event handlers bound.\n")


# ==============================================================================
# PART C: FULL 10-TEST PRE-RELEASE VERIFICATION (4 STRESS + 6 BENCHMARKS)
# ==============================================================================
print("=== PART C: 4 STRESS TESTS (INCLUDING ANP ZERO-FALSE-POSITIVE & GUANIDINIUM) ===")
sys.path.insert(0, ".")
import importlib
import unippi_cli as u
importlib.reload(u)

# Stress Test 0: Verify 5ZC9 ANP N6 false positive is ELIMINATED!
atoms_5zc9 = u.parse_structure_file("benchmarks/5ZC9.cif")[0]
ord_5zc9, _ = u.build_residues(atoms_5zc9)
for r in ord_5zc9:
    if r.entity_type == "LIGAND":
        cat_set = u.get_ligand_cationic_nitrogens(r)
        assert len(cat_set) == 0, f"False positive cationic N on {r.res_name}: {cat_set}"
print("  [PASS] Stress Test 0 (5ZC9 ANP & RCG Zero False-Positive Cationic N): ANP=[], RCG=[]")

# Stress Test 1: Aliphatic basic amine drug + RNA (SB >= 1, Cat-pi >= 1)
subprocess.run([sys.executable, "unippi_cli.py", "-i", "/tmp/st1_amine_rna.pdb", "-o", "/tmp/st1_out"], check=True, stdout=subprocess.DEVNULL)
st1_r = json.loads(Path("/tmp/st1_out/unippi_full_report.json").read_text(encoding="utf-8"))[0]
assert st1_r["salt_bridges"] >= 1 and st1_r["cation_pi"] >= 1, f"Stress Test 1 failed: {st1_r}"
print(f"  [PASS] Stress Test 1 (Aliphatic Amine Drug-RNA): SB={st1_r['salt_bridges']} Cat-pi={st1_r['cation_pi']}")

# Stress Test 2: Acyclic Guanidinium drug + Protein ASP carboxylate (true guanidinium MUST fire!)
st4_pdb = """EXPDTA    X-RAY DIFFRACTION
ATOM      1  CA  ASP A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CG  ASP A   1       1.500   0.000   0.000  1.00 20.00           C
ATOM      3  OD1 ASP A   1       2.100   1.050   0.000  1.00 20.00           O
ATOM      4  OD2 ASP A   1       2.100  -1.050   0.000  1.00 20.00           O
HETATM    5  C1  GND B   1       5.700   0.000   0.000  1.00 22.00           C
HETATM    6  N1  GND B   1       5.035   1.152   0.000  1.00 22.00           N
HETATM    7  N2  GND B   1       5.035  -1.152   0.000  1.00 22.00           N
HETATM    8  N3  GND B   1       7.030   0.000   0.000  1.00 22.00           N
HETATM    9  C2  GND B   1       7.750   1.250   0.000  1.00 22.00           C
END
"""
Path("/tmp/st4_guan_asp.pdb").write_text(st4_pdb, encoding="utf-8")
subprocess.run([sys.executable, "unippi_cli.py", "-i", "/tmp/st4_guan_asp.pdb", "-o", "/tmp/st4_out"], check=True, stdout=subprocess.DEVNULL)
st4_r = json.loads(Path("/tmp/st4_out/unippi_full_report.json").read_text(encoding="utf-8"))[0]
assert st4_r["salt_bridges"] >= 2 and st4_r["pct_disordered"] == 0.0, f"Stress Test 2 failed: {st4_r}"
print(f"  [PASS] Stress Test 2 (Acyclic Guanidinium Drug-ASP): SB={st4_r['salt_bridges']} dG={st4_r['dG_uni_kcal']}")

# Stress Test 3: Boltz-1 confidence_*.json (`complex_iptm` + `pair_chains_iptm` dict)
Path("/tmp/st5_boltz").mkdir(exist_ok=True)
Path("/tmp/st5_boltz/my_boltz_model_0.pdb").write_text(st4_pdb.replace("EXPDTA    X-RAY DIFFRACTION", "REMARK   ALPHAFOLD"), encoding="utf-8")
Path("/tmp/st5_boltz/confidence_my_boltz_model_0.json").write_text(
    json.dumps({"complex_iptm": 0.88, "pair_chains_iptm": {"0": {"0": 0.9, "1": 0.84}, "1": {"0": 0.84, "1": 0.9}}}),
    encoding="utf-8",
)
subprocess.run([sys.executable, "unippi_cli.py", "-i", "/tmp/st5_boltz/my_boltz_model_0.pdb", "-o", "/tmp/st5_out"], check=True, stdout=subprocess.DEVNULL)
st5_r = json.loads(Path("/tmp/st5_out/unippi_full_report.json").read_text(encoding="utf-8"))[0]
assert abs(st5_r["ipTM"] - 0.84) < 0.01, f"Stress Test 3 (Boltz-1 JSON) failed: {st5_r}"
print(f"  [PASS] Stress Test 3 (Boltz-1 confidence_*.json): parsed pair_chains_iptm={st5_r['ipTM']}")

print("\n=== PART D: ALL 6 LOCKED BENCHMARKS FINAL PRE-RELEASE VERIFICATION ===")
runs = [
    ("1. VIP1-bZIP61 PPI (AF3 5-model .zip)", "benchmarks/fold_vip1_bzip61.zip", "/tmp/v21_fin_vip1", None, 0.566, -10.99),
    ("2. 1URN U1A-RNA ITC Benchmark", "benchmarks/1URN.cif", "/tmp/v21_fin_1urn", None, 0.983, -13.61),
    ("3. 5ZC9 MARVAS eIF4A1+RNA+RocA+ANP+Mg", "benchmarks/5ZC9.cif", "/tmp/v21_fin_5zc9", None, 0.958, -12.88),
    ("4. 4WNP MARVAS ULK1+3RJ Kinase", "benchmarks/4WNP.cif", "/tmp/v21_fin_4wnp", None, 0.728, -8.41),
    ("5. 7DIY MARVAS Nsp14-Nsp10+Zn+Mg", "benchmarks/7DIY.cif", "/tmp/v21_fin_7diy", None, 0.847, -10.62),
    ("6. 5ZC9 Intra-Domain Split (RecA1/RecA2)", "benchmarks/5ZC9.cif", "/tmp/v21_fin_split", "A:1-236,237-406", 0.907, -11.80),
]
for label, inp, outp, sp, exp_uis, exp_dg in runs:
    cmd = [sys.executable, "unippi_cli.py", "-i", inp, "-o", outp]
    if sp:
        cmd += ["--split-chain", sp]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    rows = json.loads((Path(outp) / "unippi_full_report.json").read_text(encoding="utf-8"))
    top = rows[0]
    assert abs(top["UIS"] - exp_uis) < 0.002 and abs(top["dG_uni_kcal"] - exp_dg) < 0.02, (
        f"Regression on {label}: got UIS={top['UIS']} dG={top['dG_uni_kcal']} vs expected {exp_uis} / {exp_dg}"
    )
    print(f"  [PASS] {label:<42} | Top={top['chain_pair']:<18} ({top['modality']:<15}) UIS={top['UIS']:.3f} dG={top['dG_uni_kcal']:>6.2f} {top['pattern_id']}")

print("\n[FINAL PRE-RELEASE AUDIT PASSED 100%] ZERO UNRESOLVED FLAWS — READY FOR GITHUB PAGES DEPLOYMENT!")
