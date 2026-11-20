#!/usr/bin/env python3
"""
Phase 3 - Step 5 (Final Master Audit Upgrade — Fixes All 5 Deep-Audit Flaws in CLI & Browser):
  1. Fixes monomer-chop bug in Tab 2 (`buildResidues(atoms, autoSplitMonomer=true)`, passing
     `false` in Tab 2 so single-chain PDBs keep 100% of their sequence intact).
  2. Fixes `MODE_C_EXP_CRYSTAL` high-B ligand/ion artifact bug in BOTH `unippi_cli.py` and
     `index.html` (references single-residue ligands/ions to global polymer `bMinGlobal`)
     and routes `LIGAND-LIGAND` into `isLigandPair`.
  3. Fixes Chai-1 ionic SMILES (`[Mg+2]`, `[Zn+2]`, ...), keeps AF3 JSON strictly schema-valid
     (blocks invalid custom SMILES in AF3 JSON), and adds 1-click Boltz-1 CLI/Colab script export.
  4. Adds Live RCSB PDB ID Fetch (`e.g. 5ZC9`, `4WNP`, `7DIY`, `1URN`) in Tab 2 Slots X & Y,
     plus a "Merge Uploaded PDB/CIFs as Single X+Y Complex" engine in Tab 1 & Tab 2 (with an
     explicit >30 A un-docked crystal frame warning).
  5. Wires `splitChainInput` + exact JSON matching into Tab 3 Diff Mode, cleans dead `normKey`,
     and upgrades composite `A+B` SVG Contact Network to dual color-coded left-axis bars.
  6. Executes FULL 6-benchmark CLI audit + 34-function AST scope audit + JS handler audit.
"""

import ast
import builtins
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

cli_path = Path("unippi_cli.py")
html_path = Path("index.html")


def ensure_replaced(text: str, old: str, new: str, label: str, count: int = 1) -> str:
    if new in text:
        print(f"  [OK-ALREADY] {label}")
        return text
    if old in text:
        print(f"  [APPLIED]    {label}")
        return text.replace(old, new, count)
    raise AssertionError(f"Anchor mismatch for [{label}]:\n{old[:180]!r}")


def ensure_regex_replaced(text: str, pattern: str, new: str, label: str) -> str:
    if new in text:
        print(f"  [OK-ALREADY] {label}")
        return text
    m = re.search(pattern, text, flags=re.DOTALL)
    if m:
        print(f"  [APPLIED]    {label}")
        return text[:m.start()] + new + text[m.end():]
    raise AssertionError(f"Regex mismatch for [{label}]: {pattern}")


# ==============================================================================
# PART A: unippi_cli.py (HIGH-B LIGAND/ION GLOBAL BASELINE + LIGAND-LIGAND ROUTING)
# ==============================================================================
print("=== PART A: UPGRADING unippi_cli.py ===")
cli = cli_path.read_text(encoding="utf-8")

# A1. Route LIGAND-LIGAND into is_ligand_pair
cli = ensure_replaced(
    cli,
    '    is_ligand_pair = modality in ("PROTEIN-LIGAND", "RNA-LIGAND", "DNA-LIGAND", "TERNARY-GLUE")',
    '    is_ligand_pair = modality in ("PROTEIN-LIGAND", "RNA-LIGAND", "DNA-LIGAND", "TERNARY-GLUE", "LIGAND-LIGAND")',
    "CLI: Route LIGAND-LIGAND into MODE_LIGAND physics",
)

# A2. Fix MODE_C high-B single-residue ligand/ion B-factor baseline (use global polymer b_min)
old_bmin_py = '''    b_min_chain: Dict[str, float] = {}
    for _r in (residues_a + residues_b):
        _c = _r.chain_id
        if _c not in b_min_chain or _r.mean_bfactor < b_min_chain[_c]:
            b_min_chain[_c] = _r.mean_bfactor'''
new_bmin_py = '''    poly_bfacs = [_r.mean_bfactor for _r in (residues_a + residues_b) if _r.entity_type in ("PROTEIN", "NUCLEIC") and _r.atoms]
    b_min_global = min(poly_bfacs) if poly_bfacs else min((_r.mean_bfactor for _r in (residues_a + residues_b) if _r.atoms), default=10.0)
    b_min_chain: Dict[str, float] = {}
    for _r in (residues_a + residues_b):
        if _r.entity_type in ("LIGAND", "ION"):
            b_min_chain[_r.chain_id] = b_min_global
        elif _r.chain_id not in b_min_chain or _r.mean_bfactor < b_min_chain[_r.chain_id]:
            b_min_chain[_r.chain_id] = _r.mean_bfactor'''
if old_bmin_py in cli or new_bmin_py in cli:
    cli = ensure_replaced(cli, old_bmin_py, new_bmin_py, "CLI: Reference LIGAND/ION MODE_C B-factors to global polymer baseline")
else:
    # Fallback regex in case variable names differ slightly
    cli = ensure_regex_replaced(
        cli,
        r'    b_min_chain[^\n]*\n(?:    [^\n]*b_min_chain[^\n]*\n)+',
        new_bmin_py + '\n',
        "CLI: Reference LIGAND/ION MODE_C B-factors to global polymer baseline (regex)",
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
print(f"[PASS] unippi_cli.py upgraded & verified across ALL {fn_cnt} functions (0 undefined names).\n")

# ==============================================================================
# PART B: index.html (ALL 5 AUDIT FIXES + LIVE RCSB FETCH + TWO-FILE X+Y MERGE)
# ==============================================================================
print("=== PART B: UPGRADING index.html ===")
html = html_path.read_text(encoding="utf-8")

# B1. Add Merge-Files toggle in Tab 1 param-bar
html = ensure_replaced(
    html,
    '''        <label>Intra-Chain Domain Split (Optional)
          <input type="text" id="splitChainInput" placeholder="e.g. A:1-236,237-406" onchange="recomputeActiveSuite()">
        </label>''',
    '''        <label>Intra-Chain Domain Split (Optional)
          <input type="text" id="splitChainInput" placeholder="e.g. A:1-236,237-406" onchange="recomputeActiveSuite()">
        </label>
        <label style="justify-content:center;">Multi-File Docking Mode
          <span style="display:flex; align-items:center; gap:6px; font-weight:600; color:var(--text-main); padding:6px 0;">
            <input type="checkbox" id="mergeUploadedCheck" onchange="recomputeActiveSuite()" style="width:auto;">
            Merge uploaded PDB/CIFs as 1 X+Y complex
          </span>
        </label>''',
    "HTML: Add Merge-Files X+Y toggle in Tab 1 param-bar",
)

# B2. Add Live RCSB PDB ID Fetch inputs to Tab 2 Slot X and Slot Y + Merge X+Y Structures button
old_slot_x_top = '''            <div style="display:grid; grid-template-columns:1fr 1fr; gap:8px;">
              <select id="slotXType" onchange="updateStudioSlotUI('X')">
                <option value="protein">Protein (aa sequence / PDB)</option>
                <option value="rna">RNA (A/C/G/U / PDB)</option>
                <option value="dna">DNA (A/C/G/T / PDB)</option>
                <option value="ligand">Small Molecule (SMILES / CCD)</option>
                <option value="ion">Catalytic Metal / Ion (MG/ZN...)</option>
              </select>
              <label class="btn btn-secondary btn-sm" style="cursor:pointer; margin:0;">
                + Load Known .PDB/.CIF
                <input type="file" id="slotXFile" accept=".cif,.mmcif,.pdb,.ent" style="display:none" onchange="loadSlotStructureFile('X', this.files)">
              </label>
            </div>'''
new_slot_x_top = old_slot_x_top + '''
            <div style="display:flex; gap:6px;">
              <input type="text" id="slotXPdbId" placeholder="Or fetch RCSB PDB ID (e.g. 5ZC9, 4WNP)" style="flex:1;">
              <button class="btn btn-secondary btn-sm" onclick="fetchRcsbForSlot('X')">Fetch PDB</button>
            </div>'''
html = ensure_replaced(html, old_slot_x_top, new_slot_x_top, "HTML: Add live RCSB PDB ID fetch to Slot X")

old_slot_y_top = '''            <div style="display:grid; grid-template-columns:1fr 1fr; gap:8px;">
              <select id="slotYType" onchange="updateStudioSlotUI('Y')">
                <option value="rna">RNA (A/C/G/U / PDB)</option>
                <option value="protein">Protein (aa sequence / PDB)</option>
                <option value="dna">DNA (A/C/G/T / PDB)</option>
                <option value="ligand">Small Molecule (SMILES / CCD)</option>
                <option value="ion">Catalytic Metal / Ion (MG/ZN...)</option>
              </select>
              <label class="btn btn-secondary btn-sm" style="cursor:pointer; margin:0;">
                + Load Known .PDB/.CIF
                <input type="file" id="slotYFile" accept=".cif,.mmcif,.pdb,.ent" style="display:none" onchange="loadSlotStructureFile('Y', this.files)">
              </label>
            </div>'''
new_slot_y_top = old_slot_y_top + '''
            <div style="display:flex; gap:6px;">
              <input type="text" id="slotYPdbId" placeholder="Or fetch RCSB PDB ID (e.g. 1URN, 7DIY)" style="flex:1;">
              <button class="btn btn-secondary btn-sm" onclick="fetchRcsbForSlot('Y')">Fetch PDB</button>
            </div>'''
html = ensure_replaced(html, old_slot_y_top, new_slot_y_top, "HTML: Add live RCSB PDB ID fetch to Slot Y")

old_studio_btns = '''          <button class="btn btn-secondary btn-sm" onclick="exportStudioJob('chai')">Download Chai-1 .FASTA</button>
          <button class="btn btn-secondary btn-sm" onclick="copyStudioPreview()">Copy Preview</button>'''
new_studio_btns = '''          <button class="btn btn-secondary btn-sm" onclick="exportStudioJob('chai')">Download Chai-1 .FASTA</button>
          <button class="btn btn-secondary btn-sm" onclick="exportStudioJob('colab')">Download Boltz-1 Colab/CLI Runner (.sh)</button>
          <button id="btnMergeXYStructs" class="btn btn-secondary btn-sm" style="display:none; border-color:var(--emerald); color:var(--emerald);" onclick="mergeLoadedXYIntoTab1()">Merge Loaded X + Y PDBs &amp; Score in Tab 1</button>
          <button class="btn btn-secondary btn-sm" onclick="copyStudioPreview()">Copy Preview</button>'''
html = ensure_replaced(html, old_studio_btns, new_studio_btns, "HTML: Add Boltz-1 runner & Merge X+Y PDB buttons in Tab 2")

# B3. Fix buildResidues() monomer-chop bug via optional `autoSplitMonomer=true`
html = ensure_replaced(
    html,
    '    function buildResidues(atoms) {',
    '    function buildResidues(atoms, autoSplitMonomer=true) {',
    "JS: Add autoSplitMonomer parameter to buildResidues()",
)
html = ensure_replaced(
    html,
    '      if (Object.keys(splitChains).length === 1) {',
    '      if (autoSplitMonomer && Object.keys(splitChains).length === 1) {',
    "JS: Guard monomer auto-split with autoSplitMonomer flag",
)

# B4. Fix MODE_C high-B ligand/ion baseline & LIGAND-LIGAND routing in analyzePairJS()
old_bmin_js = '''      const bMinChain = {};
      for (const _r of [...resA, ...resB]) {
        if (bMinChain[_r.chainId] === undefined || _r.meanBfac < bMinChain[_r.chainId]) bMinChain[_r.chainId] = _r.meanBfac;
      }'''
new_bmin_js = '''      const polyB = [...resA, ...resB].filter(r => (r.entityType === "PROTEIN" || r.entityType === "NUCLEIC") && r.atoms.length).map(r => r.meanBfac);
      const bMinGlobal = polyB.length ? Math.min(...polyB) : Math.min(...[...resA, ...resB].map(r => r.meanBfac), 10.0);
      const bMinChain = {};
      for (const _r of [...resA, ...resB]) {
        if (_r.entityType === "LIGAND" || _r.entityType === "ION") bMinChain[_r.chainId] = bMinGlobal;
        else if (bMinChain[_r.chainId] === undefined || _r.meanBfac < bMinChain[_r.chainId]) bMinChain[_r.chainId] = _r.meanBfac;
      }'''
html = ensure_replaced(html, old_bmin_js, new_bmin_js, "JS: Reference LIGAND/ION MODE_C B-factors to global polymer baseline")

html = ensure_replaced(
    html,
    '      const isLigandPair = (modality === "PROTEIN-LIGAND" || modality === "RNA-LIGAND" || modality === "DNA-LIGAND" || modality === "TERNARY-GLUE");',
    '      const isLigandPair = (modality === "PROTEIN-LIGAND" || modality === "RNA-LIGAND" || modality === "DNA-LIGAND" || modality === "TERNARY-GLUE" || modality === "LIGAND-LIGAND");',
    "JS: Route LIGAND-LIGAND into MODE_LIGAND physics",
)

# B5. Support multi-file X+Y merging in Tab 1 handleFileSelect()
old_merge_hook = '''      structEntries.sort((a,b) => a.name.localeCompare(b.name));
      const queueLogs = [];
      ALL_RESULTS = [];'''
new_merge_hook = '''      structEntries.sort((a,b) => a.name.localeCompare(b.name));
      const mergeMode = document.getElementById("mergeUploadedCheck")?.checked;
      if (mergeMode && structEntries.length >= 2) {
        const mergedAtoms = [];
        const usedCids = new Set();
        const pool = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        const fileCentroids = [];
        for (const sf of structEntries) {
          const { atoms: fa } = parseStructureText(sf.text, sf.name);
          if (!fa.length) continue;
          const cx = meanArr(fa.map(a=>a.x)), cy = meanArr(fa.map(a=>a.y)), cz = meanArr(fa.map(a=>a.z));
          fileCentroids.push([sf.name, [cx, cy, cz]]);
          const localMap = {};
          for (const a of fa) {
            if (!localMap[a.chainId]) {
              const cid = !usedCids.has(a.chainId) ? a.chainId : (Array.from(pool).find(ch => !usedCids.has(ch)) || `M${usedCids.size}`);
              localMap[a.chainId] = cid;
              usedCids.add(cid);
            }
            a.chainId = localMap[a.chainId];
            mergedAtoms.push(a);
          }
        }
        if (mergedAtoms.length) {
          const { ordered: mOrd } = buildResidues(mergedAtoms, false);
          const mergedPdb = toSanitizedPDBString(mOrd);
          const mName = structEntries.map(s => s.name.replace(/\\.(cif|mmcif|pdb|ent|pqr)$/i, "")).join("_PLUS_").slice(0, 60) + ".pdb";
          structEntries.length = 0;
          structEntries.push({ name: mName, text: mergedPdb, mergedMeta: fileCentroids });
        }
      }
      const queueLogs = [];
      ALL_RESULTS = [];'''
html = ensure_replaced(html, old_merge_hook, new_merge_hook, "JS: Add Multi-File X+Y Complex Merge mode in Tab 1")

old_qlog_line = '''        queueLogs.push(
          `<div><strong>${sf.name}</strong> (${cIds.map(c=>`${c}:${chains[c].length}`).join(", ")}) &leftrightarrow; JSON: <strong>${pairedNames.length ? [...new Set(pairedNames)].join(" + ") : "None (Coordinate Mode)"}</strong> · Regime: <span style="color:var(--primary); font-weight:700;">${mode}</span></div>`
        );'''
new_qlog_line = old_qlog_line + '''
        if (sf.mergedMeta && sf.mergedMeta.length >= 2) {
          const dCent = dist3(sf.mergedMeta[0][1], sf.mergedMeta[1][1]);
          if (dCent > 35.0) {
            queueLogs.push(`<div style="color:var(--amber); font-weight:700;">[Spatial Frame Notice] Merged files (${sf.mergedMeta[0][0]} &amp; ${sf.mergedMeta[1][0]}) have centroids ${dCent.toFixed(1)} Å apart — if they came from separate crystal frames without docking/co-folding, use Tab 2 (X vs. Y Studio) to co-fold them into a single 3D complex first!</div>`);
          }
        }'''
html = ensure_replaced(html, old_qlog_line, new_qlog_line, "JS: Warn if merged separate PDBs are in un-docked crystal frames (>35 A)")

# B6. Dual color-coded left-axis bars on composite A+B SVG Contact Network
old_svg_bars = '''          <text x="120" y="22" text-anchor="middle" font-size="12" font-weight="bold" fill="#0f172a">Chain ${cA} (${minSeqA}–${maxSeqA})</text>
          <text x="340" y="22" text-anchor="middle" font-size="12" font-weight="bold" fill="#0f172a">Chain ${cB} (${minSeqB}–${maxSeqB})</text>
          <line x1="120" y1="35" x2="120" y2="${h-35}" stroke="#38bdf8" stroke-width="6" stroke-linecap="round"/>'''
new_svg_bars = '''          <text x="120" y="22" text-anchor="middle" font-size="12" font-weight="bold" fill="#0f172a">${isCompA ? `Composite ${cA}` : `Chain ${cA} (${minSeqA}–${maxSeqA})`}</text>
          <text x="340" y="22" text-anchor="middle" font-size="12" font-weight="bold" fill="#0f172a">Chain ${cB} (${minSeqB === maxSeqB ? minSeqB : `${minSeqB}–${maxSeqB}`})</text>
          ${isCompA
            ? `<line x1="120" y1="35" x2="120" y2="138" stroke="#38bdf8" stroke-width="6" stroke-linecap="round"/><text x="108" y="90" text-anchor="end" font-size="10" font-weight="bold" fill="#0284c7">${compChParts[0]}</text><line x1="120" y1="158" x2="120" y2="265" stroke="#f59e0b" stroke-width="6" stroke-linecap="round"/><text x="108" y="215" text-anchor="end" font-size="10" font-weight="bold" fill="#d97706">${compChParts[1] || ""}</text>`
            : `<line x1="120" y1="35" x2="120" y2="${h-35}" stroke="#38bdf8" stroke-width="6" stroke-linecap="round"/>`}'''
html = ensure_replaced(html, old_svg_bars, new_svg_bars, "JS: Dual color-coded left-axis bars on composite A+B SVG Contact Network")

# B7. Fix Tab 2 monomer-chop bug + add Live RCSB PDB Fetch + Merge X+Y button + valid Chai-1/AF3/Boltz-1 export
old_slot_load = '''    async function loadSlotStructureFile(slot, fileList) {
      if (!fileList || !fileList.length) return;
      const f = fileList[0];
      const text = await f.text();
      STUDIO_SLOT_FILES[slot] = { name: f.name, file: f, text };
      const { atoms } = parseStructureText(text, f.name);
      if (!atoms.length) { alert(`No valid atoms parsed from ${f.name}`); return; }
      const { chains } = buildResidues(atoms);
      STUDIO_SLOT_PARSED[slot] = chains;'''

new_slot_load = '''    function populateSlotFromText(slot, fname, text, fileObj=null) {
      STUDIO_SLOT_FILES[slot] = { name: fname, file: fileObj || new File([text], fname, { type: "text/plain" }), text };
      const { atoms } = parseStructureText(text, fname);
      if (!atoms.length) { alert(`No valid atoms parsed from ${fname}`); return; }
      const { chains } = buildResidues(atoms, false);
      STUDIO_SLOT_PARSED[slot] = chains;
      const pick = document.getElementById(`slot${slot}ChainPick`);
      const cIds = Object.keys(chains);
      if (pick) {
        pick.innerHTML = cIds.map(cid => {
          const rlist = chains[cid];
          const et = rlist[0]?.entityType || "PROTEIN";
          return `<option value="${cid}">Chain ${cid} (${rlist.length} ${et}) — from ${fname}</option>`;
        }).join("");
        pick.style.display = "block";
      }
      if (slot === "X") {
        const iw = document.getElementById("slotXIntraWrap");
        if (iw) iw.style.display = "block";
      }
      const mergeBtn = document.getElementById("btnMergeXYStructs");
      if (mergeBtn) mergeBtn.style.display = (STUDIO_SLOT_FILES.X && STUDIO_SLOT_FILES.Y) ? "inline-flex" : "none";
      applyPickedSlotChain(slot);
    }

    async function loadSlotStructureFile(slot, fileList) {
      if (!fileList || !fileList.length) return;
      const f = fileList[0];
      const text = await f.text();
      populateSlotFromText(slot, f.name, text, f);
    }

    async function fetchRcsbForSlot(slot) {
      const idIn = (document.getElementById(`slot${slot}PdbId`)?.value || "").trim().toUpperCase().replace(/[^A-Z0-9]/g, "");
      if (idIn.length !== 4) { alert("Enter a valid 4-character RCSB PDB ID (e.g. 5ZC9, 4WNP, 7DIY, 1URN)."); return; }
      const badge = document.getElementById(`slot${slot}Badge`);
      if (badge) badge.textContent = `Fetching ${idIn}...`;
      try {
        const resp = await fetch(`https://files.rcsb.org/download/${idIn}.cif`);
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
        const text = await resp.text();
        populateSlotFromText(slot, `${idIn}.cif`, text, null);
      } catch (err) {
        if (badge) badge.textContent = "Fetch failed";
        alert(`Could not fetch ${idIn}.cif from RCSB (${err.message}). If offline, download ${idIn}.cif manually and click + Load Known .PDB/.CIF.`);
      }
    }

    async function mergeLoadedXYIntoTab1() {
      if (!STUDIO_SLOT_FILES.X || !STUDIO_SLOT_FILES.Y) { alert("Load or fetch both Partner X and Partner Y structures first."); return; }
      const chk = document.getElementById("mergeUploadedCheck");
      if (chk) chk.checked = true;
      switchWorkspaceTab("auditor");
      await handleFileSelect([STUDIO_SLOT_FILES.X.file, STUDIO_SLOT_FILES.Y.file]);
    }'''

html = ensure_regex_replaced(
    html,
    r'    async function loadSlotStructureFile\(slot, fileList\) \{[\s\S]*?applyPickedSlotChain\(slot\);\n    \}',
    new_slot_load,
    "JS: Fix Tab 2 monomer-chop bug + add Live RCSB PDB Fetch & Merge X+Y into Tab 1",
)

# B8. Fix AF3 JSON (strictly schema-valid), Chai-1 ionic SMILES, & Boltz-1 Colab/CLI runner export
old_export_build = '''      // Build AF3 Server JSON, Boltz-1/2 YAML, and Chai-1 FASTA
      const jobName = entities.map(e => e.name).join("_vs_").slice(0, 64);
      const af3Seqs = entities.map(e => {
        if (e.type === "protein") return { proteinChain: { sequence: e.seq, count: 1 } };
        if (e.type === "rna") return { rnaSequence: { sequence: e.seq, count: 1 } };
        if (e.type === "dna") return { dnaSequence: { sequence: e.seq, count: 1 } };
        if (e.type === "ion") return { ion: { ion: e.seq.toUpperCase(), count: 1 } };
        const st = analyzeLigandStats(e.seq);
        return st.isCcd ? { ligand: { ligand: `CCD_${st.code}`, count: 1 } } : { ligand: { ligand: "CCD_ATP", _note_use_boltz_for_custom_smiles: e.seq, count: 1 } };
      });
      const af3Json = JSON.stringify([{ name: jobName, modelSeeds: [], sequences: af3Seqs, dialect: "alphafoldserver", version: 1 }], null, 2);

      const boltzLines = ["version: 1", "sequences:"];
      for (const e of entities) {
        if (e.type === "protein" || e.type === "rna" || e.type === "dna") {
          boltzLines.push(`  - ${e.type}:`, `      id: ${e.chainId}`, `      sequence: ${e.seq}`);
        } else if (e.type === "ion") {
          boltzLines.push(`  - ligand:`, `      id: ${e.chainId}`, `      ccd: ${e.seq.toUpperCase()}`);
        } else {
          const st = analyzeLigandStats(e.seq);
          boltzLines.push(`  - ligand:`, `      id: ${e.chainId}`, st.isCcd ? `      ccd: ${st.code}` : `      smiles: '${e.seq}'`);
        }
      }
      const boltzYaml = boltzLines.join("\\n");

      const chaiLines = [];
      for (const e of entities) {
        const cType = e.type === "ion" ? "ligand" : e.type;
        chaiLines.push(`>${cType}|name=${e.name}`, e.seq);
      }
      const chaiFasta = chaiLines.join("\\n");

      LAST_STUDIO_JOBS = { jobName, af3: af3Json, boltz: boltzYaml, chai: chaiFasta };'''

new_export_build = '''      // Build strictly schema-valid AF3 Server JSON, Boltz-1/2 YAML, Chai-1 FASTA & Boltz-1 Runner
      const jobName = entities.map(e => e.name).join("_vs_").slice(0, 64);
      const ION_TO_SMILES = { MG:"[Mg+2]", ZN:"[Zn+2]", CA:"[Ca+2]", MN:"[Mn+2]", FE:"[Fe+2]", CU:"[Cu+2]", CO:"[Co+2]", NI:"[Ni+2]", NA:"[Na+]", K:"[K+]", CL:"[Cl-]" };
      const af3Seqs = entities.filter(e => !(e.type === "ligand" && !analyzeLigandStats(e.seq).isCcd)).map(e => {
        if (e.type === "protein") return { proteinChain: { sequence: e.seq, count: 1 } };
        if (e.type === "rna") return { rnaSequence: { sequence: e.seq, count: 1 } };
        if (e.type === "dna") return { dnaSequence: { sequence: e.seq, count: 1 } };
        if (e.type === "ion") return { ion: { ion: e.seq.toUpperCase(), count: 1 } };
        return { ligand: { ligand: `CCD_${analyzeLigandStats(e.seq).code}`, count: 1 } };
      });
      const af3Json = hasCustomSmiles
        ? JSON.stringify({ error: "AlphaFold 3 Server web UI does not accept arbitrary custom SMILES ligands. Use the Boltz-1/2 YAML or Chai-1 FASTA tab instead (both natively co-fold custom SMILES + Protein + RNA + Ions)." }, null, 2)
        : JSON.stringify([{ name: jobName, modelSeeds: [], sequences: af3Seqs, dialect: "alphafoldserver", version: 1 }], null, 2);

      const boltzLines = ["version: 1", "sequences:"];
      for (const e of entities) {
        if (e.type === "protein" || e.type === "rna" || e.type === "dna") {
          boltzLines.push(`  - ${e.type}:`, `      id: ${e.chainId}`, `      sequence: ${e.seq}`);
        } else if (e.type === "ion") {
          boltzLines.push(`  - ligand:`, `      id: ${e.chainId}`, `      ccd: ${e.seq.toUpperCase()}`);
        } else {
          const st = analyzeLigandStats(e.seq);
          boltzLines.push(`  - ligand:`, `      id: ${e.chainId}`, st.isCcd ? `      ccd: ${st.code}` : `      smiles: '${e.seq}'`);
        }
      }
      const boltzYaml = boltzLines.join("\\n");

      const chaiLines = [];
      for (const e of entities) {
        if (e.type === "ion") {
          const iSm = ION_TO_SMILES[e.seq.toUpperCase()] || `[${e.seq.toUpperCase()}+2]`;
          chaiLines.push(`>ligand|name=${e.name}`, iSm);
        } else {
          const cType = e.type;
          chaiLines.push(`>${cType}|name=${e.name}`, e.seq);
        }
      }
      const chaiFasta = chaiLines.join("\\n");

      const colabSh = [
        "#!/usr/bin/env bash",
        `# UniPPI-OWL v2.0-MARVAS — 1-Click Boltz-1/2 Co-Folding Script (${jobName})`,
        "# Works on any Linux GPU or Google Colab T4 cell (!bash run_cofold.sh)",
        "pip install -q boltz",
        `cat << 'YAML_EOF' > ${jobName}_boltz.yaml`,
        boltzYaml,
        "YAML_EOF",
        `boltz predict ${jobName}_boltz.yaml --use_msa_server --diffusion_samples 5 --out_dir ./boltz_${jobName}`,
        `echo "Done! Drop ./boltz_${jobName}/predictions/${jobName}_boltz/*.cif + confidence JSONs into UniPPI-OWL Tab 1."`
      ].join("\\n");

      LAST_STUDIO_JOBS = { jobName, af3: af3Json, boltz: boltzYaml, chai: chaiFasta, colab: colabSh, hasCustomSmiles };'''
html = ensure_replaced(html, old_export_build, new_export_build, "JS: Schema-valid AF3 JSON, ionic SMILES for Chai-1, & Boltz-1 Colab runner")

html = ensure_replaced(
    html,
    '''    function exportStudioJob(fmt) {
      if (!LAST_STUDIO_JOBS) runUniversalStudioPreScreen();
      if (!LAST_STUDIO_JOBS) return;
      const extMap = { af3: "_af3_server.json", boltz: "_boltz.yaml", chai: "_chai1.fasta" };
      downloadTextFile(`${LAST_STUDIO_JOBS.jobName}${extMap[fmt]}`, LAST_STUDIO_JOBS[fmt]);
    }''',
    '''    function exportStudioJob(fmt) {
      if (!LAST_STUDIO_JOBS) runUniversalStudioPreScreen();
      if (!LAST_STUDIO_JOBS) return;
      if (fmt === "af3" && LAST_STUDIO_JOBS.hasCustomSmiles) {
        alert("AlphaFold 3 Server web UI cannot accept arbitrary custom SMILES strings. Downloading Boltz-1/2 YAML instead (supports any custom SMILES + Protein + RNA + Ions)!");
        fmt = "boltz";
      }
      const extMap = { af3: "_af3_server.json", boltz: "_boltz.yaml", chai: "_chai1.fasta", colab: "_run_boltz.sh" };
      downloadTextFile(`${LAST_STUDIO_JOBS.jobName}${extMap[fmt]}`, LAST_STUDIO_JOBS[fmt]);
    }''',
    "JS: Guard AF3 custom-SMILES export & support .sh runner download",
)

# B9. Sync Tab 3 evaluateFilesQuiet() with splitChainInput + exact prefix JSON matching + remove dead normKey
old_diff_build = '''        const { ordered, chains } = buildResidues(atoms);
        const cIds = Object.keys(chains);'''
new_diff_build = '''        let { ordered, chains } = buildResidues(atoms);
        const splitSpec = (document.getElementById("splitChainInput")?.value || "").trim();
        if (splitSpec && splitSpec.includes(":")) {
          const [scCidRaw, scRange] = splitSpec.split(":", 2);
          const scCid = scCidRaw.trim();
          if (chains[scCid]) {
            const srcList = chains[scCid], d1 = [], d2 = [];
            if (scRange.includes(",")) {
              const [p1, p2] = scRange.split(",", 2);
              const [lo1, hi1] = p1.split("-").map(x => parseInt(x.trim(), 10));
              const [lo2, hi2] = p2.split("-").map(x => parseInt(x.trim(), 10));
              for (const r of srcList) { if (r.resSeq >= lo1 && r.resSeq <= hi1) d1.push(r); else if (r.resSeq >= lo2 && r.resSeq <= hi2) d2.push(r); }
            } else {
              const cutRes = parseInt(scRange.trim(), 10);
              for (const r of srcList) (r.resSeq <= cutRes ? d1 : d2).push(r);
            }
            if (d1.length && d2.length) {
              const c1 = `${scCid}_D1`, c2 = `${scCid}_D2`;
              for (const r of d1) { r.chainId = c1; r.key = `${c1}|${r.resSeq}|${r.ins}`; for (const a of r.atoms) a.chainId = c1; }
              for (const r of d2) { r.chainId = c2; r.key = `${c2}|${r.resSeq}|${r.ins}`; for (const a of r.atoms) a.chainId = c2; }
              delete chains[scCid]; chains[c1] = d1; chains[c2] = d2;
              const kSet = new Set(Object.values(chains).flat());
              ordered = ordered.filter(r => kSet.has(r));
            }
          }
        }
        const cIds = Object.keys(chains);'''
html = ensure_replaced(html, old_diff_build, new_diff_build, "JS: Wire Intra-Chain Domain Split into Tab 3 Diff Mode")

html = ensure_replaced(
    html,
    '        const activeJsons = jsonParsed.filter(j => (!mIdx || j.lower.includes(`_${mIdx}.json`) || j.lower.includes(`model_${mIdx}`)) && (!cleanStem || j.lower.includes(cleanStem) || structEntries.length === 1));',
    '        const exactJs = jsonParsed.filter(j => (!mIdx || j.lower.includes(`_${mIdx}.json`) || j.lower.includes(`model_${mIdx}`)) && (!cleanStem || j.lower.startsWith(`${cleanStem}_`) || j.lower === `${cleanStem}.json` || structEntries.length === 1));\n'
    '        const activeJsons = exactJs.length ? exactJs : jsonParsed.filter(j => (!mIdx || j.lower.includes(`_${mIdx}.json`) || j.lower.includes(`model_${mIdx}`)) && (!cleanStem || j.lower.includes(cleanStem) || structEntries.length === 1));',
    "JS: Sync exact JSON prefix matching in Tab 3 Diff Mode",
)
html = ensure_replaced(
    html,
    '            const normKey = `${c.resA.replace(/\\d+$/, m=>m)} ↔ ${c.resB}`;\n',
    '',
    "JS: Remove dead unused normKey variable in Tab 3",
)

html_path.write_text(html, encoding="utf-8")

# Verify HTML/JS event handlers & optional node --check
m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
assert m_js, "Missing inline <script> block"
handlers = set(re.findall(r'on(?:click|change|input)="([A-Za-z_]\w*)\(', html))
defined_fns = set(re.findall(r'\bfunction\s+([A-Za-z_]\w*)\s*\(', m_js[-1]))
missing_js = sorted(handlers - defined_fns)
assert not missing_js, f"Undefined JS event handlers: {missing_js}"
print(f"[PASS] index.html upgraded cleanly: all {len(handlers)} event handlers bound.")

node_bin = shutil.which("node") or shutil.which("node.exe")
if node_bin:
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
        tf.write(m_js[-1])
        tf_path = tf.name
    res_node = subprocess.run([node_bin, "--check", tf_path], capture_output=True, text=True)
    Path(tf_path).unlink(missing_ok=True)
    if res_node.returncode != 0:
        raise RuntimeError(f"JavaScript syntax error:\n{res_node.stderr}")
    print(f"[PASS] Verified JS syntax cleanly via `{node_bin} --check`.")

# ==============================================================================
# PART C: FINAL FULL 6-BENCHMARK MASTER REGRESSION SUITE (CLI GROUND TRUTH)
# ==============================================================================
print("\n=== PART C: FINAL 6-BENCHMARK MASTER REGRESSION VERIFICATION ===")
runs = [
    ("1. VIP1-bZIP61 PPI (AF3 5-model .zip)", "benchmarks/fold_vip1_bzip61.zip", "/tmp/v2_final_vip1", None),
    ("2. 1URN U1A-RNA ITC Benchmark", "benchmarks/1URN.cif", "/tmp/v2_final_1urn", None),
    ("3. 5ZC9 MARVAS eIF4A1+RNA+RocA+ANP+Mg", "benchmarks/5ZC9.cif", "/tmp/v2_final_5zc9", None),
    ("4. 4WNP MARVAS ULK1+3RJ Kinase", "benchmarks/4WNP.cif", "/tmp/v2_final_4wnp", None),
    ("5. 7DIY MARVAS Nsp14-Nsp10+Zn+Mg", "benchmarks/7DIY.cif", "/tmp/v2_final_7diy", None),
    ("6. 5ZC9 Intra-Domain Split (RecA1/RecA2)", "benchmarks/5ZC9.cif", "/tmp/v2_final_split", "A:1-236,237-406"),
]
for label, inp, outp, sp in runs:
    cmd = [sys.executable, "unippi_cli.py", "-i", inp, "-o", outp]
    if sp:
        cmd += ["--split-chain", sp]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    rows = json.loads((Path(outp) / "unippi_full_report.json").read_text(encoding="utf-8"))
    for r in rows:
        if r["pattern_id"] == "Pattern A":
            assert r["UIS"] >= 0.50, f"Pattern A / BELOW_THRESHOLD violation in {label}: {r['chain_pair']}"
    top = rows[0]
    print(f"  [PASS] {label:<42} | Top={top['chain_pair']:<18} ({top['modality']:<15}) UIS={top['UIS']:.3f} dG={top['dG_uni_kcal']:>6.2f} {top['pattern_id']}")

# Verify key invariants explicitly
r_vip = json.loads(Path("/tmp/v2_final_vip1/unippi_full_report.json").read_text(encoding="utf-8"))
assert abs(r_vip[0]["UIS"] - 0.5661) < 0.002 and abs(r_vip[0]["dG_uni_kcal"] - (-10.99)) < 0.02
r_5zc9 = json.loads(Path("/tmp/v2_final_5zc9/unippi_full_report.json").read_text(encoding="utf-8"))
assert any(r["modality"] == "TERNARY-GLUE" and r["UIS"] >= 0.50 for r in r_5zc9)
assert any(r["modality"] == "METAL-BRIDGE" and r["UIS"] >= 0.50 and r["pattern_id"] == "Pattern A" for r in r_5zc9)
r_4wnp = json.loads(Path("/tmp/v2_final_4wnp/unippi_full_report.json").read_text(encoding="utf-8"))
assert len([r for r in r_4wnp if r["modality"] == "PROTEIN-LIGAND" and r["UIS"] >= 0.50]) == 4
r_7diy = json.loads(Path("/tmp/v2_final_7diy/unippi_full_report.json").read_text(encoding="utf-8"))
assert len([r for r in r_7diy if r["modality"] == "PROTEIN-METAL" and r["metal_coord_bonds"] == 4 and r["UIS"] >= 0.65]) == 4
print("\n[MASTER AUDIT COMPLETE] ALL 6 BENCHMARKS, AST SCOPE AUDIT & BROWSER HANDLERS PASSED 100%!")
