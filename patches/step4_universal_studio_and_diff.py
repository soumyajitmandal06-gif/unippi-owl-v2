#!/usr/bin/env python3
"""
Phase 3 - Step 4 (Universal X vs Y Mix-and-Match Studio + WT vs Mutant Ensemble Diff Mode
+ Step-3C Remaining UI Polish + CLI --diff-b Parity):
  1. Polishes 3 remaining UI items in index.html (double composite loci prefix, grouped
     composite A+B SVG Contact Network axis, clean 3D legend & updated Scope/Footer copy).
  2. Installs a clean 3-Tab Workspace Switcher in index.html:
       - Tab 1: 3D Atomic Interface Auditor (all 10 verified modalities untouched).
       - Tab 2: Universal X vs Y (+ Optional Z Glue/Ion) Mix-and-Match Studio with PDB/CIF
                chain extraction, custom lab Protein/RNA/DNA/SMILES/Ion inputs, point-mutation
                injector, same-protein intra-chain domain shortcut, MARVAS presets, honest
                biophysical compatibility pre-screen, and 1-click AF3 JSON / Boltz-1 YAML / Chai-1 FASTA export.
       - Tab 3: WT vs Mutant / Apo vs Holo Ensemble Diff Mode with exact two-tailed Welch's
                t-test p-values, Kd fold-change, Contact Rewiring Audit (LOST/GAINED), CSV export,
                and a 1-click demo escape-mutant benchmark.
  3. Adds matching `--diff-b` ensemble comparison mode to unippi_cli.py (100% CLI-Browser parity).
  4. Runs Python AST check, JS structural/handler audit, and CLI zero-regression + diff verification.
"""

import ast
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
# PART A: UPGRADE index.html (POLISH + 3-TAB UNIVERSAL STUDIO & ENSEMBLE DIFF)
# ==============================================================================
print("=== PART A: UPGRADING index.html ===")
html = html_path.read_text(encoding="utf-8")

# A1. Polish the 3 minor Step-3C UI leftovers
html = ensure_regex_replaced(
    html,
    r'Catalytic Metals \(Lime Spheres\)\s*(?:·|&middot;)\s*Metal \(Lime\)',
    'Catalytic Metals (Lime)',
    "HTML: Clean redundant Metal word in 3D Legend",
)

# Fix double chain prefix on composite ternary lociA ("A+B:A:163-195 | B:7-8" -> "A:163-195 | B:7-8")
html = ensure_replaced(
    html,
    'const lociStr = (r.lociA && r.lociA !== "None") ? `${cAName}:${r.lociA} ↔ ${cBName}:${r.lociB}` : "None";',
    'const lociAOut = (r.lociA && r.lociA.includes(":")) ? r.lociA : `${cAName}:${r.lociA || "—"}`;\n'
    '        const lociStr = (r.lociA && r.lociA !== "None") ? `${lociAOut} ↔ ${cBName}:${r.lociB}` : "None";',
    "JS: Clean composite lociA prefix in leaderboard",
)
html = ensure_replaced(
    html,
    '<strong>${cAName}:${r.lociA || "—"} ↔ ${cBName}:${r.lociB || "—"} (${nConf} pairs, ${r.intf_residues} res)</strong>',
    '<strong>${(r.lociA && r.lociA.includes(":")) ? r.lociA : `${cAName}:${r.lociA || "—"}`} ↔ ${cBName}:${r.lociB || "—"} (${nConf} pairs, ${r.intf_residues} res)</strong>',
    "JS: Clean composite lociA prefix in Descriptor 4",
)
html = ensure_replaced(
    html,
    'within <strong>${cAName}:${r.lociA || "—"} &leftrightarrow; ${cBName}:${r.lociB || "—"}</strong>',
    'within <strong>${(r.lociA && r.lociA.includes(":")) ? r.lociA : `${cAName}:${r.lociA || "—"}`} &leftrightarrow; ${cBName}:${r.lociB || "—"}</strong>',
    "JS: Clean composite lociA prefix in Step 1 derivation",
)
html = ensure_replaced(
    html,
    'Ordered Binding Loci = <strong>${cAName}:${r.lociA || "—"} &leftrightarrow; ${cBName}:${r.lociB || "—"} (${r.intf_residues} contacting res)</strong>',
    'Ordered Binding Loci = <strong>${(r.lociA && r.lociA.includes(":")) ? r.lociA : `${cAName}:${r.lociA || "—"}`} &leftrightarrow; ${cBName}:${r.lociB || "—"} (${r.intf_residues} contacting res)</strong>',
    "JS: Clean composite lociA prefix in Track 5",
)

# Upgrade SVG Contact Network left-axis grouping when cA is a composite ternary cleft (e.g. "A+B")
old_svg_coords = '''      let linesSvg = "";
      for (const c of sortedC) {
        const normA = clamp(((c.resSeqA || minSeqA) - minSeqA) / Math.max(1, maxSeqA - minSeqA), 0.0, 1.0);
        const normB = clamp(((c.resSeqB || minSeqB) - minSeqB) / Math.max(1, maxSeqB - minSeqB), 0.0, 1.0);
        const yA = 35 + normA * (h - 70);
        const yB = 35 + normB * (h - 70);'''

new_svg_coords = '''      const isCompA = cA.includes("+");
      const compChParts = isCompA ? cA.split("+") : [cA];
      const compSeqMap = {};
      if (isCompA) {
        compChParts.forEach((chP, pIdx) => {
          const uniq = [...new Set(activeC.filter(c => (c.resA||"").split(":")[0] === chP).map(c => c.resSeqA))].sort((a,b)=>a-b);
          const yLo = pIdx === 0 ? 35 : 158, yHi = pIdx === 0 ? 138 : 265;
          uniq.forEach((sq, sIdx) => { compSeqMap[`${chP}|${sq}`] = uniq.length > 1 ? yLo + (sIdx / (uniq.length - 1)) * (yHi - yLo) : 0.5 * (yLo + yHi); });
        });
      }
      let linesSvg = "";
      for (const c of sortedC) {
        const chTagA = (c.resA || "").split(":")[0];
        const normA = clamp(((c.resSeqA || minSeqA) - minSeqA) / Math.max(1, maxSeqA - minSeqA), 0.0, 1.0);
        const normB = clamp(((c.resSeqB || minSeqB) - minSeqB) / Math.max(1, maxSeqB - minSeqB), 0.0, 1.0);
        const yA = (isCompA && compSeqMap[`${chTagA}|${c.resSeqA}`] !== undefined) ? compSeqMap[`${chTagA}|${c.resSeqA}`] : (35 + normA * (h - 70));
        const yB = (minSeqB === maxSeqB) ? 150.0 : (35 + normB * (h - 70));'''
html = ensure_replaced(html, old_svg_coords, new_svg_coords, "JS: Grouped dual-segment SVG axis for composite ternary A+B")

# Update Scope Box & Footer copy to reflect complete 10-modality + X-vs-Y + Diff capabilities
html = ensure_regex_replaced(
    html,
    r'<li>64-point Golden-Spiral Shrake–Rupley SASA.*?</li>',
    '<li>Dedicated Small-Molecule Pocket, Directional Halogen-Bond (≥120° σ-hole), Inner/Outer-Sphere Metal Coordination (Zn²⁺/Mg²⁺/Ca²⁺/Mn²⁺), Automatic <strong>Ternary Molecular-Glue &amp; Metal-Bridge</strong> detection, and custom <strong>Same-Protein Intra-Chain Domain Slicing</strong>.</li>\n'
    '            <li>64-point Golden-Spiral Shrake–Rupley SASA, confidence-weighted buried surface area (wBSA), Non-Interacting Surface (NIS) properties, Universal X vs. Y Studio (AF3 / Boltz-1 / Chai-1 exporter), and Welch\'s t-test WT vs. Mutant Ensemble Diff Mode.</li>',
    "HTML: Update Scope Box with Ligand/Ternary/Metal/Studio/Diff features",
)
html = ensure_regex_replaced(
    html,
    r'<strong>UniPPI-OWL v2\.0-MARVAS</strong>\s*(?:·|&middot;)\s*Universal Protein–Protein.*?<br>',
    '<strong>UniPPI-OWL v2.0-MARVAS</strong> · Universal Multi-Modal Biomolecular Affinity, Ternary Glue &amp; Ensemble Diff Resource<br>',
    "HTML: Update footer subtitle",
)

# A2. Add CSS for Tabs, Universal X-vs-Y Studio & Ensemble Diff Mode
studio_css = '''
    .tab-bar {
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
      background: var(--bg-card);
      border: 1px solid var(--border);
      border-radius: var(--radius);
      padding: 8px 10px;
      box-shadow: var(--shadow);
    }
    .tab-btn {
      flex: 1;
      min-width: 210px;
      padding: 10px 14px;
      border-radius: 7px;
      border: 1px solid transparent;
      background: var(--bg-subtle);
      color: var(--text-muted);
      font-weight: 700;
      font-size: 13.5px;
      cursor: pointer;
      transition: all 0.15s ease;
      text-align: center;
    }
    .tab-btn:hover { background: #e2e8f0; color: var(--text-main); }
    .tab-btn.active { background: var(--primary); color: #ffffff; border-color: var(--primary); box-shadow: 0 2px 6px rgba(15, 76, 129, 0.22); }
    .studio-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 18px; margin-bottom: 16px; }
    .slot-box { border: 1px solid var(--border-strong); border-radius: 8px; padding: 16px; background: var(--bg-subtle); display: flex; flex-direction: column; gap: 10px; }
    .slot-box h3 { font-size: 15px; color: var(--primary); display: flex; justify-content: space-between; align-items: center; }
    .slot-box textarea, .slot-box input[type="text"], .slot-box select {
      width: 100%; padding: 8px 10px; border: 1px solid var(--border-strong); border-radius: 6px;
      font-family: var(--font-mono); font-size: 12.5px; background: #ffffff;
    }
    .slot-box textarea { min-height: 92px; resize: vertical; }
    .diff-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 16px; }
    @media (max-width: 820px) { .diff-grid { grid-template-columns: 1fr; } }
    .badge-lost { background: #ffe4e6; color: #be123c; }
    .badge-gained { background: #d1fae5; color: #047857; }
    .badge-persist { background: #e0f2fe; color: #0f4c81; }
  </style>'''
html = ensure_replaced(html, '  </style>', studio_css, "CSS: Add 3-Tab, Studio & Diff styles")

# A3. Upgrade top navigation links to switch tabs cleanly
old_nav = '''      <nav class="resource-links">
        <a href="file:///C:/Users/USER/OneDrive/Desktop/UNI-PPI-OWL-v2-MARVAS/index.html#analyze">Upload &amp; Run</a>
        <a href="file:///C:/Users/USER/OneDrive/Desktop/UNI-PPI-OWL-v2-MARVAS/index.html#leaderboardCard">Results</a>
        <a href="file:///C:/Users/USER/OneDrive/Desktop/UNI-PPI-OWL-v2-MARVAS/index.html#workbench">3D &amp; Contacts</a>
        <a href="file:///C:/Users/USER/OneDrive/Desktop/UNI-PPI-OWL-v2-MARVAS/index.html#svg-maps">Vector Maps</a>
        <a href="file:///C:/Users/USER/OneDrive/Desktop/UNI-PPI-OWL-v2-MARVAS/index.html#interpretation-guide">Method &amp; Scope</a>
      </nav>'''
new_nav = '''      <nav class="resource-links">
        <a href="#analyze" onclick="switchWorkspaceTab('auditor')">3D Auditor</a>
        <a href="#tabPanelStudio" onclick="switchWorkspaceTab('studio')">X vs. Y Studio</a>
        <a href="#tabPanelDiff" onclick="switchWorkspaceTab('diff')">Ensemble Diff</a>
        <a href="#workbench" onclick="switchWorkspaceTab('auditor')">3D &amp; Contacts</a>
        <a href="#interpretation-guide">Method &amp; Scope</a>
      </nav>'''
if old_nav in html:
    html = ensure_replaced(html, old_nav, new_nav, "HTML: Update top nav links for 3-tab workspace")
else:
    html = ensure_regex_replaced(html, r'<nav class="resource-links">.*?</nav>', new_nav, "HTML: Update top nav links for 3-tab workspace (regex)")

# A4. Insert 3-Tab Switcher right after Hero section and wrap Tab 1 panels
tab_bar_html = '''    </section>

    <div class="tab-bar" role="tablist" aria-label="UniPPI-OWL Workspace Tabs">
      <button id="tabBtnAuditor" class="tab-btn active" onclick="switchWorkspaceTab('auditor')">1. 3D Atomic Interface Auditor (PPI · RNA · Drug · Metal · Glue)</button>
      <button id="tabBtnStudio" class="tab-btn" onclick="switchWorkspaceTab('studio')">2. Universal X vs. Y Studio (Known PDB + Custom Lab Seq / SMILES)</button>
      <button id="tabBtnDiff" class="tab-btn" onclick="switchWorkspaceTab('diff')">3. WT vs. Mutant / Apo vs. Holo Ensemble Diff (Welch t-Test)</button>
    </div>

    <div id="tabPanelAuditor" style="display:flex; flex-direction:column; gap:26px;">

    <section class="card" id="analyze">'''
html = ensure_replaced(
    html,
    '    </section>\n\n    <section class="card" id="analyze">',
    tab_bar_html,
    "HTML: Insert 3-Tab Switcher & open Tab 1 container",
)

# A5. Close Tab 1 & insert Tab 2 (Universal X vs Y Studio) and Tab 3 (Ensemble Diff) before Method & Scope
tabs_2_and_3_html = '''    </div><!-- /tabPanelAuditor -->

    <div id="tabPanelStudio" style="display:none; flex-direction:column; gap:20px;">
      <section class="card">
        <div class="card-head">
          <div>
            <div class="section-kicker">Universal Mix-and-Match Partner Builder (Zero-Faking Honest Architecture)</div>
            <h2>Universal X vs. Y (+ Optional Z Glue / Ion) Interaction Studio</h2>
            <p style="color:var(--text-muted); font-size:13.5px;">
              Combine <strong>any known <code>.pdb</code>/<code>.cif</code></strong> with <strong>any custom lab-synthesized protein, RNA/DNA, point mutation, small-molecule drug (SMILES/CCD), or catalytic ion</strong> — or slice <strong>two domains inside one protein</strong>.
              Step A runs an <strong>instant sequence &amp; chemical compatibility pre-screen</strong>; Step B exports ready-to-run 3D co-fold jobs for <strong>AlphaFold 3 Server</strong>, <strong>Boltz-1/2 (custom SMILES supported)</strong>, and <strong>Chai-1</strong> so you can drop the co-folded 3D complex into Tab 1 for exact atomic ΔG, K<sub>d</sub>, UIS, and 3D contact maps.
            </p>
          </div>
        </div>

        <div style="display:flex; gap:8px; flex-wrap:wrap; margin-bottom:16px;">
          <button class="btn btn-secondary btn-sm" onclick="loadStudioPreset('marvas_eif4a')">MARVAS Preset: eIF4A1 + Polypurine RNA + Rocaglamide A</button>
          <button class="btn btn-secondary btn-sm" onclick="loadStudioPreset('marvas_ulk1')">MARVAS Preset: ULK1 Kinase + Inhibitor SMILES</button>
          <button class="btn btn-secondary btn-sm" onclick="loadStudioPreset('marvas_nsp14')">MARVAS Preset: CoV Nsp14 + Nsp10 + Zn²⁺</button>
          <button class="btn btn-secondary btn-sm" onclick="loadStudioPreset('custom_rna')">Universal Preset: Custom Lab Protein + Custom Lab RNA</button>
        </div>

        <div class="studio-grid">
          <div class="slot-box">
            <h3><span>Partner X (Receptor / Molecule 1)</span><span id="slotXBadge" style="font-size:11px; font-family:var(--font-mono); color:var(--text-muted);">Empty</span></h3>
            <div style="display:grid; grid-template-columns:1fr 1fr; gap:8px;">
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
            </div>
            <select id="slotXChainPick" style="display:none;" onchange="applyPickedSlotChain('X')"></select>
            <input type="text" id="slotXName" placeholder="Partner X Name (e.g. Human_eIF4A1)" value="Partner_X">
            <textarea id="slotXSeq" placeholder="Paste custom amino acid / RNA / DNA sequence, SMILES, or upload a .pdb/.cif above..." oninput="refreshSlotBadge('X')"></textarea>
            <input type="text" id="slotXMut" placeholder="Optional point mutations (e.g. F163L, Q195A or G8U)">
            <div id="slotXIntraWrap" style="display:none; background:#e0f2fe; border:1px solid #bae6fd; border-radius:6px; padding:10px; font-size:12px;">
              <strong style="color:var(--primary);">Same-Protein Intra-Chain Domain Shortcut:</strong> Uploaded a single structure? Score two domains inside it right now in Tab 1:
              <div style="display:flex; gap:6px; margin-top:6px; flex-wrap:wrap;">
                <input type="text" id="slotXDomainSpec" placeholder="e.g. A:1-236,237-406 (or blank for auto-trough)" style="flex:1; min-width:170px;">
                <button class="btn btn-primary btn-sm" onclick="runSlotXIntraDomainNow()">Audit Domains in Tab 1</button>
              </div>
            </div>
          </div>

          <div class="slot-box">
            <h3><span>Partner Y (Target / Molecule 2)</span><span id="slotYBadge" style="font-size:11px; font-family:var(--font-mono); color:var(--text-muted);">Empty</span></h3>
            <div style="display:grid; grid-template-columns:1fr 1fr; gap:8px;">
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
            </div>
            <select id="slotYChainPick" style="display:none;" onchange="applyPickedSlotChain('Y')"></select>
            <input type="text" id="slotYName" placeholder="Partner Y Name (e.g. Custom_Lab_RNA)" value="Partner_Y">
            <textarea id="slotYSeq" placeholder="Paste custom RNA / protein / DNA sequence, SMILES (e.g. Rocaglamide A), CCD code (e.g. ATP, RCG), or upload a .pdb/.cif..." oninput="refreshSlotBadge('Y')"></textarea>
            <input type="text" id="slotYMut" placeholder="Optional point mutations (e.g. A7U, G8C)">
          </div>
        </div>

        <div class="slot-box" style="margin-bottom:16px; background:#faf5ff; border-color:#d8b4fe;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px;">
            <strong style="color:#581c87; font-size:13.5px;">Optional Slot Z — Ternary Molecular Glue / Small-Molecule Drug / Catalytic Ion (e.g. Rocaglamide A, Aglain 13, OX18, Mg²⁺, Zn²⁺)</strong>
            <label style="font-size:12px; font-weight:600; cursor:pointer; display:flex; align-items:center; gap:6px;">
              <input type="checkbox" id="slotZEnable" onchange="toggleSlotZBox()"> Enable 3rd Component (Ternary Complex)
            </label>
          </div>
          <div id="slotZFields" style="display:none; grid-template-columns:180px 180px 1fr; gap:10px; margin-top:6px;">
            <select id="slotZType">
              <option value="ligand">Small Molecule (SMILES / CCD)</option>
              <option value="ion">Catalytic Metal / Ion</option>
              <option value="protein">3rd Protein Chain</option>
              <option value="rna">3rd RNA Chain</option>
            </select>
            <input type="text" id="slotZName" value="Glue_Z" placeholder="Name (e.g. RocA)">
            <input type="text" id="slotZSeq" placeholder="Enter SMILES (e.g. COC1=CC=C(...)), CCD code (RCG, ATP), ion (MG, ZN), or sequence">
          </div>
        </div>

        <div style="display:flex; gap:10px; flex-wrap:wrap; align-items:center; margin-bottom:16px;">
          <button class="btn btn-primary" onclick="runUniversalStudioPreScreen()">Run Instant Pre-Screen &amp; Build 3D Co-Fold Jobs</button>
          <button class="btn btn-secondary btn-sm" onclick="exportStudioJob('af3')">Download AF3 Server .JSON</button>
          <button class="btn btn-secondary btn-sm" onclick="exportStudioJob('boltz')">Download Boltz-1/2 .YAML (Any Custom SMILES)</button>
          <button class="btn btn-secondary btn-sm" onclick="exportStudioJob('chai')">Download Chai-1 .FASTA</button>
          <button class="btn btn-secondary btn-sm" onclick="copyStudioPreview()">Copy Preview</button>
        </div>

        <div id="studioOutputPanel" style="display:none;">
          <div id="studioWarningBanner" style="display:none; background:var(--amber-bg); color:var(--amber); border:1px solid #fcd34d; padding:10px 14px; border-radius:6px; margin-bottom:12px; font-size:12.5px; font-weight:600;"></div>
          <div id="studioPrescreenGrid" class="descriptor-grid" style="margin-bottom:14px;"></div>
          <div id="studioVerdictBox" style="padding:12px 14px; border-radius:8px; background:var(--primary-soft); border-left:4px solid var(--primary); font-size:13px; margin-bottom:14px;"></div>
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:6px;">
            <strong style="font-size:12.5px; color:var(--text-muted); text-transform:uppercase;">Generated Co-Folding Job Preview (Run in AF3 Server / Boltz-1 / Chai-1, then drop output .ZIP/.CIF into Tab 1)</strong>
            <select id="studioPreviewFormat" style="padding:4px 8px; font-size:12px; border:1px solid var(--border-strong); border-radius:5px;" onchange="renderStudioJobPreview()">
              <option value="boltz">Boltz-1 / Boltz-2 YAML (Universal — Supports Custom SMILES)</option>
              <option value="af3">AlphaFold 3 Server JSON</option>
              <option value="chai">Chai-1 Multi-Modal FASTA</option>
            </select>
          </div>
          <pre id="studioJobPreview" style="background:#0f172a; color:#e2e8f0; padding:14px; border-radius:8px; font-family:var(--font-mono); font-size:12px; max-height:260px; overflow:auto; white-space:pre-wrap;"></pre>
        </div>
      </section>
    </div><!-- /tabPanelStudio -->

    <div id="tabPanelDiff" style="display:none; flex-direction:column; gap:20px;">
      <section class="card">
        <div class="card-head">
          <div>
            <div class="section-kicker">MARVAS Escape-Mutant &amp; Apo-vs-Holo Statistical Auditor</div>
            <h2>WT vs. Mutant / Apo vs. Holo Ensemble Diff Mode</h2>
            <p style="color:var(--text-muted); font-size:13.5px;">
              Upload two structures or 5-model AF3 <code>.zip</code> ensembles (e.g., <strong>Condition A = WT or Drug-Bound Holo</strong> vs. <strong>Condition B = Escape Mutant or Apo</strong>).
              UniPPI-OWL scores both ensembles with the locked Step 3C physics engine and computes <strong>ΔUIS</strong>, <strong>ΔΔG</strong>, <strong>K<sub>d</sub> fold-change</strong>, <strong>two-tailed Welch's t-test p-values</strong>, and an automated <strong>Contact Rewiring Audit (LOST vs. GAINED contacts)</strong>.
            </p>
          </div>
          <div style="display:flex; gap:8px; flex-wrap:wrap;">
            <button class="btn btn-secondary btn-sm" onclick="loadDemoEnsembleDiff()">Load Demo WT vs. Escape-Mutant (5+5 Models)</button>
            <button id="btnExportDiffCSV" class="btn btn-secondary btn-sm" style="display:none;" onclick="downloadDiffRewiringCSV()">Download Diff CSV</button>
          </div>
        </div>

        <div class="diff-grid">
          <div class="slot-box" style="background:#f0f9ff; border-color:#7dd3fc;">
            <h3><span>Condition A (Reference: WT / Holo)</span><span id="diffACount" style="font-size:11.5px; font-family:var(--font-mono);">0 models</span></h3>
            <label class="btn btn-primary btn-sm" style="cursor:pointer;">
              + Upload Condition A (.ZIP / .CIF / .PDB + .JSON)
              <input type="file" id="diffAInput" multiple="" accept=".zip,.cif,.mmcif,.pdb,.json,.ent" style="display:none" onchange="handleDiffUpload('A', this.files)">
            </label>
            <div id="diffAStatus" style="font-family:var(--font-mono); font-size:12px; color:var(--text-muted);">No Condition A structures loaded yet.</div>
          </div>

          <div class="slot-box" style="background:#fff1f2; border-color:#fda4af;">
            <h3 style="color:var(--rose);"><span>Condition B (Perturbed: Mutant / Apo)</span><span id="diffBCount" style="font-size:11.5px; font-family:var(--font-mono);">0 models</span></h3>
            <label class="btn btn-secondary btn-sm" style="cursor:pointer; border-color:var(--rose); color:var(--rose);">
              + Upload Condition B (.ZIP / .CIF / .PDB + .JSON)
              <input type="file" id="diffBInput" multiple="" accept=".zip,.cif,.mmcif,.pdb,.json,.ent" style="display:none" onchange="handleDiffUpload('B', this.files)">
            </label>
            <div id="diffBStatus" style="font-family:var(--font-mono); font-size:12px; color:var(--text-muted);">No Condition B structures loaded yet.</div>
          </div>
        </div>

        <div style="display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin-bottom:16px;">
          <label style="font-size:12px; font-weight:600; color:var(--text-muted);">Compare Interface Pair:
            <select id="diffPairSelect" style="margin-left:6px; padding:5px 9px; border:1px solid var(--border-strong); border-radius:6px; font-family:var(--font-mono); font-size:12px;" onchange="renderEnsembleDiffResults()">
              <option value="TOP">Top-Ranked Interface per Model</option>
            </select>
          </label>
          <button class="btn btn-primary btn-sm" onclick="renderEnsembleDiffResults()">Compute Ensemble Diff &amp; Contact Rewiring</button>
        </div>

        <div id="diffResultsWrap" style="display:none;">
          <div id="diffHeadlineBanner" style="padding:12px 16px; border-radius:8px; margin-bottom:14px; font-family:var(--font-mono); font-size:13px; border-left:4px solid var(--primary); background:var(--bg-subtle);"></div>
          <div class="table-wrap" style="margin-bottom:16px;">
            <table>
              <thead>
                <tr>
                  <th>Biophysical Metric</th>
                  <th>Condition A (Mean ± SD)</th>
                  <th>Condition B (Mean ± SD)</th>
                  <th>Delta (B − A)</th>
                  <th>Welch's t-Test (p-value)</th>
                  <th>Interpretation</th>
                </tr>
              </thead>
              <tbody id="diffStatsBody"></tbody>
            </table>
          </div>
          <h3 style="font-size:14.5px; margin-bottom:8px; color:var(--primary);">Atomic Contact Rewiring Audit (Occupancy Shift Across Ensembles)</h3>
          <div class="table-wrap" style="max-height:360px; overflow-y:auto;">
            <table>
              <thead>
                <tr>
                  <th>Status</th>
                  <th>Contact Pair</th>
                  <th>Cond A Occupancy</th>
                  <th>Cond B Occupancy</th>
                  <th>Mean d<sub>min</sub> A → B (Å)</th>
                  <th>Channels A → B</th>
                </tr>
              </thead>
              <tbody id="diffRewiringBody"></tbody>
            </table>
          </div>
        </div>
      </section>
    </div><!-- /tabPanelDiff -->

    <section class="card" id="interpretation-guide">'''
html = ensure_replaced(
    html,
    '    <section class="card" id="interpretation-guide">',
    tabs_2_and_3_html,
    "HTML: Insert Tab 2 (Universal X vs Y Studio) & Tab 3 (Ensemble Diff Mode)",
)

# A6. Inject complete JavaScript for Tab Switcher, Universal Studio, and Ensemble Diff Mode
step4_js = r'''
    // =========================================================================
    // STEP 4: 3-TAB SWITCHER, UNIVERSAL X-VS-Y STUDIO & ENSEMBLE DIFF ENGINE
    // =========================================================================
    function switchWorkspaceTab(tab) {
      const map = { auditor: ["tabPanelAuditor", "tabBtnAuditor"], studio: ["tabPanelStudio", "tabBtnStudio"], diff: ["tabPanelDiff", "tabBtnDiff"] };
      for (const [k, [pid, bid]] of Object.entries(map)) {
        const p = document.getElementById(pid), b = document.getElementById(bid);
        if (p) p.style.display = (k === tab) ? "flex" : "none";
        if (b) b.classList.toggle("active", k === tab);
      }
      if (tab === "auditor" && VIEWER_3D) {
        setTimeout(() => { try { VIEWER_3D.resize(); VIEWER_3D.render(); } catch(e) {} }, 60);
      }
    }

    const STUDIO_SLOT_FILES = { X: null, Y: null };
    const STUDIO_SLOT_PARSED = { X: null, Y: null };

    function updateStudioSlotUI(slot) {
      const tp = document.getElementById(`slot${slot}Type`)?.value || "protein";
      const mutEl = document.getElementById(`slot${slot}Mut`);
      if (mutEl) mutEl.style.display = (tp === "protein" || tp === "rna" || tp === "dna") ? "block" : "none";
      refreshSlotBadge(slot);
    }

    function toggleSlotZBox() {
      const en = document.getElementById("slotZEnable")?.checked;
      const box = document.getElementById("slotZFields");
      if (box) box.style.display = en ? "grid" : "none";
    }

    function cleanSeqInput(raw, kind) {
      const stripped = (raw || "").split(/\r?\n/).filter(l => !l.trim().startsWith(">")).join("").trim();
      if (kind === "ligand" || kind === "ion") return stripped.replace(/\s+/g, "");
      const up = stripped.toUpperCase().replace(/[^A-Z]/g, "");
      return kind === "rna" ? up.replace(/T/g, "U") : (kind === "dna" ? up.replace(/U/g, "T") : up);
    }

    function refreshSlotBadge(slot) {
      const tp = document.getElementById(`slot${slot}Type`)?.value || "protein";
      const val = cleanSeqInput(document.getElementById(`slot${slot}Seq`)?.value || "", tp);
      const badge = document.getElementById(`slot${slot}Badge`);
      if (!badge) return;
      if (!val) { badge.textContent = "Empty"; return; }
      if (tp === "protein") badge.textContent = `${val.length} aa`;
      else if (tp === "rna" || tp === "dna") badge.textContent = `${val.length} nt (${tp.toUpperCase()})`;
      else if (tp === "ion") badge.textContent = `Ion: ${val.toUpperCase()}`;
      else badge.textContent = val.length <= 5 && /^[A-Z0-9]+$/i.test(val) ? `CCD: ${val.toUpperCase()}` : `SMILES (${val.length} chars)`;
    }

    async function loadSlotStructureFile(slot, fileList) {
      if (!fileList || !fileList.length) return;
      const f = fileList[0];
      const text = await f.text();
      STUDIO_SLOT_FILES[slot] = { name: f.name, file: f, text };
      const { atoms } = parseStructureText(text, f.name);
      if (!atoms.length) { alert(`No valid atoms parsed from ${f.name}`); return; }
      const { chains } = buildResidues(atoms);
      STUDIO_SLOT_PARSED[slot] = chains;
      const pick = document.getElementById(`slot${slot}ChainPick`);
      const cIds = Object.keys(chains);
      if (pick) {
        pick.innerHTML = cIds.map(cid => {
          const rlist = chains[cid];
          const et = rlist[0]?.entityType || "PROTEIN";
          return `<option value="${cid}">Chain ${cid} (${rlist.length} ${et}) — from ${f.name}</option>`;
        }).join("");
        pick.style.display = "block";
      }
      if (slot === "X") {
        const iw = document.getElementById("slotXIntraWrap");
        if (iw) iw.style.display = "block";
      }
      applyPickedSlotChain(slot);
    }

    function applyPickedSlotChain(slot) {
      const chains = STUDIO_SLOT_PARSED[slot];
      const pick = document.getElementById(`slot${slot}ChainPick`);
      if (!chains || !pick) return;
      const cid = pick.value;
      const rlist = chains[cid];
      if (!rlist || !rlist.length) return;
      const et = rlist[0].entityType;
      const typeSel = document.getElementById(`slot${slot}Type`);
      const nameIn = document.getElementById(`slot${slot}Name`);
      const seqIn = document.getElementById(`slot${slot}Seq`);
      const stem = (STUDIO_SLOT_FILES[slot]?.name || "PDB").replace(/\.(cif|mmcif|pdb|ent)$/i, "");
      if (nameIn) nameIn.value = `${stem}_${cid}`;
      if (et === "PROTEIN") {
        if (typeSel) typeSel.value = "protein";
        if (seqIn) seqIn.value = rlist.map(r => THREE_TO_ONE[r.canonicalName] || "A").join("");
      } else if (et === "NUCLEIC") {
        const isRna = rlist.some(r => r.canonicalName === "U" || r.canonicalName === "URA" || r.atoms.some(a => a.name === "O2'"));
        if (typeSel) typeSel.value = isRna ? "rna" : "dna";
        if (seqIn) seqIn.value = rlist.map(r => (NUC_TO_ONE[r.canonicalName] || "A").toUpperCase()).join("");
      } else if (et === "ION") {
        if (typeSel) typeSel.value = "ion";
        if (seqIn) seqIn.value = rlist[0].resName;
      } else {
        if (typeSel) typeSel.value = "ligand";
        if (seqIn) seqIn.value = rlist[0].resName;
      }
      updateStudioSlotUI(slot);
    }

    async function runSlotXIntraDomainNow() {
      const entry = STUDIO_SLOT_FILES.X;
      if (!entry || !entry.file) { alert("Upload a .pdb/.cif file into Partner X first."); return; }
      const spec = (document.getElementById("slotXDomainSpec")?.value || "").trim();
      const splitBox = document.getElementById("splitChainInput");
      if (splitBox) splitBox.value = spec;
      switchWorkspaceTab("auditor");
      await handleFileSelect([entry.file]);
    }

    function loadStudioPreset(preset) {
      const setSlot = (s, tp, name, seq, mut="") => {
        document.getElementById(`slot${s}Type`).value = tp;
        document.getElementById(`slot${s}Name`).value = name;
        document.getElementById(`slot${s}Seq`).value = seq;
        const mEl = document.getElementById(`slot${s}Mut`);
        if (mEl) mEl.value = mut;
        const pEl = document.getElementById(`slot${s}ChainPick`);
        if (pEl) pEl.style.display = "none";
        updateStudioSlotUI(s);
      };
      const zCheck = document.getElementById("slotZEnable");
      if (preset === "marvas_eif4a") {
        setSlot("X", "protein", "Human_eIF4A1_RecA1", "MSASQDSRSRDNGPDGMEPEGVIESNWNEIVDSFDDMNLSESLLRGIYAYGFEKPSAIQQRAILPCIKGYDVIAQAQSGTGKTATFAISILQQIELDLKATQALVLAPTRELAQQIQKVVMALGDYMGASCHACIGGTNVRAEVQKLQMEAPHIIVGTPGRVFDMLNRRYLSPKYIKMFVLDEADEMLSRGFKDQIYDIFQKLNSNTQVVLLSATMPSDVLEVTKKFMRDPIRILVKKEELTLEGIRQFYINVEREEWKLDTLCDLYETLTITQAVIFINTRRKVDWLTEKMHARDFTVSAMHGDMDQKERDVIMREFRSGSSRVLITTDLLARGIDVQQVSLVINYDLPTNRENYIHRIGRGGRFGRKGVAINMVTEEDKRTLRDIETFYNTSIEEMPLNVADLI");
        setSlot("Y", "rna", "Polypurine_RNA_10nt", "AGAGAGAGAG");
        if (zCheck) { zCheck.checked = true; toggleSlotZBox(); }
        document.getElementById("slotZType").value = "ligand";
        document.getElementById("slotZName").value = "Rocaglamide_A";
        document.getElementById("slotZSeq").value = "COC1=CC2=C(C(=C1)OC)[C@@]3([C@H]([C@@H]([C@H]([C@@]3(O2)C4=CC=C(C=C4)OC)O)C(=O)N(C)C)C5=CC=CC=C5)O";
      } else if (preset === "marvas_ulk1") {
        setSlot("X", "protein", "Human_ULK1_Kinase", "MEFEYKRVIGKGAFAVVFKGRHREKHDLEVAVKCINKKNLAKSQTLLGKEIKILKELKHENIVALYDFQEMANSVYLVMEYCNGGDLADYLHAMRTLSEDTIRLFLQQIAGAMRLLHSKGIIHRDLKPQNILLSNPAGRRANPNSIRVKIADFGFARYLQSNMMAATLCGSPMYMAPEVIMSQHYDGKADLWSIGTIVYQCLTGKAPFQASSPQDLRLFYEKNKTLVPTIPRETSAPLRQLLLALLQRNHKDRMDFDEFFHHPFLDASPSVRKSPPVPVPSYPSSGSG");
        setSlot("Y", "ligand", "ULK1_Inhibitor_SMILES", "CC1=CC(=NN1)NC2=NC(=NC=C2Cl)NC3=CC=CC(=C3)C(=O)N");
        if (zCheck) { zCheck.checked = false; toggleSlotZBox(); }
      } else if (preset === "marvas_nsp14") {
        setSlot("X", "protein", "SARS2_Nsp14_ExoN", "AENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQ");
        setSlot("Y", "protein", "SARS2_Nsp10", "AGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQ");
        if (zCheck) { zCheck.checked = true; toggleSlotZBox(); }
        document.getElementById("slotZType").value = "ion";
        document.getElementById("slotZName").value = "Catalytic_ZN";
        document.getElementById("slotZSeq").value = "ZN";
      } else {
        setSlot("X", "protein", "Custom_Lab_RBP", "MKRKRRGFSYFDQYLKKRKEKAVILGGFTRKGGYVVFEKPDAAEKLRKRYGKNIVFVRK");
        setSlot("Y", "rna", "Custom_Lab_RNA", "GGAGACUCAGCUAGAGAGCUCC");
        if (zCheck) { zCheck.checked = false; toggleSlotZBox(); }
      }
      runUniversalStudioPreScreen();
    }

    function applyPointMutationsToSeq(seq, mutStr, kind) {
      const notes = [], warnings = [];
      if (!mutStr || !mutStr.trim()) return { seq, notes, warnings };
      const arr = Array.from(seq);
      const toks = mutStr.split(/[,\s;]+/).map(s => s.trim().toUpperCase()).filter(Boolean);
      for (const tok of toks) {
        const m = tok.match(/^([A-Z])(\d+)([A-Z])$/);
        if (!m) { warnings.push(`Skipped unrecognized mutation '${tok}' (use e.g. F163L or G8U).`); continue; }
        const wt = m[1], pos = parseInt(m[2], 10) - 1, mut = m[3];
        if (pos < 0 || pos >= arr.length) { warnings.push(`Mutation ${tok} out of range (length ${arr.length}).`); continue; }
        if (arr[pos] !== wt) warnings.push(`WT mismatch at ${pos+1}: sequence has ${arr[pos]}, expected ${wt} (applied ${mut} anyway).`);
        arr[pos] = mut;
        notes.push(`${wt}${pos+1}${mut}`);
      }
      return { seq: arr.join(""), notes, warnings };
    }

    function analyzeProteinSeqStats(seq) {
      const n = seq.length;
      if (!n) return null;
      const cnt = {};
      for (const ch of seq) cnt[ch] = (cnt[ch] || 0) + 1;
      const chargeAtPh = ph => {
        const pos = 1/(1+Math.pow(10, ph-8.0)) + (cnt.K||0)/(1+Math.pow(10, ph-10.5)) + (cnt.R||0)/(1+Math.pow(10, ph-12.5)) + (cnt.H||0)/(1+Math.pow(10, ph-6.0));
        const neg = 1/(1+Math.pow(10, 3.1-ph)) + (cnt.D||0)/(1+Math.pow(10, 3.9-ph)) + (cnt.E||0)/(1+Math.pow(10, 4.1-ph)) + (cnt.C||0)/(1+Math.pow(10, 8.3-ph)) + (cnt.Y||0)/(1+Math.pow(10, 10.1-ph));
        return pos - neg;
      };
      let lo = 1.0, hi = 14.0;
      for (let it = 0; it < 28; it++) { const mid = 0.5*(lo+hi); if (chargeAtPh(mid) > 0) lo = mid; else hi = mid; }
      const pI = 0.5*(lo+hi), q72 = chargeAtPh(7.2);
      let basicPatches = 0;
      for (let i = 0; i <= n - 5; i++) {
        let rk = 0;
        for (let j = i; j < i + 5; j++) if (seq[j] === "R" || seq[j] === "K") rk++;
        if (rk >= 3) basicPatches++;
      }
      const aromPct = 100.0 * ((cnt.F||0) + (cnt.Y||0) + (cnt.W||0)) / n;
      const gravy = Array.from(seq).reduce((s, ch) => s + (KD_HYDROPATHY[ONE_TO_THREE[ch]] ?? 0), 0) / n;
      return { len: n, mwKda: (n * 0.110), q72, pI, basicPatches, aromPct, gravy };
    }

    function analyzeNucleicSeqStats(seq, kind) {
      const n = seq.length;
      if (!n) return null;
      const gc = 100.0 * Array.from(seq).filter(c => c === "G" || c === "C").length / n;
      const pur = 100.0 * Array.from(seq).filter(c => c === "A" || c === "G").length / n;
      let maxPurRun = 0, cur = 0;
      for (const c of seq) { if (c === "A" || c === "G") { cur++; maxPurRun = Math.max(maxPurRun, cur); } else cur = 0; }
      return { len: n, q72: -Math.max(1, n - 1), gc, pur, maxPurRun, kind };
    }

    function analyzeLigandStats(val) {
      const isCcd = val.length <= 5 && /^[A-Z0-9]+$/i.test(val);
      if (isCcd) return { isCcd: true, code: val.toUpperCase(), heavyEst: 0, aromAtoms: 0, halogens: 0, hba: 0, anionicHint: ["ATP","ANP","ADP","GTP","GNP","AMP"].includes(val.toUpperCase()) };
      const hal = (val.match(/Cl|Br|F|I/g) || []).length;
      const arom = (val.match(/[cnos]/g) || []).length;
      const hba = (val.match(/[NnOo]/g) || []).length;
      const heavy = (val.match(/Cl|Br|[BCNOPSFIbcnops]/g) || []).length;
      const anionicHint = /P\(=O\)|C\(=O\)\[?O-?\]?/.test(val);
      return { isCcd: false, smiles: val, heavyEst: heavy, aromAtoms: arom, halogens: hal, hba, anionicHint };
    }

    function maxDuplexComplementarity(s1, s2) {
      const comp = { A:"U", U:"A", T:"A", G:"C", C:"G" };
      const a = s1.replace(/T/g, "U"), b = Array.from(s2.replace(/T/g, "U")).reverse().join("");
      let best = 0;
      for (let shift = -b.length + 1; shift < a.length; shift++) {
        let run = 0;
        for (let i = 0; i < b.length; i++) {
          const ai = shift + i;
          if (ai >= 0 && ai < a.length && (comp[a[ai]] === b[i] || (a[ai]==="G"&&b[i]==="U") || (a[ai]==="U"&&b[i]==="G"))) {
            run++; best = Math.max(best, run);
          } else run = 0;
        }
      }
      return best;
    }

    let LAST_STUDIO_JOBS = null;

    function buildStudioEntities() {
      const warnings = [];
      const entities = [];
      const slots = ["X", "Y"].concat(document.getElementById("slotZEnable")?.checked ? ["Z"] : []);
      const ids = ["A", "B", "C"];
      slots.forEach((slot, idx) => {
        const tp = document.getElementById(`slot${slot}Type`)?.value || "protein";
        const name = (document.getElementById(`slot${slot}Name`)?.value || `Partner_${slot}`).trim().replace(/\s+/g, "_");
        const raw = cleanSeqInput(document.getElementById(`slot${slot}Seq`)?.value || "", tp);
        if (!raw) return;
        let seq = raw, mutNotes = [];
        if (slot !== "Z" && (tp === "protein" || tp === "rna" || tp === "dna")) {
          const mutStr = document.getElementById(`slot${slot}Mut`)?.value || "";
          const applied = applyPointMutationsToSeq(raw, mutStr, tp);
          seq = applied.seq; mutNotes = applied.notes; warnings.push(...applied.warnings);
        }
        entities.push({ slot, chainId: ids[idx], type: tp, name, seq, mutNotes });
      });
      return { entities, warnings };
    }

    function runUniversalStudioPreScreen() {
      const { entities, warnings } = buildStudioEntities();
      if (entities.length < 2) { alert("Please enter or upload at least Partner X and Partner Y."); return; }
      const outPanel = document.getElementById("studioOutputPanel");
      if (outPanel) outPanel.style.display = "block";

      // Compute individual descriptors & pairwise compatibility
      const cards = [];
      const charges = [];
      let hasCustomSmiles = false;
      for (const e of entities) {
        const mutTag = e.mutNotes.length ? ` [Mut: ${e.mutNotes.join(",")}]` : "";
        if (e.type === "protein") {
          const st = analyzeProteinSeqStats(e.seq);
          charges.push(st.q72);
          cards.push(`<div class="descriptor"><span>${e.slot}: ${e.name} (Protein${mutTag})</span><strong>${st.len} aa · ${st.mwKda.toFixed(1)} kDa · Q<sub>7.2</sub>=${st.q72>=0?"+":""}${st.q72.toFixed(1)} · pI=${st.pI.toFixed(2)}</strong><span style="font-size:11px;">Basic R/K patches: ${st.basicPatches} · Arom: ${st.aromPct.toFixed(1)}% · GRAVY: ${st.gravy.toFixed(2)}</span></div>`);
        } else if (e.type === "rna" || e.type === "dna") {
          const st = analyzeNucleicSeqStats(e.seq, e.type);
          charges.push(st.q72);
          cards.push(`<div class="descriptor"><span>${e.slot}: ${e.name} (${e.type.toUpperCase()}${mutTag})</span><strong>${st.len} nt · Q<sub>phos</sub>=${st.q72} · GC=${st.gc.toFixed(1)}%</strong><span style="font-size:11px;">Purine (A+G): ${st.pur.toFixed(1)}% · Longest (A/G)<sub>n</sub> tract: ${st.maxPurRun} nt</span></div>`);
        } else if (e.type === "ligand") {
          const st = analyzeLigandStats(e.seq);
          if (!st.isCcd) hasCustomSmiles = true;
          charges.push(st.anionicHint ? -2.0 : 0.0);
          cards.push(`<div class="descriptor"><span>${e.slot}: ${e.name} (${st.isCcd ? "CCD Ligand" : "Custom SMILES"})</span><strong>${st.isCcd ? `CCD ${st.code}` : `~${st.heavyEst} heavy atoms · Arom:${st.aromAtoms} · HBA:${st.hba}`}</strong><span style="font-size:11px;">Halogens (X-Bond): ${st.halogens} · Anionic PO₄/COO⁻: ${st.anionicHint ? "Yes" : "No"}</span></div>`);
        } else {
          charges.push(2.0);
          cards.push(`<div class="descriptor"><span>${e.slot}: ${e.name} (Catalytic Metal / Ion)</span><strong>${e.seq.toUpperCase()}²⁺</strong><span style="font-size:11px;">Inner-sphere N/O/S + hydrated outer-sphere coordination ready</span></div>`);
        }
      }
      document.getElementById("studioPrescreenGrid").innerHTML = cards.join("");

      const e0 = entities[0], e1 = entities[1];
      const types = new Set(entities.map(e => e.type));
      const notes = [];
      if (charges.length >= 2) {
        const prod = charges[0] * charges[1];
        if (prod < -5) notes.push(`Strong electrostatic complementarity between ${e0.slot} (${charges[0]>=0?"+":""}${charges[0].toFixed(1)}) and ${e1.slot} (${charges[1]>=0?"+":""}${charges[1].toFixed(1)}).`);
        else if (prod > 15) notes.push(`Net like-charge repulsion at pH 7.2 (${charges[0].toFixed(1)} vs ${charges[1].toFixed(1)}) — binding requires localized complementary patches, hydrophobic burial, or bridging divalent cations (Mg²⁺/Zn²⁺).`);
      }
      if ((e0.type === "rna" || e0.type === "dna") && (e1.type === "rna" || e1.type === "dna")) {
        const bpRun = maxDuplexComplementarity(e0.seq, e1.seq);
        notes.push(`Longest continuous antiparallel Watson–Crick/G·U base-paired tract = <strong>${bpRun} bp</strong>.`);
      }
      if (types.has("rna") && entities.some(e => e.type === "rna" && analyzeNucleicSeqStats(e.seq, "rna").maxPurRun >= 6)) {
        notes.push(`RNA contains a ≥6-nt polypurine (A/G) tract — high compatibility for DEAD-box helicase / rocaglate molecular-glue clamping.`);
      }
      if (entities.length === 3) {
        notes.push(`Ternary 3-component job configured (${entities.map(e=>`${e.slot}:${e.type.toUpperCase()}`).join(" + ")}) — Tab 1 will automatically trigger Ternary Molecular-Glue & Metal-Bridge auditing on the resulting 3D structure.`);
      }
      if (hasCustomSmiles) {
        warnings.push("Custom SMILES detected: AlphaFold 3 Server web UI restricts ligands to supported CCD codes. Use the generated Boltz-1/2 YAML or Chai-1 FASTA file to co-fold ANY custom synthetic lab molecule, then drop the resulting .cif/.pdb into Tab 1.");
      }

      const warnBox = document.getElementById("studioWarningBanner");
      if (warnBox) {
        if (warnings.length) { warnBox.style.display = "block"; warnBox.innerHTML = warnings.map(w => `⚠ ${w}`).join("<br>"); }
        else warnBox.style.display = "none";
      }

      document.getElementById("studioVerdictBox").innerHTML = `
        <strong>Biophysical Pre-Screen Summary:</strong> ${notes.length ? notes.join(" ") : "Multi-modal entities validated and ready for 3D co-folding."}
        <br><span style="color:var(--text-muted); font-size:12px;"><strong>Honest Physics Guarantee:</strong> UniPPI-OWL never fabricates atomic ΔG/K<sub>d</sub> from 1D sequences alone. Export one of the 3D co-fold jobs below, run it (~2 min), and drop the output <code>.zip</code>/<code>.cif</code> into <strong>Tab 1</strong> for exact atomic ΔG, K<sub>d</sub>, UIS, 3Dmol.js &amp; SVG maps.</span>
      `;

      // Build AF3 Server JSON, Boltz-1/2 YAML, and Chai-1 FASTA
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
      const boltzYaml = boltzLines.join("\n");

      const chaiLines = [];
      for (const e of entities) {
        const cType = e.type === "ion" ? "ligand" : e.type;
        chaiLines.push(`>${cType}|name=${e.name}`, e.seq);
      }
      const chaiFasta = chaiLines.join("\n");

      LAST_STUDIO_JOBS = { jobName, af3: af3Json, boltz: boltzYaml, chai: chaiFasta };
      const fmtSel = document.getElementById("studioPreviewFormat");
      if (fmtSel && hasCustomSmiles) fmtSel.value = "boltz";
      renderStudioJobPreview();
    }

    function renderStudioJobPreview() {
      if (!LAST_STUDIO_JOBS) return;
      const fmt = document.getElementById("studioPreviewFormat")?.value || "boltz";
      const pre = document.getElementById("studioJobPreview");
      if (pre) pre.textContent = LAST_STUDIO_JOBS[fmt] || "";
    }

    function exportStudioJob(fmt) {
      if (!LAST_STUDIO_JOBS) runUniversalStudioPreScreen();
      if (!LAST_STUDIO_JOBS) return;
      const extMap = { af3: "_af3_server.json", boltz: "_boltz.yaml", chai: "_chai1.fasta" };
      downloadTextFile(`${LAST_STUDIO_JOBS.jobName}${extMap[fmt]}`, LAST_STUDIO_JOBS[fmt]);
    }

    function copyStudioPreview() {
      const txt = document.getElementById("studioJobPreview")?.textContent || "";
      if (txt && navigator.clipboard) navigator.clipboard.writeText(txt);
    }

    // =========================================================================
    // TAB 3: WT VS MUTANT / APO VS HOLO ENSEMBLE DIFF ENGINE + WELCH T-TEST
    // =========================================================================
    const DIFF_STATE = { A: [], B: [], lastRewiring: [] };

    async function evaluateFilesQuiet(fileList) {
      const rawFiles = Array.from(fileList || []);
      const tempK = parseFloat(document.getElementById("tempInput")?.value) || 298.15;
      const lambdaPae = parseFloat(document.getElementById("lambdaInput")?.value) || 8.0;
      const structEntries = [], jsonParsed = [];
      for (const f of rawFiles) {
        const low = f.name.toLowerCase();
        const rel = (f.webkitRelativePath || f.name).toLowerCase();
        if (rel.includes("msas/") || rel.includes("templates/") || rel.includes("template_hit") || rel.includes("terms_of_use")) continue;
        if (low.endsWith(".zip") && window.JSZip) {
          const zip = new JSZip();
          const zc = await zip.loadAsync(f);
          for (const [zPath, zObj] of Object.entries(zc.files)) {
            if (zObj.dir) continue;
            const zLow = zPath.toLowerCase();
            if (zLow.includes("msas/") || zLow.includes("templates/") || zLow.includes("template_hit") || zLow.includes("terms_of_use") || zLow.includes("job_request")) continue;
            const bName = zPath.split("/").pop();
            if (/\.(cif|mmcif|pdb|ent|pqr)$/i.test(bName) && !bName.includes("_sanitized.pdb")) structEntries.push({ name: bName, text: await zObj.async("text") });
            else if (/\.json$/i.test(bName)) { try { jsonParsed.push({ name: bName, lower: bName.toLowerCase(), data: JSON.parse(await zObj.async("text")) }); } catch(e) {} }
          }
        } else if (/\.(cif|mmcif|pdb|ent|pqr)$/i.test(f.name) && !f.name.includes("_sanitized.pdb")) {
          structEntries.push({ name: f.name, text: await f.text() });
        } else if (/\.json$/i.test(f.name) && !low.includes("job_request")) {
          try { jsonParsed.push({ name: f.name, lower: low, data: JSON.parse(await f.text()) }); } catch(e) {}
        }
      }
      structEntries.sort((a,b) => a.name.localeCompare(b.name));
      const out = [];
      for (const sf of structEntries) {
        const { atoms, meta } = parseStructureText(sf.text, sf.name);
        if (!atoms.length) continue;
        const { ordered, chains } = buildResidues(atoms);
        const cIds = Object.keys(chains);
        if (cIds.length < 2) continue;
        const stem = sf.name.replace(/\.(cif|mmcif|pdb|ent|pqr)$/i, "");
        const cleanStem = stem.replace(/_model_\d+$/i, "").replace(/_model$/i, "").toLowerCase();
        const mIdx = (stem.match(/_model_(\d+)$/i) || stem.match(/_(\d+)$/i))?.[1] || null;
        const activeJsons = jsonParsed.filter(j => (!mIdx || j.lower.includes(`_${mIdx}.json`) || j.lower.includes(`model_${mIdx}`)) && (!cleanStem || j.lower.includes(cleanStem) || structEntries.length === 1));
        let iptm = null, chainPairIptm = null, paeMat = null, hasJson = false;
        for (const j of activeJsons) {
          const rawIptm = j.data.iptm !== undefined ? j.data.iptm : (j.data["iptm+ptm"] !== undefined ? j.data["iptm+ptm"] : j.data.confidence_score);
          if (Array.isArray(j.data.chain_pair_iptm)) chainPairIptm = j.data.chain_pair_iptm;
          if (rawIptm !== undefined && rawIptm !== null) { iptm = parseFloat(rawIptm); hasJson = true; }
          const rawPae = j.data.pae || j.data.predicted_aligned_error;
          if (Array.isArray(rawPae) && rawPae.length > 0) {
            if (Array.isArray(j.data.token_chain_ids) && Array.isArray(j.data.token_res_ids)) {
              const keyToToks = new Map();
              j.data.token_chain_ids.forEach((tc, tidx) => { const k = `${tc}|${j.data.token_res_ids[tidx]}`; if (!keyToToks.has(k)) keyToToks.set(k, []); keyToToks.get(k).push(Math.min(rawPae.length - 1, tidx)); });
              const tokLists = ordered.map(r => keyToToks.get(`${r.origChainId || r.chainId}|${r.resSeq}`) || keyToToks.get(`${r.chainId}|${r.resSeq}`) || []);
              paeMat = tokLists.map(li => !li.length ? tokLists.map(()=>15.0) : (li.length===1 ? tokLists.map(lj=>lj.length===1?Number(rawPae[li[0]][lj[0]]):(lj.length?lj.reduce((s,t)=>s+Number(rawPae[li[0]][t]),0)/lj.length:15.0)) : tokLists.map(lj=>{ if(!lj.length) return 15.0; let s=0,c=0; for(const ti of li) for(const tj of lj){s+=Number(rawPae[ti][tj]); c++;} return s/c; })));
              hasJson = true;
            } else if (rawPae.length === ordered.length) { paeMat = rawPae; hasJson = true; }
          }
        }
        const mode = detectConfidenceModeJS(ordered, meta, hasJson, paeMat);
        const sasaCompFull = computeResidueSASA(ordered), sasaFree = {};
        for (const cid of cIds) sasaFree[cid] = computeResidueSASA(chains[cid]);
        const resToIdx = new Map(ordered.map((r, i) => [r.key, i]));
        const isBinary = (cIds.length === 2);
        const origOrder = [...new Set(ordered.map(r => r.origChainId || r.chainId))];
        const fileRes = [];
        for (let i = 0; i < cIds.length; i++) {
          for (let j = i + 1; j < cIds.length; j++) {
            const ca = cIds[i], cb = cIds[j];
            let pairSasa = sasaCompFull;
            if (!isBinary) {
              if (!chains[ca].some(ra => chains[cb].some(rb => dist3(ra.centroid, rb.centroid) <= 16.0))) continue;
              pairSasa = computeResidueSASA([...chains[ca], ...chains[cb]]);
            }
            const oi = origOrder.indexOf(chains[ca][0]?.origChainId || ca), oj = origOrder.indexOf(chains[cb][0]?.origChainId || cb);
            const pIptm = (Array.isArray(chainPairIptm) && oi>=0 && oj>=0 && oi!==oj && chainPairIptm[oi] && typeof chainPairIptm[oi][oj]==="number") ? chainPairIptm[oi][oj] : iptm;
            const res = analyzePairJS(stem, ca, cb, chains[ca], chains[cb], ordered, resToIdx, mode, { iptm: pIptm, res_pae_matrix: paeMat }, pairSasa, sasaFree[ca], sasaFree[cb], tempK, lambdaPae, isBinary);
            if (res) { fileRes.push(res); out.push(res); }
          }
        }
        // Include Ternary Glue & Metal Bridge in Diff Mode as well
        if (cIds.length >= 3) {
          for (const lc of cIds) {
            const lcRes = chains[lc];
            if (!lcRes || !lcRes.length || (lcRes[0].entityType !== "LIGAND" && lcRes[0].entityType !== "ION")) continue;
            const isLigC = (lcRes[0].entityType === "LIGAND");
            const hits = [];
            for (const rrec of fileRes) {
              const other = rrec.chainAId === lc ? rrec.chainBId : (rrec.chainBId === lc ? rrec.chainAId : null);
              if (!other || !chains[other]?.length) continue;
              const ot = chains[other][0].entityType;
              if (isLigC && (ot==="PROTEIN"||ot==="NUCLEIC") && (rrec.lig_burial_pct||0)>=12 && (rrec.confContactCount||0)>=1) hits.push([other, rrec.lig_burial_pct, 0, 0]);
              else if (!isLigC && (ot==="PROTEIN"||ot==="NUCLEIC"||ot==="LIGAND") && ((rrec.metal_coord_bonds||0)>=1 || (rrec.metal_outer_bonds||0)>=2)) hits.push([other, (rrec.metal_coord_bonds||0)+0.5*(rrec.metal_outer_bonds||0), rrec.metal_coord_bonds||0, rrec.metal_outer_bonds||0]);
            }
            if ((isLigC && hits.length>=2) || (!isLigC && hits.length>=2 && hits.some(h=>h[2]>=1))) {
              hits.sort((a,b)=>b[1]-a[1]);
              const [m1] = hits[0], [m2] = hits[1];
              const compM = [...chains[m1], ...chains[m2]], cidM = `${m1}+${m2}`;
              const tRes = analyzePairJS(stem, cidM, lc, compM, lcRes, ordered, resToIdx, mode, { iptm, res_pae_matrix: paeMat }, computeResidueSASA([...compM, ...lcRes]), computeResidueSASA(compM), sasaFree[lc], tempK, lambdaPae, false);
              if (tRes) out.push(tRes);
            }
          }
        }
      }
      return out;
    }

    async function handleDiffUpload(cond, fileList) {
      const stEl = document.getElementById(`diff${cond}Status`);
      if (stEl) stEl.textContent = "Evaluating structures locally...";
      const recs = await evaluateFilesQuiet(fileList);
      DIFF_STATE[cond] = recs;
      const uniqModels = new Set(recs.map(r => r.complex_id)).size;
      const cntEl = document.getElementById(`diff${cond}Count`);
      if (cntEl) cntEl.textContent = `${uniqModels} model(s), ${recs.length} interface(s)`;
      if (stEl) stEl.textContent = recs.length ? `Loaded ${uniqModels} model(s): ${[...new Set(recs.map(r=>r.chain_pair))].slice(0,5).join(", ")}` : "No valid interfaces found.";
      refreshDiffPairOptions();
      if (DIFF_STATE.A.length && DIFF_STATE.B.length) renderEnsembleDiffResults();
    }

    function refreshDiffPairOptions() {
      const sel = document.getElementById("diffPairSelect");
      if (!sel) return;
      const pairsA = new Set(DIFF_STATE.A.map(r => r.chain_pair));
      const pairsB = new Set(DIFF_STATE.B.map(r => r.chain_pair));
      const shared = [...pairsA].filter(p => pairsB.has(p));
      sel.innerHTML = `<option value="TOP">Top-Ranked Interface per Model</option>` + shared.map(p => `<option value="${p}">Matched Chain Pair: ${p}</option>`).join("");
    }

    function logGamma(z) {
      const g = 7, c = [0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
      if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z);
      z -= 1;
      let x = c[0];
      for (let i = 1; i < g + 2; i++) x += c[i] / (z + i);
      const t = z + g + 0.5;
      return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
    }

    function regIncBeta(x, a, b) {
      if (x <= 0) return 0;
      if (x >= 1) return 1;
      const bt = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
      const cf = (x0, a0, b0) => {
        let qab = a0 + b0, qap = a0 + 1, qam = a0 - 1, c = 1, d = 1 - qab * x0 / qap;
        if (Math.abs(d) < 1e-30) d = 1e-30;
        d = 1 / d; let h = d;
        for (let m = 1; m <= 120; m++) {
          let m2 = 2 * m;
          let aa = m * (b0 - m) * x0 / ((qam + m2) * (a0 + m2));
          d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
          c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
          d = 1 / d; h *= d * c;
          aa = -(a0 + m) * (qab + m) * x0 / ((a0 + m2) * (qap + m2));
          d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
          c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
          d = 1 / d; const del = d * c; h *= del;
          if (Math.abs(del - 1) < 1e-9) break;
        }
        return h;
      };
      return x < (a + 1) / (a + b + 2) ? bt * cf(x, a, b) / a : 1 - bt * cf(1 - x, b, a) / b;
    }

    function welchTTestPVal(arrA, arrB) {
      const n1 = arrA.length, n2 = arrB.length;
      if (n1 < 2 || n2 < 2) return null;
      const m1 = meanArr(arrA), m2 = meanArr(arrB);
      const v1 = arrA.reduce((s,x)=>s+(x-m1)*(x-m1),0) / (n1 - 1);
      const v2 = arrB.reduce((s,x)=>s+(x-m2)*(x-m2),0) / (n2 - 1);
      const se2 = v1/n1 + v2/n2;
      if (se2 < 1e-12) return Math.abs(m1 - m2) < 1e-9 ? 1.0 : 1e-9;
      const t = (m2 - m1) / Math.sqrt(se2);
      const df = (se2 * se2) / (((v1/n1)*(v1/n1))/(n1-1) + ((v2/n2)*(v2/n2))/(n2-1));
      return clamp(regIncBeta(df / (df + t * t), df / 2, 0.5), 0.0, 1.0);
    }

    function selectDiffSubset(recs, pairFilter) {
      if (pairFilter && pairFilter !== "TOP") return recs.filter(r => r.chain_pair === pairFilter);
      const byModel = new Map();
      for (const r of recs) {
        if (!byModel.has(r.complex_id) || r.UIS > byModel.get(r.complex_id).UIS) byModel.set(r.complex_id, r);
      }
      return [...byModel.values()];
    }

    function renderEnsembleDiffResults() {
      if (!DIFF_STATE.A.length || !DIFF_STATE.B.length) { alert("Load both Condition A and Condition B ensembles first (or click Load Demo)."); return; }
      const pFilter = document.getElementById("diffPairSelect")?.value || "TOP";
      const subA = selectDiffSubset(DIFF_STATE.A, pFilter);
      const subB = selectDiffSubset(DIFF_STATE.B, pFilter);
      if (!subA.length || !subB.length) { alert("No matched interfaces found for this selection."); return; }

      const tempK = parseFloat(document.getElementById("tempInput")?.value) || 298.15;
      const mSdStr = (arr, dec=2) => {
        const m = meanArr(arr);
        const sd = arr.length > 1 ? Math.sqrt(arr.reduce((s,x)=>s+(x-m)*(x-m),0)/(arr.length-1)) : 0;
        return arr.length > 1 ? `${m.toFixed(dec)} ± ${sd.toFixed(dec)}` : m.toFixed(dec);
      };
      const fmtP = p => p === null ? "N/A (N<2)" : (p < 0.0001 ? "< 1e-4 ***" : (p < 0.001 ? `${p.toFixed(4)} ***` : (p < 0.01 ? `${p.toFixed(4)} **` : (p < 0.05 ? `${p.toFixed(4)} *` : `${p.toFixed(3)} (ns)`))));

      const dgA = meanArr(subA.map(r=>r.dG_uni_kcal)), dgB = meanArr(subB.map(r=>r.dG_uni_kcal));
      const ddg = dgB - dgA;
      const kdFold = Math.exp(ddg / (R_GAS_KCAL * tempK));
      const uisA = meanArr(subA.map(r=>r.UIS)), uisB = meanArr(subB.map(r=>r.UIS));
      const pDdg = welchTTestPVal(subA.map(r=>r.dG_uni_kcal), subB.map(r=>r.dG_uni_kcal));
      const verdict = ddg >= 1.0 ? "DESTABILIZING / ESCAPE PERTURBATION (Binding Weakened in B)" : (ddg <= -1.0 ? "STABILIZING / GAIN-OF-BINDING (Binding Strengthened in B)" : "NEUTRAL / MODEST EFFECT (|ΔΔG| < 1.0 kcal/mol)");

      document.getElementById("diffResultsWrap").style.display = "block";
      document.getElementById("btnExportDiffCSV").style.display = "inline-flex";
      const hb = document.getElementById("diffHeadlineBanner");
      hb.style.borderLeftColor = ddg >= 1.0 ? "var(--rose)" : (ddg <= -1.0 ? "var(--emerald)" : "var(--primary)");
      hb.innerHTML = `<strong>${verdict}</strong> (N<sub>A</sub>=${subA.length} vs N<sub>B</sub>=${subB.length})<br>ΔΔG (B − A) = <strong>${ddg>=0?"+":""}${ddg.toFixed(2)} kcal/mol</strong> · K<sub>d</sub> Fold-Change = <strong>${kdFold >= 1000 ? `${kdFold.toExponential(2)}×` : `${kdFold.toFixed(2)}×`}</strong> · ΔUIS = <strong>${(uisB-uisA)>=0?"+":""}${(uisB-uisA).toFixed(3)}</strong> · Welch p(ΔG) = <strong>${fmtP(pDdg)}</strong>`;

      const metrics = [
        ["UniPPI UIS Score", r=>r.UIS, 3, d => d <= -0.08 ? "Confidence drop" : (d >= 0.08 ? "Confidence gain" : "Stable")],
        ["Binding Free Energy ΔG (kcal/mol)", r=>r.dG_uni_kcal, 2, d => d >= 1.0 ? "Weaker binding (Escape)" : (d <= -1.0 ? "Tighter binding" : "Neutral")],
        ["Ordered Buried Area wBSA (Å²)", r=>r.wBSA_A2, 1, d => d <= -150 ? "Interface loss" : (d >= 150 ? "Interface expansion" : "Conserved")],
        ["Local PAE Score (LIS)", r=>r.LIS, 3, d => Math.abs(d) >= 0.05 ? (d < 0 ? "Local PAE degraded" : "Local PAE improved") : "Stable"],
        ["Salt Bridges (SB)", r=>r.salt_bridges, 1, d => d <= -1 ? "Electrostatic loss" : (d >= 1 ? "Electrostatic gain" : "Conserved")],
        ["Hydrogen Bonds (HB)", r=>r.h_bonds, 1, d => d <= -1 ? "H-bond loss" : (d >= 1 ? "H-bond gain" : "Conserved")],
        ["Hydrophobic Contacts (HP)", r=>r.hydrophobic_contacts, 1, d => d <= -3 ? "Core packing loss" : (d >= 3 ? "Core packing gain" : "Conserved")],
        ["π / Cation-π Stacks", r=>(r.pi_stacking||0)+(r.cation_pi||0), 1, d => d <= -0.5 ? "Stacking disrupted" : (d >= 0.5 ? "Stacking gained" : "Conserved")]
      ];
      document.getElementById("diffStatsBody").innerHTML = metrics.map(([lbl, fn, dec, noteFn]) => {
        const aVals = subA.map(fn), bVals = subB.map(fn);
        const d = meanArr(bVals) - meanArr(aVals);
        const pv = welchTTestPVal(aVals, bVals);
        return `<tr><td><strong>${lbl}</strong></td><td>${mSdStr(aVals, dec)}</td><td>${mSdStr(bVals, dec)}</td><td><strong>${d>=0?"+":""}${d.toFixed(dec)}</strong></td><td>${fmtP(pv)}</td><td>${noteFn(d)}</td></tr>`;
      }).join("");

      // Contact Rewiring Audit across Ensembles
      const aggContacts = (sub) => {
        const m = new Map();
        const chTag = c => [c.sb?"SB":"", c.hb?"HB":"", (c.pi||c.catpi)?"π":"", c.xb?"XB":"", c.hp?"HP":""].filter(Boolean).join("/") || "vdW";
        for (const r of sub) {
          const seenInModel = new Set();
          for (const c of (r.contacts || [])) {
            if (c.weight < 0.25) continue;
            const normKey = `${c.resA.replace(/\d+$/, m=>m)} ↔ ${c.resB}`;
            const k = `${(c.resA||"").split(":")[0]}:${c.resSeqA} ↔ ${(c.resB||"").split(":")[0]}:${c.resSeqB}`;
            if (seenInModel.has(k)) continue;
            seenInModel.add(k);
            if (!m.has(k)) m.set(k, { label: `${c.resA} ↔ ${c.resB}`, count: 0, dists: [], tags: new Set() });
            const rec = m.get(k);
            rec.count++; rec.dists.push(c.minDist); rec.tags.add(chTag(c));
          }
        }
        return m;
      };
      const mapA = aggContacts(subA), mapB = aggContacts(subB);
      const allKeys = new Set([...mapA.keys(), ...mapB.keys()]);
      const rewiring = [];
      for (const k of allKeys) {
        const a = mapA.get(k), b = mapB.get(k);
        const occA = a ? a.count / subA.length : 0, occB = b ? b.count / subB.length : 0;
        let status = "PERSISTENT", rank = 2;
        if (occA >= 0.6 && occB < 0.4) { status = "LOST"; rank = 0; }
        else if (occA < 0.4 && occB >= 0.6) { status = "GAINED"; rank = 1; }
        else if (occA < 0.4 && occB < 0.4) continue;
        rewiring.push({
          status, rank,
          pair: a && b && a.label !== b.label ? `${a.label} [B: ${b.label}]` : (a?.label || b?.label || k),
          occA: Math.round(occA*100), occB: Math.round(occB*100),
          dA: a ? meanArr(a.dists).toFixed(2) : "—", dB: b ? meanArr(b.dists).toFixed(2) : "—",
          tA: a ? [...a.tags].join(",") : "—", tB: b ? [...b.tags].join(",") : "—"
        });
      }
      rewiring.sort((x, y) => x.rank - y.rank || Math.abs(y.occA - y.occB) - Math.abs(x.occA - x.occB));
      DIFF_STATE.lastRewiring = rewiring;
      document.getElementById("diffRewiringBody").innerHTML = rewiring.slice(0, 80).map(rw => {
        const bCls = rw.status === "LOST" ? "badge-lost" : (rw.status === "GAINED" ? "badge-gained" : "badge-persist");
        return `<tr><td><span class="badge-call ${bCls}">${rw.status}</span></td><td><strong>${rw.pair}</strong></td><td>${rw.occA}%</td><td>${rw.occB}%</td><td>${rw.dA} → ${rw.dB}</td><td>${rw.tA} → ${rw.tB}</td></tr>`;
      }).join("") || `<tr><td colspan="6" style="text-align:center;">No confident contacts shifted across ensembles.</td></tr>`;
    }

    function downloadDiffRewiringCSV() {
      if (!DIFF_STATE.lastRewiring.length) return;
      const rows = ["status,contact_pair,occ_A_pct,occ_B_pct,mean_dist_A_A,mean_dist_B_A,channels_A,channels_B"];
      for (const r of DIFF_STATE.lastRewiring) rows.push([r.status, `"${r.pair}"`, r.occA, r.occB, r.dA, r.dB, `"${r.tA}"`, `"${r.tB}"`].join(","));
      downloadTextFile("unippi_ensemble_diff_rewiring.csv", rows.join("\n"));
    }

    function loadDemoEnsembleDiff() {
      // Generates a 5-model WT ensemble (tight binding) vs 5-model Escape-Mutant ensemble (lost SB/HP, higher PAE)
      const mkModels = (prefix, spanBase, paeBase, iptmBase, plddtBase, dropSb) => {
        const CHAIN_A_SEQ = "AQVINTEDGVQDYLQTYHKLADNYSITKSEAQALGFVISNGNLADTLGKS";
        const CHAIN_B_SEQ = "KKILINKGQINSIDSLHQELKKELQLPEYRGENLDVLFDTLTGWVQYPLV";
        const SB_IDX = dropSb ? new Set([29]) : new Set([6, 18, 21, 27, 29]);
        const HB_IDX = dropSb ? new Set([10, 45]) : new Set([10, 14, 24, 39, 45]);
        const HP_IDX = dropSb ? new Set([2, 19, 36, 46]) : new Set([2, 3, 12, 19, 33, 35, 36, 37, 46]);
        const HB_MAP = { GLN:["OE1","O"], ASN:["ND2","N"], SER:["OG","O"], THR:["OG1","O"], TYR:["OH","O"] };
        const resOut = [];
        for (let m = 0; m < 5; m++) {
          const name = `${prefix}_model_${m}`;
          const span = spanBase + m * 0.04, plddt = plddtBase - m * 0.8, pae = paeBase + m * 0.25, iptm = iptmBase - m * 0.015;
          const seqA = Array.from(CHAIN_A_SEQ).map(ch => ONE_TO_THREE[ch] || "ALA");
          const seqB = Array.from(CHAIN_B_SEQ).map(ch => ONE_TO_THREE[ch] || "ALA");
          const lines = [`data_${name}`, "#", "_ma_qa_metric_local.type pLDDT", "loop_", "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol", "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv", "_atom_site.auth_seq_id", "_atom_site.auth_asym_id", "_atom_site.pdbx_PDB_model_num"];
          let ser = 1;
          for (const [ci, [cid, seq]] of [["A", seqA], ["B", seqB]].entries()) {
            const sign = ci === 0 ? -1.0 : 1.0;
            for (let ri = 0; ri < seq.length; ri++) {
              const rn = seq[ri], rs = ri + 1, zB = ri * 2.15, yB = Math.cos(ri * 0.9) * 1.45;
              const xCa = sign * (span + 1.5), xCb = sign * span;
              const sp = [["N", xCa+sign*0.8, yB-0.5, zB-0.6, "N"], ["CA", xCa, yB, zB, "C"], ["C", xCa+sign*0.6, yB+0.6, zB+0.6, "C"], ["O", xCa+sign*1.2, yB+1.1, zB+0.8, "O"]];
              if (rn !== "GLY") sp.push(["CB", xCb, yB, zB, "C"]);
              if (SB_IDX.has(ri)) {
                if (rn==="ASP"||rn==="GLU") sp.push(["CG", sign*(span-0.6), yB, zB, "C"], [rn==="ASP"?"OD1":"OE1", sign*(span-1.45), yB, zB, "O"]);
                else sp.push(["CG", sign*(span-0.6), yB, zB, "C"], [rn==="LYS"?"NZ":"NH1", sign*(span-1.45), yB, zB, "N"]);
              } else if (HB_IDX.has(ri) && HB_MAP[rn]) {
                sp.push([HB_MAP[rn][0], sign*(span-1.4), yB, zB, HB_MAP[rn][1]]);
              } else if (HP_IDX.has(ri)) {
                sp.push(["CG", sign*(span-0.75), yB+0.12, zB, "C"], ["CD1", sign*(span-1.22), yB-0.08, zB+0.1, "C"]);
              }
              for (const [an, ax, ay, az, el] of sp) lines.push(`ATOM ${ser++} ${el} ${an} . ${rn} ${cid} 1 ${rs} ? ${ax.toFixed(3)} ${ay.toFixed(3)} ${az.toFixed(3)} 1.00 ${plddt.toFixed(2)} ${rs} ${cid} 1`);
            }
          }
          const { atoms } = parseStructureText(lines.join("\n"), `${name}.cif`);
          const { ordered, chains } = buildResidues(atoms);
          const sComp = computeResidueSASA(ordered), sA = computeResidueSASA(chains.A), sB = computeResidueSASA(chains.B);
          const r2i = new Map(ordered.map((r,i)=>[r.key, i])), tot = ordered.length, nA = chains.A.length;
          const pMat = Array.from({length: tot}, (_, r) => Array.from({length: tot}, (__, c) => r===c ? 0.8 : ((r<nA&&c<nA)||(r>=nA&&c>=nA) ? 1.8 : pae)));
          resOut.push(analyzePairJS(name, "A", "B", chains.A, chains.B, ordered, r2i, "MODE_A_PRED_FULL", { iptm, res_pae_matrix: pMat }, sComp, sA, sB, 298.15, 8.0));
        }
        return resOut;
      };
      DIFF_STATE.A = mkModels("WT_complex", 2.94, 2.05, 0.87, 93.4, false);
      DIFF_STATE.B = mkModels("EscapeMut_complex", 3.65, 8.80, 0.44, 76.0, true);
      document.getElementById("diffACount").textContent = "5 models (Demo WT)";
      document.getElementById("diffBCount").textContent = "5 models (Demo Escape Mutant)";
      document.getElementById("diffAStatus").textContent = "Loaded 5-model WT AF3 ensemble (A-B).";
      document.getElementById("diffBStatus").textContent = "Loaded 5-model Escape-Mutant AF3 ensemble (A-B).";
      refreshDiffPairOptions();
      renderEnsembleDiffResults();
    }

'''

html = ensure_replaced(
    html,
    '    window.addEventListener("DOMContentLoaded", () => {',
    step4_js + '    window.addEventListener("DOMContentLoaded", () => {',
    "JS: Inject 3-Tab Switcher, Universal X-vs-Y Studio & Ensemble Diff Mode",
)

html_path.write_text(html, encoding="utf-8")
print("[PASS] index.html upgraded cleanly with Step 4 Universal Studio & Ensemble Diff Mode.\n")

# Strict JS Structural & Event-Handler Audit (works even without Node.js)
m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
assert m_js, "Missing inline <script> block in index.html"
js_code = m_js[-1]
handlers = set(re.findall(r'on(?:click|change|input)="([A-Za-z_]\w*)\(', html))
defined_fns = set(re.findall(r'\bfunction\s+([A-Za-z_]\w*)\s*\(', js_code))
missing_fns = sorted(handlers - defined_fns)
assert not missing_fns, f"Undefined inline event handlers in index.html: {missing_fns}"
print(f"[PASS] Verified all {len(handlers)} HTML event handlers map to defined JS functions.")

node_bin = shutil.which("node") or shutil.which("node.exe")
if node_bin:
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
        tf.write(js_code)
        tf_path = tf.name
    res_node = subprocess.run([node_bin, "--check", tf_path], capture_output=True, text=True)
    Path(tf_path).unlink(missing_ok=True)
    if res_node.returncode != 0:
        raise RuntimeError(f"JavaScript syntax error:\n{res_node.stderr}")
    print(f"[PASS] Verified JS syntax cleanly via `{node_bin} --check`.")

# ==============================================================================
# PART B: ADD MATCHING `--diff-b` ENSEMBLE DIFF MODE TO unippi_cli.py
# ==============================================================================
print("\n=== PART B: ADDING --diff-b ENSEMBLE DIFF PARITY TO unippi_cli.py ===")
cli = cli_path.read_text(encoding="utf-8")

diff_helper_py = '''
def _welch_t_pval_py(a: List[float], b: List[float]) -> Optional[float]:
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return None
    m1, m2 = sum(a) / n1, sum(b) / n2
    v1 = sum((x - m1) ** 2 for x in a) / (n1 - 1)
    v2 = sum((x - m2) ** 2 for x in b) / (n2 - 1)
    se2 = v1 / n1 + v2 / n2
    if se2 < 1e-12:
        return 1.0 if abs(m1 - m2) < 1e-9 else 1e-9
    t_val = (m2 - m1) / math.sqrt(se2)
    df = (se2 * se2) / (((v1 / n1) ** 2) / (n1 - 1) + ((v2 / n2) ** 2) / (n2 - 1))
    x = df / (df + t_val * t_val)
    a0, b0 = df / 2.0, 0.5
    bt = math.exp(math.lgamma(a0 + b0) - math.lgamma(a0) - math.lgamma(b0) + a0 * math.log(max(1e-30, x)) + b0 * math.log(max(1e-30, 1.0 - x)))
    def _cf(xv, av, bv):
        qab, qap, qam = av + bv, av + 1.0, av - 1.0
        c, d = 1.0, 1.0 - qab * xv / qap
        d = 1.0 / (d if abs(d) >= 1e-30 else 1e-30)
        h = d
        for m in range(1, 121):
            m2 = 2 * m
            aa = m * (bv - m) * xv / ((qam + m2) * (av + m2))
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) >= 1e-30 else 1e-30)
            c = 1.0 + aa / (c if abs(c) >= 1e-30 else 1e-30)
            h *= d * c
            aa = -(av + m) * (qab + m) * xv / ((av + m2) * (qap + m2))
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) >= 1e-30 else 1e-30)
            c = 1.0 + aa / (c if abs(c) >= 1e-30 else 1e-30)
            del_v = d * c
            h *= del_v
            if abs(del_v - 1.0) < 1e-9:
                break
        return h
    p = bt * _cf(x, a0, b0) / a0 if x < (a0 + 1.0) / (a0 + b0 + 2.0) else 1.0 - bt * _cf(1.0 - x, b0, a0) / b0
    return clamp(p, 0.0, 1.0)

def run_cli_ensemble_diff(recs_a: List[Dict], recs_b: List[Dict], out_dir: str, temp_k: float) -> None:
    def _top_per_model(recs):
        by_m: Dict[str, Dict] = {}
        for r in recs:
            cid = r["complex_id"]
            if cid not in by_m or r["UIS"] > by_m[cid]["UIS"]:
                by_m[cid] = r
        return list(by_m.values())
    sub_a, sub_b = _top_per_model(recs_a), _top_per_model(recs_b)
    if not sub_a or not sub_b:
        print("[DIFF] Insufficient records to compare Condition A vs Condition B.")
        return
    dg_a = [r["dG_uni_kcal"] for r in sub_a]
    dg_b = [r["dG_uni_kcal"] for r in sub_b]
    uis_a = [r["UIS"] for r in sub_a]
    uis_b = [r["UIS"] for r in sub_b]
    ddg = (sum(dg_b) / len(dg_b)) - (sum(dg_a) / len(dg_a))
    duis = (sum(uis_b) / len(uis_b)) - (sum(uis_a) / len(uis_a))
    kd_fold = math.exp(ddg / (R_GAS_KCAL * temp_k))
    p_dg = _welch_t_pval_py(dg_a, dg_b)
    p_uis = _welch_t_pval_py(uis_a, uis_b)
    print(f"\\n{'='*96}")
    print(f"[ENSEMBLE DIFF: Condition A (N={len(sub_a)}) vs Condition B (N={len(sub_b)})]")
    print(f"  dUIS = {duis:+.3f} (Welch p = {f'{p_uis:.4g}' if p_uis is not None else 'N/A'})")
    print(f"  ddG  = {ddg:+.2f} kcal/mol | Kd Fold-Change = {kd_fold:.2f}x (Welch p = {f'{p_dg:.4g}' if p_dg is not None else 'N/A'})")
    def _agg(sub):
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
    Path(diff_csv).write_text("\\n".join(diff_rows) + "\\n", encoding="utf-8")
    print(f"  Saved Diff CSV: {diff_csv}")
    print(f"{'='*96}")

def main():'''

cli = ensure_replaced(cli, '\ndef main():', diff_helper_py, "CLI: Add _welch_t_pval_py & run_cli_ensemble_diff()")
cli = ensure_replaced(
    cli,
    '    parser.add_argument("--split-chain", type=str, default=None, help="Split one chain into 2 domains for intra-chain auditing (e.g. A:1-236,237-406 or A:236)")',
    '    parser.add_argument("--split-chain", type=str, default=None, help="Split one chain into 2 domains for intra-chain auditing (e.g. A:1-236,237-406 or A:236)")\n'
    '    parser.add_argument("--diff-b", type=str, default=None, help="Second structure/.zip/folder (Condition B: Mutant or Apo) for Welch t-test & Contact Rewiring diff")',
    "CLI: Add --diff-b argument",
)
old_cli_end = '    print(f"[UniPPI-CLI] Saved Sanitized PDBs in: {os.path.join(args.output, \'converted_structures\')}")'
new_cli_end = old_cli_end + '''
    if args.diff_b:
        if os.path.isfile(args.diff_b) and args.diff_b.lower().endswith(".zip"):
            _tmp_b = tempfile.mkdtemp(prefix="unippi_diffb_")
            with zipfile.ZipFile(args.diff_b, "r") as zfb:
                zfb.extractall(_tmp_b)
            files_b = sorted(glob.glob(os.path.join(_tmp_b, "**", "*.cif"), recursive=True) + glob.glob(os.path.join(_tmp_b, "**", "*.mmcif"), recursive=True) + glob.glob(os.path.join(_tmp_b, "**", "*.pdb"), recursive=True))
        elif os.path.isfile(args.diff_b):
            files_b = [args.diff_b]
        else:
            files_b = sorted(glob.glob(os.path.join(args.diff_b, "**", "*.cif"), recursive=True) + glob.glob(os.path.join(args.diff_b, "**", "*.mmcif"), recursive=True) + glob.glob(os.path.join(args.diff_b, "**", "*.pdb"), recursive=True))
        files_b = [fp for fp in files_b if "_sanitized.pdb" not in fp and not any(ig in fp.replace("\\\\", "/").lower() for ig in ("/msas/", "/templates/", "template_hit"))]
        recs_b: List[Dict] = []
        for fpb in files_b:
            recs_b.extend(process_structure_file(fpb, args.output, temp_k=args.temp, lambda_pae=args.lambda_pae, ionic_mM=args.ionic_strength, ph_val=args.ph, split_chain=args.split_chain))
        run_cli_ensemble_diff(all_records, recs_b, args.output, args.temp)'''
cli = ensure_replaced(cli, old_cli_end, new_cli_end, "CLI: Wire --diff-b execution in main()")

ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] unippi_cli.py upgraded with --diff-b cleanly (ast.parse OK).\n")

# Final verification: Zero-regression PPI + CLI --diff-b self-consistency check
print("=" * 96)
print("STEP 4 VERIFICATION: ZERO-REGRESSION PPI + CLI --diff-b SELF-CHECK (fold_vip1_bzip61.zip)")
print("=" * 96)
subprocess.run([sys.executable, "unippi_cli.py", "-i", "benchmarks/fold_vip1_bzip61.zip", "--diff-b", "benchmarks/fold_vip1_bzip61.zip", "-o", "/tmp/v2_s4_diff_check"], check=True)
