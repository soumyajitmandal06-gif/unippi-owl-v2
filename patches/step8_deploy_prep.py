#!/usr/bin/env python3
"""
Step 8 (Final Pre-Deployment Patch & New GitHub Repo Preparation):
  1. Guards Merge Mode against fusing AF3/Boltz multi-model .zip or _model_0..4 ensembles.
  2. Auto-switches to Tab 1 when Hero 'Load Sample' is clicked from Tab 2/3.
  3. Updates remaining 'UniPPI-v2.0 Engine' card label & aligns hotspots CSV headers.
  4. Rewrites canonical/OG/JSON-LD URLs + generates clean publication README.md for the new repo.
"""
import re
import sys
from pathlib import Path

repo_name = sys.argv[1].strip() if len(sys.argv) > 1 else "unippi-owl-v2"
gh_user = "soumyajitmandal06-gif"
live_url = f"https://{gh_user}.github.io/{repo_name}/"
repo_url = f"https://github.com/{gh_user}/{repo_name}"

html_path = Path("index.html")
html = html_path.read_text(encoding="utf-8")


def rep(text, old, new, lbl):
    if new in text:
        print(f"  [OK-ALREADY] {lbl}")
        return text
    assert old in text, f"Anchor not found for [{lbl}]"
    print(f"  [APPLIED]    {lbl}")
    return text.replace(old, new, 1)


# 1. Tag ZIP entries with source archive name and block accidental multi-model ensemble merging
html = rep(
    html,
    '                structEntries.push({ name: bName, text: await zObj.async("text") });',
    '                structEntries.push({ name: bName, text: await zObj.async("text"), fromZip: f.name });',
    "JS: Tag ZIP-extracted structures with source archive name",
)
html = rep(
    html,
    '''      const mergeMode = document.getElementById("mergeUploadedCheck")?.checked;
      if (mergeMode && structEntries.length >= 2) {''',
    '''      const mergeMode = document.getElementById("mergeUploadedCheck")?.checked;
      const sameZipEns = structEntries.length >= 2 && structEntries.every(s => s.fromZip && s.fromZip === structEntries[0].fromZip);
      const modelEns = structEntries.length >= 2 && structEntries.every(s => /_model_\\d+\\.(cif|mmcif|pdb|ent|pqr)$/i.test(s.name));
      let skippedMergeNote = "";
      if (mergeMode && structEntries.length >= 2 && (sameZipEns || modelEns)) {
        skippedMergeNote = `<div style="color:var(--primary); font-weight:600;">[Auto-Guard] Skipped Merge Mode because uploaded files are a multi-model ensemble (${structEntries.length} models); scoring each model independently.</div>`;
      } else if (mergeMode && structEntries.length >= 2) {''',
    "JS: Prevent Merge Mode from fusing multi-model AF3/Boltz ensembles",
)
html = rep(
    html,
    '      if (qList && queueLogs.length) {\n        qList.innerHTML = queueLogs.join("")',
    '      if (qList && queueLogs.length) {\n        qList.innerHTML = (skippedMergeNote || "") + queueLogs.join("")',
    "JS: Display Auto-Guard notice in batch manifest when merge is skipped",
)

# 2. Switch to Tab 1 when Hero 'Load Sample' is clicked from Tab 2/3
html = rep(
    html,
    '''    function loadSampleAF3Ensemble() {
      LAST_UPLOADED_FILES = null;''',
    '''    function loadSampleAF3Ensemble() {
      switchWorkspaceTab("auditor");
      LAST_UPLOADED_FILES = null;''',
    "JS: Auto-switch to Tab 1 when Load Sample is triggered",
)

# 3. Polish remaining v2.0 badge & hotspots CSV header
html = rep(
    html,
    '<div class="model-card-head"><span>UniPPI-v2.0 Engine</span><span>CLIENT-SIDE</span></div>',
    '<div class="model-card-head"><span>UniPPI-v2.1 Engine</span><span>CLIENT-SIDE</span></div>',
    "HTML: Update hero model-card header to UniPPI-v2.1 Engine",
)
html = rep(
    html,
    'const out = ["rank,residue,buried_pct,wSB,wHB,wHP,wPi_CatPi,est_ddG_trunc_kcal,tier"];',
    'const out = ["rank,residue,buried_pct,wSB,wHB,wHP,wPi_CatPi_XB,est_ddG_pert_kcal,tier"];',
    "JS: Align Hotspots CSV header with XB & perturbation column",
)

# 4. Rewrite canonical / OG / JSON-LD URLs to the brand-new GitHub Pages URL
html = re.sub(
    r"https://soumyajitmandal06-gif\.github\.io/[A-Za-z0-9._-]+/",
    live_url,
    html,
)
print(f"  [APPLIED]    SEO/OG/Canonical URLs updated to {live_url}")

html_path.write_text(html, encoding="utf-8")

# Verify all JS handlers remain 100% bound
m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
handlers = set(re.findall(r'on(?:click|change|input)="([A-Za-z_]\w*)\(', html))
defined_fns = set(re.findall(r'\bfunction\s+([A-Za-z_]\w*)\s*\(', m_js[-1]))
missing = sorted(handlers - defined_fns)
assert not missing, f"Undefined JS handlers: {missing}"
print(f"[PASS] index.html verified ({len(handlers)} event handlers bound).")

# 5. Generate clean README.md for the new repository
FENCE = chr(96) * 3
readme = "\n".join([
    "# UniPPI-OWL v2.1 (MARVAS) — Universal Multi-Modal Structure-to-Affinity & Interface Auditing Engine",
    "",
    f"**Live Zero-Server Web App:** [{live_url}]({live_url})",
    "",
    "UniPPI-OWL v2.1-MARVAS is a dual-layer (**100% client-side Browser SPA + pure-Python CLI**) post-prediction biophysical auditing engine for **AlphaFold 3**, **Boltz-1/2**, **Chai-1**, and experimental **PDB/mmCIF** complexes. All calculations run locally with zero external server uploads.",
    "",
    "## Supported Biomolecular Modalities (10 Regimes)",
    "- **Protein-Protein (`PPI`)**: Enthalpy-entropy compensated 8-channel solvation + contact-gated PRODIGY hybrid, PAE-weighted `w_ij` IDR damping, Local Interaction Score (`LIS`), and strained helical Proline kink detection.",
    "- **Protein-RNA / DNA & Nucleic Duplexes**: Dedicated 9-channel physics separating phosphate salt bridges (`pKa ~ 1.0`), base-specific vs. ribose `2'-OH` vs. backbone H-bonds, plane-normal nucleobase pi-stacking, `60 deg`-cone cation-pi, and Sequence-Specificity Index (`SSI`).",
    "- **Small Molecules, Drugs & Cofactors**: 3D aromatic ring detector, anionic phosphate/carboxylate & cationic aliphatic amine/guanidinium salt-bridge detectors, directional `>=120 deg` sigma-hole halogen bonds, fractional ligand burial `Q_pocket`, and `LIGAND-LIGAND` support.",
    "- **Catalytic & Structural Metals (`Zn2+`, `Mg2+`, `Ca2+`, `Mn2+`, `Fe`, `Cu`, `Co`, `Ni`)**: Inner-sphere (`<=2.65-2.75 A`) and hydrated outer-sphere (`2.75-4.15 A`) oxygen coordination.",
    "- **Automatic Ternary Molecular-Glue & Metal-Bridge Detection**: Detects composite `A+B-Ligand` wedges and `A+B-Metal` catalytic bridges with burial splits and cooperativity balance.",
    "- **Same-Protein Intra-Chain Domain Slicing**: Custom `--split-chain A:1-236,237-406` or automatic contact-trough monomer splitting.",
    "- **Universal X vs. Y (+ Unlimited Z1..Zn) Studio**: Mix known PDB/UniProt entries with custom lab-synthesized sequences, PDB-numbered mutations, and custom drug SMILES; exports ready-to-run **Boltz-1/2 YAML + 1-click Colab `.sh`**, **AF3 Server JSON**, and **Chai-1 FASTA**.",
    "- **WT vs. Mutant / Apo vs. Holo Ensemble Diff Mode**: Computes `dUIS`, `ddG`, `Kd` fold-change, two-tailed Welch's *t*-test *p*-values, and atomic `LOST`/`GAINED` Contact Rewiring audits.",
    "",
    "## Quick CLI Usage",
    f"{FENCE}bash",
    "# Score any AF3 .zip, Boltz-1/Chai-1 .cif + JSON, or experimental PDB/mmCIF",
    "python3 unippi_cli.py -i benchmarks/5ZC9.cif -o results_5zc9",
    "",
    "# Custom physiological pH, ionic strength, or intra-chain domain split",
    "python3 unippi_cli.py -i benchmarks/5ZC9.cif --split-chain A:1-236,237-406 --ph 7.2 --ionic-strength 150 -o results_split",
    "",
    "# WT vs. Mutant / Apo vs. Holo statistical ensemble diff",
    "python3 unippi_cli.py -i wt_ensemble.zip --diff-b mutant_ensemble.zip -o results_diff",
    f"{FENCE}",
    "",
    "## Validated Regression Benchmarks",
    "| Benchmark | Top Interface | Modality | UIS | dG_uni (kcal/mol) | Pattern |",
    "| :--- | :--- | :--- | :---: | :---: | :---: |",
    "| `fold_vip1_bzip61.zip` (AF3 5-model) | `A-B` | `PROTEIN-PROTEIN` | `0.566` | `-10.99` | `Pattern B` |",
    "| `1URN.cif` (U1A-RNA ITC) | `P-A` | `PROTEIN-RNA` | `0.983` | `-13.61` | `Pattern A` |",
    "| `5ZC9.cif` (eIF4A1+RNA+RocA+ANP+Mg) | `A-B` (+ Glue & Mg-Bridge) | `PROTEIN-RNA` | `0.958` | `-12.88` | `Pattern A` |",
    "| `4WNP.cif` (ULK1+3RJ Kinase) | `D-D_3RJ401` | `PROTEIN-LIGAND` | `0.728` | `-8.41` | `Pattern A` |",
    "| `7DIY.cif` (Nsp14-Nsp10+Zn+Mg) | `A-B` (+ 4x Zn2+) | `PROTEIN-PROTEIN` | `0.847` | `-10.62` | `Pattern A` |",
    "| `5ZC9.cif` (`--split-chain A:1-236,237-406`) | `A_D1+A_D2-A_ANP500` | `TERNARY-GLUE` | `0.907` | `-11.80` | `Pattern A` |",
    "",
    "## Citation & Author",
    "Architected by **Soumyajit Mandal** (Department of Biochemistry, Central University of Rajasthan, India).",
    f"Repository: [{repo_url}]({repo_url}) · License: MIT / BSD-3.",
    "",
])
Path("README.md").write_text(readme, encoding="utf-8")
print(f"[PASS] Wrote clean README.md pointing to {live_url}")
