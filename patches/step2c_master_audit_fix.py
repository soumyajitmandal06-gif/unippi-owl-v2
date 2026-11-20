#!/usr/bin/env python3
"""
Phase 2 - Step 2C (Entity-Safe, Idempotent Master Audit Fix for unippi_cli.py + index.html):
  1. Unique RNA/DNA SHA-256 hashing via NUC_TO_ONE (no "XXXXX" collisions).
  2. Extends is_nucleic_pair to RNA-RNA / DNA-DNA / RNA-DNA duplexes + zeroes helical_prolines.
  3. Fixes multi-chain AF3 chain_pair_iptm indexing via origChainId (immune to HETATM splitting).
  4. Replaces hardcoded "aa" & "A:"/"B:" UI labels (including Burial Profile) with dynamic chain IDs
     and units (aa/nt/lig/ion), adds Modality column, and updates RNA patternCallout & Track 3/5.
  5. Upgrades Contact Inspector filter (PI/CAT/2OH/BASE/BB), Contacts CSV export, 3Dmol.js
     multi-modal rendering, and SVG Contact Network + Burial Profile legends.
  6. Entity-safe HTML prose updates (handles both raw &pi;/&middot; and Unicode π/·) + CLI verification.
"""

import ast
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
    raise AssertionError(f"Neither old nor new state found for [{label}].\nExpected snippet: {old[:140]!r}")


def ensure_regex_replaced(text: str, pattern: str, new: str, label: str) -> str:
    if new in text:
        print(f"  [OK-ALREADY] {label}")
        return text
    m = re.search(pattern, text, flags=re.DOTALL)
    if m:
        print(f"  [APPLIED]    {label}")
        return text[:m.start()] + new + text[m.end():]
    raise AssertionError(f"Neither pattern nor target found for [{label}].\nPattern: {pattern}")


# ==============================================================================
# PART A: IDEMPOTENT VERIFICATION / PATCH OF unippi_cli.py
# ==============================================================================
cli = cli_path.read_text(encoding="utf-8")
print("=== PART A: unippi_cli.py ===")

anchor_n = 'NUCLEIC_CARBON_ATOMS = {"C1\'", "C2\'", "C3\'", "C4\'", "C5\'", "C2", "C4", "C5", "C6", "C7", "C8"}'
nuc_map_line = 'NUC_TO_ONE = {"A": "a", "C": "c", "G": "g", "U": "u", "T": "t", "DA": "a", "DC": "c", "DG": "g", "DT": "t", "ADE": "a", "CYT": "c", "GUA": "g", "URA": "u", "THY": "t"}'
cli = ensure_replaced(cli, anchor_n, anchor_n + "\n" + nuc_map_line, "CLI: NUC_TO_ONE constant")

cli_sub, n_seq = re.subn(
    r'THREE_TO_ONE\.get\((r\.canonical_name),\s*"X"\)',
    r'THREE_TO_ONE.get(\1, NUC_TO_ONE.get(\1, \1))',
    cli,
)
if n_seq > 0:
    cli = cli_sub
    print(f"  [APPLIED]    CLI: NUC_TO_ONE hash fallback ({n_seq} site(s))")
else:
    print("  [OK-ALREADY] CLI: NUC_TO_ONE hash fallback")

target_is_nuc = 'is_nucleic_pair = modality in ("PROTEIN-RNA", "PROTEIN-DNA", "RNA-RNA", "DNA-DNA", "RNA-DNA")'
if target_is_nuc in cli:
    print("  [OK-ALREADY] CLI: is_nucleic_pair duplex support")
else:
    cli_sub, n_nuc = re.subn(r'is_nucleic_pair\s*=\s*modality\s+in\s+\([^)]+\)', target_is_nuc, cli, count=1)
    assert n_nuc == 1, "Could not locate is_nucleic_pair assignment in unippi_cli.py"
    cli = cli_sub
    print("  [APPLIED]    CLI: is_nucleic_pair duplex support")

target_cli_is_nuc = '_is_nuc = _mod in ("PROTEIN-RNA", "PROTEIN-DNA", "RNA-RNA", "DNA-DNA", "RNA-DNA")'
if target_cli_is_nuc in cli:
    print("  [OK-ALREADY] CLI: _is_nuc printout duplex support")
else:
    cli_sub, n_cn = re.subn(r'_is_nuc\s*=\s*_mod\s+in\s+\([^)]+\)', target_cli_is_nuc, cli, count=1)
    if n_cn == 1:
        cli = cli_sub
        print("  [APPLIED]    CLI: _is_nuc printout duplex support")

cli = ensure_replaced(
    cli,
    '"helical_prolines": helical_prolines,',
    '"helical_prolines": 0 if is_nucleic_pair else helical_prolines,',
    "CLI: zero helical_prolines on nucleic pairs",
)

ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] unippi_cli.py verified (ast.parse OK).\n")

# ==============================================================================
# PART B: ENTITY-SAFE IDEMPOTENT PATCH OF index.html
# ==============================================================================
html = html_path.read_text(encoding="utf-8")
print("=== PART B: index.html ===")

# B1. Stale v1.8 metadata, hero H1, model card, scope boxes, footer (entity-tolerant for prose)
html = ensure_regex_replaced(
    html,
    r'content="UniPPI-OWL v1\.8\s*(?:—|&mdash;|-+)\s*Universal Structure-to-Affinity\s*(?:&amp;|&)\s*Atomic Interface Auditing Resource for AlphaFold3 and Experimental Complexes\."',
    'content="UniPPI-OWL v2.0 (MARVAS) — Universal Protein–Protein &amp; Protein–RNA/DNA Structure-to-Affinity &amp; Atomic Interface Auditing Resource."',
    "HTML: meta description",
)
html = ensure_replaced(
    html,
    'content="UniPPI-OWL v1.8 (unippiowl) | Universal Structure-to-Affinity Resource"',
    'content="UniPPI-OWL v2.0 (MARVAS) | Universal Protein–Protein &amp; Protein–RNA Affinity Resource"',
    "HTML: og:title",
)
html = ensure_replaced(
    html,
    '"alternateName": ["unippiowl", "unippi-owl", "UniPPI-OWL v1.8"],',
    '"alternateName": ["unippiowl", "unippi-owl", "UniPPI-OWL v2.0", "UniPPI-v2.0-MARVAS"],',
    "HTML: JSON-LD alternateName",
)
html = ensure_replaced(
    html,
    '<h1 style="font-size:26px;">UniPPI-OWL <span style="font-size:14px; font-weight:600; color:var(--text-muted);">v1.8</span></h1>',
    '<h1 style="font-size:26px;">UniPPI-OWL <span style="font-size:14px; font-weight:600; color:var(--primary);">v2.0-MARVAS</span></h1>',
    "HTML: hero H1 version badge",
)
html = ensure_replaced(
    html,
    '<aside class="model-card" aria-label="UniPPI-v1.8 model summary">',
    '<aside class="model-card" aria-label="UniPPI-v2.0-MARVAS model summary">',
    "HTML: model card aria-label",
)
html = ensure_replaced(
    html,
    '<span>Scoring model</span><strong>Hybrid PRODIGY + 8-Ch Solvation (tanh)</strong>',
    '<span>Scoring models</span><strong>8-Ch PPI + 9-Ch Protein–RNA/DNA (tanh)</strong>',
    "HTML: model card scoring models",
)
html = ensure_replaced(
    html,
    '<strong>What UniPPI-OWL v1.8 Evaluates:</strong>',
    '<strong>What UniPPI-OWL v2.0-MARVAS Evaluates:</strong>',
    "HTML: scope evaluates header",
)
html = ensure_regex_replaced(
    html,
    r'<li>Confidence-weighted 8-channel atomic interactions.*?</li>',
    '<li>Confidence-weighted multi-modal atomic physics: 8-channel Protein–Protein solvation + dedicated 9-channel Protein–RNA/DNA recognition (phosphate salt bridges, base-specific vs. ribose 2\'-OH vs. backbone H-bonds, plane-normal nucleobase π-stacking, 60°-cone cation–π, and Sequence-Specificity Index SSI).</li>',
    "HTML: scope bullet 1",
)
html = ensure_replaced(
    html,
    '<strong>UniPPI-OWL v1.8 Scope Boundaries (Does Not):</strong>',
    '<strong>UniPPI-OWL v2.0-MARVAS Scope Boundaries (Does Not):</strong>',
    "HTML: scope boundaries header",
)
html = ensure_regex_replaced(
    html,
    r'<strong>UniPPI-OWL v1\.8</strong>\s*(?:·|&middot;)\s*Universal Structure-to-Affinity\s*(?:&amp;|&)\s*Atomic Interface Auditing Resource<br>',
    '<strong>UniPPI-OWL v2.0-MARVAS</strong> · Universal Protein–Protein &amp; Protein–RNA/DNA Structure-to-Affinity Auditing Resource<br>',
    "HTML: footer title",
)
html = ensure_replaced(
    html,
    'Frozen Scientific Model: <code>UniPPI-v1.8</code><br>',
    'Frozen Scientific Model: <code>UniPPI-v2.0-MARVAS</code><br>',
    "HTML: footer model version",
)

# B2. NUC_TO_ONE constant
old_t2o = '    const ONE_TO_THREE = Object.fromEntries(Object.entries(THREE_TO_ONE).map(([k,v]) => [v,k]));'
new_t2o = old_t2o + '\n    const NUC_TO_ONE = { A:"a", C:"c", G:"g", U:"u", T:"t", DA:"a", DC:"c", DG:"g", DT:"t", ADE:"a", CYT:"c", GUA:"g", URA:"u", THY:"t" };'
html = ensure_replaced(html, old_t2o, new_t2o, "JS: NUC_TO_ONE constant")

# B3. Nucleotide SHA-256 hashing + dynamic chain units + duplex isNucleicPair + helical_prolines
html = ensure_replaced(
    html,
    '''      const seqA = resA.map(r => THREE_TO_ONE[r.canonicalName] || "X").join("");
      const seqB = resB.map(r => THREE_TO_ONE[r.canonicalName] || "X").join("");''',
    '''      const seqA = resA.map(r => THREE_TO_ONE[r.canonicalName] || NUC_TO_ONE[r.canonicalName] || r.canonicalName).join("");
      const seqB = resB.map(r => THREE_TO_ONE[r.canonicalName] || NUC_TO_ONE[r.canonicalName] || r.canonicalName).join("");
      const getChainUnit = (rlist) => {
        const ts = new Set(rlist.map(r => r.entityType));
        if (ts.has("PROTEIN")) return "aa";
        if (ts.has("NUCLEIC")) return "nt";
        if (ts.has("LIGAND")) return "lig";
        if (ts.has("ION")) return "ion";
        return "res";
      };
      const unitA = getChainUnit(resA), unitB = getChainUnit(resB);''',
    "JS: unique RNA/DNA seq hash & chain unit detection",
)

html = ensure_replaced(
    html,
    '      const isNucleicPair = (modality === "PROTEIN-RNA" || modality === "PROTEIN-DNA");',
    '      const isNucleicPair = (modality === "PROTEIN-RNA" || modality === "PROTEIN-DNA" || modality === "RNA-RNA" || modality === "DNA-DNA" || modality === "RNA-DNA");',
    "JS: isNucleicPair duplex support",
)

html = ensure_replaced(
    html,
    '        complex_id: complexId, chain_pair: `${caId}-${cbId}`, modality, sha256_id, mode,',
    '        complex_id: complexId, chain_pair: `${caId}-${cbId}`, chainAId: caId, chainBId: cbId, unitA, unitB, modality, sha256_id, mode,',
    "JS: return chainAId/chainBId/unitA/unitB",
)

html = ensure_replaced(
    html,
    '        helical_prolines: helPro, f_apolar: Number(apFrac.toFixed(3)),',
    '        helical_prolines: effHelPro, f_apolar: Number(apFrac.toFixed(3)),',
    "JS: helical_prolines zero on nucleic pairs",
)

# B4. AF3 chain_pair_iptm alignment by origChainId
html = ensure_replaced(
    html,
    '''        const isBinary = (cIds.length === 2);

        for (let i = 0; i < cIds.length; i++) {
          for (let j = i + 1; j < cIds.length; j++) {
            const ca = cIds[i], cb = cIds[j];
            let pairSasaComp = sasaCompFull;
            if (!isBinary) {
              let hasProx = false;
              for (const ra of chains[ca]) {
                for (const rb of chains[cb]) {
                  if (dist3(ra.centroid, rb.centroid) <= 16.0) { hasProx = true; break; }
                }
                if (hasProx) break;
              }
              if (!hasProx) continue;
              pairSasaComp = computeResidueSASA([...chains[ca], ...chains[cb]]);
            }
            const pairIptm = (Array.isArray(chainPairIptm) && chainPairIptm[i] && typeof chainPairIptm[i][j] === "number") ? chainPairIptm[i][j] : iptm;''',
    '''        const isBinary = (cIds.length === 2);
        const origChainOrder = [...new Set(ordered.map(r => r.origChainId || r.chainId))];

        for (let i = 0; i < cIds.length; i++) {
          for (let j = i + 1; j < cIds.length; j++) {
            const ca = cIds[i], cb = cIds[j];
            let pairSasaComp = sasaCompFull;
            if (!isBinary) {
              let hasProx = false;
              for (const ra of chains[ca]) {
                for (const rb of chains[cb]) {
                  if (dist3(ra.centroid, rb.centroid) <= 16.0) { hasProx = true; break; }
                }
                if (hasProx) break;
              }
              if (!hasProx) continue;
              pairSasaComp = computeResidueSASA([...chains[ca], ...chains[cb]]);
            }
            const oi = origChainOrder.indexOf(chains[ca][0]?.origChainId || ca);
            const oj = origChainOrder.indexOf(chains[cb][0]?.origChainId || cb);
            const pairIptm = (Array.isArray(chainPairIptm) && oi >= 0 && oj >= 0 && oi !== oj && chainPairIptm[oi] && typeof chainPairIptm[oi][oj] === "number") ? chainPairIptm[oi][oj] : iptm;''',
    "JS: origChainId-aligned chain_pair_iptm lookup",
)

html = ensure_replaced(
    html,
    'cIds.map(c=>`${c}:${chains[c].length}aa`).join(", ")',
    'cIds.map(c=>`${c}:${chains[c].length}`).join(", ")',
    "JS: neutral queue chain length unit",
)

# B5. Leaderboard Modality column + dynamic chain labels/units
html = ensure_replaced(
    html,
    '''              <th>Complex ID</th>
              <th>Chains (Lengths)</th>
              <th>Binding Loci</th>''',
    '''              <th>Complex ID</th>
              <th>Chains (Lengths)</th>
              <th>Modality</th>
              <th>Binding Loci</th>''',
    "HTML: Leaderboard Modality header",
)

html = ensure_replaced(
    html,
    '''        const lociStr = (r.lociA && r.lociA !== "None") ? `A:${r.lociA} ↔ B:${r.lociB}` : "None";
        tr.innerHTML = `
          <td>#${idx+1}</td>
          <td><strong>${r.complex_id}</strong></td>
          <td>${r.chain_pair} (${fLenA}+${fLenB} aa)</td>
          <td>${lociStr} (${nConf} pairs)</td>''',
    '''        const cAName = r.chainAId || (r.chain_pair || "A-B").split("-")[0];
        const cBName = r.chainBId || (r.chain_pair || "A-B").split("-").slice(1).join("-") || "B";
        const uA = r.unitA || "aa", uB = r.unitB || "aa";
        const lociStr = (r.lociA && r.lociA !== "None") ? `${cAName}:${r.lociA} ↔ ${cBName}:${r.lociB}` : "None";
        tr.innerHTML = `
          <td>#${idx+1}</td>
          <td><strong>${r.complex_id}</strong></td>
          <td>${r.chain_pair} (${fLenA} ${uA} + ${fLenB} ${uB})</td>
          <td><span class="badge-call" style="background:var(--bg-subtle);color:var(--primary);border:1px solid var(--border-strong);">${r.modality || "PROTEIN-PROTEIN"}</span></td>
          <td>${lociStr} (${nConf} pairs)</td>''',
    "JS: Leaderboard Modality badge & dynamic chain units",
)

# B6. renderActiveComplex dynamic chain labels, units, SSI descriptor, RNA callout & Track 3/5
html = ensure_replaced(
    html,
    '      document.getElementById("activeComplexTitle").textContent = `${r.complex_id} (Chains ${r.chain_pair}: ${fLenA} aa + ${fLenB} aa)`;',
    '''      const cAName = r.chainAId || (r.chain_pair || "A-B").split("-")[0];
      const cBName = r.chainBId || (r.chain_pair || "A-B").split("-").slice(1).join("-") || "B";
      const uA = r.unitA || "aa", uB = r.unitB || "aa";
      const isNucUI = ["PROTEIN-RNA","PROTEIN-DNA","RNA-RNA","DNA-DNA","RNA-DNA"].includes(r.modality);
      document.getElementById("activeComplexTitle").textContent = `${r.complex_id} (${cAName}: ${fLenA} ${uA} ↔ ${cBName}: ${fLenB} ${uB} · ${r.modality || "PROTEIN-PROTEIN"})`;''',
    "JS: activeComplexTitle dynamic units & modality",
)

html = ensure_replaced(
    html,
    '''      document.getElementById("stableDescriptorsGrid").innerHTML = `
        <div class="descriptor"><span>1. Apolar Burial Fraction (f<sub>apolar</sub>)</span><strong>${r.f_apolar.toFixed(3)} (${r.hydrophobic_contacts} HP contacts)</strong></div>
        <div class="descriptor"><span>2. Electrostatic Density / 1000 Å²</span><strong>${r.elec_density_1000A2.toFixed(2)} (${r.salt_bridges} SB, ${r.h_bonds} HB)</strong></div>
        <div class="descriptor"><span>3. Mean Joint Interface Confidence (w̄<sub>intf</sub>)</span><strong>${r.mean_w_intf.toFixed(3)} (wBSA: ${r.wBSA_A2} Å²)</strong></div>
        <div class="descriptor"><span>4. Ordered Binding Loci (w ≥ 0.25)</span><strong>A:${r.lociA || "—"} ↔ B:${r.lociB || "—"} (${nConf} pairs, ${r.intf_residues} res)</strong></div>
      `;''',
    '''      document.getElementById("stableDescriptorsGrid").innerHTML = `
        <div class="descriptor"><span>${isNucUI ? "1. RNA/DNA Specificity Index (SSI)" : "1. Apolar Burial Fraction (f<sub>apolar</sub>)"}</span><strong>${isNucUI ? `${(r.rna_specificity_idx||0).toFixed(2)} (π:${r.pi_stacking||0}, Cat-π:${r.cation_pi||0})` : `${r.f_apolar.toFixed(3)} (${r.hydrophobic_contacts} HP contacts)`}</strong></div>
        <div class="descriptor"><span>2. Electrostatic Density / 1000 Å²</span><strong>${r.elec_density_1000A2.toFixed(2)} (${r.salt_bridges} SB, ${r.h_bonds} HB${isNucUI ? ` · 2'-OH:${r.hb_2oh||0}` : ""})</strong></div>
        <div class="descriptor"><span>3. Mean Joint Interface Confidence (w̄<sub>intf</sub>)</span><strong>${r.mean_w_intf.toFixed(3)} (wBSA: ${r.wBSA_A2} Å²)</strong></div>
        <div class="descriptor"><span>4. Ordered Binding Loci (w ≥ 0.25)</span><strong>${cAName}:${r.lociA || "—"} ↔ ${cBName}:${r.lociB || "—"} (${nConf} pairs, ${r.intf_residues} res)</strong></div>
      `;''',
    "JS: stableDescriptorsGrid RNA SSI & dynamic loci",
)

html = ensure_replaced(
    html,
    '''              Parsed 3D coordinates: Chain A = <strong>${fLenA} aa</strong>, Chain B = <strong>${fLenB} aa</strong> (${r.pct_disordered ?? 0}% disordered residues at pLDDT &lt; 70).
              Detected <strong>${nConf} confident inter-chain contact pairs</strong> (${nRaw} total raw pairs) across <strong>${r.intf_residues} unique contacting interface residues</strong> within <strong>A:${r.lociA || "—"} &leftrightarrow; B:${r.lociB || "—"}</strong>''',
    '''              Parsed 3D coordinates: Chain ${cAName} = <strong>${fLenA} ${uA}</strong>, Chain ${cBName} = <strong>${fLenB} ${uB}</strong> (${r.pct_disordered ?? 0}% flexible/disordered residues).
              Detected <strong>${nConf} confident inter-chain contact pairs</strong> (${nRaw} total raw pairs) across <strong>${r.intf_residues} unique contacting interface residues</strong> within <strong>${cAName}:${r.lociA || "—"} &leftrightarrow; ${cBName}:${r.lociB || "—"}</strong>''',
    "JS: Step 1 derivation dynamic chain IDs & units",
)

html = html.replace('(r.modality==="PROTEIN-RNA"||r.modality==="PROTEIN-DNA")', 'isNucUI')
html = html.replace('(r.modality === "PROTEIN-RNA" || r.modality === "PROTEIN-DNA")', 'isNucUI')

html = ensure_replaced(
    html,
    '''      document.getElementById("patternCallout").innerHTML = `
        <strong>${r.pattern_id} Classification:</strong> ${r.pattern_desc}.
        <br><span style="color:var(--text-muted); font-size:13px;">Note: ΔG<sub>uni</sub> (${r.dG_uni_kcal.toFixed(2)} kcal/mol) weights published PRODIGY (${Math.round(r.w_prodigy_blend*100)}%: ${r.dG_prodigy_kcal.toFixed(2)} kcal/mol) and enthalpy–entropy compensated 8-channel solvation (${Math.round(r.w_phys_blend*100)}%: ${r.dG_phys_kcal.toFixed(2)} kcal/mol).</span>
      `;''',
    '''      document.getElementById("patternCallout").innerHTML = `
        <strong>${r.pattern_id} Classification:</strong> ${r.pattern_desc}.
        <br><span style="color:var(--text-muted); font-size:13px;">${isNucUI
          ? `Note: Dedicated 9-channel ${r.modality} physics active (Protein–Protein PRODIGY disabled; SSI = ${(r.rna_specificity_idx||0).toFixed(2)}, Base H-Bonds = ${r.hb_base||0}, Ribose 2'-OH = ${r.hb_2oh||0}, Base π = ${r.pi_stacking||0}, Cation-π = ${r.cation_pi||0}).`
          : `Note: ΔG<sub>uni</sub> (${r.dG_uni_kcal.toFixed(2)} kcal/mol) weights published PRODIGY (${Math.round(r.w_prodigy_blend*100)}%: ${r.dG_prodigy_kcal.toFixed(2)} kcal/mol) and enthalpy–entropy compensated 8-channel solvation (${Math.round(r.w_phys_blend*100)}%: ${r.dG_phys_kcal.toFixed(2)} kcal/mol).`
        }</span>
      `;''',
    "JS: patternCallout RNA-aware explanation",
)

html = ensure_replaced(
    html,
    '<td>ΔG<sub>uni</sub> = ${r.dG_uni_kcal} kcal/mol · SB = ${r.salt_bridges} · HB = ${r.h_bonds} · HP = ${r.hydrophobic_contacts} · π = ${r.pi_stacking || 0} · SS = ${r.disulfides || 0}</td>',
    '<td>ΔG<sub>uni</sub> = ${r.dG_uni_kcal} kcal/mol · SB = ${r.salt_bridges} · HB = ${r.h_bonds}${isNucUI ? ` (Base:${r.hb_base||0}, 2\'-OH:${r.hb_2oh||0}, BB:${r.hb_backbone||0})` : ""} · HP = ${r.hydrophobic_contacts} · π = ${r.pi_stacking || 0}${isNucUI ? ` · Cat-π = ${r.cation_pi || 0} · SSI = ${(r.rna_specificity_idx||0).toFixed(2)}` : ""} · SS = ${r.disulfides || 0}</td>',
    "JS: Track 3 RNA breakdown",
)

html = ensure_replaced(
    html,
    '<td>Modality = <strong>${r.modality || \'PROTEIN-PROTEIN\'}</strong> · Chain A = <strong>${fLenA} aa</strong>, Chain B = <strong>${fLenB} aa</strong> · Ordered Binding Loci = <strong>A:${r.lociA || "—"} &leftrightarrow; B:${r.lociB || "—"} (${r.intf_residues} contacting res)</strong> · SHA-256 = <code>${r.sha256_id || "—"}</code></td>',
    '<td>Modality = <strong>${r.modality || \'PROTEIN-PROTEIN\'}</strong> · Chain ${cAName} = <strong>${fLenA} ${uA}</strong>, Chain ${cBName} = <strong>${fLenB} ${uB}</strong> · Ordered Binding Loci = <strong>${cAName}:${r.lociA || "—"} &leftrightarrow; ${cBName}:${r.lociB || "—"} (${r.intf_residues} contacting res)</strong> · SHA-256 = <code>${r.sha256_id || "—"}</code></td>',
    "JS: Track 5 dynamic chain IDs & units",
)

# B7. Multi-modal 3Dmol.js rendering
html = ensure_replaced(
    html,
    '''      const [rawCA, rawCB] = r.chain_pair.split("-");
      const cA = r.pdbChainA || (rawCA.length > 1 ? "A" : rawCA);
      const cB = r.pdbChainB || (rawCB.length > 1 ? "B" : rawCB);
      const dimDisordered = document.getElementById("dimDisorderedCheck")?.checked;

      if (dimDisordered && r.mode !== "MODE_C_EXP_CRYSTAL") {
        VIEWER_3D.setStyle({ chain: cA, b: { gte: 70 } }, { cartoon: { color: "#38bdf8" } });
        VIEWER_3D.setStyle({ chain: cA, b: { lt: 70 } }, { cartoon: { color: "#334155", opacity: 0.22 } });
        VIEWER_3D.setStyle({ chain: cB, b: { gte: 70 } }, { cartoon: { color: "#34d399" } });
        VIEWER_3D.setStyle({ chain: cB, b: { lt: 70 } }, { cartoon: { color: "#334155", opacity: 0.22 } });
      } else {
        VIEWER_3D.setStyle({ chain: cA }, { cartoon: { color: "#38bdf8" } });
        VIEWER_3D.setStyle({ chain: cB }, { cartoon: { color: "#34d399" } });
      }

      if (r.contacts) {
        for (const c of r.contacts) {
          if (c.weight < 0.25) continue;
          VIEWER_3D.addStyle({ chain: cA, resi: c.resSeqA }, { stick: { radius: 0.15, color: "#38bdf8" } });
          VIEWER_3D.addStyle({ chain: cB, resi: c.resSeqB }, { stick: { radius: 0.15, color: "#34d399" } });
          if (c.sb > 0 && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.11, color: "#f43f5e", dashed: true });
          } else if (c.hb > 0 && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.09, color: "#facc15", dashed: true });
          }
        }
      }''',
    '''      const cA = r.pdbChainA || (r.chainAId && r.chainAId.length === 1 ? r.chainAId : "A");
      const cB = r.pdbChainB || (r.chainBId && r.chainBId.length === 1 ? r.chainBId : "B");
      const colA = r.unitA === "nt" ? "#f59e0b" : "#38bdf8";
      const colB = r.unitB === "nt" ? "#f59e0b" : "#34d399";
      const dimDisordered = document.getElementById("dimDisorderedCheck")?.checked;

      // Subtle bystander context + visible HETATM drugs (magenta) & catalytic metals (lime)
      VIEWER_3D.setStyle({}, { cartoon: { color: "#475569", opacity: 0.30 } });
      VIEWER_3D.setStyle({ hetflag: true }, { stick: { radius: 0.22, colorscheme: "magentaCarbon" }, sphere: { scale: 0.26, colorscheme: "magentaCarbon" } });
      VIEWER_3D.setStyle({ hetflag: true, elem: ["MG","ZN","MN","CA","FE","CU","CO","NI"] }, { sphere: { scale: 0.55, color: "#84cc16" } });

      const styleActiveChain = (cid, col, unit) => {
        if (unit === "lig") {
          VIEWER_3D.setStyle({ chain: cid }, { stick: { radius: 0.26, colorscheme: "magentaCarbon" }, sphere: { scale: 0.28, colorscheme: "magentaCarbon" } });
        } else if (unit === "ion") {
          VIEWER_3D.setStyle({ chain: cid }, { sphere: { scale: 0.65, color: "#84cc16" } });
        } else if (unit === "nt") {
          VIEWER_3D.setStyle({ chain: cid }, { cartoon: { color: col }, stick: { radius: 0.14, color: col } });
        } else if (dimDisordered && r.mode !== "MODE_C_EXP_CRYSTAL") {
          VIEWER_3D.setStyle({ chain: cid, b: { gte: 70 } }, { cartoon: { color: col } });
          VIEWER_3D.setStyle({ chain: cid, b: { lt: 70 } }, { cartoon: { color: "#334155", opacity: 0.22 } });
        } else {
          VIEWER_3D.setStyle({ chain: cid }, { cartoon: { color: col } });
        }
      };
      styleActiveChain(cA, colA, r.unitA || "aa");
      styleActiveChain(cB, colB, r.unitB || "aa");

      if (r.contacts) {
        for (const c of r.contacts) {
          if (c.weight < 0.25) continue;
          if ((r.unitA || "aa") !== "lig" && (r.unitA || "aa") !== "ion") VIEWER_3D.addStyle({ chain: cA, resi: c.resSeqA }, { stick: { radius: 0.15, color: colA } });
          if ((r.unitB || "aa") !== "lig" && (r.unitB || "aa") !== "ion") VIEWER_3D.addStyle({ chain: cB, resi: c.resSeqB }, { stick: { radius: 0.15, color: colB } });
          if (c.sb > 0 && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.11, color: "#f43f5e", dashed: true });
          } else if ((c.pi > 0 || c.catpi > 0) && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.10, color: "#a855f7", dashed: true });
          } else if (c.hb > 0 && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.09, color: "#facc15", dashed: true });
          }
        }
      }''',
    "JS: 3Dmol.js multi-modal rendering",
)

html = ensure_regex_replaced(
    html,
    r'<strong>3D Legend:</strong>\s*Chain A \(Sky Blue\).*?H-Bonds \(Yellow dashed\)\.',
    '<strong>3D Legend:</strong> Protein (Sky Blue/Emerald) · RNA/DNA (Amber) · Drugs/Ligands (Magenta Sticks) · Catalytic Metals (Lime Spheres) · Salt Bridges (Red) · π/Cation-π (Purple) · H-Bonds (Yellow).',
    "HTML: 3D viewer multi-modal legend",
)

# B8. SVG Contact Network & Burial Profile (purple pi/cat-pi + dynamic burial unit)
html = ensure_replaced(
    html,
    '''      const [cA, cB] = (r.chain_pair || "A-B").split("-");
      const confC = (r.contacts || []).filter(c => c.weight >= 0.25);
      const activeC = confC.length ? confC : (r.contacts || []);
      const conA = activeC.map(c => c.resSeqA || 1);
      const conB = activeC.map(c => c.resSeqB || 1);
      const minSeqA = conA.length ? Math.min(...conA) : 1, maxSeqA = conA.length ? Math.max(minSeqA + 1, ...conA) : (r.lenA || 50);
      const minSeqB = conB.length ? Math.min(...conB) : 1, maxSeqB = conB.length ? Math.max(minSeqB + 1, ...conB) : (r.lenB || 50);

      // Layer pure vdW contacts underneath strong specific contacts (SB/HB/HP/pi)
      const rankType = c => c.sb > 0 ? 3 : (c.hb > 0 ? 2 : ((c.hp > 0 || c.pi > 0) ? 1 : 0));
      const sortedC = activeC.slice().sort((a, b) => rankType(a) - rankType(b));
      let linesSvg = "";
      for (const c of sortedC) {
        const normA = clamp(((c.resSeqA || minSeqA) - minSeqA) / Math.max(1, maxSeqA - minSeqA), 0.0, 1.0);
        const normB = clamp(((c.resSeqB || minSeqB) - minSeqB) / Math.max(1, maxSeqB - minSeqB), 0.0, 1.0);
        const yA = 35 + normA * (h - 70);
        const yB = 35 + normB * (h - 70);
        const rt = rankType(c);
        const color = rt === 3 ? "#e11d48" : (rt === 2 ? "#d97706" : (rt === 1 ? "#0284c7" : "#cbd5e1"));
        const strokeW = rt === 0 ? 1.0 : Math.max(1.2, 2.8 * (c.weight || 1.0));
        const opacity = rt === 0 ? "0.45" : "0.88";
        const tag = rt === 3 ? "Salt Bridge" : (rt === 2 ? "H-Bond" : (rt === 1 ? "Hydrophobic/π" : "vdW"));''',
    '''      const cA = r.chainAId || (r.chain_pair || "A-B").split("-")[0];
      const cB = r.chainBId || (r.chain_pair || "A-B").split("-").slice(1).join("-") || "B";
      const confC = (r.contacts || []).filter(c => c.weight >= 0.25);
      const activeC = confC.length ? confC : (r.contacts || []);
      const conA = activeC.map(c => c.resSeqA || 1);
      const conB = activeC.map(c => c.resSeqB || 1);
      const minSeqA = conA.length ? Math.min(...conA) : 1, maxSeqA = conA.length ? Math.max(minSeqA + 1, ...conA) : (r.lenA || 50);
      const minSeqB = conB.length ? Math.min(...conB) : 1, maxSeqB = conB.length ? Math.max(minSeqB + 1, ...conB) : (r.lenB || 50);

      // Layer pure vdW underneath HP, HB, pi/Cat-pi, and SB
      const rankType = c => c.sb > 0 ? 4 : ((c.pi > 0 || c.catpi > 0) ? 3 : (c.hb > 0 ? 2 : (c.hp > 0 ? 1 : 0)));
      const sortedC = activeC.slice().sort((a, b) => rankType(a) - rankType(b));
      let linesSvg = "";
      for (const c of sortedC) {
        const normA = clamp(((c.resSeqA || minSeqA) - minSeqA) / Math.max(1, maxSeqA - minSeqA), 0.0, 1.0);
        const normB = clamp(((c.resSeqB || minSeqB) - minSeqB) / Math.max(1, maxSeqB - minSeqB), 0.0, 1.0);
        const yA = 35 + normA * (h - 70);
        const yB = 35 + normB * (h - 70);
        const rt = rankType(c);
        const color = rt === 4 ? "#e11d48" : (rt === 3 ? "#9333ea" : (rt === 2 ? "#d97706" : (rt === 1 ? "#0284c7" : "#cbd5e1")));
        const strokeW = rt === 0 ? 1.0 : Math.max(1.2, 2.8 * (c.weight || 1.0));
        const opacity = rt === 0 ? "0.45" : "0.90";
        const tag = rt === 4 ? "Salt Bridge" : (rt === 3 ? (c.catpi > 0 ? "Cation-π" : "π-Stack") : (rt === 2 ? "H-Bond" : (rt === 1 ? "Hydrophobic" : "vdW")));''',
    "JS: SVG Contact Network purple pi/Cat-pi edges",
)

html = ensure_replaced(
    html,
    'Red: Salt Bridge · Amber: H-Bond · Blue: Hydrophobic/π · Gray: vdW',
    'Red: Salt Bridge · Purple: π/Cat-π · Amber: H-Bond · Blue: Hydrophobic · Gray: vdW',
    "JS: SVG Contact Network legend",
)

html = ensure_replaced(
    html,
    'windowTag = `Ordered Binding Locus (${viewRes[0].resSeq}–${viewRes[viewRes.length-1].resSeq} of ${fullChainRes.length} aa)`;',
    'windowTag = `Ordered Binding Locus (${viewRes[0].resSeq}–${viewRes[viewRes.length-1].resSeq} of ${fullChainRes.length} ${BURIAL_CHAIN_IDX === 0 ? (r.unitA || "aa") : (r.unitB || "aa")})`;',
    "JS: Burial SVG dynamic unit (aa/nt)",
)

html = ensure_replaced(
    html,
    '''        if (ci.sb > 0) { col = "#e11d48"; chTag = "S"; }
        else if (ci.hb > 0) { col = "#d97706"; chTag = "H"; }
        else if (ci.hp > 0 || ci.pi > 0) { col = "#0284c7"; chTag = "A"; }''',
    '''        if (ci.sb > 0) { col = "#e11d48"; chTag = "S"; }
        else if (ci.pi > 0) { col = "#9333ea"; chTag = "π"; }
        else if (ci.hb > 0) { col = "#d97706"; chTag = "H"; }
        else if (ci.hp > 0) { col = "#0284c7"; chTag = "A"; }''',
    "JS: Burial SVG purple pi bar",
)

html = ensure_replaced(
    html,
    'Red (S): Salt Bridge · Amber (H): H-Bond · Blue (A): Apolar/π · Gray: Exposed/vdW',
    'Red (S): Salt Bridge · Purple (π): π/Cat-π · Amber (H): H-Bond · Blue (A): Apolar · Gray: Exposed/vdW',
    "JS: Burial SVG legend",
)

# B9. Contact Inspector filter (PI/CAT/2OH/BASE/BB) & Contacts CSV export
html = ensure_replaced(
    html,
    'const filt = allC.filter(c => (!hideLow || c.weight >= 0.25) && (!q || c.resA.toUpperCase().includes(q) || c.resB.toUpperCase().includes(q) || (q==="SB" && c.sb>0) || (q==="HB" && c.hb>0) || (q==="HP" && c.hp>0)));',
    'const filt = allC.filter(c => (!hideLow || c.weight >= 0.25) && (!q || c.resA.toUpperCase().includes(q) || c.resB.toUpperCase().includes(q) || (q==="SB" && c.sb>0) || (q==="HB" && c.hb>0) || (q==="HP" && c.hp>0) || ((q==="PI"||q==="π") && c.pi>0) || (q.includes("CAT") && c.catpi>0) || ((q==="2OH"||q==="2\'-OH") && c.hb2oh>0) || (q==="BASE" && c.hbBase>0) || (q==="BB" && c.hbBb>0)));',
    "JS: Contact Inspector RNA filter",
)

html = ensure_replaced(
    html,
    'const cols = ["resA","resB","minDist","pae","weight","sb","hb","hp","pi","ss","clash"];',
    'const cols = ["resA","resB","minDist","pae","weight","sb","hb","hbBase","hb2oh","hbBb","hp","pi","catpi","ss","clash"];',
    "JS: Contacts CSV RNA columns",
)

html_path.write_text(html, encoding="utf-8")
print("[PASS] index.html verified and saved cleanly.")

# Optional Node.js syntax check
m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
assert m_js, "No inline <script> block found in index.html"
if shutil.which("node"):
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
        tf.write(m_js[-1])
        tf_path = tf.name
    res_node = subprocess.run(["node", "--check", tf_path], capture_output=True, text=True)
    Path(tf_path).unlink(missing_ok=True)
    if res_node.returncode != 0:
        raise RuntimeError(f"JavaScript syntax error in index.html:\n{res_node.stderr}")
    print("[PASS] Verified index.html inline JS syntax via `node --check`.")

# ==============================================================================
# PART C: FULL CLI BENCHMARK & SHA-256 UNIQUENESS AUDIT
# ==============================================================================
for title, inp, outp in [
    ("STEP 2C ZERO-REGRESSION GUARD — fold_vip1_bzip61.zip", "benchmarks/fold_vip1_bzip61.zip", "/tmp/v2c_vip1"),
    ("STEP 2C VERIFY — 5ZC9.cif (eIF4A1 + Polypurine RNA)", "benchmarks/5ZC9.cif", "/tmp/v2c_5zc9"),
    ("STEP 2C VERIFY — 1URN.cif (U1A RRM + Hairpin RNA)", "benchmarks/1URN.cif", "/tmp/v2c_1urn"),
]:
    print("\n" + "=" * 95)
    print(title)
    print("=" * 95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", inp, "-o", outp], check=True)

r5 = json.loads(Path("/tmp/v2c_5zc9/unippi_full_report.json").read_text(encoding="utf-8"))[0]
r1 = json.loads(Path("/tmp/v2c_1urn/unippi_full_report.json").read_text(encoding="utf-8"))[0]
print("\n=== STEP 2C INVARIANT CHECKS ===")
print(f"  5ZC9 A-B: modality={r5['modality']} | hel_pro={r5['helical_prolines']} | sha256={r5['sha256_id']} | UIS={r5['UIS']} | dG={r5['dG_uni_kcal']}")
print(f"  1URN P-A: modality={r1['modality']} | hel_pro={r1['helical_prolines']} | sha256={r1['sha256_id']} | UIS={r1['UIS']} | dG={r1['dG_uni_kcal']}")
assert r5["helical_prolines"] == 0 and r1["helical_prolines"] == 0, "helical_prolines must be 0 on nucleic pairs!"
assert r5["sha256_id"] != r1["sha256_id"], "RNA SHA-256 fingerprints collided!"
print("[ALL STEP 2C INVARIANTS PASSED]")
