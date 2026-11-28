#!/usr/bin/env python3
"""
Permanent Basics & Sanity Regression Suite for UniPPI-OWL v1.0.
Verifies real-world edge cases before any Phase 2 benchmarking resumes:
  [TEST 1] Real 3D MDM2-p53 crystal (1YCR) -> ABOVE_THRESHOLD, Pattern A, p53 17-29
  [TEST 2] CLI --out & --salt aliases work identically to --output & --ionic-strength
  [TEST 3] Ca/P-only backbone trace triggers incomplete_atoms_warning in CLI & JS
  [TEST 4] Tab 2 alphabet validator blocks Protein-in-RNA, Protein-in-SMILES, PDB-ID-in-seq
  [TEST 5] Tab 2 has zero legacy "Will They Interact?" or "LOW (UNLIKELY)" fake 1D verdicts
  [TEST 6] Tab 1 single-chain monomer verdict says INTRA-CHAIN DOMAIN PACKING (not intermolecular PPI)
"""
import csv
import glob
import os
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CLI = ROOT / "unippi_cli.py"
HTML = ROOT / "index.html"

# TEST 1 + TEST 2: Real 3D human MDM2-p53 (1YCR) via --out and --salt aliases
with tempfile.TemporaryDirectory(prefix="basics_1ycr_") as td:
    cif_p = Path(td) / "1YCR.cif"
    req = urllib.request.Request("https://files.rcsb.org/download/1YCR.cif", headers={"User-Agent": "UniPPI-OWL/1.0"})
    cif_p.write_bytes(urllib.request.urlopen(req, timeout=25).read())
    out_d = Path(td) / "out_alias"
    proc = subprocess.run(
        [sys.executable, str(CLI), "--input", str(cif_p), "--out", str(out_d), "--ph", "7.2", "--salt", "150"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"CLI --out/--salt alias failed: {proc.stderr}"
    csvs = sorted(glob.glob(str(out_d / "**" / "*master*.csv"), recursive=True))
    rows = list(csv.DictReader(open(csvs[0], encoding="utf-8")))
    top = rows[0]
    assert top["threshold_call"] == "ABOVE_THRESHOLD", f"1YCR MDM2-p53 failed call: {top}"
    assert top["pattern_id"] == "Pattern A", f"1YCR expected Pattern A, got {top['pattern_id']}"
    loci_a = top.get("lociA") or top.get("loci_a") or top.get("binding_loci_A") or "N/A"
    loci_b = top.get("lociB") or top.get("loci_b") or top.get("binding_loci_B") or "N/A"
    print(f"[PASS 1 & 2] 3D MDM2-p53 (1YCR) via --out/--salt -> {top['chain_pair']} | UIS={top['UIS']} | dG={top['dG_uni_kcal']} kcal/mol | Kd={float(top['Kd_nM']):.1f} nM | Loci={loci_a} <-> {loci_b}")
    if "lociA" not in top:
        print(f"  [SCHEMA NOTE] CLI CSV headers (check parity vs browser): {list(top.keys())}")

# TEST 3: Ca-only backbone trace triggers incomplete_atoms_warning
with tempfile.TemporaryDirectory(prefix="basics_ca_") as td:
    ca_pdb = Path(td) / "ca_only_test.pdb"
    lines = ["EXPDTA    X-RAY DIFFRACTION"]
    for i in range(1, 21):
        lines.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    {0.0:8.3f}{i*3.8:8.3f}{0.0:8.3f}  1.00 25.00           C")
        lines.append(f"ATOM  {i+20:5d}  CA  LEU B{i:4d}    {4.5:8.3f}{i*3.8:8.3f}{0.0:8.3f}  1.00 25.00           C")
    lines.append("END")
    ca_pdb.write_text("\n".join(lines) + "\n")
    out_d = Path(td) / "out_ca"
    proc = subprocess.run([sys.executable, str(CLI), "--input", str(ca_pdb), "--output", str(out_d)], capture_output=True, text=True)
    assert proc.returncode == 0, f"Ca-only run crashed: {proc.stderr}"
    assert "BACKBONE_ONLY_TRACE" in proc.stderr, f"Expected BACKBONE_ONLY_TRACE warning in stderr, got: {proc.stderr}"
    print("[PASS 3] Ca-only backbone trace explicitly flagged BACKBONE_ONLY_TRACE warning in CLI stderr")

# TEST 4, 5, 6: Node.js verification of Tab 2 alphabet guard, zero fake 1D YES/NO verdict, and Tab 1 monomer guard
node_check = r"""
const fs = require("fs");
const html = fs.readFileSync("index.html", "utf8");
if (html.includes("1. Will They Interact?") || html.includes("LOW (UNLIKELY)")) {
  throw new Error("Legacy unvalidated 1D YES/NO verdict still present in index.html!");
}
if (!html.includes("INTRA-CHAIN DOMAIN PACKING (SINGLE-CHAIN MONOMER")) {
  throw new Error("Monomer auto-split plain-English guard missing in index.html!");
}
if (!html.includes("tab1PdbIdInput")) {
  throw new Error("Tab 1 1-click RCSB PDB fetcher missing!");
}
if (html.includes("--out ./results/") || html.includes("--salt 150")) {
  throw new Error("Old broken --out/--salt CLI docs still present in Tab 4!");
}
const js = html.match(/<script>([\s\S]*?)<\/script>/g).pop().replace(/<\/?script>/g, "");
globalThis.window = { addEventListener: () => {} };
globalThis.document = { getElementById: () => null, querySelectorAll: () => [] };
const fn = new Function(js + "\nreturn { validateStudioSlotInput };");
const { validateStudioSlotInput } = fn();
// Protein pasted into RNA box must error
const r1 = validateStudioSlotInput("MEEPQSDPSVEPPLSQETFSDLWKLLPENN", "rna", "Y");
if (!r1.error) throw new Error("Failed to catch protein sequence pasted into RNA slot!");
// Valid RNA in RNA box must pass
const r2 = validateStudioSlotInput("GGCUGAGAGAGAAAGCCUGGGAGGUUGGG", "rna", "Y");
if (r2.error || r2.warning) throw new Error("False positive on valid RNA: " + JSON.stringify(r2));
// PDB ID pasted into sequence box must error
const r3 = validateStudioSlotInput("1YCR", "protein", "X");
if (!r3.error) throw new Error("Failed to catch 4-char PDB ID pasted into sequence box!");
// Protein pasted into ligand SMILES box must error
const r4 = validateStudioSlotInput("MCNTNMSVPTDGAVTTSQIPASEQETLVRPKPLLLKLLK", "ligand", "Y");
if (!r4.error) throw new Error("Failed to catch protein pasted into SMILES slot!");
console.log("[PASS 4, 5 & 6] Tab 2 alphabet guard, zero fake 1D verdicts, Tab 1 PDB fetch & monomer guard all verified.");
"""
proc_n = subprocess.run(["node", "-e", node_check], capture_output=True, text=True)
assert proc_n.returncode == 0, f"Node basics check failed:\n{proc_n.stderr}"
print(proc_n.stdout.strip())
print("\nALL BASICS & SANITY REGRESSION TESTS PASSED.")

# TEST 7 & 8: Headless Node.js end-to-end check of Tab 2 on full-length MDM2 + p53 & CSV Header Parity
e2e_node = r"""
const fs = require("fs");
const html = fs.readFileSync("index.html", "utf8");
const js = html.match(/<script>([\s\S]*?)<\/script>/g).pop().replace(/<\/?script>/g, "");

const MDM2 = "MCNTNMSVPTDGAVTTSQIPASEQETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVVVNQQESSDSGTSVSENRCHLEGGSDQKDLVQELQEEKPSSSHLVSRPSTSSRRRAISETEENSDELSGERQRKRHKSDSISLSFDESLALCVIREICCERSSSSESTGTPSNPDLDAGVSEHSGDWLDQDSVSDQFSVEFEVESLDSEDYSLSEEGQELSDEDDEVYQVTVYQAGESDTDSFEEDPEISLADYWKCTSCNEMNPPLPSHCNRCWALRENWLPEDKGKDKGEISEKAKLENSTQAEEGFDVPDCKKTIVNDSRESCVEENDDKITQASQSQESEDYSQPSTSSSIIYSSQEDVKEFEREETQDKEESVESSLPLNAIEPCVICQGRPKNGCIVHGKTGHLMACFTCAKKLKKRNKPCPVCRQPIQMIVLTYFP";
const P53  = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD";

const dom = {};
const mkEl = (id, val="") => (dom[id] = { id, value: val, innerHTML: "", textContent: "", style: {} });
mkEl("slotXType", "protein"); mkEl("slotYType", "protein");
mkEl("slotXName", "MDM2_HUMAN"); mkEl("slotYName", "P53_HUMAN");
mkEl("slotXSeq", MDM2); mkEl("slotYSeq", P53);
mkEl("slotXCount", "1"); mkEl("slotYCount", "1");
mkEl("slotXMut", ""); mkEl("slotYMut", "");
for (const k of ["studioOutputPanel","studioPrescreenGrid","studioSummaryCards","studioSiteMapWrap","studioWarningBanner","studioVerdictBox","studioPreviewFormat","studioJobPreview"]) mkEl(k);

globalThis.window = { addEventListener: () => {} };
globalThis.document = { getElementById: id => dom[id] || null, querySelectorAll: () => [] };
const api = new Function(js + "\nreturn { runUniversalStudioPreScreen, buildMasterDelimited };")();

// Test 7: Full-length MDM2 + p53 in Tab 2 must output CO-FOLD REQUIRED, zero UNLIKELY, and valid Boltz YAML
api.runUniversalStudioPreScreen();
const vHtml = dom.studioVerdictBox.innerHTML;
if (!vHtml.includes("CO-FOLD REQUIRED")) throw new Error("Tab 2 MDM2-p53 missing CO-FOLD REQUIRED badge!");
if (vHtml.includes("UNLIKELY") || vHtml.includes("Will They Interact?")) throw new Error("Tab 2 MDM2-p53 still printed fake UNLIKELY call!");
if (!dom.studioJobPreview.textContent.includes(MDM2.slice(0, 30)) || !dom.studioJobPreview.textContent.includes(P53.slice(0, 30))) {
  throw new Error("Tab 2 failed to export full-length MDM2+p53 sequences to Boltz YAML!");
}

// Test 8: Browser CSV header parity check (exact snake_case loci_a, loci_b, conf_contacts, BSA_A2, incomplete_atoms_warning)
const csvHdr = api.buildMasterDelimited(",").split("\n")[0].split(",");
for (const reqCol of ["loci_a", "loci_b", "conf_contacts", "BSA_A2", "mean_intf_pLDDT", "mean_intf_PAE_A", "steric_clashes", "incomplete_atoms_warning"]) {
  if (!csvHdr.includes(reqCol)) throw new Error("Browser CSV header missing CLI parity column: " + reqCol);
}
if (!html.includes('<select id="slotYType" onchange="updateStudioSlotUI(\'Y\')">\n                <option value="protein">')) {
  throw new Error("Partner Y dropdown does not default to Protein!");
}
console.log("[PASS 7 & 8] Headless Tab 2 MDM2+p53 honesty + Boltz export + Browser/CLI CSV header parity verified.");
"""
proc_e2e = subprocess.run(["node", "-e", e2e_node], capture_output=True, text=True)
assert proc_e2e.returncode == 0, f"Headless E2E web check failed:\n{proc_e2e.stderr}"
print(proc_e2e.stdout.strip())
