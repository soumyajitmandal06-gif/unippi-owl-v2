#!/usr/bin/env python3
"""
Phase 2 - Step 2B-1.5 (Verbatim Anchors — CLI RNA Audit Polish & Ground-Truth Verification):
  1. Preserves ribose prime (') in mmCIF atom names (fixes O2'/O4'/C4' corruption & 5ZC9 RNA vs DNA).
  2. Exempts catalytic metal-ion inner-sphere coordination (<= 2.8 A) from clash penalties.
  3. Adds per-chain crystallographic B-factor plateau in MODE_C_EXP_CRYSTAL (zero change to AF3 MODE_A/B).
  4. Fixes Ensemble Summary guard (multi-model only), bumps MODEL_VERSION to UniPPI-v2.0-MARVAS,
     adds modality & RNA specificity columns to Master CSV, and prints full RNA channels in Step 1 & Step 4.
  5. Re-runs fold_vip1_bzip61.zip (zero PPI regression), 5ZC9.cif, and 1URN.cif + prints exact pi/Cat-pi/2'-OH pairs.
"""

import ast
import json
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
src = cli_path.read_text(encoding="utf-8")

# 1. Bump version string to UniPPI-v2.0-MARVAS
assert 'MODEL_VERSION = "UniPPI-v1.8"' in src, "Missing MODEL_VERSION anchor"
src = src.replace('MODEL_VERSION = "UniPPI-v1.8"', 'MODEL_VERSION = "UniPPI-v2.0-MARVAS"', 1)

# 2. Preserve ribose prime (') in mmCIF atom names (strip matching outer quotes only)
old_atom_strip = 'name = row.get("label_atom_id", row.get("auth_atom_id", "CA")).strip("\'\\"")'
new_atom_strip = (
    'name = row.get("label_atom_id", row.get("auth_atom_id", "CA")).strip()\n'
    '                if len(name) >= 2 and name[0] == name[-1] and name[0] in ("\'", \'"\'):\n'
    '                    name = name[1:-1]'
)
assert old_atom_strip in src, "Missing mmCIF label_atom_id strip anchor"
src = src.replace(old_atom_strip, new_atom_strip, 1)

# 3. Exempt catalytic metal-ion inner-sphere coordination (<= 2.8 A) from clash penalty
old_clash_blk = '''                    if d < 2.2:
                        pair_clash += 1'''
new_clash_blk = '''                    is_met_a = (ra.entity_type in ("ION", "LIGAND") and atom_a.element in CATALYTIC_METALS)
                    is_met_b = (rb.entity_type in ("ION", "LIGAND") and atom_b.element in CATALYTIC_METALS)
                    if d <= 2.8 and ((is_met_a and atom_b.element in ("O", "N", "S")) or (is_met_b and atom_a.element in ("O", "N", "S"))):
                        pair_sb += 1
                        continue

                    if d < 2.2:
                        pair_clash += 1'''
assert old_clash_blk in src, "Missing d < 2.2 clash anchor"
src = src.replace(old_clash_blk, new_clash_blk, 1)

# 4. Per-chain crystallographic B-factor plateau in MODE_C_EXP_CRYSTAL (leaves MODE_A/B 100% identical)
old_plat_blk = '''    all_bfacs = [r.mean_bfactor for r in (residues_a + residues_b)]
    b_min = min(all_bfacs) if all_bfacs else 10.0

    def residue_raw_plddt(res: Residue) -> float:
        b = res.mean_bfactor
        if mode in ("MODE_A_PRED_FULL", "MODE_B_PRED_COORD"):
            return clamp(b, 0.0, 100.0)
        return clamp(100.0 * math.exp(-max(0.0, b - b_min) / 60.0), 5.0, 100.0)

    def residue_plateau_conf(res: Residue) -> float:
        if mode in ("MODE_A_PRED_FULL", "MODE_B_PRED_COORD"):
            return clamp((res.mean_bfactor - 30.0) / 55.0, 0.02, 1.0)
        return math.exp(-max(0.0, res.mean_bfactor - b_min) / 60.0)'''

new_plat_blk = '''    b_min_a = min((r.mean_bfactor for r in residues_a), default=10.0)
    b_min_b = min((r.mean_bfactor for r in residues_b), default=10.0)
    b_min_chain = {chain_a_id: b_min_a, chain_b_id: b_min_b}

    def residue_raw_plddt(res: Residue) -> float:
        b = res.mean_bfactor
        if mode in ("MODE_A_PRED_FULL", "MODE_B_PRED_COORD"):
            return clamp(b, 0.0, 100.0)
        b_ref = max(b_min_chain.get(res.chain_id, 10.0) + 20.0, 45.0)
        return clamp(100.0 * math.exp(-max(0.0, b - b_ref) / 100.0), 15.0, 100.0)

    def residue_plateau_conf(res: Residue) -> float:
        if mode in ("MODE_A_PRED_FULL", "MODE_B_PRED_COORD"):
            return clamp((res.mean_bfactor - 30.0) / 55.0, 0.02, 1.0)
        b_ref = max(b_min_chain.get(res.chain_id, 10.0) + 20.0, 45.0)
        return clamp(math.exp(-max(0.0, res.mean_bfactor - b_ref) / 100.0), 0.15, 1.0)'''
assert old_plat_blk in src, "Missing verbatim residue_raw_plddt/residue_plateau_conf anchor"
src = src.replace(old_plat_blk, new_plat_blk, 1)

# 5. Add modality & RNA specificity columns to Master CSV
old_csv_head = '        "complex_id", "chain_pair", "sha256_id", "mode", "UIS", "UIS_threshold", "UIS_margin",'
new_csv_head = '        "complex_id", "chain_pair", "modality", "sha256_id", "mode", "UIS", "UIS_threshold", "UIS_margin",'
assert old_csv_head in src, "Missing csv_cols header anchor"
src = src.replace(old_csv_head, new_csv_head, 1)

old_csv_tail = '        "mean_w_intf", "top_hotspot", "top_hotspot_ddG_Ala"'
new_csv_tail = '        "mean_w_intf", "cation_pi", "hb_base", "hb_2oh", "hb_backbone", "rna_specificity_idx", "top_hotspot", "top_hotspot_ddG_Ala"'
assert old_csv_tail in src, "Missing csv_cols tail anchor"
src = src.replace(old_csv_tail, new_csv_tail, 1)

# 6. Fix Ensemble Summary guard + Step 1 & Step 4 modality-aware printouts via exact slice
start_tok = '    if len(all_records) > 1:\n        mean_uis = sum(r["UIS"] for r in all_records) / len(all_records)'
end_tok = '    print(f"[UniPPI-CLI] Saved Master CSV: {csv_path}")'
i0 = src.find(start_tok)
i1 = src.find(end_tok, i0)
assert i0 > 0 and i1 > i0, "Missing ensemble & step-by-step printout block"

new_ens_and_steps = '''    unique_complexes = {r["complex_id"] for r in all_records}
    if len(unique_complexes) > 1:
        mean_uis = sum(r["UIS"] for r in all_records) / len(all_records)
        mean_dg = sum(r["dG_uni_kcal"] for r in all_records) / len(all_records)
        mean_lis = sum(r["LIS"] for r in all_records) / len(all_records)
        mean_wbsa = sum(r["wBSA_A2"] for r in all_records) / len(all_records)
        print(f"[Ensemble Summary across N={len(all_records)} models] Mean UIS = {mean_uis:.3f} | Mean dG_uni = {mean_dg:.2f} kcal/mol | Mean LIS = {mean_lis:.3f} | Mean wBSA = {mean_wbsa:.1f} A^2")
    elif len(all_records) > 1:
        top0 = all_records[0]
        print(f"[Multi-Chain Complex Summary: {len(all_records)} interfaces detected in {top0['complex_id']}] Top Interface = {top0['chain_pair']} [{top0.get('modality', 'PPI')}] (UIS = {top0['UIS']:.3f}, dG = {top0['dG_uni_kcal']:.2f} kcal/mol)")
    if all_records:
        top = all_records[0]
        _mod = top.get("modality", "PROTEIN-PROTEIN")
        _c_pair = top["chain_pair"].split("-")
        _ca_lbl = _c_pair[0] if len(_c_pair) >= 2 else "A"
        _cb_lbl = _c_pair[1] if len(_c_pair) >= 2 else "B"
        print(f"\\n[Step-by-Step Mathematical Derivation for Top Interface: {top['complex_id']} ({top['chain_pair']} | {_mod})]")
        _conf_c = [c for c in top.get("contacts", []) if c.get("weight_w", 0) >= 0.25]
        def _fmt_segs(keys):
            nums = sorted({int("".join(ch for ch in k.split(":")[-1] if ch.isdigit()) or "0") for k in keys})
            if not nums: return "None"
            segs, s, p = [], nums[0], nums[0]
            for v in nums[1:]:
                if v - p > 15:
                    segs.append(f"{s}" if s == p else f"{s}-{p}")
                    s = v
                p = v
            segs.append(f"{s}" if s == p else f"{s}-{p}")
            return ", ".join(segs)
        _loc_a = _fmt_segs([c["res_a"] for c in _conf_c])
        _loc_b = _fmt_segs([c["res_b"] for c in _conf_c])
        _is_nuc = _mod in ("PROTEIN-RNA", "PROTEIN-DNA")
        _rna_str = f", pi={top.get('pi_stacking', 0)}, Cat-pi={top.get('cation_pi', 0)}, HB_base={top.get('hb_base', 0)}, HB_2OH={top.get('hb_2oh', 0)}, HB_bb={top.get('hb_backbone', 0)}, SSI={top.get('rna_specificity_idx', 0.0):.2f}" if _is_nuc else ""
        print(f"  1. Chain Topology & Contacts (d_min <= 5.5 A): {_ca_lbl} = {top.get('len_a', '?')} res, {_cb_lbl} = {top.get('len_b', '?')} res | Ordered Core Loci (w>=0.25): {_ca_lbl}:{_loc_a} <-> {_cb_lbl}:{_loc_b} | Confident Contacts: {len(_conf_c)} ({len(top.get('contacts', []))} raw) across {top['intf_residues']} interface residues (SB={top['salt_bridges']}, HB={top['h_bonds']}, HP={top['hydrophobic_contacts']}{_rna_str})")
        print(f"  2. Local Confidence Coherence: Mean Intf Score = {top['mean_intf_pLDDT']:.1f}, Mean Contact PAE = {top['mean_intf_PAE_A']:.2f} A, Mean Weight w_intf = {top['mean_w_intf']:.3f}, Cross-Chain LIS (PAE<=12A) = {top['LIS']:.3f}")
        print(f"  3. Ordered Buried Surface Area: 2-Chain dSASA = {top['BSA_A2']:.1f} A^2 (PISA 1/2 dSASA = {0.5*top['BSA_A2']:.1f} A^2) -> Confidence-Weighted wBSA = {top['wBSA_A2']:.1f} A^2")
        _kd_disp = f"{top['Kd_nM']*1000:.2f} pM" if top["Kd_nM"] < 0.05 else f"{top['Kd_nM']:.2f} nM"
        if _is_nuc:
            print(f"  4. Thermodynamics (Dedicated 9-Ch {_mod} Model, PRODIGY disabled): dG_RNA = {top['dG_uni_kcal']:.2f} kcal/mol (Kd = {_kd_disp})")
        else:
            _pct_p = round(100*(top['dG_uni_kcal']-top['dG_phys_kcal'])/(top['dG_prodigy_kcal']-top['dG_phys_kcal'])) if abs(top['dG_prodigy_kcal']-top['dG_phys_kcal'])>1e-3 else 0
            print(f"  4. Thermodynamics (Enthalpy-Entropy Compensated): dG_PRODIGY = {top['dG_prodigy_kcal']:.2f} kcal/mol, dG_phys = {top['dG_phys_kcal']:.2f} kcal/mol (blend: {_pct_p}% PRODIGY / {100-_pct_p}% phys) -> Consensus dG_uni = {top['dG_uni_kcal']:.2f} kcal/mol (Kd = {_kd_disp})")
'''
src = src[:i0] + new_ens_and_steps + src[i1:]

ast.parse(src)
cli_path.write_text(src, encoding="utf-8")
print("[PASS] Step 2B-1.5 polish applied cleanly to unippi_cli.py (ast.parse OK).\n")

for title, inp, outp in [
    ("ZERO-REGRESSION PPI GUARD (fold_vip1_bzip61.zip)", "benchmarks/fold_vip1_bzip61.zip", "/tmp/v2_p15_vip1"),
    ("POLISHED BENCHMARK 1: 5ZC9.cif (eIF4A1 + Polypurine RNA + Rocaglamide A)", "benchmarks/5ZC9.cif", "/tmp/v2_p15_5zc9"),
    ("POLISHED BENCHMARK 2: 1URN.cif (U1A RRM + Hairpin RNA — ITC Ground Truth dG ~ -13.8 kcal/mol)", "benchmarks/1URN.cif", "/tmp/v2_p15_1urn"),
]:
    print("\n" + "=" * 95)
    print(title)
    print("=" * 95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", inp, "-o", outp], check=True)

print("\n" + "=" * 95)
print("BIOLOGICAL GROUND-TRUTH CHANNEL AUDIT (5ZC9 vs 1URN — PRIMES, 2'-OH, PI & CAT-PI)")
print("=" * 95)
for label, jpath in [("5ZC9", "/tmp/v2_p15_5zc9/unippi_full_report.json"), ("1URN", "/tmp/v2_p15_1urn/unippi_full_report.json")]:
    rep = json.loads(Path(jpath).read_text(encoding="utf-8"))
    for r in rep:
        if r.get("modality") in ("PROTEIN-RNA", "PROTEIN-DNA") and r["UIS"] >= 0.5:
            print(f"[{label}] {r['chain_pair']} ({r['modality']}): UIS={r['UIS']} | dG={r['dG_uni_kcal']} | Kd={r['Kd_nM']} nM | "
                  f"SB={r['salt_bridges']} HB={r['h_bonds']} (Base={r['hb_base']}, 2'OH={r['hb_2oh']}, BB={r['hb_backbone']}) | "
                  f"pi={r['pi_stacking']} Cat-pi={r['cation_pi']} HP={r['hydrophobic_contacts']} | SSI={r['rna_specificity_idx']} | mean_w={r['mean_w_intf']}")
            cs = r.get("contacts", [])
            pi_pairs = [c for c in cs if c.get("pi", c.get("pi_stack", 0)) > 0 or c.get("cat_pi", 0) > 0]
            oh_pairs = [c for c in cs if c.get("hb_2oh", 0) > 0]
            if pi_pairs:
                print("   pi / Cat-pi :", ", ".join(f"{c['res_a']}<->{c['res_b']}(pi={c.get('pi',c.get('pi_stack',0))},catpi={c.get('cat_pi',0)})" for c in pi_pairs))
            if oh_pairs:
                print("   Ribose 2'-OH:", ", ".join(f"{c['res_a']}<->{c['res_b']}(2'OH={c.get('hb_2oh',0)})" for c in oh_pairs[:8]))
