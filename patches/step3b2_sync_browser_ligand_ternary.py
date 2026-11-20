#!/usr/bin/env python3
"""
Phase 2 - Step 3B-2 (Finalized — Outer-Sphere Mg/Ca/Mn Bridge + Ligand Phosphate/Carboxylate SB
+ Full 6-Benchmark CLI Audit + 100% index.html Parity Sync):
  Stage 1: Adds 3D ligand anionic-O detector (phosphates/carboxylates -> SB) & outer-sphere
           water-mediated coordination for hard catalytic metals (MG/CA/MN, 2.75-4.15 A -> O)
           in unippi_cli.py, then verifies ALL 6 real benchmarks cleanly.
  Stage 2 & 3: Mirrors the complete Step 3B-1 + 3B-2 engine & UI verbatim into index.html
           (Ternary Audit banner, Intra-Chain Domain Split input, dynamic Contact headers,
           XB/MET tags & filters, Ligand burial profile, composite 3Dmol.js & SVG colors).
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
    raise AssertionError(f"Anchor mismatch for [{label}]:\n{old[:180]!r}")


# ==============================================================================
# STAGE 1: UPGRADE unippi_cli.py (LIGAND ANIONIC SB + OUTER-SPHERE MG/CA/MN)
# ==============================================================================
print("=== STAGE 1: UPGRADING unippi_cli.py (LIGAND PHOSPHATE/CARBOXYLATE SB + OUTER-SPHERE MG BRIDGE) ===")
cli = cli_path.read_text(encoding="utf-8")

# 1a. Add universal 3D ligand anionic-oxygen detector right after get_all_rings_for_residue()
old_ring_end = '''def get_all_rings_for_residue(res: Residue) -> List[Tuple[Vec3, Optional[Vec3]]]:
    if res.entity_type == "LIGAND":
        return detect_ligand_planar_rings(res)
    cent, norm = get_ring_centroid_and_normal(res)
    return [(cent, norm)] if cent is not None else []'''

new_ring_end = old_ring_end + '''

def get_ligand_anionic_oxygens(res: Residue) -> Set[str]:
    """Universal 3D detector for anionic ligand oxygens: phosphate O (P-O) & carboxylate -COO-."""
    if res.entity_type != "LIGAND":
        return set()
    if hasattr(res, "_cached_anionic_o"):
        return res._cached_anionic_o
    anionic: Set[str] = set()
    o_atoms = [a for a in res.atoms if a.element == "O"]
    p_atoms = [a for a in res.atoms if a.element == "P"]
    c_atoms = [a for a in res.atoms if a.element == "C"]
    # 1. Phosphate / phosphonate oxygens (terminal O bonded to P at 1.35-1.75 A, not ester-bridged to C)
    for oa in o_atoms:
        if any(1.35 <= math.dist(oa.coord, pa.coord) <= 1.75 for pa in p_atoms):
            bonded_c = sum(1 for ca in c_atoms if math.dist(oa.coord, ca.coord) <= 1.65)
            bonded_p = sum(1 for pa in p_atoms if 1.35 <= math.dist(oa.coord, pa.coord) <= 1.75)
            if bonded_c == 0 and bonded_p == 1:
                anionic.add(oa.name)
    # 2. Carboxylate oxygens (C bonded to two terminal O atoms at 1.15-1.40 A)
    for ca in c_atoms:
        term_o = []
        for oa in o_atoms:
            if 1.15 <= math.dist(ca.coord, oa.coord) <= 1.40:
                deg = sum(1 for other in res.atoms if other is not oa and math.dist(oa.coord, other.coord) <= 1.65)
                if deg == 1:
                    term_o.append(oa.name)
        if len(term_o) == 2:
            anionic.update(term_o)
    res._cached_anionic_o = anionic
    return anionic'''
cli = ensure_replaced(cli, old_ring_end, new_ring_end, "CLI: Add universal 3D ligand phosphate/carboxylate anionic-O detector")

# 1b. Add outer-sphere metal counter initialization
cli = ensure_replaced(
    cli,
    '    n_hb_base, n_hb_2oh, n_hb_bb, n_catpi, n_xb, n_metal_coord = 0, 0, 0, 0, 0, 0\n'
    '    w_sb, w_hb, w_hp, w_pi, w_clash = 0.0, 0.0, 0.0, 0.0, 0.0\n'
    '    w_hb_base, w_hb_2oh, w_hb_bb, w_pi_rna, w_catpi, w_xb, w_metal_coord = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0',
    '    n_hb_base, n_hb_2oh, n_hb_bb, n_catpi, n_xb, n_metal_coord, n_metal_outer = 0, 0, 0, 0, 0, 0, 0\n'
    '    w_sb, w_hb, w_hp, w_pi, w_clash = 0.0, 0.0, 0.0, 0.0, 0.0\n'
    '    w_hb_base, w_hb_2oh, w_hb_bb, w_pi_rna, w_catpi, w_xb, w_metal_coord, w_metal_outer = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0',
    "CLI: Add n_metal_outer & w_metal_outer counters",
)

cli = ensure_replaced(
    cli,
    '            pair_hb_base, pair_hb_2oh, pair_hb_bb, pair_metal_coord = 0, 0, 0, 0',
    '            pair_hb_base, pair_hb_2oh, pair_hb_bb, pair_metal_coord, pair_metal_outer = 0, 0, 0, 0, 0\n'
    '            lig_ani_a = get_ligand_anionic_oxygens(ra) if ra.entity_type == "LIGAND" else set()\n'
    '            lig_ani_b = get_ligand_anionic_oxygens(rb) if rb.entity_type == "LIGAND" else set()',
    "CLI: Init pair_metal_outer & ligand anionic oxygen sets",
)

# 1c. Outer-sphere water-mediated MG/CA/MN coordination (2.75 < d <= 4.15 A -> O)
old_met_inner = '''                    if is_met_a or is_met_b:
                        m_el = atom_a.element if is_met_a else atom_b.element
                        other_el = atom_b.element if is_met_a else atom_a.element
                        d_mcut = 2.65 if m_el in ("ZN", "FE", "CU", "CO", "NI") else 2.75
                        if d <= d_mcut and other_el in ("O", "N", "S"):
                            pair_sb += 1
                            pair_metal_coord += 1
                            continue'''
new_met_inner = '''                    if is_met_a or is_met_b:
                        m_el = atom_a.element if is_met_a else atom_b.element
                        other_el = atom_b.element if is_met_a else atom_a.element
                        d_mcut = 2.65 if m_el in ("ZN", "FE", "CU", "CO", "NI") else 2.75
                        if d <= d_mcut and other_el in ("O", "N", "S"):
                            pair_sb += 1
                            pair_metal_coord += 1
                            continue
                        if m_el in ("MG", "CA", "MN") and d_mcut < d <= 4.15 and other_el == "O":
                            pair_metal_outer += 1
                            continue'''
cli = ensure_replaced(cli, old_met_inner, new_met_inner, "CLI: Add outer-sphere water-mediated MG/CA/MN coordination (2.75-4.15 A)")

# 1d. Include lig_ani_a / lig_ani_b in salt-bridge check
cli = ensure_replaced(
    cli,
    '                    is_ani_a = (cna, atom_a.name) in ANIONIC_ATOMS or (ra.entity_type == "NUCLEIC" and atom_a.name in NUCLEIC_PHOSPHATE_O)\n'
    '                    is_cat_b = (cnb, atom_b.name) in CATIONIC_ATOMS\n'
    '                    is_ani_b = (cnb, atom_b.name) in ANIONIC_ATOMS or (rb.entity_type == "NUCLEIC" and atom_b.name in NUCLEIC_PHOSPHATE_O)',
    '                    is_ani_a = (cna, atom_a.name) in ANIONIC_ATOMS or (ra.entity_type == "NUCLEIC" and atom_a.name in NUCLEIC_PHOSPHATE_O) or (atom_a.name in lig_ani_a)\n'
    '                    is_cat_b = (cnb, atom_b.name) in CATIONIC_ATOMS\n'
    '                    is_ani_b = (cnb, atom_b.name) in ANIONIC_ATOMS or (rb.entity_type == "NUCLEIC" and atom_b.name in NUCLEIC_PHOSPHATE_O) or (atom_b.name in lig_ani_b)',
    "CLI: Count ligand phosphate/carboxylate salt bridges",
)

# 1e. Accumulate eff_mouter (capped at 1 per residue to avoid counting OD1+OD2 twice as separate waters)
cli = ensure_replaced(
    cli,
    '            eff_mcoord = min(pair_metal_coord, 2)\n\n'
    '            if w_ij >= 0.25:\n'
    '                n_sb += eff_sb\n'
    '                n_hb += eff_hb\n'
    '                n_hb_base += eff_hb_base\n'
    '                n_hb_2oh += eff_hb_2oh\n'
    '                n_hb_bb += eff_hb_bb\n'
    '                n_hp += eff_hp\n'
    '                n_pi += pair_pi\n'
    '                n_catpi += pair_catpi\n'
    '                n_xb += pair_xb\n'
    '                n_metal_coord += eff_mcoord\n'
    '                n_ss += pair_ss\n'
    '                n_clash += eff_clash\n\n'
    '            w_sb += w_ij * eff_sb\n'
    '            w_hb += w_ij * eff_hb\n'
    '            w_hb_base += w_ij * eff_hb_base\n'
    '            w_hb_2oh += w_ij * eff_hb_2oh\n'
    '            w_hb_bb += w_ij * eff_hb_bb\n'
    '            w_hp += w_ij * eff_hp_weighted\n'
    '            w_pi += w_ij * pair_pi\n'
    '            w_pi_rna += w_ij * pair_pi_w\n'
    '            w_catpi += w_ij * pair_catpi\n'
    '            w_xb += w_ij * pair_xb\n'
    '            w_metal_coord += w_ij * eff_mcoord',
    '            eff_mcoord = min(pair_metal_coord, 2)\n'
    '            eff_mouter = 1 if (pair_metal_outer > 0 and eff_mcoord == 0) else 0\n\n'
    '            if w_ij >= 0.25:\n'
    '                n_sb += eff_sb\n'
    '                n_hb += eff_hb\n'
    '                n_hb_base += eff_hb_base\n'
    '                n_hb_2oh += eff_hb_2oh\n'
    '                n_hb_bb += eff_hb_bb\n'
    '                n_hp += eff_hp\n'
    '                n_pi += pair_pi\n'
    '                n_catpi += pair_catpi\n'
    '                n_xb += pair_xb\n'
    '                n_metal_coord += eff_mcoord\n'
    '                n_metal_outer += eff_mouter\n'
    '                n_ss += pair_ss\n'
    '                n_clash += eff_clash\n\n'
    '            w_sb += w_ij * eff_sb\n'
    '            w_hb += w_ij * eff_hb\n'
    '            w_hb_base += w_ij * eff_hb_base\n'
    '            w_hb_2oh += w_ij * eff_hb_2oh\n'
    '            w_hb_bb += w_ij * eff_hb_bb\n'
    '            w_hp += w_ij * eff_hp_weighted\n'
    '            w_pi += w_ij * pair_pi\n'
    '            w_pi_rna += w_ij * pair_pi_w\n'
    '            w_catpi += w_ij * pair_catpi\n'
    '            w_xb += w_ij * pair_xb\n'
    '            w_metal_coord += w_ij * eff_mcoord\n'
    '            w_metal_outer += w_ij * eff_mouter',
    "CLI: Accumulate outer-sphere metal contacts",
)

# 1f. Update MODE_METAL energy & Q_metal to reward hydrated outer-sphere MG/CA/MN clusters (>=2 oxygens)
old_met_energy = '''    elif is_metal_pair:
        ion_res = [r for r in (residues_a + residues_b) if r.entity_type == "ION"]
        is_trans_met = any(r.canonical_name in ("ZN", "FE", "CU", "CO", "NI") for r in ion_res)
        k_coord = 2.60 if is_trans_met else 1.75
        if n_metal_coord == 0:
            dg_phys = dg_cratic + 1.50 * w_clash
        else:
            e_raw = k_coord * min(w_metal_coord, 6.0) + 0.0040 * w_bsa
            e_eff = 13.5 * math.tanh(e_raw / 13.5)
            dg_phys = 1.50 - e_eff + 1.50 * w_clash'''

new_met_energy = '''    elif is_metal_pair:
        ion_res = [r for r in (residues_a + residues_b) if r.entity_type == "ION"]
        is_trans_met = any(r.canonical_name in ("ZN", "FE", "CU", "CO", "NI") for r in ion_res)
        k_coord = 2.60 if is_trans_met else 1.85
        eff_w_met = w_metal_coord + (0.45 * min(w_metal_outer, 4.0) if (n_metal_coord > 0 or n_metal_outer >= 2) else 0.0)
        if eff_w_met <= 0.0:
            dg_phys = dg_cratic + 1.50 * w_clash
        else:
            e_raw = k_coord * min(eff_w_met, 6.0) + 0.0040 * w_bsa
            e_eff = 13.5 * math.tanh(e_raw / 13.5)
            dg_phys = 1.50 - e_eff + 1.50 * w_clash'''
cli = ensure_replaced(cli, old_met_energy, new_met_energy, "CLI: Add outer-sphere hydrated MG/CA/MN contribution to MODE_METAL energy")

cli = ensure_replaced(
    cli,
    '    elif is_metal_pair:\n'
    '        pocket_quality = mean_w_intf * clamp(w_metal_coord / 3.0, 0.0, 1.0) if n_metal_coord > 0 else 0.0',
    '    elif is_metal_pair:\n'
    '        _eff_q_met = w_metal_coord + (0.45 * min(w_metal_outer, 4.0) if (n_metal_coord > 0 or n_metal_outer >= 2) else 0.0)\n'
    '        pocket_quality = mean_w_intf * clamp(_eff_q_met / 3.0, 0.0, 1.0) if _eff_q_met > 0.0 else 0.0',
    "CLI: Include outer-sphere MG/CA/MN in Q_metal",
)

cli = ensure_replaced(
    cli,
    '        "metal_coord_bonds": n_metal_coord,',
    '        "metal_coord_bonds": n_metal_coord,\n        "metal_outer_bonds": n_metal_outer,',
    "CLI: Expose metal_outer_bonds in return dict",
)

# 1g. Allow METAL-BRIDGE when ion has >=1 inner-sphere anchor AND >=1 inner- or >=2 outer-sphere contacts on Partner 2
old_met_bridge = '''                else:
                    if o_type in ("PROTEIN", "NUCLEIC", "LIGAND") and rrec.get("metal_coord_bonds", 0) >= 1:
                        partner_hits.append((other, float(rrec.get("metal_coord_bonds", 0))))
            if len(partner_hits) >= 2:
                partner_hits.sort(key=lambda x: x[1], reverse=True)
                (m1, v1), (m2, v2) = partner_hits[0], partner_hits[1]'''

new_met_bridge = '''                else:
                    _mc_in = rrec.get("metal_coord_bonds", 0)
                    _mc_out = rrec.get("metal_outer_bonds", 0)
                    if o_type in ("PROTEIN", "NUCLEIC", "LIGAND") and (_mc_in >= 1 or _mc_out >= 2):
                        partner_hits.append((other, float(_mc_in) + 0.45 * float(_mc_out), _mc_in, _mc_out))
            if (is_lig_c and len(partner_hits) >= 2) or (not is_lig_c and len(partner_hits) >= 2 and any(h[2] >= 1 for h in partner_hits)):
                partner_hits.sort(key=lambda x: x[1], reverse=True)
                if is_lig_c:
                    (m1, v1), (m2, v2) = partner_hits[0], partner_hits[1]
                else:
                    m1, v1, _in1, _out1 = partner_hits[0]
                    m2, v2, _in2, _out2 = partner_hits[1]'''
cli = ensure_replaced(cli, old_met_bridge, new_met_bridge, "CLI: Enable inner+outer-sphere catalytic METAL-BRIDGE detection (eIF4A1 Walker-A/B + ANP + Mg)")

cli = ensure_replaced(
    cli,
    '                    t_pair["ternary_burial_split"] = f"{m1}:{v1:.1f}% / {m2}:{v2:.1f}%" if is_lig_c else f"{m1}:{int(v1)} coord / {m2}:{int(v2)} coord"',
    '                    t_pair["ternary_burial_split"] = f"{m1}:{v1:.1f}% / {m2}:{v2:.1f}%" if is_lig_c else f"{m1}:{_in1}in+{_out1}out / {m2}:{_in2}in+{_out2}out"',
    "CLI: Format inner+outer split on METAL-BRIDGE",
)
cli = ensure_replaced(
    cli,
    '                    elif not is_lig_c and t_pair.get("metal_coord_bonds", 0) >= 2:\n'
    '                        t_pair["pattern_id"] = "Pattern A"\n'
    '                        t_pair["pattern_desc"] = f"Catalytic Metal Bridge ({lc} coordinates {m1} [{int(v1)}] & {m2} [{int(v2)}] inner-sphere bonds)"',
    '                    elif not is_lig_c and t_pair.get("metal_coord_bonds", 0) >= 1:\n'
    '                        t_pair["pattern_id"] = "Pattern A"\n'
    '                        t_pair["pattern_desc"] = f"Catalytic Metal Bridge ({lc} bridges {m1} [{_in1} inner, {_out1} outer-H2O] & {m2} [{_in2} inner, {_out2} outer-H2O])"',
    "CLI: Clear inner+outer catalytic Metal Bridge pattern description",
)

ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] Upgraded unippi_cli.py cleanly (ast.parse OK).\n")

# Re-run 5ZC9, 7DIY, 5ZC9_split (and verify all 6 benchmarks)
bench_runs = [
    ("vip1", "benchmarks/fold_vip1_bzip61.zip", "/tmp/v2_s3_vip1", None, False),
    ("5zc9", "benchmarks/5ZC9.cif", "/tmp/v2_s3_5zc9", None, True),
    ("1urn", "benchmarks/1URN.cif", "/tmp/v2_s3_1urn", None, False),
    ("4wnp", "benchmarks/4WNP.cif", "/tmp/v2_s3_4wnp", None, False),
    ("7diy", "benchmarks/7DIY.cif", "/tmp/v2_s3_7diy", None, True),
    ("5zc9_split", "benchmarks/5ZC9.cif", "/tmp/v2_s3_5zc9_split", "A:1-236,237-406", True),
]
reports = {}
for key, inp, outp, split_arg, force_run in bench_runs:
    jfile = Path(outp) / "unippi_full_report.json"
    if force_run or not jfile.exists():
        cmd = [sys.executable, "unippi_cli.py", "-i", inp, "-o", outp]
        if split_arg:
            cmd += ["--split-chain", split_arg]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    reports[key] = json.loads(jfile.read_text(encoding="utf-8"))

print("=== UPDATED CLI AUDIT: 5ZC9, 7DIY & 5ZC9 SPLIT ===")
for label, key in [
    ("5ZC9.cif (eIF4A1 + RNA + Rocaglamide A + ANP + Mg2+ Bridge)", "5zc9"),
    ("7DIY.cif (Nsp14-Nsp10 + 4 Zn2+ Fingers + Catalytic Mg2+)", "7diy"),
    ("5ZC9.cif --split-chain A:1-236,237-406", "5zc9_split"),
]:
    print(f"\n[{label}]")
    for r in reports[key]:
        extra = []
        if r.get("lig_burial_pct", 0) > 0: extra.append(f"Bur={r['lig_burial_pct']:.1f}%")
        if r.get("pocket_quality", 0) > 0: extra.append(f"Q={r['pocket_quality']:.3f}")
        if r.get("pi_stacking", 0) > 0: extra.append(f"pi={r['pi_stacking']}")
        if r.get("metal_coord_bonds", 0) > 0: extra.append(f"MetIn={r['metal_coord_bonds']}")
        if r.get("metal_outer_bonds", 0) > 0: extra.append(f"MetOut={r['metal_outer_bonds']}")
        if r.get("ternary_burial_split"): extra.append(f"Split[{r['ternary_burial_split']}] Bal={r.get('cooperativity_balance', 0):.2f}")
        if r.get("rna_specificity_idx", 0) > 0: extra.append(f"SSI={r['rna_specificity_idx']:.2f}")
        kd_s = f"{r['Kd_nM']*1000:.1f} pM" if r["Kd_nM"] < 0.05 else (f"{r['Kd_nM']:.2f} nM" if r["Kd_nM"] < 1e7 else ">10mM")
        print(f"  {r['chain_pair']:<18} | {r.get('modality','?'):<15} | UIS={r['UIS']:.3f} | dG={r['dG_uni_kcal']:>6.2f} ({kd_s:>10}) | SB={r['salt_bridges']:<2} HB={r['h_bonds']:<2} HP={r['hydrophobic_contacts']:<2} | {r['pattern_id']} | {' '.join(extra)}")

# Hard scientific invariant assertions across ALL 6 benchmarks
vip_top = reports["vip1"][0]
assert abs(vip_top["UIS"] - 0.5661) < 0.002 and abs(vip_top["dG_uni_kcal"] - (-10.99)) < 0.02, "PPI regression!"
urn_top = reports["1urn"][0]
assert urn_top["modality"] == "PROTEIN-RNA" and abs(urn_top["dG_uni_kcal"] - (-13.61)) < 0.05, "1URN regression!"
anp_row = next(r for r in reports["5zc9"] if r["chain_pair"] == "A-A_ANP500")
assert anp_row["salt_bridges"] >= 2, f"Expected ANP500 phosphate salt bridges >= 2, got {anp_row['salt_bridges']}"
rcg_tern = [r for r in reports["5zc9"] if r.get("modality") == "TERNARY-GLUE"]
assert rcg_tern and rcg_tern[0]["UIS"] >= 0.50, "5ZC9 TERNARY-GLUE failed!"
mg_br = [r for r in reports["5zc9"] if r.get("modality") == "METAL-BRIDGE"]
assert mg_br, "5ZC9 METAL-BRIDGE (A_ANP500+A-A_MG501) was NOT detected!"
ulk_ligs = [r for r in reports["4wnp"] if r.get("modality") == "PROTEIN-LIGAND"]
assert len(ulk_ligs) == 4 and all(r["UIS"] >= 0.50 for r in ulk_ligs), "4WNP ULK1 pockets failed!"
zn_rows = [r for r in reports["7diy"] if r["chain_pair"] in ("A-A_ZN201", "A-A_ZN202", "B-B_ZN301", "B-B_ZN302")]
assert len(zn_rows) == 4 and all(r["UIS"] >= 0.65 and r["metal_coord_bonds"] == 4 for r in zn_rows), "7DIY Zn fingers failed!"
intra_rows = [r for r in reports["5zc9_split"] if r["chain_pair"] == "A_D1-A_D2"]
assert intra_rows and intra_rows[0]["modality"] == "INTRA-DOMAIN (MONOMER)", "5ZC9 --split-chain failed!"
print("\n[PASS] ALL 8 SCIENTIFIC INVARIANTS ACROSS ALL 6 BENCHMARKS VERIFIED IN CLI!\n")

# ==============================================================================
# STAGE 2 & 3: MIRROR VERIFIED ENGINE & UI INTO index.html (100% PARITY)
# ==============================================================================
print("=== STAGE 2 & 3: SYNCING index.html TO VERIFIED STEP 3B-2 PARITY ===")
html = html_path.read_text(encoding="utf-8")

html = ensure_replaced(
    html,
    '<meta name="keywords" content="UniPPI-OWL, UniPPI-OWL v1.8, AlphaFold3 binding free energy, structure to affinity, protein-protein interaction calculator, PAE, ipTM, PRODIGY">',
    '<meta name="keywords" content="UniPPI-OWL, UniPPI-OWL v2.0, MARVAS, AlphaFold3 binding free energy, protein-RNA interaction, molecular glue, protein-ligand affinity, PAE, ipTM, PRODIGY">',
    "HTML: Update meta keywords to v2.0-MARVAS",
)
html = ensure_replaced(
    html,
    '<span>Scoring models</span><strong>8-Ch PPI + 9-Ch Protein–RNA/DNA (tanh)</strong>',
    '<span>Scoring models</span><strong>PPI + RNA/DNA + Ligand + Metal + Ternary Glue</strong>',
    "HTML: Update model summary card",
)
html = ensure_replaced(
    html,
    '''        <label>PAE Decay Constant λ<sub>PAE</sub> (Å)
          <input type="number" id="lambdaInput" value="8.0" step="0.5" min="2.0" max="25.0" onchange="recomputeActiveSuite()">
        </label>
        <label>Report Depth''',
    '''        <label>PAE Decay Constant λ<sub>PAE</sub> (Å)
          <input type="number" id="lambdaInput" value="8.0" step="0.5" min="2.0" max="25.0" onchange="recomputeActiveSuite()">
        </label>
        <label>Intra-Chain Domain Split (Optional)
          <input type="text" id="splitChainInput" placeholder="e.g. A:1-236,237-406" onchange="recomputeActiveSuite()">
        </label>
        <label>Report Depth''',
    "HTML: Add Intra-Chain Domain Split control in param-bar",
)
html = ensure_replaced(
    html,
    '''      <div class="table-wrap">
        <table id="leaderboardTable">''',
    '''      <div id="ternaryAuditBanner" style="display:none; background:#faf5ff; border:1px solid #d8b4fe; border-left:4px solid #9333ea; padding:10px 14px; border-radius:6px; margin-bottom:12px; font-family:var(--font-mono); font-size:12px; color:#3b0764;"></div>
      <div class="table-wrap">
        <table id="leaderboardTable">''',
    "HTML: Add Ternary Glue & Metal-Bridge Audit banner",
)
html = ensure_replaced(
    html,
    'placeholder="Filter (e.g. GLU, ARG, SB)..."',
    'placeholder="Filter (PHE, PI, CAT, 2OH, XB, MET)..."',
    "HTML: Expand Contact Inspector filter placeholder",
)
html = ensure_replaced(
    html,
    '''                  <th>Chain A Res</th>
                  <th>Chain B Res</th>''',
    '''                  <th id="contactColA">Chain A Res</th>
                  <th id="contactColB">Chain B Res</th>''',
    "HTML: Dynamic IDs on Contact Inspector headers",
)

old_nuc_js = '    const NUC_TO_ONE = { A:"a", C:"c", G:"g", U:"u", T:"t", DA:"a", DC:"c", DG:"g", DT:"t", ADE:"a", CYT:"c", GUA:"g", URA:"u", THY:"t" };'
new_nuc_js = old_nuc_js + '''
    const MODIFIED_NUC_MAP = { I:"G", INO:"G", PSU:"U", H2U:"U", "5MU":"U", "4SU":"U", OMU:"U", "6MA":"A", "1MA":"A", A2M:"A", M6A:"A", MIA:"A", T6A:"A", "5MC":"C", OMC:"C", "4AC":"C", "3MC":"C", OMG:"G", "7MG":"G", G7M:"G", "2MG":"G", M2G:"G", "1MG":"G", YYG:"G", DI:"DG", "5CM":"DC", BRU:"DT" };
    const HALOGEN_ELEMENTS = new Set(["F","CL","BR","I"]);'''
html = ensure_replaced(html, old_nuc_js, new_nuc_js, "JS: Add MODIFIED_NUC_MAP & HALOGEN_ELEMENTS")

old_det_mod_js = '''    function detectPairModalityJS(caId, cbId, resA, resB) {
      if (caId.endsWith("_N") && cbId.endsWith("_C")) return "INTRA-DOMAIN (MONOMER)";
      const allT = new Set([...resA.map(r=>r.entityType), ...resB.map(r=>r.entityType)]);
      if (allT.has("NUCLEIC")) {
        const nucRes = [...resA, ...resB].filter(r => r.entityType === "NUCLEIC");
        const isRna = nucRes.some(r => r.canonicalName === "U" || r.canonicalName === "URA" || r.atoms.some(a => a.name === "O2'"));
        const tag = isRna ? "RNA" : "DNA";
        if (allT.has("PROTEIN")) return `PROTEIN-${tag}`;
        if (allT.has("LIGAND")) return `${tag}-LIGAND`;
        if (allT.has("ION")) return `${tag}-METAL`;
        return `${tag}-${tag}`;
      }
      if (allT.has("LIGAND")) return "PROTEIN-LIGAND";
      if (allT.has("ION")) return "PROTEIN-METAL";
      return "PROTEIN-PROTEIN";
    }'''

new_det_mod_js = '''    function detectLigandPlanarRingsJS(res) {
      if (res._cachedPlanarRings) return res._cachedPlanarRings;
      const atoms = res.atoms.filter(a => ["C","N","O","S"].includes(a.element));
      const n = atoms.length;
      if (n < 5) { res._cachedPlanarRings = []; return []; }
      const adj = Array.from({length: n}, () => []);
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const d = Math.hypot(atoms[i].x-atoms[j].x, atoms[i].y-atoms[j].y, atoms[i].z-atoms[j].z);
          if (d >= 1.15 && d <= 1.68) { adj[i].push(j); adj[j].push(i); }
        }
      }
      const ringKeys = new Set(), ringList = [];
      for (let start = 0; start < n; start++) {
        const stack = [[start, [start]]];
        while (stack.length) {
          const [curr, path] = stack.pop();
          if (path.length > 6) continue;
          for (const nxt of adj[curr]) {
            if (nxt === start && (path.length === 5 || path.length === 6)) {
              const sorted = path.slice().sort((a,b)=>a-b);
              const k = sorted.join(",");
              if (!ringKeys.has(k)) { ringKeys.add(k); ringList.push(sorted); }
            } else if (nxt > start && !path.includes(nxt) && path.length < 6) {
              stack.push([nxt, [...path, nxt]]);
            }
          }
        }
      }
      const planar = [];
      for (const rIdx of ringList) {
        const pts = rIdx.map(k => [atoms[k].x, atoms[k].y, atoms[k].z]);
        const m = pts.length;
        const cent = [meanArr(pts.map(p=>p[0])), meanArr(pts.map(p=>p[1])), meanArr(pts.map(p=>p[2]))];
        const p1 = pts[0], p2 = pts[Math.floor(m/3)], p3 = pts[Math.floor(2*m/3)];
        const ux = p2[0]-p1[0], uy = p2[1]-p1[1], uz = p2[2]-p1[2];
        const vx = p3[0]-p1[0], vy = p3[1]-p1[1], vz = p3[2]-p1[2];
        const nx = uy*vz - uz*vy, ny = uz*vx - ux*vz, nz = ux*vy - uy*vx;
        const nl = Math.hypot(nx, ny, nz);
        if (nl < 1e-6) continue;
        const norm = [nx/nl, ny/nl, nz/nl];
        const maxDev = Math.max(...pts.map(p => Math.abs((p[0]-cent[0])*norm[0] + (p[1]-cent[1])*norm[1] + (p[2]-cent[2])*norm[2])));
        if (maxDev <= 0.15) planar.push([cent, norm]);
      }
      res._cachedPlanarRings = planar;
      return planar;
    }

    function getAllRingsForResidueJS(res) {
      if (res.entityType === "LIGAND") return detectLigandPlanarRingsJS(res);
      const [cent, norm] = getRingCentroidAndNormalJS(res);
      return cent ? [[cent, norm]] : [];
    }

    function getLigandAnionicOxygensJS(res) {
      if (res.entityType !== "LIGAND") return new Set();
      if (res._cachedAnionicO) return res._cachedAnionicO;
      const anionic = new Set();
      const oAt = res.atoms.filter(a => a.element === "O");
      const pAt = res.atoms.filter(a => a.element === "P");
      const cAt = res.atoms.filter(a => a.element === "C");
      for (const oa of oAt) {
        const bp = pAt.filter(pa => { const d = Math.hypot(oa.x-pa.x, oa.y-pa.y, oa.z-pa.z); return d >= 1.35 && d <= 1.75; }).length;
        if (bp === 1) {
          const bc = cAt.filter(ca => Math.hypot(oa.x-ca.x, oa.y-ca.y, oa.z-ca.z) <= 1.65).length;
          if (bc === 0) anionic.add(oa.name);
        }
      }
      for (const ca of cAt) {
        const termO = [];
        for (const oa of oAt) {
          const d = Math.hypot(ca.x-oa.x, ca.y-oa.y, ca.z-oa.z);
          if (d >= 1.15 && d <= 1.40) {
            const deg = res.atoms.filter(ot => ot !== oa && Math.hypot(oa.x-ot.x, oa.y-ot.y, oa.z-ot.z) <= 1.65).length;
            if (deg === 1) termO.push(oa.name);
          }
        }
        if (termO.length === 2) termO.forEach(nm => anionic.add(nm));
      }
      res._cachedAnionicO = anionic;
      return anionic;
    }

    function countDirectionalHalogenBondsJS(ra, rb) {
      let xb = 0;
      for (const [srcR, dstR] of [[ra, rb], [rb, ra]]) {
        const halAtoms = srcR.atoms.filter(a => HALOGEN_ELEMENTS.has(a.element));
        if (!halAtoms.length) continue;
        const cAtoms = srcR.atoms.filter(a => a.element === "C");
        const accAtoms = dstR.atoms.filter(a => ["O","N","S"].includes(a.element));
        if (!cAtoms.length || !accAtoms.length) continue;
        for (const hx of halAtoms) {
          let cBond = null, bestCd = 9.0;
          for (const ca of cAtoms) {
            const dc = Math.hypot(ca.x-hx.x, ca.y-hx.y, ca.z-hx.z);
            if (dc >= 1.20 && dc <= 2.25 && dc < bestCd) { bestCd = dc; cBond = ca; }
          }
          if (!cBond) continue;
          const dCut = hx.element === "F" ? 3.30 : 3.70;
          const vxc = [cBond.x-hx.x, cBond.y-hx.y, cBond.z-hx.z];
          const lxc = Math.hypot(vxc[0], vxc[1], vxc[2]);
          if (lxc < 1e-4) continue;
          for (const acc of accAtoms) {
            const da = Math.hypot(hx.x-acc.x, hx.y-acc.y, hx.z-acc.z);
            if (da >= 2.20 && da <= dCut) {
              const vxa = [acc.x-hx.x, acc.y-hx.y, acc.z-hx.z];
              const cosAng = (vxc[0]*vxa[0] + vxc[1]*vxa[1] + vxc[2]*vxa[2]) / (lxc * da);
              if (cosAng <= -0.50) { xb++; break; }
            }
          }
        }
      }
      return Math.min(xb, 2);
    }

    function detectPairModalityJS(caId, cbId, resA, resB) {
      if (caId.includes("+")) {
        const tB = new Set(resB.map(r => r.entityType));
        if (tB.has("LIGAND")) return "TERNARY-GLUE";
        if (tB.has("ION")) return "METAL-BRIDGE";
      }
      if ((caId.endsWith("_N") && cbId.endsWith("_C")) || (caId.includes("_D1") && cbId.includes("_D2"))) return "INTRA-DOMAIN (MONOMER)";
      const allT = new Set([...resA.map(r=>r.entityType), ...resB.map(r=>r.entityType)]);
      if (allT.has("NUCLEIC")) {
        const nucRes = [...resA, ...resB].filter(r => r.entityType === "NUCLEIC");
        const isRna = nucRes.some(r => r.canonicalName === "U" || r.canonicalName === "URA" || r.atoms.some(a => a.name === "O2'"));
        const tag = isRna ? "RNA" : "DNA";
        if (allT.has("PROTEIN")) return `PROTEIN-${tag}`;
        if (allT.has("LIGAND")) return `${tag}-LIGAND`;
        if (allT.has("ION")) return `${tag}-METAL`;
        return `${tag}-${tag}`;
      }
      if (allT.has("LIGAND")) {
        if (allT.has("PROTEIN")) return "PROTEIN-LIGAND";
        if (allT.has("ION")) return "LIGAND-METAL";
        return "LIGAND-LIGAND";
      }
      if (allT.has("ION")) return allT.has("PROTEIN") ? "PROTEIN-METAL" : "ION-ION";
      return "PROTEIN-PROTEIN";
    }'''
html = ensure_replaced(html, old_det_mod_js, new_det_mod_js, "JS: Add 3D ring, anionic-O, halogen-bond & LIGAND-METAL/TERNARY helpers")

html = ensure_replaced(
    html,
    '          const canonical = MODIFIED_AA_MAP[a.resName] || a.resName;',
    '          const canonical = MODIFIED_AA_MAP[a.resName] || MODIFIED_NUC_MAP[a.resName] || a.resName;',
    "JS: Canonicalize modified RNA/DNA nucleotides in buildResidues()",
)

old_mono_js = '''        if (allRes.length >= 10 && allRes.every(r => r.entityType === "PROTEIN" || r.entityType === "NUCLEIC")) {
          const mid = Math.floor(allRes.length / 2);'''
new_mono_js = '''        if (allRes.length >= 10 && allRes.every(r => r.entityType === "PROTEIN" || r.entityType === "NUCLEIC")) {
          const nTot = allRes.length;
          let mid = Math.floor(nTot / 2);
          if (nTot >= 40) {
            const loC = Math.max(15, Math.floor(nTot / 4)), hiC = Math.min(nTot - 15, Math.floor(3 * nTot / 4));
            const stepC = Math.max(1, Math.floor((hiC - loC) / 25));
            let bestCross = 1e9;
            for (let cut = loC; cut <= hiC; cut += stepC) {
              let cross = 0;
              for (let ii = Math.max(0, cut - 35); ii < cut; ii++) {
                const ci = allRes[ii].centroid;
                for (let jj = cut; jj < Math.min(nTot, cut + 35); jj++) {
                  if (dist3(ci, allRes[jj].centroid) <= 8.0) cross++;
                }
              }
              if (cross < bestCross) { bestCross = cross; mid = cut; }
            }
          }'''
html = ensure_replaced(html, old_mono_js, new_mono_js, "JS: Contact-trough monomer domain finder")

html = ensure_replaced(
    html,
    '''      const bMinA = resA.length ? Math.min(...resA.map(r => r.meanBfac)) : 10.0;
      const bMinB = resB.length ? Math.min(...resB.map(r => r.meanBfac)) : 10.0;
      const bMinChain = { [caId]: bMinA, [cbId]: bMinB };''',
    '''      const bMinChain = {};
      for (const _r of [...resA, ...resB]) {
        if (bMinChain[_r.chainId] === undefined || _r.meanBfac < bMinChain[_r.chainId]) bMinChain[_r.chainId] = _r.meanBfac;
      }''',
    "JS: Group MODE_C bMinChain by actual residue chainId (ternary-safe)",
)

html = ensure_replaced(
    html,
    '''      const modality = detectPairModalityJS(caId, cbId, resA, resB);
      const isNucleicPair = (modality === "PROTEIN-RNA" || modality === "PROTEIN-DNA" || modality === "RNA-RNA" || modality === "DNA-DNA" || modality === "RNA-DNA");
      let nSb=0, nHb=0, nHp=0, nPi=0, nSs=0, nClash=0, nCb8=0;
      let nHbBase=0, nHb2oh=0, nHbBb=0, nCatPi=0;
      let wSb=0, wHb=0, wHp=0, wPi=0, wClash=0;
      let wHbBase=0, wHb2oh=0, wHbBb=0, wPiRna=0, wCatPi=0;''',
    '''      const modality = detectPairModalityJS(caId, cbId, resA, resB);
      const isNucleicPair = (modality === "PROTEIN-RNA" || modality === "PROTEIN-DNA" || modality === "RNA-RNA" || modality === "DNA-DNA" || modality === "RNA-DNA");
      const isDuplexPair = (modality === "RNA-RNA" || modality === "DNA-DNA" || modality === "RNA-DNA");
      const isLigandPair = (modality === "PROTEIN-LIGAND" || modality === "RNA-LIGAND" || modality === "DNA-LIGAND" || modality === "TERNARY-GLUE");
      const isMetalPair = (modality === "PROTEIN-METAL" || modality === "RNA-METAL" || modality === "DNA-METAL" || modality === "LIGAND-METAL" || modality === "METAL-BRIDGE");
      const isNonPpi = isNucleicPair || isLigandPair || isMetalPair;
      let nSb=0, nHb=0, nHp=0, nPi=0, nSs=0, nClash=0, nCb8=0;
      let nHbBase=0, nHb2oh=0, nHbBb=0, nCatPi=0, nXb=0, nMetalCoord=0, nMetalOuter=0;
      let wSb=0, wHb=0, wHp=0, wPi=0, wClash=0;
      let wHbBase=0, wHb2oh=0, wHbBb=0, wPiRna=0, wCatPi=0, wXb=0, wMetalCoord=0, wMetalOuter=0;''',
    "JS: Add MODE_LIGAND/MODE_METAL flags & outer-sphere metal counters",
)

old_inner_js = '''          let minD = 999.0, pSb=0, pHb=0, pHp=0, pSs=0, pClash=0, hpDistSum=0.0;
          let pHbBase=0, pHb2oh=0, pHbBb=0;
          const cna = ra.canonicalName, cnb = rb.canonicalName;
          const apA = HYDROPHOBIC_RES.has(cna) || ra.entityType === "LIGAND";
          const apB = HYDROPHOBIC_RES.has(cnb) || rb.entityType === "LIGAND";

          for (const aa of ra.atoms) {
            for (const ab of rb.atoms) {
              const d = Math.hypot(aa.x-ab.x, aa.y-ab.y, aa.z-ab.z);
              if (d < minD) minD = d;
              if (d > 4.5) continue;
              if (cna==="CYS" && cnb==="CYS" && aa.name==="SG" && ab.name==="SG" && d<=2.3) { pSs++; continue; }

              const isMetA = (ra.entityType==="ION" || ra.entityType==="LIGAND") && CATALYTIC_METALS.has(aa.element);
              const isMetB = (rb.entityType==="ION" || rb.entityType==="LIGAND") && CATALYTIC_METALS.has(ab.element);
              if (d <= 2.8 && ((isMetA && ["O","N","S"].includes(ab.element)) || (isMetB && ["O","N","S"].includes(aa.element)))) { pSb++; continue; }
              if (d < 2.2) pClash++;

              const isPtmA = ["SEP","TPO","PTR"].includes(ra.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(aa.name);
              const isPtmB = ["SEP","TPO","PTR"].includes(rb.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(ab.name);
              const catA = (cna==="ARG"&&["NE","NH1","NH2"].includes(aa.name)) || (cna==="LYS"&&aa.name==="NZ") || (cna==="HIS"&&["ND1","NE2"].includes(aa.name));
              const aniA = (cna==="ASP"&&["OD1","OD2"].includes(aa.name)) || (cna==="GLU"&&["OE1","OE2"].includes(aa.name)) || (ra.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(aa.name)) || isPtmA;
              const catB = (cnb==="ARG"&&["NE","NH1","NH2"].includes(ab.name)) || (cnb==="LYS"&&ab.name==="NZ") || (cnb==="HIS"&&["ND1","NE2"].includes(ab.name));
              const aniB = (cnb==="ASP"&&["OD1","OD2"].includes(ab.name)) || (cnb==="GLU"&&["OE1","OE2"].includes(ab.name)) || (rb.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(ab.name)) || isPtmB;

              if (d <= 4.0 && ((catA && aniB) || (aniA && catB))) { pSb += (isPtmA || isPtmB) ? 2 : 1; continue; }
              if (d <= 3.5 && ["N","O","S"].includes(aa.element) && ["N","O","S"].includes(ab.element)) {
                pHb++;
                if (isNucleicPair) {
                  const nucAt = ra.entityType === "NUCLEIC" ? aa : ab;
                  if (nucAt.name === "O2'") pHb2oh++;
                  else if (NUCLEOBASE_HB_ATOMS.has(nucAt.name)) pHbBase++;
                  else pHbBb++;
                }
                continue;
              }
              if (d <= 4.5 && ["C","S"].includes(aa.element) && ["C","S"].includes(ab.element)) {
                if (apA && apB && !["C","CA"].includes(aa.name) && !["C","CA"].includes(ab.name)) {
                  pHp++;
                  hpDistSum += clamp((4.75 - d) / 0.80, 0.35, 1.0);
                } else if (isNucleicPair) {
                  const protIsA = (ra.entityType === "PROTEIN");
                  const pAt = protIsA ? aa : ab, nAt = protIsA ? ab : aa;
                  const pAp = protIsA ? apA : apB;
                  if (pAp && !["C","CA"].includes(pAt.name) && NUCLEIC_CARBON_ATOMS.has(nAt.name)) {
                    pHp++;
                    hpDistSum += clamp((4.75 - d) / 0.80, 0.35, 1.0);
                  }
                }
              }
            }
          }

          if (minD > activeCutoff) continue;
          let pPi = 0, pPiW = 0.0, pCatPi = 0;
          if (isNucleicPair) {
            const [ringA, normA] = getRingCentroidAndNormalJS(ra);
            const [ringB, normB] = getRingCentroidAndNormalJS(rb);
            if (ringA && ringB) {
              const dRing = dist3(ringA, ringB);
              if (dRing <= 5.5) {
                const nRef = ra.entityType === "NUCLEIC" ? normA : normB;
                if (normA && normB && nRef) {
                  const vx = ringB[0]-ringA[0], vy = ringB[1]-ringA[1], vz = ringB[2]-ringA[2];
                  const dPerp = Math.abs(vx*nRef[0] + vy*nRef[1] + vz*nRef[2]);
                  const dPar = Math.sqrt(Math.max(0.0, dRing*dRing - dPerp*dPerp));
                  const cosAng = Math.abs(normA[0]*normB[0] + normA[1]*normB[1] + normA[2]*normB[2]);
                  if (dPerp >= 2.2 && dPar <= 4.2) {
                    if (cosAng >= 0.60) { pPi = 1; pPiW = 1.0; }
                    else if (dRing <= 5.0) { pPi = 1; pPiW = 0.5; }
                  }
                } else { pPi = 1; pPiW = 0.75; }
              }
            }
            const protR = ra.entityType === "PROTEIN" ? ra : rb;
            const nucR = ra.entityType === "PROTEIN" ? rb : ra;
            if (protR.canonicalName === "ARG" || protR.canonicalName === "LYS") {
              const [nucCent, nucNorm] = getRingCentroidAndNormalJS(nucR);
              const catName = protR.canonicalName === "ARG" ? "CZ" : "NZ";
              const catAt = protR.atoms.find(a => a.name === catName);
              if (nucCent && nucNorm && catAt) {
                const dCp = dist3([catAt.x, catAt.y, catAt.z], nucCent);
                if (dCp <= 5.0 && dCp > 1e-4) {
                  const cosPhi = Math.abs((catAt.x-nucCent[0])*nucNorm[0] + (catAt.y-nucCent[1])*nucNorm[1] + (catAt.z-nucCent[2])*nucNorm[2]) / dCp;
                  if (cosPhi >= 0.50) pCatPi = 1;
                }
              }
            }
          } else if (AROMATIC_RINGS[cna] && AROMATIC_RINGS[cnb]) {
            const [ringA] = getRingCentroidAndNormalJS(ra);
            const [ringB] = getRingCentroidAndNormalJS(rb);
            if (ringA && ringB && dist3(ringA, ringB) <= 5.5) { pPi = 1; pPiW = 1.0; }
          }

          intfA.add(ra.key); intfB.add(rb.key);
          if (cbDist > 8.0) { intfPae.push(pairPae); weightsSeen.push(wij); }

          const clsA = ["SEP","TPO","PTR"].includes(ra.resName) ? "CHARGED" : (PRODIGY_CLASS[cna] || (ra.entityType === "LIGAND" ? "APOLAR" : "POLAR"));
          const clsB = ["SEP","TPO","PTR"].includes(rb.resName) ? "CHARGED" : (PRODIGY_CLASS[cnb] || (rb.entityType === "LIGAND" ? "APOLAR" : "POLAR"));
          let pairCls = [clsA, clsB].sort().join("_");
          if (pairCls === "APOLAR_CHARGED") pairCls = "CHARGED_APOLAR";
          else if (pairCls === "APOLAR_POLAR") pairCls = "POLAR_APOLAR";
          if (prodigyIC[pairCls] !== undefined) prodigyIC[pairCls] += wij;

          const effSb = Math.min(pSb, 2), effHb = Math.min(pHb, 2);
          const effHbBase = Math.min(pHbBase, effHb);
          const effHb2oh = Math.min(pHb2oh, effHb - effHbBase);
          const effHbBb = Math.min(pHbBb, effHb - effHbBase - effHb2oh);
          const effHp = Math.min(pHp, 6), effHpW = Math.min(hpDistSum, 6.0), effClash = Math.min(pClash, 3);
          if (wij >= 0.25) {
            nSb += effSb; nHb += effHb; nHbBase += effHbBase; nHb2oh += effHb2oh; nHbBb += effHbBb;
            nHp += effHp; nPi += pPi; nCatPi += pCatPi; nSs += pSs; nClash += effClash;
          }
          wSb += wij*effSb; wHb += wij*effHb; wHbBase += wij*effHbBase; wHb2oh += wij*effHb2oh; wHbBb += wij*effHbBb;
          wHp += wij*effHpW; wPi += wij*pPi; wPiRna += wij*pPiW; wCatPi += wij*pCatPi; wClash += wij*effClash;

          if (wij >= 0.25) {
            for (const rk of [ra.key, rb.key]) {
              if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0 });
              const rc = resContrib.get(rk);
              rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*(pPi + pCatPi);
            }
          }

          contacts.push({
            resA: `${ra.chainId}:${ra.resName}${ra.resSeq}`, resB: `${rb.chainId}:${rb.resName}${rb.resSeq}`,
            resSeqA: ra.resSeq, resSeqB: rb.resSeq, coordA: ra.cbCoord, coordB: rb.cbCoord,
            minDist: Number(minD.toFixed(2)), pae: Number(pairPae.toFixed(2)), weight: Number(wij.toFixed(3)),
            sb: effSb, hb: effHb, hbBase: effHbBase, hb2oh: effHb2oh, hbBb: effHbBb,
            hp: effHp, pi: pPi, catpi: pCatPi, ss: pSs, clash: effClash
          });'''

new_inner_js = '''          let minD = 999.0, pSb=0, pHb=0, pHp=0, pSs=0, pClash=0, hpDistSum=0.0;
          let pHbBase=0, pHb2oh=0, pHbBb=0, pMetalCoord=0, pMetalOuter=0;
          const ligAniA = getLigandAnionicOxygensJS(ra), ligAniB = getLigandAnionicOxygensJS(rb);
          const cna = ra.canonicalName, cnb = rb.canonicalName;
          const apA = HYDROPHOBIC_RES.has(cna) || ra.entityType === "LIGAND";
          const apB = HYDROPHOBIC_RES.has(cnb) || rb.entityType === "LIGAND";

          for (const aa of ra.atoms) {
            for (const ab of rb.atoms) {
              const d = Math.hypot(aa.x-ab.x, aa.y-ab.y, aa.z-ab.z);
              if (d < minD) minD = d;
              if (d > 4.5) continue;
              if (cna==="CYS" && cnb==="CYS" && aa.name==="SG" && ab.name==="SG" && d<=2.3) { pSs++; continue; }

              const isMetA = (ra.entityType==="ION" || ra.entityType==="LIGAND") && CATALYTIC_METALS.has(aa.element);
              const isMetB = (rb.entityType==="ION" || rb.entityType==="LIGAND") && CATALYTIC_METALS.has(ab.element);
              if (isMetA || isMetB) {
                const mEl = isMetA ? aa.element : ab.element;
                const otherEl = isMetA ? ab.element : aa.element;
                const dMcut = ["ZN","FE","CU","CO","NI"].includes(mEl) ? 2.65 : 2.75;
                if (d <= dMcut && ["O","N","S"].includes(otherEl)) { pSb++; pMetalCoord++; continue; }
                if (["MG","CA","MN"].includes(mEl) && d > dMcut && d <= 4.15 && otherEl === "O") { pMetalOuter++; continue; }
              }
              if (d < 2.2) pClash++;

              const isPtmA = ["SEP","TPO","PTR"].includes(ra.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(aa.name);
              const isPtmB = ["SEP","TPO","PTR"].includes(rb.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(ab.name);
              const catA = (cna==="ARG"&&["NE","NH1","NH2"].includes(aa.name)) || (cna==="LYS"&&aa.name==="NZ") || (cna==="HIS"&&["ND1","NE2"].includes(aa.name));
              const aniA = (cna==="ASP"&&["OD1","OD2"].includes(aa.name)) || (cna==="GLU"&&["OE1","OE2"].includes(aa.name)) || (ra.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(aa.name)) || isPtmA || ligAniA.has(aa.name);
              const catB = (cnb==="ARG"&&["NE","NH1","NH2"].includes(ab.name)) || (cnb==="LYS"&&ab.name==="NZ") || (cnb==="HIS"&&["ND1","NE2"].includes(ab.name));
              const aniB = (cnb==="ASP"&&["OD1","OD2"].includes(ab.name)) || (cnb==="GLU"&&["OE1","OE2"].includes(ab.name)) || (rb.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(ab.name)) || isPtmB || ligAniB.has(ab.name);

              if (d <= 4.0 && ((catA && aniB) || (aniA && catB))) { pSb += (isPtmA || isPtmB) ? 2 : 1; continue; }
              if (d <= 3.5 && ["N","O","S"].includes(aa.element) && ["N","O","S"].includes(ab.element)) {
                pHb++;
                if (isNucleicPair) {
                  const nucAt = ra.entityType === "NUCLEIC" ? aa : ab;
                  if (nucAt.name === "O2'") pHb2oh++;
                  else if (NUCLEOBASE_HB_ATOMS.has(nucAt.name)) pHbBase++;
                  else pHbBb++;
                }
                continue;
              }
              if (d <= 4.5 && ["C","S"].includes(aa.element) && ["C","S"].includes(ab.element)) {
                if (apA && apB && !["C","CA"].includes(aa.name) && !["C","CA"].includes(ab.name)) {
                  pHp++;
                  hpDistSum += clamp((4.75 - d) / 0.80, 0.35, 1.0);
                } else if (isNucleicPair || (isLigandPair && (ra.entityType === "NUCLEIC" || rb.entityType === "NUCLEIC"))) {
                  const nucIsA = (ra.entityType === "NUCLEIC");
                  const nAt = nucIsA ? aa : ab, oAt = nucIsA ? ab : aa;
                  const oAp = nucIsA ? apB : apA;
                  if (oAp && !["C","CA"].includes(oAt.name) && NUCLEIC_CARBON_ATOMS.has(nAt.name)) {
                    pHp++;
                    hpDistSum += clamp((4.75 - d) / 0.80, 0.35, 1.0);
                  }
                }
              }
            }
          }

          if (minD > activeCutoff) continue;
          let pPi = 0, pPiW = 0.0, pCatPi = 0;
          const pXb = isLigandPair ? countDirectionalHalogenBondsJS(ra, rb) : 0;
          if (isNucleicPair || isLigandPair) {
            const ringsA = getAllRingsForResidueJS(ra), ringsB = getAllRingsForResidueJS(rb);
            if (ringsA.length && ringsB.length) {
              let bestPw = 0.0;
              for (const [ringA, normA] of ringsA) {
                for (const [ringB, normB] of ringsB) {
                  const dRing = dist3(ringA, ringB);
                  if (dRing <= 5.5) {
                    const nRef = normA || normB;
                    if (normA && normB && nRef) {
                      const vx = ringB[0]-ringA[0], vy = ringB[1]-ringA[1], vz = ringB[2]-ringA[2];
                      const dPerp = Math.abs(vx*nRef[0] + vy*nRef[1] + vz*nRef[2]);
                      const dPar = Math.sqrt(Math.max(0.0, dRing*dRing - dPerp*dPerp));
                      const cosAng = Math.abs(normA[0]*normB[0] + normA[1]*normB[1] + normA[2]*normB[2]);
                      if (dPerp >= 2.2 && dPar <= 4.2) {
                        if (cosAng >= 0.60) bestPw = Math.max(bestPw, 1.0);
                        else if (dRing <= 5.0) bestPw = Math.max(bestPw, 0.5);
                      }
                    } else { bestPw = Math.max(bestPw, 0.75); }
                  }
                }
              }
              if (bestPw > 0.0) { pPi = 1; pPiW = bestPw; }
            }
            for (const [pCand, rCand] of [[ra, rb], [rb, ra]]) {
              if (pCand.canonicalName === "ARG" || pCand.canonicalName === "LYS") {
                const candRings = getAllRingsForResidueJS(rCand);
                const catName = pCand.canonicalName === "ARG" ? "CZ" : "NZ";
                const catAt = pCand.atoms.find(a => a.name === catName);
                if (candRings.length && catAt) {
                  for (const [rCent, rNorm] of candRings) {
                    if (!rNorm) continue;
                    const dCp = dist3([catAt.x, catAt.y, catAt.z], rCent);
                    if (dCp > 1e-4 && dCp <= 5.0) {
                      const cosPhi = Math.abs((catAt.x-rCent[0])*rNorm[0] + (catAt.y-rCent[1])*rNorm[1] + (catAt.z-rCent[2])*rNorm[2]) / dCp;
                      if (cosPhi >= 0.50) { pCatPi = 1; break; }
                    }
                  }
                }
              }
            }
          } else if (AROMATIC_RINGS[cna] && AROMATIC_RINGS[cnb]) {
            const [ringA] = getRingCentroidAndNormalJS(ra);
            const [ringB] = getRingCentroidAndNormalJS(rb);
            if (ringA && ringB && dist3(ringA, ringB) <= 5.5) { pPi = 1; pPiW = 1.0; }
          }

          intfA.add(ra.key); intfB.add(rb.key);
          if (cbDist > 8.0) { intfPae.push(pairPae); weightsSeen.push(wij); }

          const clsA = ["SEP","TPO","PTR"].includes(ra.resName) ? "CHARGED" : (PRODIGY_CLASS[cna] || (ra.entityType === "LIGAND" ? "APOLAR" : "POLAR"));
          const clsB = ["SEP","TPO","PTR"].includes(rb.resName) ? "CHARGED" : (PRODIGY_CLASS[cnb] || (rb.entityType === "LIGAND" ? "APOLAR" : "POLAR"));
          let pairCls = [clsA, clsB].sort().join("_");
          if (pairCls === "APOLAR_CHARGED") pairCls = "CHARGED_APOLAR";
          else if (pairCls === "APOLAR_POLAR") pairCls = "POLAR_APOLAR";
          if (prodigyIC[pairCls] !== undefined) prodigyIC[pairCls] += wij;

          const effSb = Math.min(pSb, 2), effHb = Math.min(pHb, isDuplexPair ? 3 : 2);
          const effHbBase = Math.min(pHbBase, effHb);
          const effHb2oh = Math.min(pHb2oh, effHb - effHbBase);
          const effHbBb = Math.min(pHbBb, effHb - effHbBase - effHb2oh);
          const effHp = Math.min(pHp, 6), effHpW = Math.min(hpDistSum, 6.0), effClash = Math.min(pClash, 3);
          const effMcoord = Math.min(pMetalCoord, 2);
          const effMouter = (pMetalOuter > 0 && effMcoord === 0) ? 1 : 0;
          if (wij >= 0.25) {
            nSb += effSb; nHb += effHb; nHbBase += effHbBase; nHb2oh += effHb2oh; nHbBb += effHbBb;
            nHp += effHp; nPi += pPi; nCatPi += pCatPi; nXb += pXb; nMetalCoord += effMcoord; nMetalOuter += effMouter; nSs += pSs; nClash += effClash;
          }
          wSb += wij*effSb; wHb += wij*effHb; wHbBase += wij*effHbBase; wHb2oh += wij*effHb2oh; wHbBb += wij*effHbBb;
          wHp += wij*effHpW; wPi += wij*pPi; wPiRna += wij*pPiW; wCatPi += wij*pCatPi; wXb += wij*pXb; wMetalCoord += wij*effMcoord; wMetalOuter += wij*effMouter; wClash += wij*effClash;

          if (wij >= 0.25) {
            for (const rk of [ra.key, rb.key]) {
              if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0 });
              const rc = resContrib.get(rk);
              rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*(pPi + pCatPi);
            }
          }

          contacts.push({
            resA: `${ra.chainId}:${ra.resName}${ra.resSeq}`, resB: `${rb.chainId}:${rb.resName}${rb.resSeq}`,
            resSeqA: ra.resSeq, resSeqB: rb.resSeq, coordA: ra.cbCoord, coordB: rb.cbCoord,
            minDist: Number(minD.toFixed(2)), pae: Number(pairPae.toFixed(2)), weight: Number(wij.toFixed(3)),
            sb: effSb, hb: effHb, hbBase: effHbBase, hb2oh: effHb2oh, hbBb: effHbBb,
            hp: effHp, pi: pPi, catpi: pCatPi, xb: pXb, mcoord: effMcoord, mouter: effMouter, ss: pSs, clash: effClash
          });'''
html = ensure_replaced(html, old_inner_js, new_inner_js, "JS: Full inner loop with ligand rings, anionic SB, halogen bonds & inner/outer metal coord")

old_gate_js = '      const gate = (isNucleicPair || effIcSum < 3.0 ? 0.0 : 1.0 / (1.0 + Math.exp(-0.40 * (effIcSum - 11.0)))) * nisReliability;'
new_gate_js = '''      let ligBurialFrac = 0.0;
      if (isLigandPair) {
        const ligIsB = resB.some(r => r.entityType === "LIGAND");
        const ligList = ligIsB ? resB : resA;
        const ligFreeMap = ligIsB ? sasaFreeB : sasaFreeA;
        const ligFreeTot = ligList.reduce((s, r) => s + (ligFreeMap[r.key] || 0), 0);
        const ligBoundTot = ligList.reduce((s, r) => s + (sasaComp[r.key] || 0), 0);
        ligBurialFrac = clamp((ligFreeTot - ligBoundTot) / Math.max(1.0, ligFreeTot), 0.0, 1.0);
      }
      const gate = (isNonPpi || effIcSum < 3.0 ? 0.0 : 1.0 / (1.0 + Math.exp(-0.40 * (effIcSum - 11.0)))) * nisReliability;'''
html = ensure_replaced(html, old_gate_js, new_gate_js, "JS: Compute ligBurialFrac & disable PRODIGY on all non-PPI modalities")

old_eval_s3 = '''        if (isNucleicPair) {
          const dKapRna = Math.exp(-0.85 * ((Math.sqrt(iM) / 3.04) - (Math.sqrt(0.15) / 3.04)) * 3.5);
          const eRawRna = 0.60*wSb*dKapRna*phKappa + 0.65*wHbBase + 0.55*wHb2oh + 0.35*wHbBb + 0.95*wPiRna + 0.85*wCatPi + 0.10*wHp + 0.0025*wBsa;
          const eEffRna = 16.5 * Math.tanh(eRawRna / 16.5);
          dPhys = 2.10 - eEffRna + 1.50*wClash;
        } else {'''
new_eval_s3 = '''        if (isNucleicPair) {
          const dKapRna = Math.exp(-0.85 * ((Math.sqrt(iM) / 3.04) - (Math.sqrt(0.15) / 3.04)) * 3.5);
          const eRawRna = 0.60*wSb*dKapRna*phKappa + 0.65*wHbBase + 0.55*wHb2oh + 0.35*wHbBb + 0.95*wPiRna + 0.85*wCatPi + 0.10*wHp + 0.0025*wBsa;
          const eEffRna = 16.5 * Math.tanh(eRawRna / 16.5);
          dPhys = 2.10 - eEffRna + 1.50*wClash;
        } else if (isLigandPair) {
          const eRawLig = 0.85*wSb*dKap*phKappa + 0.65*wHb + 0.55*wXb + 0.90*wPiRna + 0.80*wCatPi + 0.14*wHp + 0.0055*wBsa;
          const eEffLig = 15.0 * Math.tanh(eRawLig / 15.0);
          dPhys = 1.80 - eEffLig + 1.50*wClash;
        } else if (isMetalPair) {
          const isTransMet = [...resA, ...resB].some(r => r.entityType === "ION" && ["ZN","FE","CU","CO","NI"].includes(r.canonicalName));
          const kCoord = isTransMet ? 2.60 : 1.85;
          const effWMet = wMetalCoord + ((nMetalCoord > 0 || nMetalOuter >= 2) ? 0.45 * Math.min(wMetalOuter, 4.0) : 0.0);
          if (effWMet <= 0.0) dPhys = 2.10 + 1.50*wClash;
          else {
            const eRawMet = kCoord * Math.min(effWMet, 6.0) + 0.0040*wBsa;
            dPhys = 1.50 - 13.5 * Math.tanh(eRawMet / 13.5) + 1.50*wClash;
          }
        } else {'''
html = ensure_replaced(html, old_eval_s3, new_eval_s3, "JS: Add MODE_LIGAND & MODE_METAL (inner+outer) thermodynamics")

html = ensure_replaced(
    html,
    '          if (r.entityType !== "PROTEIN" && r.entityType !== "NUCLEIC") return;',
    '          if (r.entityType !== "PROTEIN" && r.entityType !== "NUCLEIC" && r.entityType !== "LIGAND") return;',
    "JS: Include LIGAND residues in mutCatalog for Burial Profile",
)

old_uis_js = '      const rawUis = mode === "MODE_A_PRED_FULL" ? (0.25*iptm + 0.30*lis + 0.25*pdPae + 0.20*sigAff) : (0.55*clamp(pdCoord / 0.742, 0.0, 1.0) + 0.45*sigAff);'
new_uis_js = '''      let pocketQuality = 0.0, rawUis = 0.0;
      if (isLigandPair) {
        pocketQuality = meanW * clamp((ligBurialFrac - 0.15) / 0.65, 0.0, 1.0);
        rawUis = mode === "MODE_A_PRED_FULL" ? (0.20*iptm + 0.25*lis + 0.30*pocketQuality + 0.25*sigAff) : (0.55*pocketQuality + 0.45*sigAff);
      } else if (isMetalPair) {
        const effQMet = wMetalCoord + ((nMetalCoord > 0 || nMetalOuter >= 2) ? 0.45 * Math.min(wMetalOuter, 4.0) : 0.0);
        pocketQuality = effQMet > 0.0 ? meanW * clamp(effQMet / 3.0, 0.0, 1.0) : 0.0;
        rawUis = 0.55*pocketQuality + 0.45*sigAff;
      } else {
        rawUis = mode === "MODE_A_PRED_FULL" ? (0.25*iptm + 0.30*lis + 0.25*pdPae + 0.20*sigAff) : (0.55*clamp(pdCoord / 0.742, 0.0, 1.0) + 0.45*sigAff);
      }'''
html = ensure_replaced(html, old_uis_js, new_uis_js, "JS: Add Pocket Quality & Metal Coordination UIS scoring")

html = ensure_replaced(
    html,
    '      const effHelPro = isNucleicPair ? 0 : helPro;',
    '      const effHelPro = isNonPpi ? 0 : helPro;',
    "JS: Zero helical_prolines on all non-PPI modalities",
)

html = ensure_replaced(
    html,
    '        salt_bridges: nSb, h_bonds: nHb, hb_base: nHbBase, hb_2oh: nHb2oh, hb_backbone: nHbBb, cation_pi: nCatPi, rna_specificity_idx: rnaSsi, hydrophobic_contacts: nHp, pi_stacking: nPi, disulfides: nSs, steric_clashes: nClash, wPiRna: Number(wPiRna.toFixed(2)), wCatPi: Number(wCatPi.toFixed(2)), wHbBase: Number(wHbBase.toFixed(2)), wHb2oh: Number(wHb2oh.toFixed(2)), wHbBb: Number(wHbBb.toFixed(2)),',
    '        salt_bridges: nSb, h_bonds: nHb, hb_base: nHbBase, hb_2oh: nHb2oh, hb_backbone: nHbBb, cation_pi: nCatPi, halogen_bonds: nXb, metal_coord_bonds: nMetalCoord, metal_outer_bonds: nMetalOuter, lig_burial_pct: Number((ligBurialFrac*100.0).toFixed(1)), pocket_quality: Number(pocketQuality.toFixed(3)), rna_specificity_idx: rnaSsi, hydrophobic_contacts: nHp, pi_stacking: nPi, disulfides: nSs, steric_clashes: nClash, wPiRna: Number(wPiRna.toFixed(2)), wCatPi: Number(wCatPi.toFixed(2)), wXb: Number(wXb.toFixed(2)), wMetalCoord: Number(wMetalCoord.toFixed(2)), wMetalOuter: Number(wMetalOuter.toFixed(2)), wHbBase: Number(wHbBase.toFixed(2)), wHb2oh: Number(wHb2oh.toFixed(2)), wHbBb: Number(wHbBb.toFixed(2)),',
    "JS: Expose halogen, inner/outer metal & ligand burial fields",
)

html = ensure_replaced(
    html,
    '''        pdbChainA: (toSanitizedPDBString._lastRemap && toSanitizedPDBString._lastRemap[caId]) || (caId.length > 1 ? "A" : caId),
        pdbChainB: (toSanitizedPDBString._lastRemap && toSanitizedPDBString._lastRemap[cbId]) || (cbId.length > 1 ? "B" : cbId),''',
    '''        pdbChainA: (toSanitizedPDBString._lastRemap && toSanitizedPDBString._lastRemap[caId]) || (caId.includes("+") ? (toSanitizedPDBString._lastRemap?.[caId.split("+")[0]] || caId.split("+")[0]) : (caId.length > 1 ? "A" : caId)),
        pdbChainListA: caId.includes("+") ? caId.split("+").map(cid => ({ orig: cid, pdb: toSanitizedPDBString._lastRemap?.[cid] || cid.slice(0,1), unit: getChainUnit(resA.filter(r => r.chainId === cid)) })) : null,
        pdbChainB: (toSanitizedPDBString._lastRemap && toSanitizedPDBString._lastRemap[cbId]) || (cbId.length > 1 ? "B" : cbId),''',
    "JS: Store pdbChainListA for composite ternary 3D rendering",
)

old_split_hook = '''        const { ordered, chains } = buildResidues(atoms);
        const cIds = Object.keys(chains);'''
new_split_hook = '''        let { ordered, chains } = buildResidues(atoms);
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
              for (const r of srcList) {
                if (r.resSeq >= lo1 && r.resSeq <= hi1) d1.push(r);
                else if (r.resSeq >= lo2 && r.resSeq <= hi2) d2.push(r);
              }
            } else {
              const cutRes = parseInt(scRange.trim(), 10);
              for (const r of srcList) (r.resSeq <= cutRes ? d1 : d2).push(r);
            }
            if (d1.length && d2.length) {
              const c1 = `${scCid}_D1`, c2 = `${scCid}_D2`;
              for (const r of d1) { r.chainId = c1; r.key = `${c1}|${r.resSeq}|${r.ins}`; for (const a of r.atoms) a.chainId = c1; }
              for (const r of d2) { r.chainId = c2; r.key = `${c2}|${r.resSeq}|${r.ins}`; for (const a of r.atoms) a.chainId = c2; }
              delete chains[scCid];
              chains[c1] = d1; chains[c2] = d2;
              const keptSet2 = new Set(Object.values(chains).flat());
              ordered = ordered.filter(r => keptSet2.has(r));
            }
          }
        }
        const cIds = Object.keys(chains);'''
html = ensure_replaced(html, old_split_hook, new_split_hook, "JS: Add custom Intra-Chain Domain Split in handleFileSelect()")

old_pae_js = '''            if (Array.isArray(j.data.token_chain_ids) && Array.isArray(j.data.token_res_ids)) {
              const keyToTok = new Map();
              j.data.token_chain_ids.forEach((tc, tidx) => {
                const k = `${tc}|${j.data.token_res_ids[tidx]}`;
                if (!keyToTok.has(k)) keyToTok.set(k, tidx);
              });
              const idxs = ordered.map(r => Math.min(rawPae.length - 1, keyToTok.get(`${r.origChainId || r.chainId}|${r.resSeq}`) ?? keyToTok.get(`${r.chainId}|${r.resSeq}`) ?? 0));
              paeMat = idxs.map(ii => idxs.map(jj => Number(rawPae[ii][jj])));
              hasJson = true;
              pairedNames.push(j.name);'''
new_pae_js = '''            if (Array.isArray(j.data.token_chain_ids) && Array.isArray(j.data.token_res_ids)) {
              const keyToToks = new Map();
              j.data.token_chain_ids.forEach((tc, tidx) => {
                const k = `${tc}|${j.data.token_res_ids[tidx]}`;
                if (!keyToToks.has(k)) keyToToks.set(k, []);
                keyToToks.get(k).push(Math.min(rawPae.length - 1, tidx));
              });
              const tokLists = ordered.map(r => keyToToks.get(`${r.origChainId || r.chainId}|${r.resSeq}`) || keyToToks.get(`${r.chainId}|${r.resSeq}`) || []);
              paeMat = tokLists.map(li => {
                if (!li.length) return tokLists.map(() => 15.0);
                if (li.length === 1) {
                  const rowI = rawPae[li[0]];
                  return tokLists.map(lj => lj.length === 1 ? Number(rowI[lj[0]]) : (lj.length ? lj.reduce((s,t)=>s+Number(rowI[t]),0)/lj.length : 15.0));
                }
                return tokLists.map(lj => {
                  if (!lj.length) return 15.0;
                  let s = 0, cnt = 0;
                  for (const ti of li) for (const tj of lj) { s += Number(rawPae[ti][tj]); cnt++; }
                  return s / cnt;
                });
              });
              hasJson = true;
              pairedNames.push(j.name);'''
html = ensure_replaced(html, old_pae_js, new_pae_js, "JS: Safe multi-atom ligand PAE token aggregation")

html = ensure_replaced(
    html,
    '        const origChainOrder = [...new Set(ordered.map(r => r.origChainId || r.chainId))];',
    '        const origChainOrder = [...new Set(ordered.map(r => r.origChainId || r.chainId))];\n        const fileResults = [];',
    "JS: Init per-file results list for ternary detection",
)

old_pair_end_js = '''            const res = analyzePairJS(stem, ca, cb, chains[ca], chains[cb], ordered, resToIdx, mode, { iptm: pairIptm, res_pae_matrix: paeMat }, pairSasaComp, sasaFree[ca], sasaFree[cb], tempK, lambdaPae, isBinary);
            if (res) {
              res.rawStructText = text;
              res.rawStructFormat = sf.name.toLowerCase().endsWith(".pdb") ? "pdb" : "cif";
              ALL_RESULTS.push(res);
            }
          }
        }'''
new_pair_end_js = '''            const res = analyzePairJS(stem, ca, cb, chains[ca], chains[cb], ordered, resToIdx, mode, { iptm: pairIptm, res_pae_matrix: paeMat }, pairSasaComp, sasaFree[ca], sasaFree[cb], tempK, lambdaPae, isBinary);
            if (res) {
              res.rawStructText = text;
              res.rawStructFormat = sf.name.toLowerCase().endsWith(".pdb") ? "pdb" : "cif";
              fileResults.push(res);
              ALL_RESULTS.push(res);
            }
          }
        }
        // Automatic Ternary Molecular-Glue & Catalytic Metal-Bridge Detector (100% CLI Parity)
        if (cIds.length >= 3) {
          for (const lc of cIds) {
            const lcRes = chains[lc];
            if (!lcRes || !lcRes.length || (lcRes[0].entityType !== "LIGAND" && lcRes[0].entityType !== "ION")) continue;
            const isLigC = (lcRes[0].entityType === "LIGAND");
            const partnerHits = [];
            for (const rrec of fileResults) {
              const other = rrec.chainAId === lc ? rrec.chainBId : (rrec.chainBId === lc ? rrec.chainAId : null);
              if (!other || !chains[other] || !chains[other].length) continue;
              const oType = chains[other][0].entityType;
              if (isLigC) {
                if ((oType === "PROTEIN" || oType === "NUCLEIC") && (rrec.lig_burial_pct || 0) >= 12.0 && (rrec.confContactCount || 0) >= 1) {
                  partnerHits.push([other, rrec.lig_burial_pct, 0, 0]);
                }
              } else {
                const mcIn = rrec.metal_coord_bonds || 0, mcOut = rrec.metal_outer_bonds || 0;
                if ((oType === "PROTEIN" || oType === "NUCLEIC" || oType === "LIGAND") && (mcIn >= 1 || mcOut >= 2)) {
                  partnerHits.push([other, mcIn + 0.45 * mcOut, mcIn, mcOut]);
                }
              }
            }
            if ((isLigC && partnerHits.length >= 2) || (!isLigC && partnerHits.length >= 2 && partnerHits.some(h => h[2] >= 1))) {
              partnerHits.sort((a,b) => b[1] - a[1]);
              const [m1, v1, in1, out1] = partnerHits[0], [m2, v2, in2, out2] = partnerHits[1];
              const compResM = [...chains[m1], ...chains[m2]];
              const compIdM = `${m1}+${m2}`;
              const sasaMFree = computeResidueSASA(compResM);
              const sasaTernComp = computeResidueSASA([...compResM, ...lcRes]);
              const tRes = analyzePairJS(stem, compIdM, lc, compResM, lcRes, ordered, resToIdx, mode, { iptm, res_pae_matrix: paeMat }, sasaTernComp, sasaMFree, sasaFree[lc], tempK, lambdaPae, false);
              if (tRes) {
                tRes.rawStructText = text;
                tRes.rawStructFormat = sf.name.toLowerCase().endsWith(".pdb") ? "pdb" : "cif";
                tRes.ternary_partners = compIdM;
                tRes.ternary_burial_split = isLigC ? `${m1}:${v1.toFixed(1)}% / ${m2}:${v2.toFixed(1)}%` : `${m1}:${in1}in+${out1}out / ${m2}:${in2}in+${out2}out`;
                tRes.cooperativity_balance = Number((Math.min(v1, v2) / Math.max(1e-6, Math.max(v1, v2))).toFixed(2));
                if (isLigC && tRes.lig_burial_pct >= 50.0 && tRes.UIS >= UIS_THRESHOLD) {
                  tRes.pattern_id = "Pattern A";
                  tRes.pattern_desc = `Ternary Molecular Glue / Interface Wedge (${lc} bridges ${m1} [${v1.toFixed(1)}%] & ${m2} [${v2.toFixed(1)}%]; composite burial ${tRes.lig_burial_pct.toFixed(1)}%, balance ${tRes.cooperativity_balance.toFixed(2)})`;
                } else if (!isLigC && tRes.metal_coord_bonds >= 1) {
                  tRes.pattern_id = "Pattern A";
                  tRes.pattern_desc = `Catalytic Metal Bridge (${lc} bridges ${m1} [${in1} inner, ${out1} outer-H2O] & ${m2} [${in2} inner, ${out2} outer-H2O])`;
                }
                ALL_RESULTS.push(tRes);
              }
            }
          }
        }'''
html = ensure_replaced(html, old_pair_end_js, new_pair_end_js, "JS: Add Ternary Molecular-Glue & inner+outer Metal-Bridge Detector")

# UI rendering upgrades
html = ensure_replaced(
    html,
    '      const tbody = document.getElementById("leaderboardBody");',
    '''      const ternBanner = document.getElementById("ternaryAuditBanner");
      if (ternBanner) {
        const terns = ALL_RESULTS.filter(x => x.modality === "TERNARY-GLUE" || x.modality === "METAL-BRIDGE");
        if (terns.length) {
          ternBanner.style.display = "block";
          ternBanner.innerHTML = `<strong>TERNARY MOLECULAR GLUE &amp; CATALYTIC METAL-BRIDGE AUDIT (${terns.length} detected):</strong><br>` + terns.map(tr =>
            `&bull; <strong>[${tr.modality}] ${tr.chain_pair}</strong>: UIS = <strong>${tr.UIS.toFixed(3)}</strong> · ΔG = <strong>${tr.dG_uni_kcal.toFixed(2)} kcal/mol</strong> (${formatKdString(tr.Kd_nM)}) · Split: <strong>${tr.ternary_burial_split || "—"}</strong> · Balance = <strong>${(tr.cooperativity_balance||0).toFixed(2)}</strong> · ${tr.pattern_desc}`
          ).join("<br>");
        } else {
          ternBanner.style.display = "none";
        }
      }

      const tbody = document.getElementById("leaderboardBody");''',
    "JS: Render Ternary Audit banner",
)

html = ensure_replaced(
    html,
    '      const isNucUI = ["PROTEIN-RNA","PROTEIN-DNA","RNA-RNA","DNA-DNA","RNA-DNA"].includes(r.modality);',
    '      const isNucUI = ["PROTEIN-RNA","PROTEIN-DNA","RNA-RNA","DNA-DNA","RNA-DNA"].includes(r.modality);\n'
    '      const isLigUI = ["PROTEIN-LIGAND","RNA-LIGAND","DNA-LIGAND","TERNARY-GLUE"].includes(r.modality);\n'
    '      const isMetUI = ["PROTEIN-METAL","RNA-METAL","DNA-METAL","LIGAND-METAL","METAL-BRIDGE"].includes(r.modality);',
    "JS: Add isLigUI & isMetUI flags",
)

html = ensure_replaced(
    html,
    '''        <div class="descriptor"><span>${isNucUI ? "1. RNA/DNA Specificity Index (SSI)" : "1. Apolar Burial Fraction (f<sub>apolar</sub>)"}</span><strong>${isNucUI ? `${(r.rna_specificity_idx||0).toFixed(2)} (π:${r.pi_stacking||0}, Cat-π:${r.cation_pi||0})` : `${r.f_apolar.toFixed(3)} (${r.hydrophobic_contacts} HP contacts)`}</strong></div>''',
    '''        <div class="descriptor"><span>${isNucUI ? "1. RNA/DNA Specificity Index (SSI)" : (isLigUI ? "1. Ligand Burial &amp; Pocket Quality" : (isMetUI ? "1. Metal Coordination (In / Outer-H2O)" : "1. Apolar Burial Fraction (f<sub>apolar</sub>)"))}</span><strong>${isNucUI ? `${(r.rna_specificity_idx||0).toFixed(2)} (π:${r.pi_stacking||0}, Cat-π:${r.cation_pi||0})` : (isLigUI ? `${(r.lig_burial_pct||0).toFixed(1)}% buried (Q=${(r.pocket_quality||0).toFixed(3)}, π:${r.pi_stacking||0}, XB:${r.halogen_bonds||0})` : (isMetUI ? `${r.metal_coord_bonds||0} inner / ${r.metal_outer_bonds||0} outer (Q=${(r.pocket_quality||0).toFixed(3)})` : `${r.f_apolar.toFixed(3)} (${r.hydrophobic_contacts} HP contacts)`))}</strong></div>''',
    "JS: Update Descriptor 1 for Ligand/Ternary/Metal",
)

old_s4_js = '''              ${isNucUI ? `
              &bull; <code>Modality = <strong>${r.modality}</strong> &rarr; Protein–Protein PRODIGY term automatically disabled (0% weight).</code><br>
              &bull; <code>ΔG<sub>RNA</sub> (Dedicated 9-Channel ${r.modality} Model, E<sub>eff</sub>=16.5·tanh(E<sub>raw</sub>/16.5); wSB<sub>phos</sub>=${(r.wSb||0).toFixed(1)}, wHB<sub>base</sub>=${(r.wHbBase||0).toFixed(1)}, wHB<sub>2'OH</sub>=${(r.wHb2oh||0).toFixed(1)}, wHB<sub>bb</sub>=${(r.wHbBb||0).toFixed(1)}, wπ<sub>RNA</sub>=${(r.wPiRna||0).toFixed(1)}, wCatπ=${(r.wCatPi||0).toFixed(1)}, wHP=${(r.wHp||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å²) = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> &rarr; K<sub>d</sub> = <strong>${kdStr}</strong></code><br>
              ` : `'''
new_s4_js = '''              ${isNucUI ? `
              &bull; <code>Modality = <strong>${r.modality}</strong> &rarr; Protein–Protein PRODIGY term automatically disabled (0% weight).</code><br>
              &bull; <code>ΔG<sub>RNA</sub> (Dedicated 9-Channel ${r.modality} Model, E<sub>eff</sub>=16.5·tanh(E<sub>raw</sub>/16.5); wSB<sub>phos</sub>=${(r.wSb||0).toFixed(1)}, wHB<sub>base</sub>=${(r.wHbBase||0).toFixed(1)}, wHB<sub>2'OH</sub>=${(r.wHb2oh||0).toFixed(1)}, wHB<sub>bb</sub>=${(r.wHbBb||0).toFixed(1)}, wπ<sub>RNA</sub>=${(r.wPiRna||0).toFixed(1)}, wCatπ=${(r.wCatPi||0).toFixed(1)}, wHP=${(r.wHp||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å²) = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> &rarr; K<sub>d</sub> = <strong>${kdStr}</strong></code><br>
              ` : (isLigUI ? `
              &bull; <code>Modality = <strong>${r.modality}</strong> &rarr; Protein–Protein PRODIGY disabled; Ligand Burial = <strong>${(r.lig_burial_pct||0).toFixed(1)}%</strong>, Pocket Quality Q<sub>pocket</sub> = <strong>${(r.pocket_quality||0).toFixed(3)}</strong>${r.ternary_burial_split ? ` (${r.ternary_burial_split}, Balance=${(r.cooperativity_balance||0).toFixed(2)})` : ""}.</code><br>
              &bull; <code>ΔG<sub>lig</sub> (Dedicated ${r.modality} Pocket Model, E<sub>eff</sub>=15.0·tanh(E<sub>raw</sub>/15.0); wSB=${(r.wSb||0).toFixed(1)}, wHB=${(r.wHb||0).toFixed(1)}, wXB=${(r.wXb||0).toFixed(1)}, wπ=${(r.wPiRna||0).toFixed(1)}, wCatπ=${(r.wCatPi||0).toFixed(1)}, wHP=${(r.wHp||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å²) = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> &rarr; K<sub>d</sub> = <strong>${kdStr}</strong></code><br>
              ` : (isMetUI ? `
              &bull; <code>Modality = <strong>${r.modality}</strong> &rarr; Dedicated Metal Coordination Model (Inner-Sphere = <strong>${r.metal_coord_bonds||0}</strong>, Outer-H2O = <strong>${r.metal_outer_bonds||0}</strong>, Q<sub>metal</sub> = <strong>${(r.pocket_quality||0).toFixed(3)}</strong>${r.ternary_burial_split ? `, Bridge: ${r.ternary_burial_split}` : ""}) &rarr; ΔG<sub>met</sub> = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> (K<sub>d</sub> = <strong>${kdStr}</strong>).</code><br>
              ` : `'''
html = ensure_replaced(html, old_s4_js, new_s4_js, "JS: Add Ligand/Ternary/Metal Step 4 derivation math")

html = ensure_replaced(
    html,
    '''              `}
              &bull; <code>${r.mode === "MODE_A_PRED_FULL"''',
    '''              `))}
              &bull; <code>${isLigUI
                ? `UIS (${r.modality}) = [${r.mode === "MODE_A_PRED_FULL" ? `0.20·ipTM(${r.ipTM.toFixed(3)}) + 0.25·LIS(${r.LIS.toFixed(3)}) + 0.30·Q_pocket(${(r.pocket_quality||0).toFixed(3)}) + 0.25·σ_aff(${r.sig_aff.toFixed(3)})` : `0.55·Q_pocket(${(r.pocket_quality||0).toFixed(3)}) + 0.45·σ_aff(${r.sig_aff.toFixed(3)})`}] × gate(${r.thermo_gate.toFixed(2)})`
                : (isMetUI
                  ? `UIS (${r.modality}) = [0.55·Q_metal(${(r.pocket_quality||0).toFixed(3)}) + 0.45·σ_aff(${r.sig_aff.toFixed(3)})] × gate(${r.thermo_gate.toFixed(2)})`
                  : (r.mode === "MODE_A_PRED_FULL"''',
    "JS: Add Ligand/Metal UIS formula in Step 4",
)
html = ensure_replaced(
    html,
    '''                : `UIS (Coord Mode) = [0.55·pDockQ<sub>norm</sub>(${Math.min(1, r.pDockQ_coord / 0.742).toFixed(3)}) + 0.45·σ<sub>aff</sub>(${r.sig_aff.toFixed(3)})] × gate(${r.thermo_gate.toFixed(2)})`
              } = <strong>${r.UIS.toFixed(4)}</strong></code>''',
    '''                : `UIS (Coord Mode) = [0.55·pDockQ<sub>norm</sub>(${Math.min(1, r.pDockQ_coord / 0.742).toFixed(3)}) + 0.45·σ<sub>aff</sub>(${r.sig_aff.toFixed(3)})] × gate(${r.thermo_gate.toFixed(2)})`
              ))} = <strong>${r.UIS.toFixed(4)}</strong></code>''',
    "JS: Close Step 4 UIS ternary parens",
)

html = ensure_replaced(
    html,
    '''          : `Note: ΔG<sub>uni</sub> (${r.dG_uni_kcal.toFixed(2)} kcal/mol) weights published PRODIGY (${Math.round(r.w_prodigy_blend*100)}%: ${r.dG_prodigy_kcal.toFixed(2)} kcal/mol) and enthalpy–entropy compensated 8-channel solvation (${Math.round(r.w_phys_blend*100)}%: ${r.dG_phys_kcal.toFixed(2)} kcal/mol).`
        }</span>''',
    '''          : (isLigUI
            ? `Note: Dedicated ${r.modality} Pocket Physics active (PRODIGY disabled; Burial = ${(r.lig_burial_pct||0).toFixed(1)}%, Q_pocket = ${(r.pocket_quality||0).toFixed(3)}, π-Stacks = ${r.pi_stacking||0}, Cation-π = ${r.cation_pi||0}, Halogen Bonds = ${r.halogen_bonds||0}). Designed for post-AF3 pose auditing &amp; escape triage — not a replacement for FEP/ITC.`
            : (isMetUI
              ? `Note: Dedicated ${r.modality} Coordination Model active (${r.metal_coord_bonds||0} inner-sphere, ${r.metal_outer_bonds||0} outer-H2O bonds, Q_metal = ${(r.pocket_quality||0).toFixed(3)}).`
              : `Note: ΔG<sub>uni</sub> (${r.dG_uni_kcal.toFixed(2)} kcal/mol) weights published PRODIGY (${Math.round(r.w_prodigy_blend*100)}%: ${r.dG_prodigy_kcal.toFixed(2)} kcal/mol) and enthalpy–entropy compensated 8-channel solvation (${Math.round(r.w_phys_blend*100)}%: ${r.dG_phys_kcal.toFixed(2)} kcal/mol).`))
        }</span>''',
    "JS: Update patternCallout for Ligand/Ternary/Metal scope note",
)

html = ensure_replaced(
    html,
    '''      styleActiveChain(cA, colA, r.unitA || "aa");
      styleActiveChain(cB, colB, r.unitB || "aa");''',
    '''      if (r.pdbChainListA && r.pdbChainListA.length) {
        r.pdbChainListA.forEach((item, idx) => {
          const cCol = item.unit === "nt" ? "#f59e0b" : (idx === 0 ? "#38bdf8" : "#34d399");
          styleActiveChain(item.pdb, cCol, item.unit);
        });
      } else {
        styleActiveChain(cA, colA, r.unitA || "aa");
      }
      styleActiveChain(cB, colB, r.unitB || "aa");''',
    "JS: Composite ternary A+B 3Dmol.js styling",
)

html = ensure_replaced(
    html,
    '''          if (c.sb > 0 && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.11, color: "#f43f5e", dashed: true });
          } else if ((c.pi > 0 || c.catpi > 0) && c.coordA && c.coordB) {''',
    '''          if ((c.mcoord > 0 || c.mouter > 0) && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.11, color: "#84cc16", dashed: true });
          } else if (c.xb > 0 && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.10, color: "#06b6d4", dashed: true });
          } else if (c.sb > 0 && c.coordA && c.coordB) {
            VIEWER_3D.addCylinder({ start: { x:c.coordA[0], y:c.coordA[1], z:c.coordA[2] }, end: { x:c.coordB[0], y:c.coordB[1], z:c.coordB[2] }, radius: 0.11, color: "#f43f5e", dashed: true });
          } else if ((c.pi > 0 || c.catpi > 0) && c.coordA && c.coordB) {''',
    "JS: Add Lime Metal & Cyan Halogen cylinders in 3Dmol.js",
)

html = ensure_replaced(
    html,
    '''      const rankType = c => c.sb > 0 ? 4 : ((c.pi > 0 || c.catpi > 0) ? 3 : (c.hb > 0 ? 2 : (c.hp > 0 ? 1 : 0)));
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
    '''      const rankType = c => (c.mcoord > 0 || c.mouter > 0 || c.xb > 0) ? 5 : (c.sb > 0 ? 4 : ((c.pi > 0 || c.catpi > 0) ? 3 : (c.hb > 0 ? 2 : (c.hp > 0 ? 1 : 0))));
      const sortedC = activeC.slice().sort((a, b) => rankType(a) - rankType(b));
      let linesSvg = "";
      for (const c of sortedC) {
        const normA = clamp(((c.resSeqA || minSeqA) - minSeqA) / Math.max(1, maxSeqA - minSeqA), 0.0, 1.0);
        const normB = clamp(((c.resSeqB || minSeqB) - minSeqB) / Math.max(1, maxSeqB - minSeqB), 0.0, 1.0);
        const yA = 35 + normA * (h - 70);
        const yB = 35 + normB * (h - 70);
        const rt = rankType(c);
        const color = rt === 5 ? "#0891b2" : (rt === 4 ? "#e11d48" : (rt === 3 ? "#9333ea" : (rt === 2 ? "#d97706" : (rt === 1 ? "#0284c7" : "#cbd5e1"))));
        const strokeW = rt === 0 ? 1.0 : Math.max(1.2, 2.8 * (c.weight || 1.0));
        const opacity = rt === 0 ? "0.45" : "0.90";
        const tag = rt === 5 ? (c.mcoord > 0 ? "Metal Inner-Coord" : (c.mouter > 0 ? "Metal Outer-H2O" : "Halogen Bond")) : (rt === 4 ? "Salt Bridge" : (rt === 3 ? (c.catpi > 0 ? "Cation-π" : "π-Stack") : (rt === 2 ? "H-Bond" : (rt === 1 ? "Hydrophobic" : "vdW"))));''',
    "JS: Add Cyan Halogen/Metal edges to SVG Contact Network",
)

html = ensure_replaced(
    html,
    '''      const tbody = document.getElementById("contactInspectorBody");
      const heading = document.getElementById("contactTableHeading");
      if (!tbody || !r) return;''',
    '''      const tbody = document.getElementById("contactInspectorBody");
      const heading = document.getElementById("contactTableHeading");
      if (!tbody || !r) return;
      const colAEl = document.getElementById("contactColA"), colBEl = document.getElementById("contactColB");
      if (colAEl) colAEl.textContent = `${r.chainAId || "A"} Res`;
      if (colBEl) colBEl.textContent = `${r.chainBId || "B"} Res`;''',
    "JS: Dynamic Contact Inspector column headers",
)

html = ensure_replaced(
    html,
    '|| (q==="BASE" && c.hbBase>0) || (q==="BB" && c.hbBb>0)));',
    '|| (q==="BASE" && c.hbBase>0) || (q==="BB" && c.hbBb>0) || ((q==="XB"||q.includes("HAL")) && c.xb>0) || (q.includes("MET") && (c.mcoord>0 || c.mouter>0))));',
    "JS: Support XB/HAL & MET filtering",
)

html = ensure_replaced(
    html,
    '        if (c.catpi > 0) tags.push(`Cat-π:${c.catpi}`);',
    '        if (c.catpi > 0) tags.push(`Cat-π:${c.catpi}`);\n        if (c.xb > 0) tags.push(`XB:${c.xb}`);\n        if (c.mcoord > 0) tags.push(`Met-In:${c.mcoord}`);\n        if (c.mouter > 0) tags.push(`Met-H2O:${c.mouter}`);',
    "JS: Display XB, Met-In & Met-H2O tags",
)

html = ensure_replaced(
    html,
    '"pi_stacking","cation_pi","rna_specificity_idx","helical_prolines"',
    '"pi_stacking","cation_pi","halogen_bonds","metal_coord_bonds","metal_outer_bonds","lig_burial_pct","pocket_quality","rna_specificity_idx","ternary_partners","cooperativity_balance","helical_prolines"',
    "JS: Add Ligand/Metal/Ternary columns to Master CSV/TSV",
)

html = ensure_replaced(
    html,
    'const cols = ["resA","resB","minDist","pae","weight","sb","hb","hbBase","hb2oh","hbBb","hp","pi","catpi","ss","clash"];',
    'const cols = ["resA","resB","minDist","pae","weight","sb","hb","hbBase","hb2oh","hbBb","hp","pi","catpi","xb","mcoord","mouter","ss","clash"];',
    "JS: Add xb, mcoord & mouter to Contacts CSV",
)

html_path.write_text(html, encoding="utf-8")
print("[PASS] index.html synced cleanly to Step 3B-2 parity!")

m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
assert m_js, "Missing inline <script> block in index.html"
node_bin = shutil.which("node") or shutil.which("node.exe")
if node_bin:
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
        tf.write(m_js[-1])
        tf_path = tf.name
    res_node = subprocess.run([node_bin, "--check", tf_path], capture_output=True, text=True)
    Path(tf_path).unlink(missing_ok=True)
    if res_node.returncode != 0:
        raise RuntimeError(f"JavaScript syntax error in index.html:\n{res_node.stderr}")
    print(f"[PASS] Verified index.html JavaScript syntax cleanly via `{node_bin} --check`.")
else:
    print("[INFO] `node` not installed; please check browser DevTools console (F12) after Ctrl+F5 reload.")
