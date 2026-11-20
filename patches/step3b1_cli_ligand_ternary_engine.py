#!/usr/bin/env python3
"""
Phase 2 - Step 3B-1 (AST-Anchored & Idempotent — Universal Ligand, Halogen-Bond,
Metal-Coordination, Ternary Molecular-Glue & Custom Intra-Chain Domain Slicer):
  1. Adds MODIFIED_NUC_MAP (I, PSU, 6MA, 5MC, 5MU, OMC, OMG, A2M, G7M, etc.) so modified viral RNAs
     stay inside RNA chains instead of splitting out as fake ligands.
  2. Fixes _fmt_segs regex bug (3RJ401 -> 401) and safe multi-chain composite loci formatting.
  3. Fixes AF3 multi-atom ligand PAE token aggregation (averages all heavy-atom tokens; unmapped -> 15.0 A).
  4. Allows 3 H-bonds per base pair in RNA-RNA / DNA-DNA / RNA-DNA duplexes (G=C pairs).
  5. Adds universal 3D covalent-graph planar ring detector + directional halogen-bond detector (>= 120 deg).
  6. Adds dedicated MODE_LIGAND (pocket burial Q_pocket + calibrated dG_lig) & MODE_METAL (inner-sphere
     coordination chemistry for Zn/Mg/Mn/Ca/Fe/Cu/Co/Ni), disabling PRODIGY on all non-PPI modalities.
  7. Adds automatic TERNARY-GLUE (5ZC9 A+B <-> B_RCG600) & METAL-BRIDGE (A+A_ANP500 <-> A_MG501)
     composite cleft detection and scoring.
  8. Adds --split-chain CLI flag (same-protein domain-domain analysis) + contact-trough monomer splitter.
  9. Runs all 6 real benchmarks (fold_vip1_bzip61, 5ZC9, 1URN, 4WNP, 7DIY, and 5ZC9 --split-chain).
"""

import ast
import re
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
src = cli_path.read_text(encoding="utf-8")


def ensure_replaced(text: str, old: str, new: str, label: str, count: int = 1) -> str:
    if new in text:
        print(f"  [OK-ALREADY] {label}")
        return text
    if old in text:
        print(f"  [APPLIED]    {label}")
        return text.replace(old, new, count)
    raise AssertionError(f"Anchor mismatch for [{label}]:\n{old[:180]!r}")


def get_func_node(source: str, func_name: str) -> ast.FunctionDef:
    tree_cur = ast.parse(source)
    for node in tree_cur.body:
        if isinstance(node, ast.FunctionDef) and node.name == func_name:
            return node
    raise AssertionError(f"Function `{func_name}` not found in unippi_cli.py")


# Inspect exact sasa parameter names of analyze_interface_pair via AST
aip_node = get_func_node(src, "analyze_interface_pair")
aip_args = [a.arg for a in aip_node.args.args]
print(f"[INFO] analyze_interface_pair signature args: {aip_args[:11]}")
sasa_comp_arg, sasa_free_a_arg, sasa_free_b_arg = aip_args[7], aip_args[8], aip_args[9]

# ------------------------------------------------------------------------------
# 1. MODIFIED_NUC_MAP & HALOGEN CONSTANTS
# ------------------------------------------------------------------------------
old_nuc_map = 'NUC_TO_ONE = {"A": "a", "C": "c", "G": "g", "U": "u", "T": "t", "DA": "a", "DC": "c", "DG": "g", "DT": "t", "ADE": "a", "CYT": "c", "GUA": "g", "URA": "u", "THY": "t"}'
new_nuc_map = old_nuc_map + '''
MODIFIED_NUC_MAP = {
    "I": "G", "INO": "G", "PSU": "U", "H2U": "U", "5MU": "U", "4SU": "U", "OMU": "U",
    "6MA": "A", "1MA": "A", "A2M": "A", "M6A": "A", "MIA": "A", "T6A": "A",
    "5MC": "C", "OMC": "C", "4AC": "C", "3MC": "C",
    "OMG": "G", "7MG": "G", "G7M": "G", "2MG": "G", "M2G": "G", "1MG": "G", "YYG": "G",
    "DI": "DG", "5CM": "DC", "BRU": "DT"
}
HALOGEN_ELEMENTS = {"F", "CL", "BR", "I"}'''
src = ensure_replaced(src, old_nuc_map, new_nuc_map, "Add MODIFIED_NUC_MAP & HALOGEN_ELEMENTS")

old_canon = 'self.canonical_name = MODIFIED_AA_MAP.get(res_name, res_name)'
new_canon = 'self.canonical_name = MODIFIED_AA_MAP.get(res_name, MODIFIED_NUC_MAP.get(res_name, res_name))'
if new_canon not in src and old_canon in src:
    src = ensure_replaced(src, old_canon, new_canon, "Residue canonical_name supports MODIFIED_NUC_MAP")
else:
    src_sub, n_rep = re.subn(
        r'MODIFIED_AA_MAP\.get\((\w+),\s*\1\)',
        r'MODIFIED_AA_MAP.get(\1, MODIFIED_NUC_MAP.get(\1, \1))',
        src,
        count=1,
    )
    if n_rep == 1 and new_canon not in src:
        src = src_sub
        print("  [APPLIED]    Residue canonical_name supports MODIFIED_NUC_MAP (regex)")
    else:
        print("  [OK-ALREADY] Residue canonical_name supports MODIFIED_NUC_MAP")

# ------------------------------------------------------------------------------
# 2. SAFE MULTI-ATOM LIGAND PAE AGGREGATION IN parse_companion_json()
# ------------------------------------------------------------------------------
old_pae_map = '''                key_to_tok: Dict[Tuple[str, int], int] = {}
                for tidx, (tc, tr) in enumerate(zip(t_chains, t_res)):
                    k = (str(tc), int(tr))
                    if k not in key_to_tok:
                        key_to_tok[k] = tidx

                res_tok_indices = [
                    min(n_tok - 1, key_to_tok.get((getattr(r, "orig_chain_id", r.chain_id), r.res_seq), key_to_tok.get((r.chain_id, r.res_seq), 0)))
                    for r in ordered_residues
                ]
                aligned_pae: List[List[float]] = []
                for idx_i in res_tok_indices:
                    row_i = raw_pae[idx_i]
                    aligned_pae.append([float(row_i[idx_j]) for idx_j in res_tok_indices])'''

new_pae_map = '''                key_to_toks: Dict[Tuple[str, int], List[int]] = {}
                for tidx, (tc, tr) in enumerate(zip(t_chains, t_res)):
                    k = (str(tc), int(tr))
                    key_to_toks.setdefault(k, []).append(min(n_tok - 1, tidx))

                res_tok_lists: List[List[int]] = [
                    key_to_toks.get(
                        (getattr(r, "orig_chain_id", r.chain_id), r.res_seq),
                        key_to_toks.get((r.chain_id, r.res_seq), []),
                    )
                    for r in ordered_residues
                ]
                aligned_pae: List[List[float]] = []
                for li in res_tok_lists:
                    if not li:
                        aligned_pae.append([15.0] * len(res_tok_lists))
                    elif len(li) == 1:
                        row_i = raw_pae[li[0]]
                        aligned_pae.append([
                            (float(row_i[lj[0]]) if len(lj) == 1 else (sum(float(row_i[t]) for t in lj) / len(lj) if lj else 15.0))
                            for lj in res_tok_lists
                        ])
                    else:
                        row_out = []
                        for lj in res_tok_lists:
                            if not lj:
                                row_out.append(15.0)
                            else:
                                vals = [float(raw_pae[ti][tj]) for ti in li for tj in lj]
                                row_out.append(sum(vals) / len(vals))
                        aligned_pae.append(row_out)'''
src = ensure_replaced(src, old_pae_map, new_pae_map, "Safe multi-atom ligand PAE token aggregation (no token-0 fallback)")

# ------------------------------------------------------------------------------
# 3. UNIVERSAL 3D LIGAND PLANAR RING DETECTOR + DIRECTIONAL HALOGEN-BOND DETECTOR
# ------------------------------------------------------------------------------
old_det_mod = '''def detect_pair_modality(ca_id: str, cb_id: str, res_a: List[Residue], res_b: List[Residue]) -> str:
    if ca_id.endswith("_N") and cb_id.endswith("_C"):
        return "INTRA-DOMAIN (MONOMER)"'''

new_det_mod = '''def detect_ligand_planar_rings(res: Residue) -> List[Tuple[Vec3, Vec3]]:
    """Universal 5-/6-membered planar aromatic ring detector directly from 3D heavy-atom coordinates."""
    if hasattr(res, "_cached_planar_rings"):
        return res._cached_planar_rings
    atoms = [a for a in res.atoms if a.element in ("C", "N", "O", "S")]
    n = len(atoms)
    if n < 5:
        res._cached_planar_rings = []
        return []
    adj: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i in range(n):
        ci = atoms[i].coord
        for j in range(i + 1, n):
            d = math.dist(ci, atoms[j].coord)
            if 1.15 <= d <= 1.68:
                adj[i].append(j)
                adj[j].append(i)
    rings: Set[Tuple[int, ...]] = set()
    for start in range(n):
        stack = [(start, [start])]
        while stack:
            curr, path = stack.pop()
            if len(path) > 6:
                continue
            for nxt in adj[curr]:
                if nxt == start and len(path) in (5, 6):
                    rings.add(tuple(sorted(path)))
                elif nxt > start and nxt not in path and len(path) < 6:
                    stack.append((nxt, path + [nxt]))
    planar: List[Tuple[Vec3, Vec3]] = []
    for r_idx in sorted(rings):
        pts = [atoms[k].coord for k in r_idx]
        m = len(pts)
        cent = (sum(p[0] for p in pts) / m, sum(p[1] for p in pts) / m, sum(p[2] for p in pts) / m)
        p1, p2, p3 = pts[0], pts[m // 3], pts[(2 * m) // 3]
        ux, uy, uz = p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]
        vx, vy, vz = p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        nl = math.sqrt(nx * nx + ny * ny + nz * nz)
        if nl < 1e-6:
            continue
        norm = (nx / nl, ny / nl, nz / nl)
        max_dev = max(abs((p[0] - cent[0]) * norm[0] + (p[1] - cent[1]) * norm[1] + (p[2] - cent[2]) * norm[2]) for p in pts)
        if max_dev <= 0.15:
            planar.append((cent, norm))
    res._cached_planar_rings = planar
    return planar

def get_all_rings_for_residue(res: Residue) -> List[Tuple[Vec3, Optional[Vec3]]]:
    if res.entity_type == "LIGAND":
        return detect_ligand_planar_rings(res)
    cent, norm = get_ring_centroid_and_normal(res)
    return [(cent, norm)] if cent is not None else []

def count_directional_halogen_bonds(ra: Residue, rb: Residue) -> int:
    """Directional C-X...O/N/S halogen bond detector (angle >= 120 deg, d <= 3.70 A; 3.30 A for F)."""
    xb = 0
    for src_r, dst_r in ((ra, rb), (rb, ra)):
        hal_atoms = [a for a in src_r.atoms if a.element in HALOGEN_ELEMENTS]
        if not hal_atoms:
            continue
        c_atoms = [a for a in src_r.atoms if a.element == "C"]
        acc_atoms = [a for a in dst_r.atoms if a.element in ("O", "N", "S")]
        if not c_atoms or not acc_atoms:
            continue
        for hx in hal_atoms:
            c_bond = None
            best_cd = 9.0
            for ca in c_atoms:
                dc = math.dist(ca.coord, hx.coord)
                if 1.20 <= dc <= 2.25 and dc < best_cd:
                    best_cd = dc
                    c_bond = ca
            if c_bond is None:
                continue
            d_cut = 3.30 if hx.element == "F" else 3.70
            vxc = (c_bond.coord[0] - hx.coord[0], c_bond.coord[1] - hx.coord[1], c_bond.coord[2] - hx.coord[2])
            l_xc = math.sqrt(vxc[0] ** 2 + vxc[1] ** 2 + vxc[2] ** 2)
            if l_xc < 1e-4:
                continue
            for acc in acc_atoms:
                da = math.dist(hx.coord, acc.coord)
                if 2.20 <= da <= d_cut:
                    vxa = (acc.coord[0] - hx.coord[0], acc.coord[1] - hx.coord[1], acc.coord[2] - hx.coord[2])
                    cos_ang = (vxc[0] * vxa[0] + vxc[1] * vxa[1] + vxc[2] * vxa[2]) / (l_xc * da)
                    if cos_ang <= -0.50:
                        xb += 1
                        break
    return min(xb, 2)

def detect_pair_modality(ca_id: str, cb_id: str, res_a: List[Residue], res_b: List[Residue]) -> str:
    if "+" in ca_id:
        types_b = {r.entity_type for r in res_b}
        if "LIGAND" in types_b:
            return "TERNARY-GLUE"
        if "ION" in types_b:
            return "METAL-BRIDGE"
    if (ca_id.endswith("_N") and cb_id.endswith("_C")) or ("_D1" in ca_id and "_D2" in cb_id):
        return "INTRA-DOMAIN (MONOMER)"'''
src = ensure_replaced(src, old_det_mod, new_det_mod, "Add universal 3D ligand ring detector, directional X-bond detector & ternary modality")

# ------------------------------------------------------------------------------
# 4. UPGRADE analyze_interface_pair()
# ------------------------------------------------------------------------------
old_bmin = '''    b_min_a = min((r.mean_bfactor for r in residues_a), default=10.0)
    b_min_b = min((r.mean_bfactor for r in residues_b), default=10.0)
    b_min_chain = {chain_a_id: b_min_a, chain_b_id: b_min_b}'''
new_bmin = '''    b_min_chain: Dict[str, float] = {}
    for _r in (residues_a + residues_b):
        _c = _r.chain_id
        if _c not in b_min_chain or _r.mean_bfactor < b_min_chain[_c]:
            b_min_chain[_c] = _r.mean_bfactor'''
src = ensure_replaced(src, old_bmin, new_bmin, "Group MODE_C b_min_chain by actual residue chain_id (ternary-safe)")

old_mod_flags = '''    modality = detect_pair_modality(chain_a_id, chain_b_id, residues_a, residues_b)
    is_nucleic_pair = modality in ("PROTEIN-RNA", "PROTEIN-DNA", "RNA-RNA", "DNA-DNA", "RNA-DNA")
    n_sb, n_hb, n_hp, n_pi, n_ss, n_clash, n_cb8 = 0, 0, 0, 0, 0, 0, 0
    n_hb_base, n_hb_2oh, n_hb_bb, n_catpi = 0, 0, 0, 0
    w_sb, w_hb, w_hp, w_pi, w_clash = 0.0, 0.0, 0.0, 0.0, 0.0
    w_hb_base, w_hb_2oh, w_hb_bb, w_pi_rna, w_catpi = 0.0, 0.0, 0.0, 0.0, 0.0'''

new_mod_flags = '''    modality = detect_pair_modality(chain_a_id, chain_b_id, residues_a, residues_b)
    is_nucleic_pair = modality in ("PROTEIN-RNA", "PROTEIN-DNA", "RNA-RNA", "DNA-DNA", "RNA-DNA")
    is_duplex_pair = modality in ("RNA-RNA", "DNA-DNA", "RNA-DNA")
    is_ligand_pair = modality in ("PROTEIN-LIGAND", "RNA-LIGAND", "DNA-LIGAND", "TERNARY-GLUE")
    is_metal_pair = modality in ("PROTEIN-METAL", "RNA-METAL", "DNA-METAL", "LIGAND-METAL", "METAL-BRIDGE")
    is_non_ppi = is_nucleic_pair or is_ligand_pair or is_metal_pair
    n_sb, n_hb, n_hp, n_pi, n_ss, n_clash, n_cb8 = 0, 0, 0, 0, 0, 0, 0
    n_hb_base, n_hb_2oh, n_hb_bb, n_catpi, n_xb, n_metal_coord = 0, 0, 0, 0, 0, 0
    w_sb, w_hb, w_hp, w_pi, w_clash = 0.0, 0.0, 0.0, 0.0, 0.0
    w_hb_base, w_hb_2oh, w_hb_bb, w_pi_rna, w_catpi, w_xb, w_metal_coord = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0'''
src = ensure_replaced(src, old_mod_flags, new_mod_flags, "Add modality flags & halogen/metal counters")

old_met_check = '''                    is_met_a = (ra.entity_type in ("ION", "LIGAND") and atom_a.element in CATALYTIC_METALS)
                    is_met_b = (rb.entity_type in ("ION", "LIGAND") and atom_b.element in CATALYTIC_METALS)
                    if d <= 2.8 and ((is_met_a and atom_b.element in ("O", "N", "S")) or (is_met_b and atom_a.element in ("O", "N", "S"))):
                        pair_sb += 1
                        continue'''

new_met_check = '''                    is_met_a = (ra.entity_type in ("ION", "LIGAND") and atom_a.element in CATALYTIC_METALS)
                    is_met_b = (rb.entity_type in ("ION", "LIGAND") and atom_b.element in CATALYTIC_METALS)
                    if is_met_a or is_met_b:
                        m_el = atom_a.element if is_met_a else atom_b.element
                        other_el = atom_b.element if is_met_a else atom_a.element
                        d_mcut = 2.65 if m_el in ("ZN", "FE", "CU", "CO", "NI") else 2.75
                        if d <= d_mcut and other_el in ("O", "N", "S"):
                            pair_sb += 1
                            pair_metal_coord += 1
                            continue'''
src = ensure_replaced(src, old_met_check, new_met_check, "Track inner-sphere metal coordination count")

src = ensure_replaced(
    src,
    '            pair_hb_base, pair_hb_2oh, pair_hb_bb = 0, 0, 0',
    '            pair_hb_base, pair_hb_2oh, pair_hb_bb, pair_metal_coord = 0, 0, 0, 0',
    "Initialize pair_metal_coord per residue pair",
)

old_hp_nuc = '''                        elif is_nucleic_pair:
                            prot_is_a = (ra.entity_type == "PROTEIN")
                            p_at, n_at = (atom_a, atom_b) if prot_is_a else (atom_b, atom_a)
                            p_ap = is_apolar_a if prot_is_a else is_apolar_b
                            if p_ap and p_at.name not in ("C", "CA") and n_at.name in NUCLEIC_CARBON_ATOMS:
                                pair_hp += 1
                                hp_dist_sum += clamp((4.75 - d) / 0.80, 0.35, 1.0)'''

new_hp_nuc = '''                        elif is_nucleic_pair or (is_ligand_pair and (ra.entity_type == "NUCLEIC" or rb.entity_type == "NUCLEIC")):
                            nuc_is_a = (ra.entity_type == "NUCLEIC")
                            n_at, o_at = (atom_a, atom_b) if nuc_is_a else (atom_b, atom_a)
                            o_ap = is_apolar_b if nuc_is_a else is_apolar_a
                            if o_ap and o_at.name not in ("C", "CA") and n_at.name in NUCLEIC_CARBON_ATOMS:
                                pair_hp += 1
                                hp_dist_sum += clamp((4.75 - d) / 0.80, 0.35, 1.0)'''
src = ensure_replaced(src, old_hp_nuc, new_hp_nuc, "Enable RNA/DNA-Ligand hydrophobic groove packing")

old_pi_block = '''            pair_pi = 0
            pair_pi_w = 0.0
            pair_catpi = 0
            if is_nucleic_pair:
                ring_a, norm_a = get_ring_centroid_and_normal(ra)
                ring_b, norm_b = get_ring_centroid_and_normal(rb)
                if ring_a is not None and ring_b is not None:
                    d_ring = math.dist(ring_a, ring_b)
                    if d_ring <= 5.5:
                        n_ref = norm_a if ra.entity_type == "NUCLEIC" else norm_b
                        if norm_a is not None and norm_b is not None and n_ref is not None:
                            vec = (ring_b[0] - ring_a[0], ring_b[1] - ring_a[1], ring_b[2] - ring_a[2])
                            d_perp = abs(vec[0] * n_ref[0] + vec[1] * n_ref[1] + vec[2] * n_ref[2])
                            d_par = math.sqrt(max(0.0, d_ring * d_ring - d_perp * d_perp))
                            cos_ang = abs(norm_a[0] * norm_b[0] + norm_a[1] * norm_b[1] + norm_a[2] * norm_b[2])
                            if d_perp >= 2.2 and d_par <= 4.2:
                                if cos_ang >= 0.60:
                                    pair_pi = 1
                                    pair_pi_w = 1.0
                                elif d_ring <= 5.0:
                                    pair_pi = 1
                                    pair_pi_w = 0.5
                        else:
                            pair_pi = 1
                            pair_pi_w = 0.75
                # Geometry-aware Cation-pi (ARG:CZ / LYS:NZ above nucleobase face, |cos phi| >= 0.50)
                prot_r, nuc_r = (ra, rb) if ra.entity_type == "PROTEIN" else (rb, ra)
                if prot_r.canonical_name in ("ARG", "LYS"):
                    nuc_cent, nuc_norm = get_ring_centroid_and_normal(nuc_r)
                    cat_name = "CZ" if prot_r.canonical_name == "ARG" else "NZ"
                    cat_at = next((a for a in prot_r.atoms if a.name == cat_name), None)
                    if nuc_cent is not None and nuc_norm is not None and cat_at is not None:
                        d_cp = math.dist(cat_at.coord, nuc_cent)
                        if d_cp <= 5.0 and d_cp > 1e-4:
                            vx = cat_at.coord[0] - nuc_cent[0]
                            vy = cat_at.coord[1] - nuc_cent[1]
                            vz = cat_at.coord[2] - nuc_cent[2]
                            cos_phi = abs(vx * nuc_norm[0] + vy * nuc_norm[1] + vz * nuc_norm[2]) / d_cp
                            if cos_phi >= 0.50:
                                pair_catpi = 1
            else:
                ring_a = get_aromatic_centroid(ra)
                ring_b = get_aromatic_centroid(rb)
                if ring_a is not None and ring_b is not None:
                    if math.dist(ring_a, ring_b) <= 5.5:
                        pair_pi += 1
                        pair_pi_w = 1.0

            eff_sb = min(pair_sb, 2)
            eff_hb = min(pair_hb, 2)
            eff_hb_base = min(pair_hb_base, eff_hb)
            eff_hb_2oh = min(pair_hb_2oh, eff_hb - eff_hb_base)
            eff_hb_bb = min(pair_hb_bb, eff_hb - eff_hb_base - eff_hb_2oh)
            eff_hp = min(pair_hp, 6)
            eff_hp_weighted = min(hp_dist_sum, 6.0)
            eff_clash = min(pair_clash, 3)

            if w_ij >= 0.25:
                n_sb += eff_sb
                n_hb += eff_hb
                n_hb_base += eff_hb_base
                n_hb_2oh += eff_hb_2oh
                n_hb_bb += eff_hb_bb
                n_hp += eff_hp
                n_pi += pair_pi
                n_catpi += pair_catpi
                n_ss += pair_ss
                n_clash += eff_clash

            w_sb += w_ij * eff_sb
            w_hb += w_ij * eff_hb
            w_hb_base += w_ij * eff_hb_base
            w_hb_2oh += w_ij * eff_hb_2oh
            w_hb_bb += w_ij * eff_hb_bb
            w_hp += w_ij * eff_hp_weighted
            w_pi += w_ij * pair_pi
            w_pi_rna += w_ij * pair_pi_w
            w_catpi += w_ij * pair_catpi
            w_clash += w_ij * eff_clash'''

new_pi_block = '''            pair_pi = 0
            pair_pi_w = 0.0
            pair_catpi = 0
            pair_xb = count_directional_halogen_bonds(ra, rb) if is_ligand_pair else 0
            if is_nucleic_pair or is_ligand_pair:
                rings_a = get_all_rings_for_residue(ra)
                rings_b = get_all_rings_for_residue(rb)
                if rings_a and rings_b:
                    best_pw = 0.0
                    for ring_a, norm_a in rings_a:
                        for ring_b, norm_b in rings_b:
                            d_ring = math.dist(ring_a, ring_b)
                            if d_ring <= 5.5:
                                n_ref = norm_a if norm_a is not None else norm_b
                                if norm_a is not None and norm_b is not None and n_ref is not None:
                                    vec = (ring_b[0] - ring_a[0], ring_b[1] - ring_a[1], ring_b[2] - ring_a[2])
                                    d_perp = abs(vec[0] * n_ref[0] + vec[1] * n_ref[1] + vec[2] * n_ref[2])
                                    d_par = math.sqrt(max(0.0, d_ring * d_ring - d_perp * d_perp))
                                    cos_ang = abs(norm_a[0] * norm_b[0] + norm_a[1] * norm_b[1] + norm_a[2] * norm_b[2])
                                    if d_perp >= 2.2 and d_par <= 4.2:
                                        if cos_ang >= 0.60:
                                            best_pw = max(best_pw, 1.0)
                                        elif d_ring <= 5.0:
                                            best_pw = max(best_pw, 0.5)
                                else:
                                    best_pw = max(best_pw, 0.75)
                    if best_pw > 0.0:
                        pair_pi = 1
                        pair_pi_w = best_pw
                for p_cand, r_cand in ((ra, rb), (rb, ra)):
                    if p_cand.canonical_name in ("ARG", "LYS"):
                        cand_rings = get_all_rings_for_residue(r_cand)
                        cat_name = "CZ" if p_cand.canonical_name == "ARG" else "NZ"
                        cat_at = next((a for a in p_cand.atoms if a.name == cat_name), None)
                        if cand_rings and cat_at is not None:
                            for r_cent, r_norm in cand_rings:
                                if r_norm is None:
                                    continue
                                d_cp = math.dist(cat_at.coord, r_cent)
                                if 1e-4 < d_cp <= 5.0:
                                    vx = cat_at.coord[0] - r_cent[0]
                                    vy = cat_at.coord[1] - r_cent[1]
                                    vz = cat_at.coord[2] - r_cent[2]
                                    cos_phi = abs(vx * r_norm[0] + vy * r_norm[1] + vz * r_norm[2]) / d_cp
                                    if cos_phi >= 0.50:
                                        pair_catpi = 1
                                        break
            else:
                ring_a = get_aromatic_centroid(ra)
                ring_b = get_aromatic_centroid(rb)
                if ring_a is not None and ring_b is not None:
                    if math.dist(ring_a, ring_b) <= 5.5:
                        pair_pi += 1
                        pair_pi_w = 1.0

            eff_sb = min(pair_sb, 2)
            eff_hb = min(pair_hb, 3 if is_duplex_pair else 2)
            eff_hb_base = min(pair_hb_base, eff_hb)
            eff_hb_2oh = min(pair_hb_2oh, eff_hb - eff_hb_base)
            eff_hb_bb = min(pair_hb_bb, eff_hb - eff_hb_base - eff_hb_2oh)
            eff_hp = min(pair_hp, 6)
            eff_hp_weighted = min(hp_dist_sum, 6.0)
            eff_clash = min(pair_clash, 3)
            eff_mcoord = min(pair_metal_coord, 2)

            if w_ij >= 0.25:
                n_sb += eff_sb
                n_hb += eff_hb
                n_hb_base += eff_hb_base
                n_hb_2oh += eff_hb_2oh
                n_hb_bb += eff_hb_bb
                n_hp += eff_hp
                n_pi += pair_pi
                n_catpi += pair_catpi
                n_xb += pair_xb
                n_metal_coord += eff_mcoord
                n_ss += pair_ss
                n_clash += eff_clash

            w_sb += w_ij * eff_sb
            w_hb += w_ij * eff_hb
            w_hb_base += w_ij * eff_hb_base
            w_hb_2oh += w_ij * eff_hb_2oh
            w_hb_bb += w_ij * eff_hb_bb
            w_hp += w_ij * eff_hp_weighted
            w_pi += w_ij * pair_pi
            w_pi_rna += w_ij * pair_pi_w
            w_catpi += w_ij * pair_catpi
            w_xb += w_ij * pair_xb
            w_metal_coord += w_ij * eff_mcoord
            w_clash += w_ij * eff_clash'''
src = ensure_replaced(src, old_pi_block, new_pi_block, "Upgrade pi-pi, cation-pi, halogen bonds & duplex 3-HB cap")

src = ensure_replaced(
    src,
    '                "cat_pi": pair_catpi,\n                "hb_base": eff_hb_base,',
    '                "cat_pi": pair_catpi,\n                "x_bond": pair_xb,\n                "metal_coord": eff_mcoord,\n                "hb_base": eff_hb_base,',
    "Add x_bond & metal_coord to contact log",
)

old_thermo_s3 = '''    if is_nucleic_pair:
        debye_rna = math.exp(-0.85 * ((math.sqrt(i_m) / 3.04) - (math.sqrt(0.15) / 3.04)) * 3.5)
        e_raw = (
            0.60 * w_sb * debye_rna * ph_k
            + 0.65 * w_hb_base
            + 0.55 * w_hb_2oh
            + 0.35 * w_hb_bb
            + 0.95 * w_pi_rna
            + 0.85 * w_catpi
            + 0.10 * w_hp
            + 0.0025 * w_bsa
        )
        e_eff = 16.5 * math.tanh(e_raw / 16.5)
        dg_polar_desolv = 0.0
        dg_steric_clash = +1.50 * w_clash
        dg_proline_kink = 0.0
        dg_phys = dg_cratic - e_eff + dg_steric_clash
    else:
        e_raw = s_out * (0.110 * w_hp + 0.195 * prodigy_ic["APOLAR_APOLAR"] + 0.0032 * w_bsa) + (0.78 * w_sb * debye_k * ph_k + 0.48 * w_hb + 0.62 * w_pi)
        e_eff = 15.5 * math.tanh(e_raw / 15.5)
        dg_polar_desolv = +0.12 * prodigy_ic["POLAR_POLAR"]
        dg_steric_clash = +1.50 * w_clash
        dg_proline_kink = +2.20 * helical_prolines
        dg_phys = dg_cratic - e_eff + dg_polar_desolv + dg_steric_clash + dg_proline_kink

    eff_ic_sum = sum(prodigy_ic.values())
    if is_nucleic_pair or eff_ic_sum < 3.0:
        contact_gate = 0.0'''

new_thermo_s3 = f'''    lig_burial_frac = 0.0
    if is_ligand_pair:
        lig_is_b = any(r.entity_type == "LIGAND" for r in residues_b)
        lig_res_list = residues_b if lig_is_b else residues_a
        lig_sasa_free_map = {sasa_free_b_arg} if lig_is_b else {sasa_free_a_arg}
        lig_free_tot = sum(lig_sasa_free_map.get(r.key, 0.0) for r in lig_res_list)
        lig_bound_tot = sum({sasa_comp_arg}.get(r.key, 0.0) for r in lig_res_list)
        lig_burial_frac = clamp((lig_free_tot - lig_bound_tot) / max(1.0, lig_free_tot), 0.0, 1.0)

    if is_nucleic_pair:
        debye_rna = math.exp(-0.85 * ((math.sqrt(i_m) / 3.04) - (math.sqrt(0.15) / 3.04)) * 3.5)
        e_raw = (
            0.60 * w_sb * debye_rna * ph_k
            + 0.65 * w_hb_base
            + 0.55 * w_hb_2oh
            + 0.35 * w_hb_bb
            + 0.95 * w_pi_rna
            + 0.85 * w_catpi
            + 0.10 * w_hp
            + 0.0025 * w_bsa
        )
        e_eff = 16.5 * math.tanh(e_raw / 16.5)
        dg_polar_desolv = 0.0
        dg_steric_clash = +1.50 * w_clash
        dg_proline_kink = 0.0
        dg_phys = dg_cratic - e_eff + dg_steric_clash
    elif is_ligand_pair:
        e_raw = (
            0.85 * w_sb * debye_k * ph_k
            + 0.65 * w_hb
            + 0.55 * w_xb
            + 0.90 * w_pi_rna
            + 0.80 * w_catpi
            + 0.14 * w_hp
            + 0.0055 * w_bsa
        )
        e_eff = 15.0 * math.tanh(e_raw / 15.0)
        dg_polar_desolv = 0.0
        dg_steric_clash = +1.50 * w_clash
        dg_proline_kink = 0.0
        dg_phys = 1.80 - e_eff + dg_steric_clash
    elif is_metal_pair:
        ion_res = [r for r in (residues_a + residues_b) if r.entity_type == "ION"]
        is_trans_met = any(r.canonical_name in ("ZN", "FE", "CU", "CO", "NI") for r in ion_res)
        k_coord = 2.60 if is_trans_met else 1.75
        if n_metal_coord == 0:
            dg_phys = dg_cratic + 1.50 * w_clash
        else:
            e_raw = k_coord * min(w_metal_coord, 6.0) + 0.0040 * w_bsa
            e_eff = 13.5 * math.tanh(e_raw / 13.5)
            dg_phys = 1.50 - e_eff + 1.50 * w_clash
        dg_polar_desolv = 0.0
        dg_steric_clash = +1.50 * w_clash
        dg_proline_kink = 0.0
    else:
        e_raw = s_out * (0.110 * w_hp + 0.195 * prodigy_ic["APOLAR_APOLAR"] + 0.0032 * w_bsa) + (0.78 * w_sb * debye_k * ph_k + 0.48 * w_hb + 0.62 * w_pi)
        e_eff = 15.5 * math.tanh(e_raw / 15.5)
        dg_polar_desolv = +0.12 * prodigy_ic["POLAR_POLAR"]
        dg_steric_clash = +1.50 * w_clash
        dg_proline_kink = +2.20 * helical_prolines
        dg_phys = dg_cratic - e_eff + dg_polar_desolv + dg_steric_clash + dg_proline_kink

    eff_ic_sum = sum(prodigy_ic.values())
    if is_non_ppi or eff_ic_sum < 3.0:
        contact_gate = 0.0'''
src = ensure_replaced(src, old_thermo_s3, new_thermo_s3, "Add MODE_LIGAND & MODE_METAL thermodynamics and disable PRODIGY on non-PPI")

old_uis_calc = '''    if mode == "MODE_A_PRED_FULL":
        raw_uis = 0.25 * float(iptm_val) + 0.30 * lis_score + 0.25 * pdockq_pae + 0.20 * sigma_aff
    else:
        raw_uis = 0.55 * clamp(pdockq_coord / 0.742, 0.0, 1.0) + 0.45 * sigma_aff
    thermo_gate = clamp((-dg_uni + 0.5) / 6.0, 0.10, 1.0)
    uis = clamp(raw_uis * thermo_gate, 0.0, 1.0)'''

new_uis_calc = '''    pocket_quality = 0.0
    if is_ligand_pair:
        pocket_quality = mean_w * clamp((lig_burial_frac - 0.15) / 0.65, 0.0, 1.0)
        if mode == "MODE_A_PRED_FULL":
            raw_uis = 0.20 * float(iptm_val) + 0.25 * lis_score + 0.30 * pocket_quality + 0.25 * sigma_aff
        else:
            raw_uis = 0.55 * pocket_quality + 0.45 * sigma_aff
    elif is_metal_pair:
        pocket_quality = mean_w * clamp(w_metal_coord / 3.0, 0.0, 1.0) if n_metal_coord > 0 else 0.0
        raw_uis = 0.55 * pocket_quality + 0.45 * sigma_aff
    elif mode == "MODE_A_PRED_FULL":
        raw_uis = 0.25 * float(iptm_val) + 0.30 * lis_score + 0.25 * pdockq_pae + 0.20 * sigma_aff
    else:
        raw_uis = 0.55 * clamp(pdockq_coord / 0.742, 0.0, 1.0) + 0.45 * sigma_aff
    thermo_gate = clamp((-dg_uni + 0.5) / 6.0, 0.10, 1.0)
    uis = clamp(raw_uis * thermo_gate, 0.0, 1.0)'''
src = ensure_replaced(src, old_uis_calc, new_uis_calc, "Add Pocket Quality (MODE_LIGAND) & Coordination Quality (MODE_METAL) UIS scoring")

src = ensure_replaced(
    src,
    '        uis_rounded, float(iptm_val), lis_score, dg_uni, w_bsa, 0 if is_nucleic_pair else helical_prolines, n_clash',
    '        uis_rounded, float(iptm_val), lis_score, dg_uni, w_bsa, 0 if is_non_ppi else helical_prolines, n_clash',
    "Zero helical_prolines penalty in pattern classifier for all non-PPI modalities",
)
src = ensure_replaced(
    src,
    '"helical_prolines": 0 if is_nucleic_pair else helical_prolines,',
    '"helical_prolines": 0 if is_non_ppi else helical_prolines,',
    "Zero helical_prolines output field for all non-PPI modalities",
)

old_ret_fields = '''        "cation_pi": n_catpi,
        "rna_specificity_idx": rna_ssi,'''
new_ret_fields = '''        "cation_pi": n_catpi,
        "halogen_bonds": n_xb,
        "metal_coord_bonds": n_metal_coord,
        "lig_burial_pct": round(lig_burial_frac * 100.0, 1),
        "pocket_quality": round(pocket_quality, 3),
        "rna_specificity_idx": rna_ssi,'''
src = ensure_replaced(src, old_ret_fields, new_ret_fields, "Expose halogen_bonds, metal_coord_bonds, lig_burial_pct & pocket_quality in CLI report")

# ------------------------------------------------------------------------------
# 5. FIX _fmt_segs DIGIT BUG (3RJ401 -> 401) & MULTI-CHAIN COMPOSITE LOCI
# ------------------------------------------------------------------------------
old_fmt_segs = '''        def _fmt_segs(keys):
            nums = sorted({int("".join(ch for ch in k.split(":")[-1] if ch.isdigit()) or "0") for k in keys})
            if not nums: return "None"
            segs, s, p = [], nums[0], nums[0]
            for v in nums[1:]:
                if v - p > 15:
                    segs.append(f"{s}" if s == p else f"{s}-{p}")
                    s = v
                p = v
            segs.append(f"{s}" if s == p else f"{s}-{p}")
            return ", ".join(segs)'''

new_fmt_segs = '''        def _fmt_segs(keys):
            by_ch: Dict[str, Set[int]] = {}
            for k in keys:
                parts_k = k.split(":", 1)
                ch_k = parts_k[0] if len(parts_k) == 2 else ""
                m_num = re.search(r"(\\d+)$", parts_k[-1])
                if m_num:
                    by_ch.setdefault(ch_k, set()).add(int(m_num.group(1)))
            if not by_ch:
                return "None"
            out_chunks = []
            for ch_k, num_set in by_ch.items():
                nums = sorted(num_set)
                segs, s, p = [], nums[0], nums[0]
                for v in nums[1:]:
                    if v - p > 15:
                        segs.append(f"{s}" if s == p else f"{s}-{p}")
                        s = v
                    p = v
                segs.append(f"{s}" if s == p else f"{s}-{p}")
                out_chunks.append(f"{ch_k}:{', '.join(segs)}" if len(by_ch) > 1 else ", ".join(segs))
            return " | ".join(out_chunks)'''
src = ensure_replaced(src, old_fmt_segs, new_fmt_segs, "Fix _fmt_segs trailing-digit extraction (3RJ401 -> 401) & composite ternary loci")

# ------------------------------------------------------------------------------
# 6. AST-ANCHORED UPGRADE OF process_structure_file():
#    (a) Add split_chain param (immune to -> List[Dict]: return hints)
#    (b) Insert --split-chain slicer right after build_residues(atoms)
#    (c) Insert Automatic Ternary Glue & Metal-Bridge Detector right before return
# ------------------------------------------------------------------------------
psf_node = get_func_node(src, "process_structure_file")
psf_arg_names = [a.arg for a in psf_node.args.args]
src_lines = src.splitlines(keepends=True)

# 6a. Signature injection via AST bounds
if "split_chain" in psf_arg_names:
    print("  [OK-ALREADY] AST: split_chain parameter in process_structure_file()")
else:
    first_body_ln = psf_node.body[0].lineno
    hdr_text = "".join(src_lines[psf_node.lineno - 1 : first_body_ln - 1])
    # Find closing ')' of parameter list (before optional '-> ...:' or ':')
    m_close = re.search(r'\)\s*(?:->[^:]+)?\s*:', hdr_text)
    assert m_close, f"Could not parse header of process_structure_file:\n{hdr_text}"
    pos = m_close.start()
    before = hdr_text[:pos].rstrip()
    sep = " " if before.endswith(",") else ", "
    new_hdr = before + sep + "split_chain: Optional[str] = None" + hdr_text[pos:]
    src_lines[psf_node.lineno - 1 : first_body_ln - 1] = [new_hdr]
    src = "".join(src_lines)
    print("  [APPLIED]    AST: split_chain parameter in process_structure_file()")

# 6b. Insert --split-chain slicer right after `build_residues(atoms)` inside process_structure_file
psf_node = get_func_node(src, "process_structure_file")
src_lines = src.splitlines(keepends=True)
psf_body_txt = "".join(src_lines[psf_node.lineno - 1 : psf_node.end_lineno])
if "if split_chain:" in psf_body_txt:
    print("  [OK-ALREADY] AST: --split-chain domain slicer in process_structure_file()")
else:
    split_code = '''    if split_chain:
        # Format: "A:1-236,237-406" or "A:236" (splits chain A into A_D1 and A_D2 for intra-chain domain auditing)
        try:
            sc_cid, sc_spec = split_chain.split(":", 1)
            sc_cid = sc_cid.strip()
            if sc_cid in chains:
                src_list = chains[sc_cid]
                d1, d2 = [], []
                if "," in sc_spec:
                    p1, p2 = sc_spec.split(",", 1)
                    lo1, hi1 = (int(x) for x in p1.split("-"))
                    lo2, hi2 = (int(x) for x in p2.split("-"))
                    for r in src_list:
                        if lo1 <= r.res_seq <= hi1:
                            d1.append(r)
                        elif lo2 <= r.res_seq <= hi2:
                            d2.append(r)
                else:
                    cut_res = int(sc_spec.strip())
                    for r in src_list:
                        (d1 if r.res_seq <= cut_res else d2).append(r)
                if d1 and d2:
                    c1, c2 = f"{sc_cid}_D1", f"{sc_cid}_D2"
                    for r in d1:
                        _universal_set_chain(r, c1)
                    for r in d2:
                        _universal_set_chain(r, c2)
                    del chains[sc_cid]
                    chains[c1] = d1
                    chains[c2] = d2
                    kept_ids = {id(r) for cl in chains.values() for r in cl}
                    ordered_residues = [r for r in ordered_residues if id(r) in kept_ids]
        except Exception as _e:
            print(f"[WARN] Could not parse --split-chain '{split_chain}': {_e}")
'''
    ins_idx = None
    for idx in range(psf_node.lineno - 1, psf_node.end_lineno):
        if "build_residues(atoms)" in src_lines[idx]:
            ins_idx = idx + 1
            break
    assert ins_idx is not None, "Could not locate build_residues(atoms) call inside process_structure_file()"
    src_lines.insert(ins_idx, split_code)
    src = "".join(src_lines)
    print(f"  [APPLIED]    AST: --split-chain domain slicer after L{ins_idx}")

# 6c. Insert Automatic Ternary Glue & Catalytic Metal-Bridge Detector right before final `return results`
psf_node = get_func_node(src, "process_structure_file")
src_lines = src.splitlines(keepends=True)
psf_body_txt = "".join(src_lines[psf_node.lineno - 1 : psf_node.end_lineno])
if "Automatic Ternary Molecular-Glue & Catalytic Metal-Bridge Detector" in psf_body_txt:
    print("  [OK-ALREADY] AST: Ternary Glue & Metal-Bridge Detector in process_structure_file()")
else:
    tern_code = '''
    # Automatic Ternary Molecular-Glue & Catalytic Metal-Bridge Detector (v2.0-MARVAS):
    # Detects when a small-molecule drug bridges >= 2 macromolecules (e.g. eIF4A1 + RNA + Rocaglamide A)
    # or a catalytic metal bridges >= 2 chains (e.g. eIF4A1 + ANP500 + Mg501) and scores the composite cleft.
    if len(chain_ids) >= 3:
        for lc in chain_ids:
            lc_res = chains[lc]
            if not lc_res or lc_res[0].entity_type not in ("LIGAND", "ION"):
                continue
            is_lig_c = (lc_res[0].entity_type == "LIGAND")
            partner_hits = []
            for rrec in results:
                cp = rrec["chain_pair"].split("-")
                if len(cp) != 2 or lc not in cp:
                    continue
                other = cp[1] if cp[0] == lc else cp[0]
                if other not in chains or not chains[other]:
                    continue
                o_type = chains[other][0].entity_type
                if is_lig_c:
                    if o_type in ("PROTEIN", "NUCLEIC") and rrec.get("lig_burial_pct", 0.0) >= 12.0 and rrec.get("conf_contacts", 0) >= 1:
                        partner_hits.append((other, rrec.get("lig_burial_pct", 0.0)))
                else:
                    if o_type in ("PROTEIN", "NUCLEIC", "LIGAND") and rrec.get("metal_coord_bonds", 0) >= 1:
                        partner_hits.append((other, float(rrec.get("metal_coord_bonds", 0))))
            if len(partner_hits) >= 2:
                partner_hits.sort(key=lambda x: x[1], reverse=True)
                (m1, v1), (m2, v2) = partner_hits[0], partner_hits[1]
                comp_res_m = chains[m1] + chains[m2]
                comp_id_m = f"{m1}+{m2}"
                sasa_m_free = compute_residue_sasa(comp_res_m)
                sasa_tern_comp = compute_residue_sasa(comp_res_m + lc_res)
                t_pair = analyze_interface_pair(
                    comp_id_m, lc, comp_res_m, lc_res, res_to_global_idx,
                    mode, conf_data, sasa_tern_comp,
                    sasa_m_free, sasa_free_by_chain[lc],
                    chain_remap=chain_remap,
                    temp_k=temp_k, lambda_pae=lambda_pae, ionic_mM=ionic_mM, ph_val=ph_val,
                    force_report_if_empty=False,
                    custom_mutations=None
                )
                if t_pair:
                    t_pair["complex_id"] = os.path.splitext(os.path.basename(filepath))[0]
                    t_pair["source_file"] = os.path.basename(filepath)
                    t_pair["ternary_partners"] = f"{m1}+{m2}"
                    t_pair["ternary_burial_split"] = f"{m1}:{v1:.1f}% / {m2}:{v2:.1f}%" if is_lig_c else f"{m1}:{int(v1)} coord / {m2}:{int(v2)} coord"
                    t_pair["cooperativity_balance"] = round(min(v1, v2) / max(1e-6, max(v1, v2)), 2)
                    if is_lig_c and t_pair.get("lig_burial_pct", 0.0) >= 50.0 and t_pair["UIS"] >= UIS_DECISION_THRESHOLD:
                        t_pair["pattern_id"] = "Pattern A"
                        t_pair["pattern_desc"] = f"Ternary Molecular Glue / Interface Wedge ({lc} bridges {m1} [{v1:.1f}%] & {m2} [{v2:.1f}%]; composite burial {t_pair['lig_burial_pct']:.1f}%, balance {t_pair['cooperativity_balance']:.2f})"
                    elif not is_lig_c and t_pair.get("metal_coord_bonds", 0) >= 2:
                        t_pair["pattern_id"] = "Pattern A"
                        t_pair["pattern_desc"] = f"Catalytic Metal Bridge ({lc} coordinates {m1} [{int(v1)}] & {m2} [{int(v2)}] inner-sphere bonds)"
                    results.append(t_pair)
'''
    ret_idx = None
    for idx in range(psf_node.end_lineno - 1, psf_node.lineno - 1, -1):
        if src_lines[idx].strip() == "return results":
            ret_idx = idx
            break
    assert ret_idx is not None, "Could not locate `return results` inside process_structure_file()"
    src_lines.insert(ret_idx, tern_code)
    src = "".join(src_lines)
    print(f"  [APPLIED]    AST: Ternary Glue & Metal-Bridge Detector before L{ret_idx + 1}")

# ------------------------------------------------------------------------------
# 7. UPGRADE MONOMER CONTACT-TROUGH DOMAIN SPLITTER IN build_residues()
# ------------------------------------------------------------------------------
old_mono_cut = '''        if len(all_res) >= 10 and all(_universal_res_etype(r) in ("PROTEIN", "NUCLEIC") for r in all_res):
            mid = len(all_res) // 2'''
new_mono_cut = '''        if len(all_res) >= 10 and all(_universal_res_etype(r) in ("PROTEIN", "NUCLEIC") for r in all_res):
            # Automatic contact-trough domain boundary finder (25%..75% of chain length)
            n_tot = len(all_res)
            if n_tot >= 40:
                lo_c, hi_c = max(15, n_tot // 4), min(n_tot - 15, (3 * n_tot) // 4)
                step_c = max(1, (hi_c - lo_c) // 25)
                best_cut, best_cross = n_tot // 2, 10**9
                cents = [r.centroid for r in all_res]
                for cut in range(lo_c, hi_c + 1, step_c):
                    cross = 0
                    for ii in range(max(0, cut - 35), cut):
                        ci = cents[ii]
                        for jj in range(cut, min(n_tot, cut + 35)):
                            if math.dist(ci, cents[jj]) <= 8.0:
                                cross += 1
                    if cross < best_cross:
                        best_cross, best_cut = cross, cut
                mid = best_cut
            else:
                mid = n_tot // 2'''
src = ensure_replaced(src, old_mono_cut, new_mono_cut, "Upgrade single-chain monomer splitter to contact-trough domain finder")

# ------------------------------------------------------------------------------
# 8. ADD --split-chain CLI ARG, CSV COLUMNS, & TERNARY/LIGAND CLI OUTPUT
# ------------------------------------------------------------------------------
src = ensure_replaced(
    src,
    '    parser.add_argument("--mutate", type=str, default="", help="Comma-separated custom mutations (e.g. A:LEU5PRO,A:GLU1LYS)")',
    '    parser.add_argument("--mutate", type=str, default="", help="Comma-separated custom mutations (e.g. A:LEU5PRO,A:GLU1LYS)")\n'
    '    parser.add_argument("--split-chain", type=str, default=None, help="Split one chain into 2 domains for intra-chain auditing (e.g. A:1-236,237-406 or A:236)")',
    "Add --split-chain CLI flag",
)

src = ensure_replaced(
    src,
    '            fpath, args.output, temp_k=args.temp, lambda_pae=args.lambda_pae, ionic_mM=args.ionic_strength, ph_val=args.ph, custom_mutations=custom_mut_list',
    '            fpath, args.output, temp_k=args.temp, lambda_pae=args.lambda_pae, ionic_mM=args.ionic_strength, ph_val=args.ph, custom_mutations=custom_mut_list, split_chain=args.split_chain',
    "Pass args.split_chain into process_structure_file()",
)

src = ensure_replaced(
    src,
    '"mean_w_intf", "cation_pi", "hb_base", "hb_2oh", "hb_backbone", "rna_specificity_idx", "top_hotspot", "top_hotspot_ddG_Ala"',
    '"mean_w_intf", "cation_pi", "halogen_bonds", "metal_coord_bonds", "lig_burial_pct", "pocket_quality", "hb_base", "hb_2oh", "hb_backbone", "rna_specificity_idx", "ternary_partners", "cooperativity_balance", "top_hotspot", "top_hotspot_ddG_Ala"',
    "Add Ligand/Metal/Ternary columns to Master CSV",
)

src = ensure_replaced(
    src,
    'print(f"{\'Complex ID\':<38} | {\'Pair\':<4} | {\'UIS\':<6} | {\'Margin\':<7} | {\'dG (kcal)\':<9} | {\'Kd (nM)\':<11} | {\'LIS\':<5} | {\'Pattern\':<9} | {\'Class\'}")',
    'print(f"{\'Complex ID\':<34} | {\'Pair\':<16} | {\'Modality\':<15} | {\'UIS\':<6} | {\'dG (kcal)\':<9} | {\'Kd (nM)\':<11} | {\'LIS\':<5} | {\'Pattern\':<9} | {\'Class\'}")',
    "Add Modality column & wider Pair column to CLI table header",
)

src = ensure_replaced(
    src,
    '''        print(
            f"{r['complex_id'][:38]:<38} | {r['chain_pair']:<4} | {r['UIS']:<6.3f} | {margin_str:<7} | "
            f"{r['dG_uni_kcal']:<9.2f} | {kd_str} | {r['LIS']:<5.3f} | "
            f"{r['pattern_id']:<9} | {r['binding_class']}"
        )''',
    '''        print(
            f"{r['complex_id'][:34]:<34} | {r['chain_pair'][:16]:<16} | {r.get('modality','PPI')[:15]:<15} | {r['UIS']:<6.3f} | "
            f"{r['dG_uni_kcal']:<9.2f} | {kd_str} | {r['LIS']:<5.3f} | "
            f"{r['pattern_id']:<9} | {r['binding_class']}"
        )''',
    "Print Modality & full Pair ID in CLI leaderboard rows",
)

old_step4_cli = '''        if _is_nuc:
            print(f"  4. Thermodynamics (Dedicated 9-Ch {_mod} Model, PRODIGY disabled): dG_RNA = {top['dG_uni_kcal']:.2f} kcal/mol (Kd = {_kd_disp})")
        else:'''

new_step4_cli = '''        if _is_nuc:
            print(f"  4. Thermodynamics (Dedicated 9-Ch {_mod} Model, PRODIGY disabled): dG_RNA = {top['dG_uni_kcal']:.2f} kcal/mol (Kd = {_kd_disp})")
        elif _mod in ("PROTEIN-LIGAND", "RNA-LIGAND", "DNA-LIGAND", "TERNARY-GLUE"):
            print(f"  4. Thermodynamics (Dedicated {_mod} Pocket Model, PRODIGY disabled): Burial = {top.get('lig_burial_pct',0.0):.1f}%, Q_pocket = {top.get('pocket_quality',0.0):.3f}, pi = {top.get('pi_stacking',0)}, Cat-pi = {top.get('cation_pi',0)}, X-Bond = {top.get('halogen_bonds',0)} -> dG_lig = {top['dG_uni_kcal']:.2f} kcal/mol (Kd = {_kd_disp})")
        elif _mod in ("PROTEIN-METAL", "RNA-METAL", "DNA-METAL", "LIGAND-METAL", "METAL-BRIDGE"):
            print(f"  4. Thermodynamics (Dedicated {_mod} Coordination Model): Inner-Sphere Coord = {top.get('metal_coord_bonds',0)}, Q_metal = {top.get('pocket_quality',0.0):.3f} -> dG_met = {top['dG_uni_kcal']:.2f} kcal/mol (Kd = {_kd_disp})")
        else:'''
src = ensure_replaced(src, old_step4_cli, new_step4_cli, "Add Ligand/Ternary/Metal Step 4 CLI derivation output")

old_save_msg = '    print(f"[UniPPI-CLI] Saved Master CSV: {csv_path}")'
new_save_msg = '''    _terns = [r for r in all_records if r.get("modality") in ("TERNARY-GLUE", "METAL-BRIDGE")]
    if _terns:
        print("\\n[TERNARY MOLECULAR GLUE & CATALYTIC METAL-BRIDGE AUDIT]")
        for tr in _terns:
            _tkd = f"{tr['Kd_nM']*1000:.2f} pM" if tr["Kd_nM"] < 0.05 else (f"{tr['Kd_nM']:.2f} nM" if tr["Kd_nM"] < 1e7 else ">10mM (NB)")
            print(f"  * [{tr['modality']}] {tr['chain_pair']}: UIS = {tr['UIS']:.3f} | dG = {tr['dG_uni_kcal']:.2f} kcal/mol ({_tkd}) | Split: {tr.get('ternary_burial_split','')} | Balance = {tr.get('cooperativity_balance',0.0):.2f} | SB={tr['salt_bridges']} HB={tr['h_bonds']} pi={tr['pi_stacking']} Cat-pi={tr['cation_pi']} XB={tr.get('halogen_bonds',0)} HP={tr['hydrophobic_contacts']} | {tr['pattern_desc']}")
    print(f"[UniPPI-CLI] Saved Master CSV: {csv_path}")'''
src = ensure_replaced(src, old_save_msg, new_save_msg, "Print dedicated Ternary Glue & Metal Bridge Audit block")

ast.parse(src)
cli_path.write_text(src, encoding="utf-8")
print("\n[PASS] Step 3B-1 applied cleanly to unippi_cli.py (ast.parse OK).")

# ------------------------------------------------------------------------------
# 9. MASTER 6-TEST REAL BENCHMARK SUITE
# ------------------------------------------------------------------------------
tests = [
    ("TEST 1 (ZERO-REGRESSION PPI): fold_vip1_bzip61.zip", ["unippi_cli.py", "-i", "benchmarks/fold_vip1_bzip61.zip", "-o", "/tmp/v2_s3_vip1"]),
    ("TEST 2 (RNA + DRUG + TERNARY GLUE + MG BRIDGE): 5ZC9.cif", ["unippi_cli.py", "-i", "benchmarks/5ZC9.cif", "-o", "/tmp/v2_s3_5zc9"]),
    ("TEST 3 (RNA ITC GROUND TRUTH): 1URN.cif", ["unippi_cli.py", "-i", "benchmarks/1URN.cif", "-o", "/tmp/v2_s3_1urn"]),
    ("TEST 4 (MARVAS ULK1 KINASE + INHIBITOR 3RJ): 4WNP.cif", ["unippi_cli.py", "-i", "benchmarks/4WNP.cif", "-o", "/tmp/v2_s3_4wnp"]),
    ("TEST 5 (MARVAS CoV Nsp14-Nsp10 + ZN FINGERS + MG): 7DIY.cif", ["unippi_cli.py", "-i", "benchmarks/7DIY.cif", "-o", "/tmp/v2_s3_7diy"]),
    ("TEST 6 (CUSTOM SAME-PROTEIN INTRA-CHAIN SPLIT): 5ZC9 eIF4A1 NTD(1-236) vs CTD(237-406)", ["unippi_cli.py", "-i", "benchmarks/5ZC9.cif", "--split-chain", "A:1-236,237-406", "-o", "/tmp/v2_s3_5zc9_split"]),
]
for title, cmd in tests:
    print("\n" + "=" * 98)
    print(title)
    print("=" * 98)
    subprocess.run([sys.executable] + cmd, check=True)
