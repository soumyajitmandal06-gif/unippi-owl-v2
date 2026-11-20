#!/usr/bin/env python3
"""
Phase 2 - Step 2B-1 (Verified Anchors):
  Upgrades unippi_cli.py with:
    1. Exact _exptl.method parsing + crystallographic buffer/cryoprotectant filter
       (GOL/CL/SO4/EDO/PEG/etc.) + multi-chain biological HETATM drug/metal splitting.
    2. Removal of legacy main() post-overrides (restoring thermo_gate in MODE_C) and
       crystallographic Proline B-factor fix (residue_raw_plddt < 80.0).
    3. Geometry-aware 9-channel Protein-RNA / Protein-DNA engine (nucleobase pi-pi with
       plane-normal + lateral-slip guard, 60-deg-cone cation-pi, ribose 2'-OH vs
       Base-Specific vs Backbone H-bonds, SSI, PRODIGY auto-disable, calibrated dG_RNA).
    4. Runs Zero-Regression Guard (fold_vip1_bzip61.zip) + Real RNA Benchmarks (5ZC9, 1URN).
"""

import ast
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
src = cli_path.read_text(encoding="utf-8")

# ------------------------------------------------------------------------------
# 1. CONSTANTS: Crystal Buffer Junk, Catalytic Metals, Nucleobase Rings & RNA Sets
# ------------------------------------------------------------------------------
old_nuc_const = 'NUCLEIC_PHOSPHATE_O = {"OP1", "OP2", "O1P", "O2P"}'
new_nuc_const = '''NUCLEIC_PHOSPHATE_O = {"OP1", "OP2", "O1P", "O2P"}
CRYSTAL_BUFFER_JUNK = {
    "HOH", "WAT", "DOD", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G",
    "MPD", "SO4", "PO4", "CL", "NA", "K", "ACT", "ACY", "DMS", "MES",
    "TRS", "CIT", "FMT", "NO3", "BR", "IOD", "IPA", "EOH", "BU1", "BME",
    "NH4", "CAC", "IMD", "EPE", "HEP", "TAR", "TLA", "SUC", "MLI"
}
CATALYTIC_METALS = {"MG", "ZN", "MN", "CA", "FE", "CU", "CO", "NI"}
NUCLEOBASE_RINGS = {
    "A": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "G": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "DA": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "DG": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "ADE": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "GUA": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "C": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "U": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "T": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "DC": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "DT": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "CYT": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "URA": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "THY": {"N1", "C2", "N3", "C4", "C5", "C6"},
}
NUCLEOBASE_HB_ATOMS = {"N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"}
NUCLEIC_BB_ATOMS = {"OP1", "OP2", "O1P", "O2P", "O3'", "O4'", "O5'", "P"}
NUCLEIC_CARBON_ATOMS = {"C1'", "C2'", "C3'", "C4'", "C5'", "C2", "C4", "C5", "C6", "C7", "C8"}'''
assert old_nuc_const in src, "Missing NUCLEIC_PHOSPHATE_O anchor"
src = src.replace(old_nuc_const, new_nuc_const, 1)

# ------------------------------------------------------------------------------
# 2. PARSER FIXES: Exact _exptl.method & Crystallographic Buffer Filtering
# ------------------------------------------------------------------------------
old_exptl = '''        if s.startswith("_exptl.method"):
            meta["method"] = s.split(None, 1)[-1].strip("'\\"")'''
new_exptl = '''        if s.startswith("_exptl.method") and s.split()[0] == "_exptl.method":
            parts_ex = s.split(None, 1)
            if len(parts_ex) == 2 and parts_ex[1].strip("'\\"") not in ("?", "."):
                meta["method"] = parts_ex[1].strip("'\\"")'''
assert old_exptl in src, "Missing _exptl.method anchor"
src = src.replace(old_exptl, new_exptl, 1)

# Replace BOTH PDB (L236) and mmCIF (L304) water filters cleanly without substring collision
water_tok = 'if res_name in ("HOH", "WAT", "DOD"):'
assert src.count(water_tok) == 2, f"Expected 2 water filter sites, found {src.count(water_tok)}"
src = src.replace(water_tok, 'if res_name in CRYSTAL_BUFFER_JUNK:')

# ------------------------------------------------------------------------------
# 3. 3D RING CENTROID + UNIT NORMAL HELPER + MODALITY DETECTOR
# ------------------------------------------------------------------------------
old_arom_fn = '''def get_aromatic_centroid(res: Residue) -> Optional[Vec3]:
    cname = res.canonical_name
    if cname not in AROMATIC_RINGS:
        return None
    ring_names = AROMATIC_RINGS[cname]
    pts = [a.coord for a in res.atoms if a.name in ring_names]
    if len(pts) >= 3:
        n = len(pts)
        return (
            sum(p[0] for p in pts) / n,
            sum(p[1] for p in pts) / n,
            sum(p[2] for p in pts) / n
        )
    return None'''
new_arom_fn = '''def get_aromatic_centroid(res: Residue) -> Optional[Vec3]:
    cname = res.canonical_name
    ring_names = AROMATIC_RINGS.get(cname) or NUCLEOBASE_RINGS.get(cname)
    if not ring_names:
        return None
    pts = [a.coord for a in res.atoms if a.name in ring_names]
    if len(pts) >= 3:
        n = len(pts)
        return (
            sum(p[0] for p in pts) / n,
            sum(p[1] for p in pts) / n,
            sum(p[2] for p in pts) / n
        )
    return None

def get_ring_centroid_and_normal(res: Residue) -> Tuple[Optional[Vec3], Optional[Vec3]]:
    cname = res.canonical_name
    ring_names = AROMATIC_RINGS.get(cname) or NUCLEOBASE_RINGS.get(cname)
    if not ring_names:
        return None, None
    atom_map = {a.name: a.coord for a in res.atoms if a.name in ring_names}
    if len(atom_map) < 3:
        return None, None
    pts = list(atom_map.values())
    n = len(pts)
    cent = (sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n, sum(p[2] for p in pts) / n)
    triad_keys = ("C2", "C4", "C6") if cname in NUCLEOBASE_RINGS else ("CG", "CE1", "CE2")
    if all(k in atom_map for k in triad_keys):
        p1, p2, p3 = atom_map[triad_keys[0]], atom_map[triad_keys[1]], atom_map[triad_keys[2]]
    else:
        p1, p2, p3 = pts[0], pts[len(pts) // 3], pts[(2 * len(pts)) // 3]
    u = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    v = (p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2])
    nx = u[1] * v[2] - u[2] * v[1]
    ny = u[2] * v[0] - u[0] * v[2]
    nz = u[0] * v[1] - u[1] * v[0]
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm < 1e-6:
        return cent, None
    return cent, (nx / norm, ny / norm, nz / norm)

def detect_pair_modality(ca_id: str, cb_id: str, res_a: List[Residue], res_b: List[Residue]) -> str:
    if ca_id.endswith("_N") and cb_id.endswith("_C"):
        return "INTRA-DOMAIN (MONOMER)"
    types_a = {r.entity_type for r in res_a}
    types_b = {r.entity_type for r in res_b}
    all_t = types_a | types_b
    if "NUCLEIC" in all_t:
        nuc_res = [r for r in (res_a + res_b) if r.entity_type == "NUCLEIC"]
        is_rna = any(
            r.canonical_name in ("U", "URA") or any(a.name == "O2'" for a in r.atoms)
            for r in nuc_res
        )
        nuc_tag = "RNA" if is_rna else "DNA"
        if "PROTEIN" in all_t:
            return f"PROTEIN-{nuc_tag}"
        if "LIGAND" in all_t:
            return f"{nuc_tag}-LIGAND"
        if "ION" in all_t:
            return f"{nuc_tag}-METAL"
        return f"{nuc_tag}-{nuc_tag}"
    if "LIGAND" in all_t:
        return "PROTEIN-LIGAND"
    if "ION" in all_t:
        return "PROTEIN-METAL"
    return "PROTEIN-PROTEIN"'''
assert old_arom_fn in src, "Missing get_aromatic_centroid anchor"
src = src.replace(old_arom_fn, new_arom_fn, 1)

# ------------------------------------------------------------------------------
# 4. INTEGRATE 9-CHANNEL PROTEIN-RNA PHYSICS INSIDE analyze_interface_pair()
# ------------------------------------------------------------------------------
old_cnt_init = '''    n_sb, n_hb, n_hp, n_pi, n_ss, n_clash, n_cb8 = 0, 0, 0, 0, 0, 0, 0
    w_sb, w_hb, w_hp, w_pi, w_clash = 0.0, 0.0, 0.0, 0.0, 0.0'''
new_cnt_init = '''    modality = detect_pair_modality(chain_a_id, chain_b_id, residues_a, residues_b)
    is_nucleic_pair = modality in ("PROTEIN-RNA", "PROTEIN-DNA")
    n_sb, n_hb, n_hp, n_pi, n_ss, n_clash, n_cb8 = 0, 0, 0, 0, 0, 0, 0
    n_hb_base, n_hb_2oh, n_hb_bb, n_catpi = 0, 0, 0, 0
    w_sb, w_hb, w_hp, w_pi, w_clash = 0.0, 0.0, 0.0, 0.0, 0.0
    w_hb_base, w_hb_2oh, w_hb_bb, w_pi_rna, w_catpi = 0.0, 0.0, 0.0, 0.0, 0.0'''
assert old_cnt_init in src, "Missing n_sb/w_sb counter init anchor"
src = src.replace(old_cnt_init, new_cnt_init, 1)

old_inner_loop = '''            min_heavy_dist = 999.0
            pair_sb, pair_hb, pair_hp, pair_ss, pair_clash = 0, 0, 0, 0, 0
            hp_dist_sum = 0.0
            cna, cnb = ra.canonical_name, rb.canonical_name
            is_apolar_a = (cna in HYDROPHOBIC_RES) or (ra.entity_type == "LIGAND")
            is_apolar_b = (cnb in HYDROPHOBIC_RES) or (rb.entity_type == "LIGAND")

            for atom_a in ra.atoms:
                ca_xyz = atom_a.coord
                for atom_b in rb.atoms:
                    d = math.dist(ca_xyz, atom_b.coord)
                    if d < min_heavy_dist:
                        min_heavy_dist = d
                    if d > 4.5:
                        continue

                    if cna == "CYS" and cnb == "CYS" and atom_a.name == "SG" and atom_b.name == "SG" and d <= 2.3:
                        pair_ss += 1
                        continue

                    if d < 2.2:
                        pair_clash += 1

                    is_cat_a = (cna, atom_a.name) in CATIONIC_ATOMS
                    is_ani_a = (cna, atom_a.name) in ANIONIC_ATOMS or (ra.entity_type == "NUCLEIC" and atom_a.name in NUCLEIC_PHOSPHATE_O)
                    is_cat_b = (cnb, atom_b.name) in CATIONIC_ATOMS
                    is_ani_b = (cnb, atom_b.name) in ANIONIC_ATOMS or (rb.entity_type == "NUCLEIC" and atom_b.name in NUCLEIC_PHOSPHATE_O)

                    if d <= 4.0 and ((is_cat_a and is_ani_b) or (is_ani_a and is_cat_b)):
                        pair_sb += 1
                        continue

                    if d <= 3.5 and atom_a.element in ("N", "O", "S") and atom_b.element in ("N", "O", "S"):
                        pair_hb += 1
                        continue

                    if d <= 4.5 and is_apolar_a and is_apolar_b:
                        if atom_a.element in ("C", "S") and atom_b.element in ("C", "S"):
                            if atom_a.name not in ("C", "CA") and atom_b.name not in ("C", "CA"):
                                pair_hp += 1
                                hp_dist_sum += clamp((4.75 - d) / 0.80, 0.35, 1.0)

            if min_heavy_dist > 5.5:
                continue

            intf_res_a.add(ra.key)
            intf_res_b.add(rb.key)
            if cb_dist > 8.0:
                intf_pae_values.append(pair_pae)
                pair_weights_seen.append(w_ij)

            ca_cls = PRODIGY_CLASS.get(cna, "APOLAR" if ra.entity_type == "LIGAND" else "POLAR")
            cb_cls = PRODIGY_CLASS.get(cnb, "APOLAR" if rb.entity_type == "LIGAND" else "POLAR")
            pair_cls = "_".join(sorted([ca_cls, cb_cls]))
            if pair_cls == "APOLAR_CHARGED":
                pair_cls = "CHARGED_APOLAR"
            elif pair_cls == "APOLAR_POLAR":
                pair_cls = "POLAR_APOLAR"
            if pair_cls in prodigy_ic:
                prodigy_ic[pair_cls] += w_ij

            pair_pi = 0
            ring_a = get_aromatic_centroid(ra)
            ring_b = get_aromatic_centroid(rb)
            if ring_a is not None and ring_b is not None:
                if math.dist(ring_a, ring_b) <= 5.5:
                    pair_pi += 1

            eff_sb = min(pair_sb, 2)
            eff_hb = min(pair_hb, 2)
            eff_hp = min(pair_hp, 6)
            eff_hp_weighted = min(hp_dist_sum, 6.0)
            eff_clash = min(pair_clash, 3)

            if w_ij >= 0.25:
                n_sb += eff_sb
                n_hb += eff_hb
                n_hp += eff_hp
                n_pi += pair_pi
                n_ss += pair_ss
                n_clash += eff_clash

            w_sb += w_ij * eff_sb
            w_hb += w_ij * eff_hb
            w_hp += w_ij * eff_hp_weighted
            w_pi += w_ij * pair_pi
            w_clash += w_ij * eff_clash'''

new_inner_loop = '''            min_heavy_dist = 999.0
            pair_sb, pair_hb, pair_hp, pair_ss, pair_clash = 0, 0, 0, 0, 0
            pair_hb_base, pair_hb_2oh, pair_hb_bb = 0, 0, 0
            hp_dist_sum = 0.0
            cna, cnb = ra.canonical_name, rb.canonical_name
            is_apolar_a = (cna in HYDROPHOBIC_RES) or (ra.entity_type == "LIGAND")
            is_apolar_b = (cnb in HYDROPHOBIC_RES) or (rb.entity_type == "LIGAND")

            for atom_a in ra.atoms:
                ca_xyz = atom_a.coord
                for atom_b in rb.atoms:
                    d = math.dist(ca_xyz, atom_b.coord)
                    if d < min_heavy_dist:
                        min_heavy_dist = d
                    if d > 4.5:
                        continue

                    if cna == "CYS" and cnb == "CYS" and atom_a.name == "SG" and atom_b.name == "SG" and d <= 2.3:
                        pair_ss += 1
                        continue

                    if d < 2.2:
                        pair_clash += 1

                    is_cat_a = (cna, atom_a.name) in CATIONIC_ATOMS
                    is_ani_a = (cna, atom_a.name) in ANIONIC_ATOMS or (ra.entity_type == "NUCLEIC" and atom_a.name in NUCLEIC_PHOSPHATE_O)
                    is_cat_b = (cnb, atom_b.name) in CATIONIC_ATOMS
                    is_ani_b = (cnb, atom_b.name) in ANIONIC_ATOMS or (rb.entity_type == "NUCLEIC" and atom_b.name in NUCLEIC_PHOSPHATE_O)

                    if d <= 4.0 and ((is_cat_a and is_ani_b) or (is_ani_a and is_cat_b)):
                        pair_sb += 1
                        continue

                    if d <= 3.5 and atom_a.element in ("N", "O", "S") and atom_b.element in ("N", "O", "S"):
                        pair_hb += 1
                        if is_nucleic_pair:
                            nuc_atom = atom_a if ra.entity_type == "NUCLEIC" else atom_b
                            if nuc_atom.name == "O2'":
                                pair_hb_2oh += 1
                            elif nuc_atom.name in NUCLEOBASE_HB_ATOMS:
                                pair_hb_base += 1
                            else:
                                pair_hb_bb += 1
                        continue

                    if d <= 4.5 and atom_a.element in ("C", "S") and atom_b.element in ("C", "S"):
                        if is_apolar_a and is_apolar_b and atom_a.name not in ("C", "CA") and atom_b.name not in ("C", "CA"):
                            pair_hp += 1
                            hp_dist_sum += clamp((4.75 - d) / 0.80, 0.35, 1.0)
                        elif is_nucleic_pair:
                            prot_is_a = (ra.entity_type == "PROTEIN")
                            p_at, n_at = (atom_a, atom_b) if prot_is_a else (atom_b, atom_a)
                            p_ap = is_apolar_a if prot_is_a else is_apolar_b
                            if p_ap and p_at.name not in ("C", "CA") and n_at.name in NUCLEIC_CARBON_ATOMS:
                                pair_hp += 1
                                hp_dist_sum += clamp((4.75 - d) / 0.80, 0.35, 1.0)

            if min_heavy_dist > 5.5:
                continue

            intf_res_a.add(ra.key)
            intf_res_b.add(rb.key)
            if cb_dist > 8.0:
                intf_pae_values.append(pair_pae)
                pair_weights_seen.append(w_ij)

            ca_cls = PRODIGY_CLASS.get(cna, "APOLAR" if ra.entity_type == "LIGAND" else "POLAR")
            cb_cls = PRODIGY_CLASS.get(cnb, "APOLAR" if rb.entity_type == "LIGAND" else "POLAR")
            pair_cls = "_".join(sorted([ca_cls, cb_cls]))
            if pair_cls == "APOLAR_CHARGED":
                pair_cls = "CHARGED_APOLAR"
            elif pair_cls == "APOLAR_POLAR":
                pair_cls = "POLAR_APOLAR"
            if pair_cls in prodigy_ic:
                prodigy_ic[pair_cls] += w_ij

            pair_pi = 0
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
            eff_hb_base = min(pair_hb_base, 2)
            eff_hb_2oh = min(pair_hb_2oh, 2)
            eff_hb_bb = min(pair_hb_bb, 2)
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
assert old_inner_loop in src, "Missing inner residue-pair loop anchor"
src = src.replace(old_inner_loop, new_inner_loop, 1)

# Log RNA-specific contact fields into residue_contact_log
old_log_head = '''            residue_contact_log.append({
                "res_a": f"{ra.chain_id}:{ra.res_name}{ra.res_seq}",
                "res_b": f"{rb.chain_id}:{rb.res_name}{rb.res_seq}",'''
new_log_head = '''            residue_contact_log.append({
                "res_a": f"{ra.chain_id}:{ra.res_name}{ra.res_seq}",
                "res_b": f"{rb.chain_id}:{rb.res_name}{rb.res_seq}",
                "cat_pi": pair_catpi,
                "hb_base": eff_hb_base,
                "hb_2oh": eff_hb_2oh,
                "hb_bb": eff_hb_bb,'''
assert old_log_head in src, "Missing residue_contact_log anchor"
src = src.replace(old_log_head, new_log_head, 1)

# Fix crystallographic Proline inversion
old_pro_strain = 'is_strained = (res.mean_bfactor < 80.0) or (_mean_pae_tmp > 4.0)'
new_pro_strain = 'is_strained = (residue_raw_plddt(res) < 80.0) or (_mean_pae_tmp > 4.0)'
assert old_pro_strain in src, "Missing is_strained Proline anchor"
src = src.replace(old_pro_strain, new_pro_strain, 1)

# Add calibrated Protein-RNA / Protein-DNA thermodynamics & PRODIGY auto-disable
old_thermo = '''    e_raw = s_out * (0.110 * w_hp + 0.195 * prodigy_ic["APOLAR_APOLAR"] + 0.0032 * w_bsa) + (0.78 * w_sb * debye_k * ph_k + 0.48 * w_hb + 0.62 * w_pi)
    e_eff = 15.5 * math.tanh(e_raw / 15.5)
    dg_polar_desolv = +0.12 * prodigy_ic["POLAR_POLAR"]
    dg_steric_clash = +1.50 * w_clash
    dg_proline_kink = +2.20 * helical_prolines
    dg_phys = dg_cratic - e_eff + dg_polar_desolv + dg_steric_clash + dg_proline_kink

    eff_ic_sum = sum(prodigy_ic.values())
    if eff_ic_sum < 3.0:
        contact_gate = 0.0
    else:
        contact_gate = 1.0 / (1.0 + math.exp(-0.40 * (eff_ic_sum - 11.0)))
        nis_reliability = clamp(nis_total / 20.0, 0.0, 1.0)
        contact_gate = contact_gate * nis_reliability'''

new_thermo = '''    if is_nucleic_pair:
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
        contact_gate = 0.0
    else:
        contact_gate = 1.0 / (1.0 + math.exp(-0.40 * (eff_ic_sum - 11.0)))
        nis_reliability = clamp(nis_total / 20.0, 0.0, 1.0)
        contact_gate = contact_gate * nis_reliability'''
assert old_thermo in src, "Missing thermodynamics block anchor"
src = src.replace(old_thermo, new_thermo, 1)

old_pat_call = '''    pattern_id, pattern_desc = classify_mechanistic_pattern(
        uis_rounded, float(iptm_val), lis_score, dg_uni, w_bsa, helical_prolines, n_clash
    )'''
new_pat_call = '''    pattern_id, pattern_desc = classify_mechanistic_pattern(
        uis_rounded, float(iptm_val), lis_score, dg_uni, w_bsa, 0 if is_nucleic_pair else helical_prolines, n_clash
    )
    if mode != "MODE_A_PRED_FULL" and pattern_id in ("Pattern A", "Pattern B"):
        pattern_id = "Pattern A"
        pattern_desc = f"Coordinate-Mode Supported {modality} Interactor (Interface packing and dG_uni support binding; no PAE JSON supplied)"'''
assert old_pat_call in src, "Missing classify_mechanistic_pattern call anchor"
src = src.replace(old_pat_call, new_pat_call, 1)

old_ret_top = '''    return {
        "chain_pair": f"{chain_a_id}-{chain_b_id}",'''
new_ret_top = '''    rna_tot_spec = n_sb + n_hb_base + n_hb_2oh + n_hb_bb + n_pi + n_catpi
    rna_ssi = round((n_hb_base + n_pi + n_catpi) / max(1, rna_tot_spec), 3) if is_nucleic_pair else 0.0
    return {
        "chain_pair": f"{chain_a_id}-{chain_b_id}",
        "modality": modality,
        "hb_base": n_hb_base,
        "hb_2oh": n_hb_2oh,
        "hb_backbone": n_hb_bb,
        "cation_pi": n_catpi,
        "rna_specificity_idx": rna_ssi,'''
assert old_ret_top in src, "Missing return dict anchor"
src = src.replace(old_ret_top, new_ret_top, 1)

old_empty_ret = '''        return {
            "chain_pair": f"{chain_a_id}-{chain_b_id}",'''
if old_empty_ret in src:
    src = src.replace(old_empty_ret, '''        return {
            "chain_pair": f"{chain_a_id}-{chain_b_id}",
            "modality": modality,''', 1)

# ------------------------------------------------------------------------------
# 5. REMOVE BROKEN LEGACY MAIN() POST-OVERRIDES (L1522..L1644)
# ------------------------------------------------------------------------------
idx_override_start = src.find('                if res_pair["wBSA_A2"] >= 800.0 or res_pair["mode"] != "MODE_A_PRED_FULL":')
idx_override_end = src.find('                results.append(res_pair)', idx_override_start)
assert idx_override_start > 0 and idx_override_end > idx_override_start, "Could not locate legacy main() override block"
src = src[:idx_override_start] + src[idx_override_end:]

# ------------------------------------------------------------------------------
# 6. UPGRADE SHADOW build_residues() FOR MULTI-CHAIN HETATM DRUG/METAL SPLITTING
# ------------------------------------------------------------------------------
old_build_wrap = '''def build_residues(atoms):
    # Universal Auto-Rescue 1: Same-chain Receptor + HETATM Ligand/Ion
    ordered, chains = _orig_build_residues_fn(atoms)
    if len(chains) == 1:
        only_cid = next(iter(chains.keys()))
        all_res = chains[only_cid]
        poly_res = [r for r in all_res if _universal_res_etype(r) in ("PROTEIN", "NUCLEIC")]
        het_res = [r for r in all_res if _universal_res_etype(r) in ("LIGAND", "ION")]
        if poly_res and het_res:
            lig_cid = f"{only_cid}_LIG"
            chains[only_cid] = poly_res
            chains[lig_cid] = []
            for hr in het_res:
                _universal_set_chain(hr, lig_cid)
                chains[lig_cid].append(hr)
        elif len(poly_res) >= 10:
            # Universal Auto-Rescue 2: Single-chain monomer -> N-term vs C-term intra-domain packing
            mid = len(all_res) // 2
            n_cid, c_cid = f"{only_cid}_N", f"{only_cid}_C"
            chains = {n_cid: [], c_cid: []}
            for idx_r, r in enumerate(all_res):
                tcid = n_cid if idx_r < mid else c_cid
                _universal_set_chain(r, tcid)
                chains[tcid].append(r)
    return ordered, chains'''

new_build_wrap = '''def build_residues(atoms):
    # Universal Multi-Chain HETATM Splitter (v2.0):
    #   Splits biological small-molecule drugs/cofactors & catalytic metals out of ANY
    #   polymer chain (e.g. 5ZC9 B:RCG600 -> 'B_RCG600', A:ANP500 -> 'A_ANP500')
    #   while preserving orig_chain_id for AF3 PAE indexing.
    ordered, raw_chains = _orig_build_residues_fn(atoms)
    for r in ordered:
        if not hasattr(r, "orig_chain_id"):
            r.orig_chain_id = r.chain_id
    chains: Dict[str, List[Residue]] = {}
    for cid, rlist in raw_chains.items():
        poly_res = [r for r in rlist if _universal_res_etype(r) in ("PROTEIN", "NUCLEIC")]
        het_res = [r for r in rlist if _universal_res_etype(r) in ("LIGAND", "ION")]
        if poly_res and het_res:
            chains[cid] = poly_res
            for hr in het_res:
                if hr.res_name in CRYSTAL_BUFFER_JUNK:
                    continue
                if hr.entity_type == "ION" and hr.res_name not in CATALYTIC_METALS:
                    continue
                new_cid = f"{cid}_{hr.res_name}{hr.res_seq}"
                _universal_set_chain(hr, new_cid)
                chains.setdefault(new_cid, []).append(hr)
        else:
            kept = [
                r for r in rlist
                if r.res_name not in CRYSTAL_BUFFER_JUNK
                and not (r.entity_type == "ION" and r.res_name not in CATALYTIC_METALS)
            ]
            if kept:
                chains[cid] = kept
    kept_set = {id(r) for cl in chains.values() for r in cl}
    ordered = [r for r in ordered if id(r) in kept_set]
    if len(chains) == 1:
        only_cid = next(iter(chains.keys()))
        all_res = chains[only_cid]
        if len(all_res) >= 10 and all(_universal_res_etype(r) in ("PROTEIN", "NUCLEIC") for r in all_res):
            mid = len(all_res) // 2
            n_cid, c_cid = f"{only_cid}_N", f"{only_cid}_C"
            chains = {n_cid: [], c_cid: []}
            for idx_r, r in enumerate(all_res):
                tcid = n_cid if idx_r < mid else c_cid
                _universal_set_chain(r, tcid)
                chains[tcid].append(r)
    return ordered, chains'''
assert old_build_wrap in src, "Missing shadow build_residues wrapper anchor"
src = src.replace(old_build_wrap, new_build_wrap, 1)

src = src.replace(
    'key_to_tok.get((r.chain_id, r.res_seq), 0)',
    'key_to_tok.get((getattr(r, "orig_chain_id", r.chain_id), r.res_seq), key_to_tok.get((r.chain_id, r.res_seq), 0))'
)

# Enhance CLI Step-1 derivation printout for RNA/DNA modality
old_cli_step1 = 'print(f"1. Chain Lengths & Contacts (d_min <= 5.5 A): Chain A = {top.get(\'len_a\', 0)} aa, Chain B = {top.get(\'len_b\', 0)} aa | Ordered Core Loci (w>=0.25): A:{top.get(\'loci_a\', \'None\')} <-> B:{top.get(\'loci_b\', \'None\')} | Confident Contacts: {top.get(\'conf_contacts\', 0)} ({top.get(\'raw_contacts\', 0)} raw) across {top.get(\'intf_residues\', 0)} interface residues (SB={top.get(\'salt_bridges\', 0)}, HB={top.get(\'h_bonds\', 0)}, HP={top.get(\'hydrophobic_contacts\', 0)})")'
new_cli_step1 = '''_mod = top.get("modality", "PROTEIN-PROTEIN")
        _rna_extra = f", pi={top.get('pi_stacking', 0)}, Cat-pi={top.get('cation_pi', 0)}, HB_base={top.get('hb_base', 0)}, HB_2OH={top.get('hb_2oh', 0)}, HB_bb={top.get('hb_backbone', 0)}, SSI={top.get('rna_specificity_idx', 0.0):.2f}" if "RNA" in _mod or "DNA" in _mod else ""
        print(f"1. Modality [{_mod}] & Contacts (d_min <= 5.5 A): Chain A = {top.get('len_a', 0)}, Chain B = {top.get('len_b', 0)} | Ordered Core Loci (w>=0.25): {top.get('loci_a', 'None')} <-> {top.get('loci_b', 'None')} | Confident Contacts: {top.get('conf_contacts', 0)} ({top.get('raw_contacts', 0)} raw) across {top.get('intf_residues', 0)} interface residues (SB={top.get('salt_bridges', 0)}, HB={top.get('h_bonds', 0)}, HP={top.get('hydrophobic_contacts', 0)}{_rna_extra})")'''
if old_cli_step1 in src:
    src = src.replace(old_cli_step1, new_cli_step1, 1)

ast.parse(src)
cli_path.write_text(src, encoding="utf-8")
print("[PASS] Upgraded unippi_cli.py cleanly (ast.parse OK).\n")

tests = [
    ("ZERO-REGRESSION PPI GUARD (fold_vip1_bzip61.zip)", "benchmarks/fold_vip1_bzip61.zip", "/tmp/v2_after_vip1"),
    ("AFTER BENCHMARK 1: 5ZC9.cif (eIF4A1 + Polypurine RNA + Rocaglamide A)", "benchmarks/5ZC9.cif", "/tmp/v2_after_5zc9"),
    ("AFTER BENCHMARK 2: 1URN.cif (U1A RRM + Hairpin RNA — ITC Ground Truth dG ~ -13.8 kcal/mol)", "benchmarks/1URN.cif", "/tmp/v2_after_1urn"),
]
for title, inp, outp in tests:
    print("\n" + "=" * 95)
    print(title)
    print("=" * 95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", inp, "-o", outp], check=True)
