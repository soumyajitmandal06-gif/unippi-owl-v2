#!/usr/bin/env python3
"""Automated physics, edge-case & regression verification suite for UniPPI-OWL v2.2."""
import os
import tempfile
import unippi_cli as u

def _mk_atom(serial, name, res_name, chain, res_seq, x, y, z, el, bfac=90.0, het=False, alt=""):
    return u.Atom(serial, name, alt, res_name, chain, res_seq, "", x, y, z, 1.0, bfac, el, het)

def run_tests():
    # 1. Version & Halogen Fluorine check
    assert u.MODEL_VERSION == "UniPPI-v2.2", f"Bad version: {u.MODEL_VERSION}"
    assert "F" not in u.HALOGEN_ELEMENTS and {"CL", "BR", "I"} == u.HALOGEN_ELEMENTS
    print("[PASS 1] Model version (UniPPI-v2.2) & Halogens (Cl/Br/I only, F excluded)")

    # 2. Modified RNA nucleotide stays NUCLEIC (not LIGAND) & Calcium stays CA (not C)
    assert u.infer_element("CA", "", "CA") == "CA", "Calcium ion misclassified as Carbon!"
    assert u.infer_element("CA", "", "ALA") == "C", "Alpha carbon misclassified!"
    r_psu = u.Residue("B", 1, "", "PSU", False, [_mk_atom(1, "C4'", "PSU", "B", 1, 0, 0, 0, "C")])
    assert r_psu.canonical_name == "U" and r_psu.entity_type == "NUCLEIC", f"PSU broken: {r_psu.entity_type}"
    print("[PASS 2] Modified RNA (PSU -> NUCLEIC) & Calcium (CA -> CA)")

    # 3. altLoc 'B'-only residue is preserved (not silently dropped)
    alt_atoms = [
        _mk_atom(1, "CA", "ALA", "A", 1, 0, 0, 0, "C", alt="B"),
        _mk_atom(2, "CA", "ALA", "A", 2, 3, 0, 0, "C", alt="A"),
        _mk_atom(3, "CA", "ALA", "A", 2, 3.5, 0, 0, "C", alt="B"),
    ]
    dedup = u.deduplicate_altlocs(alt_atoms)
    assert len(dedup) == 2 and dedup[0].alt_loc == "B" and dedup[1].alt_loc == "A"
    print("[PASS 3] altLoc deduplication preserves 'B'-only residues & prefers 'A' when both exist")

    # 4. H-bond donor-acceptor typing rejects O...O & ASP...GLU, accepts SER...ASP
    r_asp = u.Residue("A", 1, "", "ASP", False, [_mk_atom(1, "OD1", "ASP", "A", 1, 0, 0, 0, "O")])
    r_glu = u.Residue("B", 1, "", "GLU", False, [_mk_atom(2, "OE1", "GLU", "B", 1, 2.9, 0, 0, "O")])
    r_ser = u.Residue("B", 2, "", "SER", False, [_mk_atom(3, "OG", "SER", "B", 2, 2.9, 0, 0, "O")])
    r_arg = u.Residue("A", 2, "", "ARG", False, [_mk_atom(4, "NH1", "ARG", "A", 2, 0, 0, 0, "N")])
    r_lys = u.Residue("B", 3, "", "LYS", False, [_mk_atom(5, "NZ", "LYS", "B", 3, 3.0, 0, 0, "N")])
    assert not u.is_valid_hbond_pair(r_asp, r_asp.atoms[0], r_glu, r_glu.atoms[0], set(), set(), set(), set())
    assert not u.is_valid_hbond_pair(r_arg, r_arg.atoms[0], r_lys, r_lys.atoms[0], set(), set(), set(), set())
    assert u.is_valid_hbond_pair(r_asp, r_asp.atoms[0], r_ser, r_ser.atoms[0], set(), set(), set(), set())
    print("[PASS 4] H-bond typing rejects Acceptor-Acceptor & Donor-Donor, accepts valid D-A")

    # 5. Endosomal Histidine pH switch: HIS-ASP is H-bond at pH 7.2, Salt Bridge at pH 5.5
    r_his = u.Residue("A", 1, "", "HIS", False, [
        _mk_atom(1, "CA", "HIS", "A", 1, -3.0, 0, 0, "C"),
        _mk_atom(2, "NE2", "HIS", "A", 1, -1.4, 0, 0, "N"),
    ])
    r_asp2 = u.Residue("B", 1, "", "ASP", False, [
        _mk_atom(3, "CA", "ASP", "B", 1, 3.0, 0, 0, "C"),
        _mk_atom(4, "OD1", "ASP", "B", 1, 1.4, 0, 0, "O"),
    ])
    ordered = [r_his, r_asp2]
    r2i = {r.key: i for i, r in enumerate(ordered)}
    sc = u.compute_residue_sasa(ordered)
    sa, sb = u.compute_residue_sasa([r_his]), u.compute_residue_sasa([r_asp2])
    res_72 = u.analyze_interface_pair("A", "B", [r_his], [r_asp2], r2i, "MODE_C_EXP_CRYSTAL", {}, sc, sa, sb, ph_val=7.2)
    res_55 = u.analyze_interface_pair("A", "B", [r_his], [r_asp2], r2i, "MODE_C_EXP_CRYSTAL", {}, sc, sa, sb, ph_val=5.5)
    assert res_72["salt_bridges"] == 0 and res_72["h_bonds"] == 1, f"pH 7.2 HIS wrong: {res_72}"
    assert res_55["salt_bridges"] == 1 and res_55["dG_uni_kcal"] < res_72["dG_uni_kcal"], f"pH 5.5 HIS switch failed: {res_55}"
    print(f"[PASS 5] Histidine pH switch verified: pH 7.2 (SB=0, HB=1, dG={res_72['dG_uni_kcal']}) -> pH 5.5 (SB=1, dG={res_55['dG_uni_kcal']})")

    # 6. Empty non-interacting decoy no longer crashes with NameError
    r_far_a = u.Residue("A", 1, "", "ALA", False, [_mk_atom(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    r_far_b = u.Residue("B", 1, "", "ALA", False, [_mk_atom(2, "CA", "ALA", "B", 1, 80, 0, 0, "C")])
    res_empty = u.analyze_interface_pair(
        "A", "B", [r_far_a], [r_far_b], {r_far_a.key: 0, r_far_b.key: 1},
        "MODE_C_EXP_CRYSTAL", {}, {r_far_a.key: 100.0, r_far_b.key: 100.0},
        {r_far_a.key: 100.0}, {r_far_b.key: 100.0}, force_report_if_empty=True
    )
    assert res_empty["UIS"] == 0.0 and res_empty["conf_contacts"] == 0 and res_empty["dG_uni_kcal"] == 2.10
    print("[PASS 6] Zero-contact binary decoy returns cleanly without NameError")

    # 7. Covalent bond (1.8 A C-S) flagged as covalent_bonds, NOT steric_clashes
    r_cys = u.Residue("A", 1, "", "CYS", False, [_mk_atom(1, "SG", "CYS", "A", 1, 0, 0, 0, "S")])
    r_lig = u.Residue("L", 1, "", "LIG", True, [
        _mk_atom(2, "C1", "LIG", "L", 1, 1.80, 0, 0, "C", het=True),
        _mk_atom(3, "C2", "LIG", "L", 1, 3.20, 0, 0, "C", het=True),
    ])
    sc_cov = u.compute_residue_sasa([r_cys, r_lig])
    res_cov = u.analyze_interface_pair(
        "A", "L", [r_cys], [r_lig], {r_cys.key: 0, r_lig.key: 1},
        "MODE_C_EXP_CRYSTAL", {}, sc_cov, u.compute_residue_sasa([r_cys]), u.compute_residue_sasa([r_lig])
    )
    assert res_cov["covalent_bonds"] == 1 and res_cov["steric_clashes"] == 0, f"Covalent misclassified: {res_cov}"
    print("[PASS 7] Covalent inhibitor bond (1.8 A C-S) counted as covalent_bonds=1, steric_clashes=0")

    # 8. Mode A vs Mode B scale check: deleting PAE JSON never inflates UIS
    with tempfile.TemporaryDirectory() as td:
        u.generate_benchmark_suite(td)
        recs_full = u.process_structure_file(os.path.join(td, "1_VIP1_Homodimer_Positive", "1_VIP1_Homodimer_Positive_model_0.cif"), td)
        os.remove(os.path.join(td, "1_VIP1_Homodimer_Positive", "1_VIP1_Homodimer_Positive_confidences_0.json"))
        recs_nojson = u.process_structure_file(os.path.join(td, "1_VIP1_Homodimer_Positive", "1_VIP1_Homodimer_Positive_model_0.cif"), td)
        assert recs_full[0]["mode"] == "MODE_A_PRED_FULL" and recs_nojson[0]["mode"] == "MODE_B_PRED_COORD"
        assert recs_nojson[0]["UIS"] <= recs_full[0]["UIS"], f"Mode B ({recs_nojson[0]['UIS']}) > Mode A ({recs_full[0]['UIS']})!"
        print(f"[PASS 8] Mode A ({recs_full[0]['UIS']}) >= Mode B without JSON ({recs_nojson[0]['UIS']})")

    # 9. Real benchmark smoke test (if benchmarks/ exists in repo)
    if os.path.isdir("benchmarks"):
        print("\n--- Running v2.2 CLI on local benchmarks/ ---")
        with tempfile.TemporaryDirectory() as td:
            for fn in sorted(os.listdir("benchmarks")):
                fp = os.path.join("benchmarks", fn)
                if not fn.lower().endswith((".cif", ".pdb", ".zip")):
                    continue
                files, tmpz = u._collect_struct_files(fp)
                recs = []
                for sf in files:
                    recs.extend(u.process_structure_file(sf, td))
                if tmpz:
                    import shutil; shutil.rmtree(tmpz, ignore_errors=True)
                recs.sort(key=lambda x: x["UIS"], reverse=True)
                if recs:
                    t = recs[0]
                    print(f"  {fn:<24} | Top={t['chain_pair']:<12} [{t['modality']:<15}] | UIS={t['UIS']:.3f} | dG={t['dG_uni_kcal']:>6.2f} | NIS={t['nis_total']:<3} | SB={t['salt_bridges']} HB={t['h_bonds']} Met={t['metal_coord_bonds']} | {t['pattern_id']}")

    print("\nALL PHASE 1A TESTS PASSED.")

if __name__ == "__main__":
    run_tests()
