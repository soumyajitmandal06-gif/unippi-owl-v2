#!/usr/bin/env python3
"""
UniPPI-OWL CLI v2.2 (Layer 1): Universal Multi-Modal Structure-to-Affinity & Interface Auditing Engine
100% Pure Python Standard Library (Zero External Dependencies - Runs on any Python 3.8+)
"""
import re

import os
import sys
import json
import math
import glob
import shlex
import shutil
import hashlib
import argparse
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set

MODEL_VERSION = "UniPPI-v1.0"

UIS_DECISION_THRESHOLD = 0.5000

UIS_HIGH_CONF_THRESHOLD = 0.6500

R_GAS_KCAL = 1.98720425864083e-3  # kcal / (mol * K)

Vec3 = Tuple[float, float, float]

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "MG": 1.73, "ZN": 1.39, "CA": 2.31, "FE": 1.40, "CO": 1.40, "MN": 1.40, "CU": 1.40, "NI": 1.63, "SE": 1.90, "NA": 2.27, "K": 2.75
}

WATER_PROBE_RADIUS = 1.40

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}

MODIFIED_AA_MAP = {
    "SEP": "SER", "TPO": "THR", "PTR": "TYR", "MSE": "MET",
    "CSO": "CYS", "HIP": "HIS", "HID": "HIS", "HIE": "HIS",
    "MLY": "LYS", "M3L": "LYS", "ALY": "LYS", "HYP": "PRO"
}

NUCLEOTIDES = {
    "A", "C", "G", "U", "DA", "DC", "DG", "DT", "ADE", "CYT", "GUA", "THY", "URA"
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

PRODIGY_CLASS = {
    "ARG": "CHARGED", "LYS": "CHARGED", "ASP": "CHARGED", "GLU": "CHARGED",
    "SEP": "CHARGED", "TPO": "CHARGED", "PTR": "CHARGED",
    "CYS": "POLAR", "HIS": "POLAR", "ASN": "POLAR", "GLN": "POLAR",
    "SER": "POLAR", "THR": "POLAR", "TRP": "POLAR",
    "ALA": "APOLAR", "PHE": "APOLAR", "GLY": "APOLAR", "ILE": "APOLAR",
    "LEU": "APOLAR", "VAL": "APOLAR", "MET": "APOLAR", "PRO": "APOLAR", "TYR": "APOLAR"
}

# Tien et al. (2013) PLoS ONE theoretical MaxASA (A^2) for PRODIGY RSA >= 5% NIS definition
MAX_ASA_TIEN2013 = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0
}

COVALENT_PAIR_ELEMENTS = {
    frozenset(("C", "S")), frozenset(("C", "N")), frozenset(("C", "O")),
    frozenset(("S", "S")), frozenset(("O", "P")), frozenset(("C", "C"))
}
PHOSPHO_PTM_RES = {"SEP", "TPO", "PTR"}
PHOSPHO_PTM_O = {"O1P", "O2P", "O3P", "OP1", "OP2", "OP3"}

KD_HYDROPATHY = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5
}

CATIONIC_ATOMS = {
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
    ("LYS", "NZ"), ("HIS", "ND1"), ("HIS", "NE2")
}

ANIONIC_ATOMS = {
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")
}

NUCLEIC_PHOSPHATE_O = {"OP1", "OP2", "O1P", "O2P"}
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
NUCLEIC_CARBON_ATOMS = {"C1'", "C2'", "C3'", "C4'", "C5'", "C2", "C4", "C5", "C6", "C7", "C8"}
NUC_TO_ONE = {"A": "a", "C": "c", "G": "g", "U": "u", "T": "t", "DA": "a", "DC": "c", "DG": "g", "DT": "t", "ADE": "a", "CYT": "c", "GUA": "g", "URA": "u", "THY": "t"}
MODIFIED_NUC_MAP = {
    "I": "G", "INO": "G", "PSU": "U", "H2U": "U", "5MU": "U", "4SU": "U", "OMU": "U",
    "6MA": "A", "1MA": "A", "A2M": "A", "M6A": "A", "MIA": "A", "T6A": "A",
    "5MC": "C", "OMC": "C", "4AC": "C", "3MC": "C",
    "OMG": "G", "7MG": "G", "G7M": "G", "2MG": "G", "M2G": "G", "1MG": "G", "YYG": "G",
    "DI": "DG", "5CM": "DC", "BRU": "DT"
}
HALOGEN_ELEMENTS = {"CL", "BR", "I"}

HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}

AROMATIC_RINGS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"}
}

CHAIN_CHAR_POOL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

def clamp(val: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, val))

def mean_list(vals: List[float], default: float = 0.0) -> float:
    return (sum(vals) / len(vals)) if vals else default

def voxel_key(coord: Vec3, cell_size: float) -> Tuple[int, int, int]:
    return (
        int(math.floor(coord[0] / cell_size)),
        int(math.floor(coord[1] / cell_size)),
        int(math.floor(coord[2] / cell_size))
    )

@dataclass
class Atom:
    serial: int
    name: str
    alt_loc: str
    res_name: str
    chain_id: str
    res_seq: int
    ins_code: str
    x: float
    y: float
    z: float
    occupancy: float
    b_factor: float
    element: str
    is_hetatm: bool

    @property
    def coord(self) -> Vec3:
        return (self.x, self.y, self.z)

@dataclass
class Residue:
    chain_id: str
    res_seq: int
    ins_code: str
    res_name: str
    is_hetatm: bool
    atoms: List[Atom] = field(default_factory=list)
    _centroid_cache: Optional[Vec3] = field(default=None, repr=False)
    _cb_cache: Optional[Vec3] = field(default=None, repr=False)

    @property
    def key(self) -> Tuple[str, int, str]:
        return (self.chain_id, self.res_seq, self.ins_code)

    @property
    def canonical_name(self) -> str:
        return MODIFIED_AA_MAP.get(self.res_name, MODIFIED_NUC_MAP.get(self.res_name, self.res_name))

    @property
    def entity_type(self) -> str:
        cname = self.canonical_name
        if cname in STANDARD_AA:
            return "PROTEIN"
        if cname in NUCLEOTIDES:
            return "NUCLEIC"
        if len(self.atoms) == 1 and self.atoms[0].element in VDW_RADII:
            return "ION"
        return "LIGAND"

    @property
    def mean_bfactor(self) -> float:
        if not self.atoms:
            return 0.0
        return sum(a.b_factor for a in self.atoms) / len(self.atoms)

    @property
    def centroid(self) -> Vec3:
        if self._centroid_cache is not None:
            return self._centroid_cache
        for a in self.atoms:
            if a.name in ("CA", "C4'"):
                self._centroid_cache = a.coord
                return self._centroid_cache
        n = len(self.atoms)
        if n == 0:
            return (0.0, 0.0, 0.0)
        self._centroid_cache = (
            sum(a.x for a in self.atoms) / n,
            sum(a.y for a in self.atoms) / n,
            sum(a.z for a in self.atoms) / n
        )
        return self._centroid_cache

    @property
    def cb_coord(self) -> Vec3:
        if self._cb_cache is not None:
            return self._cb_cache
        for a in self.atoms:
            if a.name == "CB":
                self._cb_cache = a.coord
                return self._cb_cache
        self._cb_cache = self.centroid
        return self._cb_cache

def compute_pair_sha256(residues_a: List[Residue], residues_b: List[Residue]) -> str:
    seq_a = "".join(THREE_TO_ONE.get(r.canonical_name, NUC_TO_ONE.get(r.canonical_name, r.canonical_name)) for r in residues_a)
    seq_b = "".join(THREE_TO_ONE.get(r.canonical_name, NUC_TO_ONE.get(r.canonical_name, r.canonical_name)) for r in residues_b)
    norm_pair = "::".join(sorted([seq_a, seq_b]))
    return hashlib.sha256(norm_pair.encode("utf-8")).hexdigest()[:16]

def infer_element(atom_name: str, raw_el: str = "", res_name: str = "") -> str:
    el = (raw_el or "").strip().upper()
    r_clean = (res_name or "").strip().upper()
    if el and el.isalpha():
        if el == "CA" and r_clean and r_clean != "CA":
            return "C"
        return el
    clean = "".join(ch for ch in (atom_name or "") if ch.isalpha()).upper()
    if clean == "CA" and r_clean != "CA":
        return "C"
    if len(clean) >= 2 and clean[:2] in {"CL", "BR", "MG", "ZN", "FE", "CA", "NA", "CO", "MN", "CU", "NI", "SE"}:
        if clean[:2] == "CA" and r_clean != "CA":
            return "C"
        if clean[:2] in {"CO", "MN", "CU", "NI", "FE", "ZN", "MG", "NA"} and r_clean in STANDARD_AA:
            return clean[:1]
        return clean[:2]
    return clean[:1] if clean else "C"

def deduplicate_altlocs(atoms: List[Atom]) -> List[Atom]:
    """Keep one conformer per (chain, resSeq, ins, isHet, resName, atomName), preferring '', 'A', '1' else first available."""
    rank_map = {"": 0, "A": 1, "1": 2}
    best: Dict[Tuple[str, int, str, bool, str, str], Tuple[int, int, Atom]] = {}
    for idx, a in enumerate(atoms):
        k = (a.chain_id, a.res_seq, a.ins_code, a.is_hetatm, a.res_name, a.name)
        r = rank_map.get(a.alt_loc, 10)
        if k not in best or r < best[k][0]:
            best[k] = (r, idx, a)
    return [item[2] for item in sorted(best.values(), key=lambda x: x[1])]

def parse_pdb_file(filepath: str) -> Tuple[List[Atom], Dict[str, str]]:
    atoms: List[Atom] = []
    meta = {"method": "UNKNOWN", "has_qa_header": "FALSE"}
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ENDMDL") and atoms:
                break
            if line.startswith("EXPDTA"):
                meta["method"] = line[6:].strip()
            elif line.startswith("REMARK") and ("ALPHAFOLD" in line.upper() or "PLDDT" in line.upper()):
                meta["has_qa_header"] = "TRUE"
            elif line.startswith(("ATOM  ", "HETATM")):
                is_het = line.startswith("HETATM")
                res_name = line[17:20].strip().upper()
                if res_name in CRYSTAL_BUFFER_JUNK:
                    continue
                alt_loc = line[16:17].strip()
                try:
                    serial = int(line[6:11].strip() or 0)
                    name = line[12:16].strip()
                    if res_name == "GLY" and name == "CB":
                        continue
                    chain_id = line[21:22].strip() or "A"
                    res_seq = int(line[22:26].strip() or 0)
                    ins_code = line[26:27].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    occ = float(line[54:60].strip() or 1.0)
                    bfac = float(line[60:66].strip() or 0.0)
                    raw_el = line[76:78].strip() if len(line) >= 78 else ""
                    element = infer_element(name, raw_el, res_name)
                    if element == "H":
                        continue
                    atoms.append(Atom(
                        serial, name, alt_loc, res_name, chain_id, res_seq,
                        ins_code, x, y, z, occ, bfac, element, is_het
                    ))
                except ValueError:
                    continue
    return deduplicate_altlocs(atoms), meta

def parse_mmcif_file(filepath: str) -> Tuple[List[Atom], Dict[str, str]]:
    atoms: List[Atom] = []
    meta = {"method": "UNKNOWN", "has_qa_header": "FALSE"}
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    in_atom_loop = False
    atom_headers: List[str] = []

    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            if in_atom_loop and atoms:
                in_atom_loop = False
            continue
        if s.startswith("_exptl.method") and s.split()[0] == "_exptl.method":
            parts_ex = s.split(None, 1)
            if len(parts_ex) == 2 and parts_ex[1].strip("'\"") not in ("?", "."):
                meta["method"] = parts_ex[1].strip("'\"")
        if s.startswith("_ma_qa_metric"):
            meta["has_qa_header"] = "TRUE"
        if s == "loop_":
            in_atom_loop = False
            atom_headers = []
            continue
        if s.startswith("_atom_site."):
            in_atom_loop = True
            atom_headers.append(s.split(".")[1].strip())
            continue
        if in_atom_loop and atom_headers:
            if s.startswith("_"):
                in_atom_loop = False
                continue
            tokens = shlex.split(s)
            if len(tokens) != len(atom_headers):
                continue
            row = dict(zip(atom_headers, tokens))
            group_pdb = row.get("group_PDB", "ATOM")
            is_het = (group_pdb == "HETATM")
            res_name = row.get("label_comp_id", row.get("auth_comp_id", "UNK")).upper()
            if res_name in CRYSTAL_BUFFER_JUNK:
                continue
            alt_loc = row.get("label_alt_id", ".")
            if alt_loc in (".", "?"):
                alt_loc = ""
            model_num = row.get("pdbx_PDB_model_num", "1")
            if model_num not in ("1", ".", "?"):
                continue
            try:
                serial = int(row.get("id", "0"))
                name = row.get("label_atom_id", row.get("auth_atom_id", "CA")).strip()
                if len(name) >= 2 and name[0] == name[-1] and name[0] in ("'", '"'):
                    name = name[1:-1]
                if res_name == "GLY" and name == "CB":
                    continue
                chain_id = row.get("auth_asym_id", row.get("label_asym_id", "A"))
                res_seq_str = row.get("auth_seq_id", row.get("label_seq_id", "1"))
                res_seq = int(res_seq_str) if res_seq_str not in (".", "?") else 1
                ins_code = row.get("pdbx_PDB_ins_code", "")
                if ins_code in (".", "?"):
                    ins_code = ""
                x = float(row["Cartn_x"])
                y = float(row["Cartn_y"])
                z = float(row["Cartn_z"])
                occ_str = row.get("occupancy", "1.0")
                occ = float(occ_str) if occ_str not in (".", "?") else 1.0
                bfac_str = row.get("B_iso_or_equiv", "0.0")
                bfac = float(bfac_str) if bfac_str not in (".", "?") else 0.0
                raw_el = row.get("type_symbol", "")
                element = infer_element(name, raw_el, res_name)
                if element == "H":
                    continue
                atoms.append(Atom(
                    serial, name, alt_loc, res_name, chain_id, res_seq,
                    ins_code, x, y, z, occ, bfac, element, is_het
                ))
            except (ValueError, KeyError):
                continue

    return deduplicate_altlocs(atoms), meta

def _set_residue_chain(r: Residue, new_cid: str) -> None:
    r.chain_id = new_cid
    for a in r.atoms:
        a.chain_id = new_cid

def build_residues(atoms: List[Atom], auto_split_monomer: bool = True) -> Tuple[List[Residue], Dict[str, List[Residue]]]:
    res_map: Dict[Tuple[str, int, str, bool, str], Residue] = {}
    ordered_raw: List[Residue] = []
    for atom in atoms:
        key = (atom.chain_id, atom.res_seq, atom.ins_code, atom.is_hetatm, atom.res_name if atom.is_hetatm else "")
        if key not in res_map:
            r = Residue(atom.chain_id, atom.res_seq, atom.ins_code, atom.res_name, atom.is_hetatm)
            r.orig_chain_id = atom.chain_id
            res_map[key] = r
            ordered_raw.append(r)
        res_map[key].atoms.append(atom)

    raw_chains: Dict[str, List[Residue]] = {}
    for r in ordered_raw:
        raw_chains.setdefault(r.chain_id, []).append(r)

    chains: Dict[str, List[Residue]] = {}
    for cid, rlist in raw_chains.items():
        poly_res = [r for r in rlist if r.entity_type in ("PROTEIN", "NUCLEIC")]
        het_res = [r for r in rlist if r.entity_type in ("LIGAND", "ION")]
        if poly_res and het_res:
            chains[cid] = poly_res
            for hr in het_res:
                if hr.res_name in CRYSTAL_BUFFER_JUNK:
                    continue
                if hr.entity_type == "ION" and hr.res_name not in CATALYTIC_METALS:
                    continue
                new_cid = f"{cid}_{hr.res_name}{hr.res_seq}"
                _set_residue_chain(hr, new_cid)
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
    ordered = [r for r in ordered_raw if id(r) in kept_set]

    if auto_split_monomer and len(chains) == 1:
        only_cid = next(iter(chains.keys()))
        all_res = chains[only_cid]
        if len(all_res) >= 10 and all(r.entity_type in ("PROTEIN", "NUCLEIC") for r in all_res):
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
                mid = n_tot // 2
            n_cid, c_cid = f"{only_cid}_N", f"{only_cid}_C"
            chains = {n_cid: [], c_cid: []}
            for idx_r, r in enumerate(all_res):
                tcid = n_cid if idx_r < mid else c_cid
                _set_residue_chain(r, tcid)
                chains[tcid].append(r)

    return ordered, chains

def parse_companion_json(struct_path: str, ordered_residues: List[Residue]) -> Dict:
    folder = os.path.dirname(os.path.abspath(struct_path))
    stem = os.path.splitext(os.path.basename(struct_path))[0]
    import re as _re
    m_idx_match = _re.search(r"_model_(\d+)$", stem, flags=_re.IGNORECASE) or _re.search(r"_(\d+)$", stem)
    m_idx_str = m_idx_match.group(1) if m_idx_match else None
    base_prefix = _re.sub(r"_model_\d+$", "", stem, flags=_re.IGNORECASE)

    all_jsons = [
        p for p in sorted(glob.glob(os.path.join(folder, "*.json")))
        if "job_request" not in os.path.basename(p).lower()
        and "unippi_full_report" not in os.path.basename(p).lower()
        and "terms_of_use" not in os.path.basename(p).lower()
    ]
    prefix_matched = [
        p for p in all_jsons
        if base_prefix.lower() in os.path.basename(p).lower()
        and (m_idx_str is None or f"_{m_idx_str}.json" in os.path.basename(p).lower() or f"model_{m_idx_str}" in os.path.basename(p).lower())
    ]
    json_candidates = prefix_matched if prefix_matched else all_jsons

    conf_data = {
        "has_json": False,
        "iptm": None,
        "ptm": None,
        "ranking_score": None,
        "chain_pair_iptm": {},
        "orig_chain_order": list(dict.fromkeys(getattr(r, "orig_chain_id", r.chain_id) for r in ordered_residues)),
        "res_pae_matrix": None
    }

    for jpath in json_candidates:
        try:
            with open(jpath, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            continue

        if not isinstance(data, dict):
            continue

        if "iptm" in data and data["iptm"] is not None:
            conf_data["iptm"] = float(data["iptm"])
            conf_data["has_json"] = True
        if "ptm" in data and data["ptm"] is not None:
            conf_data["ptm"] = float(data["ptm"])
            conf_data["has_json"] = True
        if "ranking_score" in data and data["ranking_score"] is not None:
            conf_data["ranking_score"] = float(data["ranking_score"])
        if "complex_iptm" in data and data["complex_iptm"] is not None and conf_data["iptm"] is None:
            conf_data["iptm"] = float(data["complex_iptm"])
            conf_data["has_json"] = True
        if "chain_pair_iptm" in data and isinstance(data["chain_pair_iptm"], list):
            conf_data["chain_pair_iptm"] = data["chain_pair_iptm"]
        elif "pair_chains_iptm" in data and isinstance(data["pair_chains_iptm"], dict):
            _pc = data["pair_chains_iptm"]
            _ks = sorted(_pc.keys(), key=lambda x: int(x) if str(x).isdigit() else str(x))
            conf_data["chain_pair_iptm"] = [[float(_pc[r].get(c, 0.0)) for c in _ks] for r in _ks]

        raw_pae = data.get("pae", data.get("predicted_aligned_error", None))
        if raw_pae is not None and isinstance(raw_pae, list) and len(raw_pae) > 0:
            n_res = len(ordered_residues)
            n_tok = len(raw_pae)

            if n_tok == n_res and len(raw_pae[0]) == n_res:
                conf_data["res_pae_matrix"] = raw_pae
                conf_data["has_json"] = True
            elif "token_chain_ids" in data and "token_res_ids" in data:
                t_chains = data["token_chain_ids"]
                t_res = data["token_res_ids"]
                key_to_toks: Dict[Tuple[str, int], List[int]] = {}
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
                        aligned_pae.append(row_out)
                conf_data["res_pae_matrix"] = aligned_pae
                conf_data["has_json"] = True

    return conf_data

def detect_confidence_mode(ordered_residues: List[Residue], meta: Dict[str, str], conf_data: Dict) -> str:
    if conf_data["has_json"] and conf_data["res_pae_matrix"] is not None:
        return "MODE_A_PRED_FULL"

    method = meta.get("method", "").upper()
    if any(k in method for k in ("X-RAY", "DIFFRACTION", "ELECTRON MICROSCOPY", "NMR")):
        return "MODE_C_EXP_CRYSTAL"

    bfacs = [r.mean_bfactor for r in ordered_residues if r.atoms]
    if not bfacs:
        return "MODE_B_PRED_COORD"

    max_b, min_b, mean_b = max(bfacs), min(bfacs), mean_list(bfacs)
    if 0.0 < max_b <= 1.0 and mean_b > 0.35:
        for r in ordered_residues:
            for a in r.atoms:
                a.b_factor *= 100.0
        return "MODE_B_PRED_COORD"

    if meta.get("has_qa_header") == "TRUE" or (min_b >= 0.0 and max_b <= 100.0 and mean_b >= 45.0):
        return "MODE_B_PRED_COORD"

    return "MODE_C_EXP_CRYSTAL"

def generate_golden_spiral_sphere(n_points: int = 64) -> List[Vec3]:
    pts: List[Vec3] = []
    golden_angle = math.pi * (1.0 + math.sqrt(5.0))
    for i in range(n_points):
        idx = i + 0.5
        phi = math.acos(1.0 - 2.0 * idx / n_points)
        theta = golden_angle * idx
        sin_p = math.sin(phi)
        pts.append((math.cos(theta) * sin_p, math.sin(theta) * sin_p, math.cos(phi)))
    return pts

UNIT_SPHERE_64 = generate_golden_spiral_sphere(64)

def compute_residue_sasa(residues: List[Residue]) -> Dict[Tuple[str, int, str], float]:
    all_coords: List[Vec3] = []
    all_radii: List[float] = []
    atom_to_res: List[Tuple[str, int, str]] = []

    for r in residues:
        rk = r.key
        for a in r.atoms:
            all_coords.append((a.x, a.y, a.z))
            all_radii.append(VDW_RADII.get(a.element, 1.70) + WATER_PROBE_RADIUS)
            atom_to_res.append(rk)

    n_atoms = len(all_coords)
    res_sasa: Dict[Tuple[str, int, str], float] = {r.key: 0.0 for r in residues}
    if n_atoms == 0:
        return res_sasa

    cell_size = 8.0
    grid: Dict[Tuple[int, int, int], List[int]] = {}
    voxels: List[Tuple[int, int, int]] = []
    for idx, c in enumerate(all_coords):
        vk = voxel_key(c, cell_size)
        voxels.append(vk)
        grid.setdefault(vk, []).append(idx)

    n_sphere = len(UNIT_SPHERE_64)
    four_pi_over_n = (4.0 * math.pi) / n_sphere

    for i in range(n_atoms):
        cx, cy, cz = all_coords[i]
        ri = all_radii[i]
        vx, vy, vz = voxels[i]

        close_nbrs: List[Tuple[float, float, float, float]] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((vx + dx, vy + dy, vz + dz), ()):
                        if j == i:
                            continue
                        rj = all_radii[j]
                        if math.dist(all_coords[i], all_coords[j]) < (ri + rj):
                            nx, ny, nz = all_coords[j]
                            close_nbrs.append((nx, ny, nz, rj * rj))

        if not close_nbrs:
            res_sasa[atom_to_res[i]] += 4.0 * math.pi * (ri * ri)
            continue

        exposed_count = 0
        for ux, uy, uz in UNIT_SPHERE_64:
            px = cx + ri * ux
            py = cy + ri * uy
            pz = cz + ri * uz
            buried = False
            for nx, ny, nz, rj_sq in close_nbrs:
                ddx = px - nx
                ddy = py - ny
                ddz = pz - nz
                if (ddx * ddx + ddy * ddy + ddz * ddz) < rj_sq:
                    buried = True
                    break
            if not buried:
                exposed_count += 1

        if exposed_count > 0:
            res_sasa[atom_to_res[i]] += four_pi_over_n * exposed_count * (ri * ri)

    return res_sasa

def get_aromatic_centroid(res: Residue) -> Optional[Vec3]:
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

def detect_ligand_planar_rings(res: Residue) -> List[Tuple[Vec3, Vec3]]:
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
    return anionic

def get_ligand_cationic_nitrogens(res: Residue) -> Set[str]:
    """Universal 3D detector for protonated cationic ligand nitrogens (aliphatic amines & guanidiniums)."""
    if res.entity_type != "LIGAND":
        return set()
    if hasattr(res, "_cached_cationic_n"):
        return res._cached_cationic_n
    cationic: Set[str] = set()
    atoms = res.atoms
    if not atoms:
        res._cached_cationic_n = cationic
        return cationic
    # Identify atoms lying inside planar 5/6-membered aromatic rings (to exclude neutral anilines & heterocycles)
    rings = detect_ligand_planar_rings(res)
    arom_c_names: Set[str] = set()
    if rings:
        for a in atoms:
            if a.element != "C":
                continue
            for r_cent, r_norm in rings:
                if math.dist(a.coord, r_cent) <= 1.75:
                    d_pl = abs((a.coord[0]-r_cent[0])*r_norm[0] + (a.coord[1]-r_cent[1])*r_norm[1] + (a.coord[2]-r_cent[2])*r_norm[2])
                    if d_pl <= 0.25:
                        arom_c_names.add(a.name)
                        break
    n_atoms = [a for a in atoms if a.element == "N"]
    c_atoms = [a for a in atoms if a.element == "C"]
    # 1. Aliphatic protonated amines (bonded strictly to non-carbonyl, non-aromatic sp3 carbons)
    for na in n_atoms:
        nbrs = [ot for ot in atoms if ot is not na and 1.10 <= math.dist(na.coord, ot.coord) <= 1.65]
        if not (1 <= len(nbrs) <= 4):
            continue
        if not all(ot.element == "C" and math.dist(na.coord, ot.coord) >= 1.38 for ot in nbrs):
            continue
        is_amide_or_aniline = False
        for cb in nbrs:
            if cb.name in arom_c_names:
                is_amide_or_aniline = True
                break
            for ot in atoms:
                if ot is cb or ot is na:
                    continue
                d_cb = math.dist(cb.coord, ot.coord)
                if (ot.element == "O" and d_cb <= 1.32) or (ot.element == "S" and d_cb <= 1.75):
                    is_amide_or_aniline = True
                    break
            if is_amide_or_aniline:
                break
        if not is_amide_or_aniline:
            cationic.add(na.name)
    # 2. Guanidinium / amidinium groups (C bonded to >=2 N at 1.25-1.42 A and 0 O)
    for ca in c_atoms:
        if ca.name in arom_c_names:
            continue
        bn = [na for na in n_atoms if 1.25 <= math.dist(ca.coord, na.coord) <= 1.42]
        bo = [ot for ot in atoms if ot.element == "O" and math.dist(ca.coord, ot.coord) <= 1.50]
        if len(bn) >= 2 and not bo:
            for na in bn:
                deg = sum(1 for ot in atoms if ot is not na and 1.10 <= math.dist(na.coord, ot.coord) <= 1.65)
                if deg == 1:
                    cationic.add(na.name)
    res._cached_cationic_n = cationic
    return cationic

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
            d_cut = 3.70
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
        if "PROTEIN" in all_t:
            return "PROTEIN-LIGAND"
        if "ION" in all_t:
            return "LIGAND-METAL"
        return "LIGAND-LIGAND"
    if "ION" in all_t:
        return "PROTEIN-METAL" if "PROTEIN" in all_t else "ION-ION"
    return "PROTEIN-PROTEIN"

def get_atom_charge_at_ph(
    res: Residue, atom: Atom, ph_val: float, lig_ani: Set[str], lig_cat: Set[str]
) -> Tuple[float, float]:
    """Per-residue Henderson-Hasselbalch positive/negative charge magnitudes at solution pH."""
    cname = res.canonical_name
    an = atom.name
    if res.res_name in PHOSPHO_PTM_RES and an in PHOSPHO_PTM_O:
        return 0.0, 1.0 + 1.0 / (1.0 + 10.0 ** (6.1 - ph_val))
    if res.entity_type == "NUCLEIC" and an in NUCLEIC_PHOSPHATE_O:
        return 0.0, 1.0 / (1.0 + 10.0 ** (1.0 - ph_val))
    if cname == "ARG" and an in ("NE", "NH1", "NH2"):
        return 1.0 / (1.0 + 10.0 ** (ph_val - 12.5)), 0.0
    if cname == "LYS" and an == "NZ":
        return 1.0 / (1.0 + 10.0 ** (ph_val - 10.5)), 0.0
    if cname == "HIS" and an in ("ND1", "NE2"):
        return 1.0 / (1.0 + 10.0 ** (ph_val - 6.0)), 0.0
    if cname == "ASP" and an in ("OD1", "OD2"):
        return 0.0, 1.0 / (1.0 + 10.0 ** (3.9 - ph_val))
    if cname == "GLU" and an in ("OE1", "OE2"):
        return 0.0, 1.0 / (1.0 + 10.0 ** (4.1 - ph_val))
    if res.entity_type == "LIGAND":
        if an in lig_cat:
            return 1.0 / (1.0 + 10.0 ** (ph_val - 9.8)), 0.0
        if an in lig_ani:
            return 0.0, 1.0 / (1.0 + 10.0 ** (3.5 - ph_val))
    return 0.0, 0.0

def get_hbond_role(res: Residue, atom: Atom, lig_ani: Set[str], lig_cat: Set[str]) -> str:
    """Return 'D' (donor-only), 'A' (acceptor-only), 'DA' (ambivalent), or '' (non-H-bonding)."""
    if atom.element not in ("N", "O", "S"):
        return ""
    an = atom.name
    cname = res.canonical_name
    if res.entity_type == "PROTEIN":
        if res.res_name in PHOSPHO_PTM_RES and an in PHOSPHO_PTM_O:
            return "A"
        if an == "N":
            return "" if cname == "PRO" else "D"
        if an in ("O", "OXT"):
            return "A"
        if (cname == "ARG" and an in ("NE", "NH1", "NH2")) or (cname == "LYS" and an == "NZ") or (cname == "TRP" and an == "NE1"):
            return "D"
        if (cname == "ASN" and an == "ND2") or (cname == "GLN" and an == "NE2"):
            return "D"
        if (cname == "ASP" and an in ("OD1", "OD2")) or (cname == "GLU" and an in ("OE1", "OE2")):
            return "A"
        if (cname == "ASN" and an == "OD1") or (cname == "GLN" and an == "OE1") or (cname == "MET" and an == "SD"):
            return "A"
        if (cname == "SER" and an == "OG") or (cname == "THR" and an == "OG1") or (cname == "TYR" and an == "OH"):
            return "DA"
        if (cname == "HIS" and an in ("ND1", "NE2")) or (cname == "CYS" and an == "SG"):
            return "DA"
        return "DA"
    if res.entity_type == "NUCLEIC":
        if an == "O2'":
            return "DA"
        if an in NUCLEIC_BB_ATOMS:
            return "A"
        if an == "N9" or (an == "N1" and cname in ("C", "DC", "CYT", "U", "URA", "T", "DT", "THY")):
            return ""
        if (cname in ("A", "DA", "ADE") and an == "N6") or (cname in ("C", "DC", "CYT") and an == "N4"):
            return "D"
        if (cname in ("G", "DG", "GUA") and an in ("N1", "N2")) or (cname in ("U", "URA", "T", "DT", "THY") and an == "N3"):
            return "D"
        if an in NUCLEOBASE_HB_ATOMS:
            return "A"
        return "DA"
    if res.entity_type == "LIGAND":
        if an in lig_ani:
            return "A"
        if an in lig_cat:
            return "D"
        return "DA"
    return ""

def is_valid_hbond_pair(
    ra: Residue, aa: Atom, rb: Residue, ab: Atom,
    lig_ani_a: Set[str], lig_cat_a: Set[str], lig_ani_b: Set[str], lig_cat_b: Set[str]
) -> bool:
    role_a = get_hbond_role(ra, aa, lig_ani_a, lig_cat_a)
    role_b = get_hbond_role(rb, ab, lig_ani_b, lig_cat_b)
    if not role_a or not role_b:
        return False
    if role_a == "D" and role_b == "D":
        return False
    if role_a == "A" and role_b == "A":
        return False
    return True

def is_in_helical_segment(res_idx: int, chain_residues: List[Residue], mode: str = "MODE_C_EXP_CRYSTAL") -> bool:
    if res_idx < 4 or res_idx >= len(chain_residues) - 4:
        return False
    r_cur = chain_residues[res_idx]
    if mode in ("MODE_A_PRED_FULL", "MODE_B_PRED_COORD") and r_cur.mean_bfactor < 50.0:
        return False
    r_prev = chain_residues[res_idx - 4]
    r_next = chain_residues[res_idx + 4]
    if (r_next.res_seq - r_cur.res_seq) != 4 or (r_cur.res_seq - r_prev.res_seq) != 4:
        return False
    d_next = math.dist(r_next.centroid, r_cur.centroid)
    d_prev = math.dist(r_prev.centroid, r_cur.centroid)
    d_span = math.dist(r_prev.centroid, r_next.centroid)
    return (5.0 <= d_next <= 7.0) and (5.0 <= d_prev <= 7.0) and (d_span >= 10.2)

def evaluate_custom_mutation_ddg(
    wt_name: str,
    mut_name: str,
    c_info: Dict[str, float],
    f_buried: float,
    in_helix: bool
) -> float:
    if wt_name == mut_name:
        return 0.0
    if wt_name == "PRO" and in_helix and mut_name != "PRO":
        return -1.85

    hb_capable = {"SER", "THR", "ASN", "GLN", "TYR", "TRP", "HIS", "ARG", "LYS", "ASP", "GLU", "CYS"}
    kd_wt = KD_HYDROPATHY.get(wt_name, 0.0)
    kd_mut = KD_HYDROPATHY.get(mut_name, 0.0)

    # Truncating to ALA/GLY or replacing an H-bonding/salt-bridging residue cannot gain fake apolar burial
    if mut_name in ("ALA", "GLY") or c_info["sb"] > 0 or c_info["hb"] > 0:
        ddg_desolv = 0.52 * max(0.0, kd_wt - kd_mut) * f_buried
    else:
        ddg_desolv = 0.52 * (kd_wt - kd_mut) * f_buried

    cls_wt = PRODIGY_CLASS.get(wt_name, "POLAR")
    cls_mut = PRODIGY_CLASS.get(mut_name, "POLAR")

    ddg_polar = 0.0
    if c_info["sb"] > 0:
        if cls_mut != "CHARGED":
            ddg_polar += 0.95 * c_info["sb"]
        elif (wt_name in ("ARG", "LYS") and mut_name in ("ASP", "GLU")) or (wt_name in ("ASP", "GLU") and mut_name in ("ARG", "LYS")):
            ddg_polar += 1.85 * c_info["sb"]
    if c_info["hb"] > 0 and mut_name not in hb_capable:
        ddg_polar += 0.85 * c_info["hb"]

    ddg_apolar = 0.0
    if wt_name in HYDROPHOBIC_RES and mut_name not in HYDROPHOBIC_RES:
        ddg_apolar += 0.28 * c_info["hp"]
    elif wt_name in ("LEU", "ILE", "PHE", "TRP", "MET", "TYR") and mut_name == "ALA":
        ddg_apolar += 0.28 * c_info["hp"] + 0.45 * f_buried
    if wt_name in AROMATIC_RINGS and mut_name not in AROMATIC_RINGS:
        ddg_apolar += 0.75 * c_info["pi"] + 0.35 * f_buried

    ddg_helix = 0.0
    if mut_name == "PRO":
        ddg_helix = 3.10 if in_helix else 0.80
    elif mut_name == "GLY" and in_helix:
        ddg_helix = 0.65

    ddg_pmim = -1.15 if (wt_name in {"SER", "THR"} and mut_name in {"ASP", "GLU"} and c_info.get("hb", 0) > 0 and c_info.get("sb", 0) == 0) else 0.0
    return round(ddg_desolv + ddg_polar + ddg_apolar + ddg_helix + ddg_pmim, 2)

def classify_mechanistic_pattern(
    uis: float,
    iptm_val: float,
    lis_score: float,
    dg_uni: float,
    w_bsa: float,
    helical_prolines: int,
    n_clash: int,
    proline_penalty: float = 0.0
) -> Tuple[str, str]:
    """Mechanistic Pattern Classification (Patterns A, B, C, D) with reachable clash/kink audit."""
    kink_penalized = (proline_penalty > 0.0 and helical_prolines > 0)
    if kink_penalized or n_clash >= 5:
        return (
            "Pattern C",
            "Thermodynamically Penalized / Steric-Disrupted Interface (Geometric proximity offset by steric clash or kink penalty)"
        )
    if uis >= UIS_DECISION_THRESHOLD and iptm_val >= 0.60 and dg_uni <= -5.8:
        return (
            "Pattern A",
            "Concordant Rigid/Coiled-Coil or Globular Interactor (Global ipTM, Local LIS, and dG_uni all support binding)"
        )
    if (uis >= UIS_DECISION_THRESHOLD or (lis_score >= 0.50 and dg_uni <= -5.8)) and iptm_val < 0.60 and dg_uni <= -5.8:
        return (
            "Pattern B",
            "IDR-Diluted Local Interactor (Strong local interface LIS & dG_uni despite low global ipTM from disordered tails)"
        )
    if uis >= UIS_DECISION_THRESHOLD:
        return (
            "Pattern A" if iptm_val >= 0.60 else "Pattern B",
            "Above-Threshold Interactor (Moderate affinity; inspect contact density & PAE)"
        )
    if w_bsa >= 250.0 and dg_uni > -5.0:
        return (
            "Pattern C",
            "Weak / Sub-Threshold Packing (Geometric burial without sufficient stabilizing contacts)"
        )
    return (
        "Pattern D",
        "Below-Threshold Decoy or Disordered Hallucination (No physical binding support)"
    )

def analyze_interface_pair(
    chain_a_id: str,
    chain_b_id: str,
    residues_a: List[Residue],
    residues_b: List[Residue],
    res_to_global_idx: Dict[Tuple[str, int, str], int],
    mode: str,
    conf_data: Dict,
    sasa_complex: Dict[Tuple[str, int, str], float],
    sasa_free_a: Dict[Tuple[str, int, str], float],
    sasa_free_b: Dict[Tuple[str, int, str], float],
    chain_remap: Optional[Dict[str, str]] = None,
    temp_k: float = 298.15,
    lambda_pae: float = 8.0,
    ionic_mM: float = 150.0,
    ph_val: float = 7.2,
    cutoff_a: float = 5.5,
    proline_penalty: float = 0.0,
    force_report_if_empty: bool = False,
    custom_mutations: Optional[List[str]] = None
) -> Optional[Dict]:
    if not residues_a or not residues_b:
        return None

    active_cutoff = clamp(cutoff_a, 3.5, 8.0)
    pair_sha256 = compute_pair_sha256(residues_a, residues_b)

    voxel_size = 14.0
    grid_b: Dict[Tuple[int, int, int], List[int]] = {}
    for j, rb in enumerate(residues_b):
        vk = voxel_key(rb.centroid, voxel_size)
        grid_b.setdefault(vk, []).append(j)

    poly_bfacs = [_r.mean_bfactor for _r in (residues_a + residues_b) if _r.entity_type in ("PROTEIN", "NUCLEIC") and _r.atoms]
    b_min_global = min(poly_bfacs) if poly_bfacs else min((_r.mean_bfactor for _r in (residues_a + residues_b) if _r.atoms), default=10.0)
    b_min_chain: Dict[str, float] = {}
    for _r in (residues_a + residues_b):
        if _r.entity_type in ("LIGAND", "ION"):
            b_min_chain[_r.chain_id] = b_min_global
        elif _r.chain_id not in b_min_chain or _r.mean_bfactor < b_min_chain[_r.chain_id]:
            b_min_chain[_r.chain_id] = _r.mean_bfactor

    def residue_raw_plddt(res: Residue) -> float:
        b = 76.75 if (mode == "MODE_B_PRED_COORD" and res.entity_type in ("LIGAND", "ION") and res.mean_bfactor <= 0.01) else res.mean_bfactor
        if mode in ("MODE_A_PRED_FULL", "MODE_B_PRED_COORD"):
            return clamp(b, 0.0, 100.0)
        b_ref = max(b_min_chain.get(res.chain_id, 10.0) + 20.0, 60.0)
        return clamp(100.0 * math.exp(-max(0.0, b - b_ref) / 250.0), 70.0, 100.0)

    def residue_plateau_conf(res: Residue) -> float:
        if mode in ("MODE_A_PRED_FULL", "MODE_B_PRED_COORD"):
            return 0.85 if (mode == "MODE_B_PRED_COORD" and res.entity_type in ("LIGAND", "ION") and res.mean_bfactor <= 0.01) else clamp((res.mean_bfactor - 30.0) / 55.0, 0.02, 1.0)
        b_ref = max(b_min_chain.get(res.chain_id, 10.0) + 20.0, 45.0)
        return 1.0

    pae_mat = conf_data.get("res_pae_matrix", None)

    def pair_weight(ra: Residue, rb: Residue) -> Tuple[float, float]:
        ca = residue_plateau_conf(ra)
        cb = residue_plateau_conf(rb)
        if pae_mat is not None:
            ia = res_to_global_idx[ra.key]
            ib = res_to_global_idx[rb.key]
            mean_pae = 0.5 * (float(pae_mat[ia][ib]) + float(pae_mat[ib][ia]))
            excess_pae = max(0.0, mean_pae - 3.5)
            w = ca * cb * math.exp(-excess_pae / lambda_pae)
            return w, mean_pae
        return ca * cb, 3.5

    modality = detect_pair_modality(chain_a_id, chain_b_id, residues_a, residues_b)
    is_intra_domain = (modality == "INTRA-DOMAIN (MONOMER)")
    is_nucleic_pair = modality in ("PROTEIN-RNA", "PROTEIN-DNA", "RNA-RNA", "DNA-DNA", "RNA-DNA")
    is_duplex_pair = modality in ("RNA-RNA", "DNA-DNA", "RNA-DNA")
    is_ligand_pair = modality in ("PROTEIN-LIGAND", "RNA-LIGAND", "DNA-LIGAND", "TERNARY-GLUE", "LIGAND-LIGAND")
    is_metal_pair = modality in ("PROTEIN-METAL", "RNA-METAL", "DNA-METAL", "LIGAND-METAL", "METAL-BRIDGE")
    is_non_ppi = is_nucleic_pair or is_ligand_pair or is_metal_pair

    n_sb, n_hb, n_hp, n_pi, n_ss, n_cov, n_clash, n_cb8 = 0, 0, 0, 0, 0, 0, 0, 0
    n_hb_base, n_hb_2oh, n_hb_bb, n_catpi, n_xb, n_metal_coord, n_metal_outer = 0, 0, 0, 0, 0, 0, 0
    w_sb, w_hb, w_hp, w_pi, w_clash = 0.0, 0.0, 0.0, 0.0, 0.0
    w_hb_base, w_hb_2oh, w_hb_bb, w_pi_rna, w_catpi, w_xb, w_metal_coord, w_metal_outer = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    pair_weights_seen: List[float] = []

    prodigy_ic = {
        "CHARGED_CHARGED": 0.0, "CHARGED_APOLAR": 0.0,
        "POLAR_POLAR": 0.0, "POLAR_APOLAR": 0.0,
        "APOLAR_APOLAR": 0.0, "CHARGED_POLAR": 0.0
    }

    cb8_intf_a: Set[Tuple[str, int, str]] = set()
    cb8_intf_b: Set[Tuple[str, int, str]] = set()
    heavy55_intf_a: Set[Tuple[str, int, str]] = set()
    heavy55_intf_b: Set[Tuple[str, int, str]] = set()
    intf_pae_values: List[float] = []
    residue_contact_log: List[Dict] = []
    res_energy_contrib: Dict[Tuple[str, int, str], Dict[str, float]] = {}

    for ra in residues_a:
        vx, vy, vz = voxel_key(ra.centroid, voxel_size)
        cand_b: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    cand_b.extend(grid_b.get((vx + dx, vy + dy, vz + dz), ()))
        if not cand_b:
            continue

        lig_ani_a = get_ligand_anionic_oxygens(ra) if ra.entity_type == "LIGAND" else set()
        lig_cat_a = get_ligand_cationic_nitrogens(ra) if ra.entity_type == "LIGAND" else set()
        cna = ra.canonical_name
        is_apolar_a = (cna in HYDROPHOBIC_RES) or (ra.entity_type == "LIGAND")

        for j in cand_b:
            rb = residues_b[j]
            if is_intra_domain and getattr(ra, "orig_chain_id", ra.chain_id) == getattr(rb, "orig_chain_id", rb.chain_id) and abs(ra.res_seq - rb.res_seq) <= 2:
                continue
            if math.dist(ra.centroid, rb.centroid) > 16.0:
                continue

            cb_dist = math.dist(ra.cb_coord, rb.cb_coord)
            w_ij, pair_pae = pair_weight(ra, rb)
            if cb_dist <= 8.0:
                n_cb8 += 1
                cb8_intf_a.add(ra.key)
                cb8_intf_b.add(rb.key)
                intf_pae_values.append(pair_pae)
                pair_weights_seen.append(w_ij)

            min_heavy_dist = 999.0
            pair_sb_cnt, pair_sb_w = 0, 0.0
            pair_hb, pair_hp, pair_ss, pair_cov, pair_clash = 0, 0, 0, 0, 0
            pair_hb_base, pair_hb_2oh, pair_hb_bb, pair_metal_coord, pair_metal_outer = 0, 0, 0, 0, 0
            lig_ani_b = get_ligand_anionic_oxygens(rb) if rb.entity_type == "LIGAND" else set()
            lig_cat_b = get_ligand_cationic_nitrogens(rb) if rb.entity_type == "LIGAND" else set()
            hp_dist_sum = 0.0
            cnb = rb.canonical_name
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

                    is_met_a = (ra.entity_type in ("ION", "LIGAND") and atom_a.element in CATALYTIC_METALS)
                    is_met_b = (rb.entity_type in ("ION", "LIGAND") and atom_b.element in CATALYTIC_METALS)
                    if is_met_a or is_met_b:
                        m_el = atom_a.element if is_met_a else atom_b.element
                        other_el = atom_b.element if is_met_a else atom_a.element
                        d_mcut = 2.65 if m_el in ("ZN", "FE", "CU", "CO", "NI") else 2.75
                        if d <= d_mcut and other_el in ("O", "N", "S"):
                            pair_metal_coord += 1
                            continue
                        if m_el in ("MG", "CA", "MN") and d_mcut < d <= 4.15 and other_el == "O":
                            pair_metal_outer += 1
                            continue

                    if d <= 1.95 and frozenset((atom_a.element, atom_b.element)) in COVALENT_PAIR_ELEMENTS:
                        pair_cov += 1
                        continue

                    if d < 2.2:
                        pair_clash += 1

                    qpa, qna = get_atom_charge_at_ph(ra, atom_a, ph_val, lig_ani_a, lig_cat_a)
                    qpb, qnb = get_atom_charge_at_ph(rb, atom_b, ph_val, lig_ani_b, lig_cat_b)
                    sb_eff = (qpa * qnb) + (qna * qpb)
                    if d <= 4.0 and sb_eff >= 0.25:
                        pair_sb_w += sb_eff
                        pair_sb_cnt += (2 if sb_eff >= 1.40 else 1)
                        continue

                    if d <= 3.5 and not ((qpa >= 0.5 and qpb >= 0.5) or (qna >= 0.5 and qnb >= 0.5)):
                        if is_valid_hbond_pair(ra, atom_a, rb, atom_b, lig_ani_a, lig_cat_a, lig_ani_b, lig_cat_b):
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
                        elif is_nucleic_pair or (is_ligand_pair and (ra.entity_type == "NUCLEIC" or rb.entity_type == "NUCLEIC")):
                            nuc_is_a = (ra.entity_type == "NUCLEIC")
                            n_at, o_at = (atom_a, atom_b) if nuc_is_a else (atom_b, atom_a)
                            o_ap = is_apolar_b if nuc_is_a else is_apolar_a
                            if o_ap and o_at.name not in ("C", "CA") and n_at.name in NUCLEIC_CARBON_ATOMS:
                                pair_hp += 1
                                hp_dist_sum += clamp((4.75 - d) / 0.80, 0.35, 1.0)

            # PRODIGY interface contacts & NIS membership are locked strictly at 5.5 A heavy-atom distance
            if min_heavy_dist <= 5.5:
                heavy55_intf_a.add(ra.key)
                heavy55_intf_b.add(rb.key)
                ca_cls = "CHARGED" if ra.res_name in PHOSPHO_PTM_RES else ("CHARGED" if ((cna) == "HIS" and ph_val <= 6.0) else PRODIGY_CLASS.get(cna, "APOLAR" if ra.entity_type == "LIGAND" else "POLAR"))
                cb_cls = "CHARGED" if rb.res_name in PHOSPHO_PTM_RES else ("CHARGED" if ((cnb) == "HIS" and ph_val <= 6.0) else PRODIGY_CLASS.get(cnb, "APOLAR" if rb.entity_type == "LIGAND" else "POLAR"))
                pair_cls = "_".join(sorted([ca_cls, cb_cls]))
                if pair_cls == "APOLAR_CHARGED":
                    pair_cls = "CHARGED_APOLAR"
                elif pair_cls == "APOLAR_POLAR":
                    pair_cls = "POLAR_APOLAR"
                if pair_cls in prodigy_ic:
                    prodigy_ic[pair_cls] += w_ij

            eval_cutoff = max(5.5, active_cutoff)
            if min_heavy_dist > eval_cutoff:
                continue

            if cb_dist > 8.0 and min_heavy_dist <= 5.5:
                intf_pae_values.append(pair_pae)
                pair_weights_seen.append(w_ij)

            pair_pi, pair_pi_w, pair_catpi = 0, 0.0, 0
            pair_xb = count_directional_halogen_bonds(ra, rb) if is_ligand_pair else 0
            if not is_metal_pair and min_heavy_dist <= 5.5:
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
                    cat_atoms_list = []
                    if p_cand.canonical_name in ("ARG", "LYS"):
                        cat_name = "CZ" if p_cand.canonical_name == "ARG" else "NZ"
                        cat_at = next((a for a in p_cand.atoms if a.name == cat_name), None)
                        if cat_at is not None:
                            cat_atoms_list.append(cat_at)
                    elif p_cand.entity_type == "LIGAND":
                        l_cats = get_ligand_cationic_nitrogens(p_cand)
                        if l_cats:
                            cat_atoms_list.extend([a for a in p_cand.atoms if a.name in l_cats])
                    if cat_atoms_list:
                        cand_rings = get_all_rings_for_residue(r_cand)
                        if cand_rings:
                            for cat_at in cat_atoms_list:
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
                                if pair_catpi:
                                    break

            eff_sb_cnt = min(pair_sb_cnt, 2)
            eff_sb_w = min(pair_sb_w, 2.0)
            eff_hb = min(pair_hb, 3 if is_duplex_pair else 2)
            eff_hb_base = min(pair_hb_base, eff_hb)
            eff_hb_2oh = min(pair_hb_2oh, eff_hb - eff_hb_base)
            eff_hb_bb = min(pair_hb_bb, eff_hb - eff_hb_base - eff_hb_2oh)
            eff_hp = min(pair_hp, 6)
            eff_hp_weighted = min(hp_dist_sum, 6.0)
            eff_clash = min(pair_clash, 3)
            eff_mcoord = min(pair_metal_coord, 2)
            eff_mouter = 1 if (pair_metal_outer > 0 and eff_mcoord == 0) else 0

            if w_ij >= 0.25:
                n_sb += eff_sb_cnt
                n_hb += eff_hb
                n_hb_base += eff_hb_base
                n_hb_2oh += eff_hb_2oh
                n_hb_bb += eff_hb_bb
                n_hp += eff_hp
                n_pi += pair_pi
                n_catpi += pair_catpi
                n_xb += pair_xb
                n_metal_coord += eff_mcoord
                n_metal_outer += eff_mouter
                n_ss += pair_ss
                n_cov += pair_cov
                n_clash += eff_clash
                for r_key in (ra.key, rb.key):
                    rec = res_energy_contrib.setdefault(r_key, {"sb": 0.0, "hb": 0.0, "hp": 0.0, "pi": 0.0, "xb": 0.0})
                    rec["sb"] += w_ij * eff_sb_w
                    rec["hb"] += w_ij * eff_hb
                    rec["hp"] += w_ij * eff_hp_weighted
                    rec["pi"] += w_ij * (pair_pi + pair_catpi)
                    rec["xb"] += w_ij * pair_xb

            w_sb += w_ij * eff_sb_w
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
            w_metal_outer += w_ij * eff_mouter
            w_clash += w_ij * eff_clash

            residue_contact_log.append({
                "res_a": f"{ra.chain_id}:{ra.res_name}{ra.res_seq}",
                "res_b": f"{rb.chain_id}:{rb.res_name}{rb.res_seq}",
                "cat_pi": pair_catpi,
                "x_bond": pair_xb,
                "metal_coord": eff_mcoord,
                "metal_outer": eff_mouter,
                "covalent": pair_cov,
                "hb_base": eff_hb_base,
                "hb_2oh": eff_hb_2oh,
                "hb_bb": eff_hb_bb,
                "min_dist_A": round(min_heavy_dist, 2),
                "pae_A": round(pair_pae, 2),
                "weight_w": round(w_ij, 3),
                "sb": eff_sb_cnt, "hb": eff_hb, "hp": eff_hp, "pi": pair_pi, "clash": eff_clash
            })

    intf_any_a = heavy55_intf_a | cb8_intf_a
    intf_any_b = heavy55_intf_b | cb8_intf_b
    if len(intf_any_a) == 0 and len(intf_any_b) == 0:
        if not force_report_if_empty:
            return None
        kd_empty = math.exp(2.10 / (R_GAS_KCAL * temp_k))
        return {
            "chain_pair": f"{chain_a_id}-{chain_b_id}",
            "modality": modality,
            "len_a": len(residues_a),
            "len_b": len(residues_b),
            "conf_contacts": 0,
            "raw_contacts": 0,
            "loci_a": "None",
            "loci_b": "None",
            "sha256_id": pair_sha256,
            "mode": mode,
            "UIS": 0.0,
            "UIS_threshold": UIS_DECISION_THRESHOLD,
            "UIS_margin": round(-UIS_DECISION_THRESHOLD, 4),
            "threshold_call": "BELOW_THRESHOLD",
            "binding_class": "NON_INTERACTING",
            "pattern_id": "Pattern D",
            "pattern_desc": "Below-Threshold Decoy or Disordered Hallucination (No physical binding support)",
            "dG_uni_kcal": 2.10,
            "dG_prodigy_kcal": 0.0,
            "dG_phys_kcal": 2.10,
            "Kd_M": float(f"{kd_empty:.3e}"),
            "Kd_nM": round(kd_empty * 1e9, 1),
            "ipTM": round(float(conf_data.get("iptm") or 0.0), 3),
            "LIS": 0.0,
            "pDockQ_PAE": 0.0,
            "pDockQ_coord": 0.0,
            "mean_intf_pLDDT": 0.0,
            "mean_intf_PAE_A": 30.0,
            "BSA_A2": 0.0,
            "wBSA_A2": 0.0,
            "intf_residues": 0,
            "pct_disordered": round(100.0 * sum(1 for r in (residues_a + residues_b) if residue_raw_plddt(r) < 70.0) / max(1, len(residues_a) + len(residues_b)), 1),
            "nis_total": 0,
            "n_cb8_contacts": 0,
            "salt_bridges": 0,
            "h_bonds": 0,
            "hb_base": 0,
            "hb_2oh": 0,
            "hb_backbone": 0,
            "cation_pi": 0,
            "halogen_bonds": 0,
            "metal_coord_bonds": 0,
            "metal_outer_bonds": 0,
            "covalent_bonds": 0,
            "lig_burial_pct": 0.0,
            "pocket_quality": 0.0,
            "rna_specificity_idx": 0.0,
            "hydrophobic_contacts": 0,
            "pi_stacking": 0,
            "disulfides": 0,
            "steric_clashes": 0,
            "helical_prolines": 0,
            "proline_warning": False,
            "f_apolar": 0.0,
            "elec_density_1000A2": 0.0,
            "mean_w_intf": 0.0,
            "top_hotspot": "NONE",
            "top_hotspot_score": 0.0,
            "hotspots": [],
            "custom_mutations": [],
            "foldx_mutations": [],
            "contacts": []
        }

    conf_contact_log = [c for c in residue_contact_log if c["weight_w"] >= 0.25]
    core_conf_contacts = [c for c in conf_contact_log if c['min_dist_A'] <= 5.5]
    has_ordered_core = len(core_conf_contacts) > 0

    helical_prolines = 0
    _mean_pae_tmp = mean_list([c["pae_A"] for c in core_conf_contacts], mean_list(intf_pae_values, 3.5))
    _ord_intf_a = {r.res_seq for r in residues_a if r.key in intf_any_a and residue_plateau_conf(r) >= 0.50}
    _ord_intf_b = {r.res_seq for r in residues_b if r.key in intf_any_b and residue_plateau_conf(r) >= 0.50}
    for chain_list, ord_seqs in ((residues_a, _ord_intf_a), (residues_b, _ord_intf_b)):
        for idx, res in enumerate(chain_list):
            if res.canonical_name == "PRO" and residue_plateau_conf(res) >= 0.50 and is_in_helical_segment(idx, chain_list, mode):
                near = any(abs(res.res_seq - s_num) <= 2 for s_num in ord_seqs)
                is_strained = (residue_raw_plddt(res) < 80.0) or (_mean_pae_tmp > 4.0)
                if near and is_strained:
                    helical_prolines += 1

    # Consistent ordered-core SASA: prevents hallucinated disordered loops from shadowing ordered atoms
    if has_ordered_core:
        ord_res_a = [r for r in residues_a if residue_plateau_conf(r) >= 0.35] or list(residues_a)
        ord_res_b = [r for r in residues_b if residue_plateau_conf(r) >= 0.35] or list(residues_b)
        if len(ord_res_a) < len(residues_a) or len(ord_res_b) < len(residues_b):
            sasa_free_a_eff = compute_residue_sasa(ord_res_a)
            sasa_free_b_eff = compute_residue_sasa(ord_res_b)
            sasa_comp_eff = compute_residue_sasa(ord_res_a + ord_res_b)
        else:
            sasa_free_a_eff, sasa_free_b_eff, sasa_comp_eff = sasa_free_a, sasa_free_b, sasa_complex
        mean_w_intf = mean_list([c["weight_w"] for c in core_conf_contacts], 0.0)
    else:
        ord_res_a, ord_res_b = residues_a, residues_b
        sasa_free_a_eff, sasa_free_b_eff, sasa_comp_eff = sasa_free_a, sasa_free_b, sasa_complex
        mean_w_intf = mean_list(pair_weights_seen, 0.0)

    total_free_a = sum(sasa_free_a_eff.get(r.key, 0.0) for r in ord_res_a)
    total_free_b = sum(sasa_free_b_eff.get(r.key, 0.0) for r in ord_res_b)
    total_comp_ab = sum(sasa_comp_eff.get(r.key, 0.0) for r in (ord_res_a + ord_res_b))
    bsa = max(0.0, (total_free_a + total_free_b) - total_comp_ab)
    w_bsa = mean_w_intf * bsa

    # True PRODIGY NIS: non-5.5A-interface protein residues with RSA >= 5% (Tien et al. 2013 MaxASA)
    nis_apolar_count, nis_charged_count, nis_total = 0, 0, 0
    for r in (ord_res_a + ord_res_b):
        if r.entity_type != "PROTEIN":
            continue
        if r.key in heavy55_intf_a or r.key in heavy55_intf_b:
            continue
        max_asa = MAX_ASA_TIEN2013.get(r.canonical_name, 180.0)
        if (sasa_comp_eff.get(r.key, 0.0) / max_asa) >= 0.05:
            nis_total += 1
            p_cls = "CHARGED" if r.res_name in PHOSPHO_PTM_RES else ("CHARGED" if ((r.canonical_name) == "HIS" and ph_val <= 6.0) else PRODIGY_CLASS.get(r.canonical_name, "POLAR"))
            if p_cls == "APOLAR":
                nis_apolar_count += 1
            elif p_cls == "CHARGED":
                nis_charged_count += 1

    pct_nis_apolar = (100.0 * nis_apolar_count / nis_total) if nis_total > 0 else 35.0
    pct_nis_charged = (100.0 * nis_charged_count / nis_total) if nis_total > 0 else 25.0

    dg_prodigy_raw = (
        - 0.09459 * prodigy_ic["CHARGED_CHARGED"]
        - 0.10007 * prodigy_ic["CHARGED_APOLAR"]
        + 0.19577 * prodigy_ic["POLAR_POLAR"]
        - 0.22671 * prodigy_ic["POLAR_APOLAR"]
        + 0.18681 * pct_nis_apolar
        + 0.13810 * pct_nis_charged
        - 15.9433
    )

    eff_ic_sum = sum(prodigy_ic.values())
    apolar_frac = prodigy_ic["APOLAR_APOLAR"] / max(1.0, eff_ic_sum)
    i_m = clamp(ionic_mM / 1000.0, 0.01, 2.0)
    shield_f = 1.0 - 0.65 * apolar_frac
    debye_k = math.exp(-0.85 * shield_f * ((math.sqrt(i_m) / 3.04) - (math.sqrt(0.15) / 3.04)) * 3.5)
    s_out = 1.0 + (0.20 + 0.25 * apolar_frac) * (i_m - 0.15)
    dg_cratic = +2.10

    lig_burial_frac = 0.0
    if is_ligand_pair:
        lig_is_b = any(r.entity_type == "LIGAND" for r in residues_b)
        lig_res_list = ord_res_b if lig_is_b else ord_res_a
        lig_sasa_free_map = sasa_free_b_eff if lig_is_b else sasa_free_a_eff
        lig_free_tot = sum(lig_sasa_free_map.get(r.key, 0.0) for r in lig_res_list)
        lig_bound_tot = sum(sasa_comp_eff.get(r.key, 0.0) for r in lig_res_list)
        lig_burial_frac = clamp((lig_free_tot - lig_bound_tot) / max(1.0, lig_free_tot), 0.0, 1.0)

    eff_hel_pro = 0 if is_non_ppi else helical_prolines
    # Note: w_sb already incorporates exact per-residue Henderson-Hasselbalch pH ionization
    if is_nucleic_pair:
        debye_rna = math.exp(-0.85 * ((math.sqrt(i_m) / 3.04) - (math.sqrt(0.15) / 3.04)) * 3.5)
        e_raw = (
            0.60 * w_sb * debye_rna
            + 0.65 * w_hb_base
            + 0.55 * w_hb_2oh
            + 0.35 * w_hb_bb
            + 0.95 * w_pi_rna
            + 0.85 * w_catpi
            + 0.10 * w_hp
            + 0.0025 * w_bsa
        )
        e_eff = 16.5 * math.tanh(e_raw / 16.5)
        dg_phys = dg_cratic - e_eff + 1.50 * w_clash
    elif is_ligand_pair:
        e_raw = (
            0.85 * w_sb * debye_k
            + 0.65 * w_hb
            + 0.55 * w_xb
            + 0.90 * w_pi_rna
            + 0.80 * w_catpi
            + 0.14 * w_hp
            + 0.0055 * w_bsa
        )
        e_eff = 15.0 * math.tanh(e_raw / 15.0)
        dg_phys = 1.80 - e_eff + 1.50 * w_clash
    elif is_metal_pair:
        ion_res = [r for r in (residues_a + residues_b) if r.entity_type == "ION"]
        is_trans_met = any(r.canonical_name in ("ZN", "FE", "CU", "CO", "NI") for r in ion_res)
        k_coord = 2.60 if is_trans_met else 2.15
        eff_w_met = w_metal_coord + (0.50 * min(w_metal_outer, 4.0) if (n_metal_coord > 0 or n_metal_outer >= 2) else 0.0)
        if eff_w_met <= 0.0:
            dg_phys = dg_cratic + 1.50 * w_clash
        else:
            e_raw = k_coord * min(eff_w_met, 6.0) + 0.0040 * w_bsa
            dg_phys = 1.50 - 13.5 * math.tanh(e_raw / 13.5) + 1.50 * w_clash
    else:
        e_raw = s_out * (0.110 * w_hp + 0.195 * prodigy_ic["APOLAR_APOLAR"] + 0.0032 * w_bsa) + (0.78 * w_sb * debye_k + 0.48 * w_hb + 0.62 * w_pi_rna + 0.75 * w_catpi)
        e_eff = 15.5 * math.tanh(e_raw / 15.5)
        dg_phys = dg_cratic - e_eff + 0.12 * prodigy_ic["POLAR_POLAR"] + 1.50 * w_clash + proline_penalty * eff_hel_pro

    if is_non_ppi or eff_ic_sum < 3.0:
        contact_gate = 0.0
    else:
        nis_reliability = clamp(nis_total / 20.0, 0.0, 1.0)
        contact_gate = (1.0 / (1.0 + math.exp(-0.40 * (eff_ic_sum - 11.0)))) * nis_reliability

    phys_weight = clamp(0.50 + 0.20 * apolar_frac + (0.10 if (proline_penalty > 0.0 and eff_hel_pro > 0) else 0.0), 0.50, 0.70)
    dg_uni = clamp(
        contact_gate * ((1.0 - phys_weight) * dg_prodigy_raw + phys_weight * dg_phys)
        + (1.0 - contact_gate) * dg_phys,
        -22.0, +5.0
    )
    kd_molar = math.exp(dg_uni / (R_GAS_KCAL * temp_k))

    # Unified Hotspot Priority Score (100% CLI-Browser parity, gated at w_ij >= 0.25)
    hotspots: List[Dict] = []
    custom_mut_results: List[Dict] = []
    foldx_mutations: List[str] = []

    parsed_custom: Dict[Tuple[str, int], List[Tuple[str, str]]] = {}
    if custom_mutations:
        for m_str in custom_mutations:
            s_clean = m_str.strip().upper()
            if not s_clean or ":" not in s_clean:
                continue
            parts = s_clean.split(":")
            c_target, rest = parts[0], parts[1]
            digits = "".join(ch for ch in rest if ch.isdigit())
            if not digits:
                continue
            r_num = int(digits)
            wt_raw = rest.split(digits)[0].strip(":")
            wt_r = ONE_TO_THREE.get(wt_raw, wt_raw)
            tail = rest.split(digits)[-1].strip(":")
            mut_3 = ONE_TO_THREE.get(tail, tail if tail in STANDARD_AA else "ALA")
            parsed_custom.setdefault((c_target, r_num), []).append((wt_r if wt_r in STANDARD_AA else "", mut_3))

    for chain_list in (residues_a, residues_b):
        for r_idx, r in enumerate(chain_list):
            if r.entity_type not in ("PROTEIN", "NUCLEIC", "LIGAND"):
                continue
            cname = r.canonical_name
            in_helix = (r.entity_type == "PROTEIN") and is_in_helical_segment(r_idx, chain_list, mode)
            c_info = res_energy_contrib.get(r.key, {"sb": 0.0, "hb": 0.0, "hp": 0.0, "pi": 0.0, "xb": 0.0})
            free_s = sasa_free_a_eff.get(r.key, sasa_free_b_eff.get(r.key, 0.0))
            f_buried = max(0.0, (free_s - sasa_comp_eff.get(r.key, 0.0)) / free_s) if free_s > 1.0 else 0.0
            fx_chain = chain_remap.get(r.chain_id, r.chain_id[:1]) if chain_remap else r.chain_id[:1]

            if r.entity_type == "PROTEIN" and (r.chain_id.upper(), r.res_seq) in parsed_custom:
                for wt_req, target_aa in [x for x in parsed_custom[(r.chain_id.upper(), r.res_seq)] if not x[0] or x[0] == cname]:
                    ddg_val = round(evaluate_custom_mutation_ddg(cname, target_aa, c_info, f_buried, in_helix) * (1.0 if cname == "PRO" else residue_plateau_conf(r)), 2)
                    dg_mut_est = round(clamp(dg_uni + ddg_val, -22.0, 5.0), 2)
                    kd_mut_nm = round(math.exp(dg_mut_est / (R_GAS_KCAL * temp_k)) * 1e9, 2)
                    custom_mut_results.append({
                        "mutation": f"{r.chain_id}:{cname}{r.res_seq}->{target_aa}",
                        "ddG_kcal": ddg_val,
                        "dG_mut_kcal": dg_mut_est,
                        "Kd_mut_nM": kd_mut_nm,
                        "in_helix": in_helix
                    })

            priority_score = round(
                0.85 * c_info["sb"] + 0.55 * c_info["hb"] + 0.85 * c_info["pi"]
                + 0.55 * c_info["xb"] + 0.14 * c_info["hp"] + 1.15 * clamp(f_buried, 0.0, 1.0),
                2
            )
            if priority_score >= 0.35 or f_buried >= 0.10:
                tier = "CRITICAL_HOTSPOT" if priority_score >= 2.0 else ("STRONG_HOTSPOT" if priority_score >= 1.0 else "SUPPORTING")
                hotspots.append({
                    "residue": f"{r.chain_id}:{cname}{r.res_seq}",
                    "priority_score": priority_score,
                    "buried_frac": round(f_buried, 3),
                    "wSB": round(c_info["sb"], 2),
                    "wHB": round(c_info["hb"], 2),
                    "wHP": round(c_info["hp"], 2),
                    "wPi_XB": round(c_info["pi"] + c_info["xb"], 2),
                    "tier": tier
                })
                if r.entity_type == "PROTEIN" and cname not in ("ALA", "GLY") and priority_score >= 1.0:
                    foldx_mutations.append(f"{THREE_TO_ONE.get(cname, 'X')}{fx_chain}{r.res_seq}A;")

    hotspots.sort(key=lambda x: x["priority_score"], reverse=True)

    if has_ordered_core:
        conf_res_labels = {c["res_a"] for c in core_conf_contacts} | {c["res_b"] for c in core_conf_contacts}
        intf_plddts = [
            residue_raw_plddt(r)
            for r in (residues_a + residues_b)
            if f"{r.chain_id}:{r.res_name}{r.res_seq}" in conf_res_labels
        ]
        mean_intf_plddt = mean_list(intf_plddts, 50.0)
        mean_intf_pae = mean_list([c["pae_A"] for c in core_conf_contacts], 3.5)
        n_cb8_eff = max(len(core_conf_contacts), sum(1 for w in pair_weights_seen if w >= 0.25))
    else:
        intf_plddts = [
            residue_raw_plddt(r)
            for r in (residues_a + residues_b)
            if (r.key in intf_any_a or r.key in intf_any_b)
        ]
        mean_intf_plddt = mean_list(intf_plddts, 50.0)
        mean_intf_pae = mean_list(intf_pae_values, 8.0)
        n_cb8_eff = n_cb8

    # Pure Bryant et al. (2022) pDockQ (size_scale = 1.0, zero 2.5x inflation)
    x_pdockq = mean_intf_plddt * math.log10(n_cb8_eff + 1.0)
    pdockq_coord = ((0.724 / (1.0 + math.exp(-0.052 * (x_pdockq - 152.611)))) + 0.018) if n_cb8_eff > 0 else 0.0
    excess_intf_pae = max(0.0, mean_intf_pae - 3.5)
    pdockq_pae = pdockq_coord * (1.0 / (1.0 + (excess_intf_pae / 8.0) ** 2)) if mode == "MODE_A_PRED_FULL" else pdockq_coord
    pdockq_norm = clamp(pdockq_pae / 0.742, 0.0, 1.0)

    lis_score = 0.0
    if pae_mat is not None:
        idx_a = [res_to_global_idx[r.key] for r in residues_a]
        idx_b = [res_to_global_idx[r.key] for r in residues_b]
        valid_terms: List[float] = []
        for ia in idx_a:
            row_ia = pae_mat[ia]
            for ib in idx_b:
                p_ab = float(row_ia[ib])
                if p_ab <= 12.0:
                    valid_terms.append(1.0 - (p_ab / 12.0))
                p_ba = float(pae_mat[ib][ia])
                if p_ba <= 12.0:
                    valid_terms.append(1.0 - (p_ba / 12.0))
        if valid_terms:
            lis_score = sum(valid_terms) / len(valid_terms)

    iptm_val = conf_data.get("iptm", None)
    _cp_iptm = conf_data.get("chain_pair_iptm", None)
    _orig_ord = conf_data.get("orig_chain_order", None)
    if isinstance(_cp_iptm, list) and isinstance(_orig_ord, list) and residues_a and residues_b:
        _a_set = sorted({getattr(r, "orig_chain_id", r.chain_id) for r in residues_a})
        _b_set = sorted({getattr(r, "orig_chain_id", r.chain_id) for r in residues_b})
        _pair_vals = []
        for _ca_o in _a_set:
            for _cb_o in _b_set:
                if _ca_o != _cb_o and _ca_o in _orig_ord and _cb_o in _orig_ord:
                    _oi, _oj = _orig_ord.index(_ca_o), _orig_ord.index(_cb_o)
                    if _oi < len(_cp_iptm) and isinstance(_cp_iptm[_oi], list) and _oj < len(_cp_iptm[_oi]):
                        _v = _cp_iptm[_oi][_oj]
                        if isinstance(_v, (int, float)):
                            _pair_vals.append(float(_v))
        if _pair_vals:
            iptm_val = sum(_pair_vals) / len(_pair_vals)
    if is_intra_domain:
        iptm_val = None
    if iptm_val is None:
        iptm_val = pdockq_norm

    sigma_aff = clamp((-dg_uni - 5.0) / 8.5, 0.0, 1.0)
    mode_b_factor = 0.85 if mode == "MODE_B_PRED_COORD" else 1.0

    pocket_quality = 0.0
    if is_ligand_pair:
        pocket_quality = mean_w_intf * clamp((lig_burial_frac - 0.15) / 0.65, 0.0, 1.0)
        if mode == "MODE_A_PRED_FULL":
            raw_uis = 0.20 * float(iptm_val) + 0.25 * lis_score + 0.30 * pocket_quality + 0.25 * sigma_aff
        else:
            raw_uis = mode_b_factor * (0.55 * pocket_quality + 0.45 * sigma_aff)
    elif is_metal_pair:
        _eff_q_met = w_metal_coord + (0.50 * min(w_metal_outer, 4.0) if (n_metal_coord > 0 or n_metal_outer >= 2) else 0.0)
        pocket_quality = mean_w_intf * clamp(_eff_q_met / 3.0, 0.0, 1.0) if _eff_q_met > 0.0 else 0.0
        raw_uis = mode_b_factor * (0.55 * pocket_quality + 0.45 * sigma_aff)
    elif mode == "MODE_A_PRED_FULL":
        raw_uis = 0.25 * float(iptm_val) + 0.30 * lis_score + 0.25 * pdockq_norm + 0.20 * sigma_aff
    else:
        raw_uis = mode_b_factor * (0.55 * pdockq_norm + 0.45 * sigma_aff)

    thermo_gate = clamp((-dg_uni + 0.5) / 6.0, 0.10, 1.0)
    uis = clamp(raw_uis * thermo_gate, 0.0, 1.0)

    uis_rounded = round(uis, 4)
    uis_margin = round(uis_rounded - UIS_DECISION_THRESHOLD, 4)
    threshold_call = "ABOVE_THRESHOLD" if uis_rounded >= UIS_DECISION_THRESHOLD else "BELOW_THRESHOLD"
    binding_class = "HIGH_CONFIDENCE" if uis_rounded >= UIS_HIGH_CONF_THRESHOLD else ("MODERATE" if uis_rounded >= 0.42 else "LOW_OR_DECOY")

    pattern_id, pattern_desc = classify_mechanistic_pattern(
        uis_rounded, float(iptm_val), lis_score, dg_uni, w_bsa, eff_hel_pro, n_clash, proline_penalty=proline_penalty
    )
    if mode == "MODE_C_EXP_CRYSTAL" and pattern_id in ("Pattern A", "Pattern B"):
        pattern_id = "Pattern A"
        pattern_desc = f"Experimental Structure Supported {modality} Interactor (Interface packing and dG_uni support binding)"
    elif mode == "MODE_B_PRED_COORD" and pattern_id in ("Pattern A", "Pattern B"):
        pattern_id = "Pattern A"
        pattern_desc = f"Unverified Coordinate-Only Prediction ({modality}; supply PAE JSON for full Mode A verification)"

    elec_density = round(1000.0 * (n_sb + n_hb) / max(100.0, bsa), 2)
    rna_tot_spec = n_sb + n_hb_base + n_hb_2oh + n_hb_bb + n_pi + n_catpi
    rna_ssi = round((n_hb_base + n_pi + n_catpi) / max(1, rna_tot_spec), 3) if is_nucleic_pair else 0.0

    def _loci_span(keys: List[str]) -> str:
        nums = [int(m.group(1)) for k in keys for m in [re.search(r"(\d+)$", k)] if m]
        return f"{min(nums)}-{max(nums)}" if nums else "None"

    return {
        "chain_pair": f"{chain_a_id}-{chain_b_id}",
        "modality": modality,
        "len_a": len(residues_a),
        "len_b": len(residues_b),
        "conf_contacts": len(core_conf_contacts),
        "raw_contacts": len([c for c in residue_contact_log if c['min_dist_A'] <= active_cutoff]),
        "loci_a": _loci_span([c["res_a"] for c in core_conf_contacts]),
        "loci_b": _loci_span([c["res_b"] for c in core_conf_contacts]),
        "sha256_id": pair_sha256,
        "mode": mode,
        "UIS": uis_rounded,
        "UIS_threshold": UIS_DECISION_THRESHOLD,
        "UIS_margin": uis_margin,
        "threshold_call": threshold_call,
        "binding_class": binding_class,
        "pattern_id": pattern_id,
        "pattern_desc": pattern_desc,
        "dG_uni_kcal": round(dg_uni, 2),
        "dG_prodigy_kcal": round(dg_prodigy_raw, 2),
        "dG_phys_kcal": round(dg_phys, 2),
        "Kd_M": float(f"{kd_molar:.3e}"),
        "Kd_nM": round(kd_molar * 1e9, 3),
        "ipTM": round(float(iptm_val), 3),
        "LIS": round(lis_score, 3),
        "pDockQ_PAE": round(pdockq_norm, 3),
        "pDockQ_coord": round(pdockq_coord, 3),
        "mean_intf_pLDDT": round(mean_intf_plddt, 1),
        "mean_intf_PAE_A": round(mean_intf_pae, 2),
        "BSA_A2": round(bsa, 1),
        "wBSA_A2": round(w_bsa, 1),
        "intf_residues": len({c["res_a"] for c in core_conf_contacts}) + len({c["res_b"] for c in core_conf_contacts}),
        "pct_disordered": round(100.0 * sum(1 for r in (residues_a + residues_b) if residue_raw_plddt(r) < 70.0) / max(1, len(residues_a) + len(residues_b)), 1),
        "nis_total": nis_total,
        "n_cb8_contacts": n_cb8_eff,
        "salt_bridges": n_sb,
        "h_bonds": n_hb,
        "hb_base": n_hb_base,
        "hb_2oh": n_hb_2oh,
        "hb_backbone": n_hb_bb,
        "cation_pi": n_catpi,
        "halogen_bonds": n_xb,
        "metal_coord_bonds": n_metal_coord,
        "metal_outer_bonds": n_metal_outer,
        "covalent_bonds": n_cov,
        "lig_burial_pct": round(lig_burial_frac * 100.0, 1),
        "pocket_quality": round(pocket_quality, 3),
        "rna_specificity_idx": rna_ssi,
        "hydrophobic_contacts": n_hp,
        "pi_stacking": n_pi,
        "disulfides": n_ss,
        "steric_clashes": n_clash,
        "helical_prolines": eff_hel_pro,
        "proline_warning": eff_hel_pro > 0,
        "f_apolar": round(apolar_frac, 3),
        "elec_density_1000A2": elec_density,
        "mean_w_intf": round(mean_w_intf, 3),
        "top_hotspot": hotspots[0]["residue"] if hotspots else "NONE",
        "top_hotspot_score": hotspots[0]["priority_score"] if hotspots else 0.0,
        "hotspots": hotspots[:15],
        "custom_mutations": custom_mut_results,
        "foldx_mutations": foldx_mutations,
        "contacts": [c for c in residue_contact_log if c['min_dist_A'] <= active_cutoff]
    }

def export_sanitized_pdb(residues: List[Residue], out_pdb_path: str, canonicalize_ptms: bool = True) -> Dict[str, str]:
    chain_remap: Dict[str, str] = {}
    ptm_logs: List[str] = []

    for r in residues:
        if r.chain_id not in chain_remap:
            idx = len(chain_remap) % len(CHAIN_CHAR_POOL)
            chain_remap[r.chain_id] = CHAIN_CHAR_POOL[idx]
        if canonicalize_ptms and r.res_name in MODIFIED_AA_MAP:
            ptm_logs.append(f"REMARK 990 PTM CANONICALIZED: {r.chain_id}:{r.res_name}{r.res_seq} -> {r.canonical_name}")

    serial = 1
    prev_chain = None
    prev_r_name = "ALA"
    prev_r_seq = 1

    with open(out_pdb_path, "w", encoding="utf-8") as f:
        f.write("HEADER    UNIPPI-OWL SANITIZED STRUCTURE FOR FOLDX/ROSETTA\n")
        for log in ptm_logs:
            f.write(log + "\n")

        for r in residues:
            record_type = "HETATM" if r.entity_type in ("LIGAND", "ION") else "ATOM  "
            c_char = chain_remap[r.chain_id]
            if prev_chain is not None and c_char != prev_chain:
                f.write(f"TER   {serial:5d}      {prev_r_name:3s} {prev_chain:1s}{prev_r_seq:4d}\n")
                serial += 1
            prev_chain = c_char

            r_name = (r.canonical_name if canonicalize_ptms else r.res_name)[:3]
            r_seq = max(-999, min(9999, r.res_seq))
            prev_r_name = r_name
            prev_r_seq = r_seq
            ins = r.ins_code[:1] if r.ins_code else " "

            for a in r.atoms:
                if r_name == "GLY" and a.name == "CB":
                    continue
                a_name = a.name[:4]
                fmt_name = f" {a_name:<3s}" if len(a_name) < 4 else a_name
                el = a.element[:2].rjust(2)
                line = (
                    f"{record_type:6s}{serial:5d} {fmt_name:4s} {r_name:3s} {c_char:1s}"
                    f"{r_seq:4d}{ins:1s}   "
                    f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occupancy:6.2f}{a.b_factor:6.2f}"
                    f"          {el:2s}  \n"
                )
                f.write(line)
                serial = (serial + 1) % 100000
        f.write("END\n")
    return chain_remap

def export_ligands_sdf(residues: List[Residue], out_sdf_path: str) -> int:
    ligands = [r for r in residues if r.entity_type == "LIGAND" and len(r.atoms) > 1]
    if not ligands:
        return 0

    with open(out_sdf_path, "w", encoding="utf-8") as f:
        for lig in ligands:
            n_atoms = len(lig.atoms)
            bonds: List[Tuple[int, int]] = []
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    d = math.dist(lig.atoms[i].coord, lig.atoms[j].coord)
                    if 0.8 <= d <= 1.95:
                        bonds.append((i + 1, j + 1))

            f.write(f"{lig.chain_id}_{lig.res_name}_{lig.res_seq}\n")
            f.write("  UniPPI-OWL 3D Ligand Extractor\n\n")
            f.write(f"{n_atoms:3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000\n")
            for a in lig.atoms:
                f.write(f"{a.x:10.4f}{a.y:10.4f}{a.z:10.4f} {a.element:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n")
            for b1, b2 in bonds:
                f.write(f"{b1:3d}{b2:3d}  1  0  0  0  0\n")
            f.write("M  END\n$$$$\n")
    return len(ligands)

def generate_benchmark_suite(demo_dir: str) -> None:
    shutil.rmtree(demo_dir, ignore_errors=True)
    os.makedirs(demo_dir, exist_ok=True)
    specs = [
        ("1_VIP1_Homodimer_Positive", "WT_VIP1", 94.0, 2.3, 0.90, 2.95),
        ("2_VIP1_bZIP61_GroupE_Positive", "HET_GROUP_E", 89.0, 3.1, 0.83, 3.10),
        ("3_Barnase_Barstar_Electrostatic_Pos", "ELECTRO_GLOBULAR", 94.0, 2.2, 0.89, 2.90),
        ("4_VIP1_GRF6_14-3-3_IDR_Diluted_Pos", "IDR_DILUTED_PAT_B", 88.5, 2.6, 0.52, 2.95),
        ("5_bZIP34_Proline_Homodimer_Negative", "HOMO_PRO_KINK", 61.0, 14.8, 0.38, 3.45),
        ("6_VIP1_IDR_Hallucination_Decoy", "IDR_DECOY", 33.0, 27.5, 0.24, 3.20),
    ]

    heptad_template = [
        "GLU", "LEU", "ALA", "ALA", "LEU", "GLU", "LYS",
        "GLU", "VAL", "ALA", "ALA", "LEU", "GLU", "LYS",
        "GLU", "ILE", "ALA", "ALA", "LEU", "GLU", "LYS",
        "GLU", "LEU", "ALA", "ALA", "LEU", "GLU", "LYS"
    ]

    for name, mode_tag, plddt_val, pae_val, iptm_val, cb_half_span in specs:
        c_dir = os.path.join(demo_dir, name)
        os.makedirs(c_dir, exist_ok=True)
        cif_path = os.path.join(c_dir, f"{name}_model_0.cif")
        json_path = os.path.join(c_dir, f"{name}_confidences_0.json")

        seq_a = list(heptad_template)
        seq_b = list(heptad_template)
        for r_i in range(0, 28, 7):
            seq_a[r_i] = "GLU"
            seq_b[r_i] = "LYS"

        if mode_tag == "HET_GROUP_E":
            seq_b[8] = "ILE"
            seq_b[15] = "VAL"
        elif mode_tag == "ELECTRO_GLOBULAR":
            for idx_e in (1, 8, 15, 22):
                seq_a[idx_e] = "ASP"
                seq_b[idx_e] = "ARG"
            for idx_h in (4, 11, 18, 25):
                seq_a[idx_h] = "TYR"
                seq_b[idx_h] = "ASN"
        elif mode_tag == "IDR_DILUTED_PAT_B":
            for idx_h in (4, 11):
                seq_a[idx_h] = "SER"
                seq_b[idx_h] = "GLN"
        elif mode_tag == "HOMO_PRO_KINK":
            for idx_p in (4, 11, 18):
                seq_a[idx_p] = "PRO"
                seq_b[idx_p] = "PRO"
        elif mode_tag == "IDR_DECOY":
            seq_a = ["SER", "GLY", "ASN", "GLN", "SER", "THR", "GLY"] * 4
            seq_b = ["SER", "GLY", "ASN", "GLN", "SER", "THR", "GLY"] * 4

        n_res_chain = len(seq_a)
        atom_lines: List[str] = []
        serial = 1

        for chain_idx, (chain_id, seq) in enumerate((("A", seq_a), ("B", seq_b))):
            sign = -1.0 if chain_idx == 0 else 1.0
            for r_i, r_name in enumerate(seq):
                r_seq = r_i + 1
                z_base = r_i * 1.50
                y_base = math.sin(r_i * 1.745) * 1.10
                kink_shift = 0.95 if (mode_tag == "HOMO_PRO_KINK" and r_name == "PRO") else 0.0
                local_span = cb_half_span + kink_shift
                x_ca = sign * (local_span + 1.50)
                x_cb = sign * local_span

                atoms_spec = [
                    ("N", x_ca + sign * 0.8, y_base - 0.5, z_base - 0.6, "N"),
                    ("CA", x_ca, y_base, z_base, "C"),
                    ("C", x_ca + sign * 0.6, y_base + 0.6, z_base + 0.6, "C"),
                    ("O", x_ca + sign * 1.2, y_base + 1.1, z_base + 0.8, "O"),
                ]
                if r_name != "GLY":
                    atoms_spec.append(("CB", x_cb, y_base, z_base, "C"))

                if r_name in ("LEU", "ILE", "VAL"):
                    atoms_spec.append(("CG", sign * (local_span - 0.70), y_base + 0.15, z_base, "C"))
                    atoms_spec.append(("CD1", sign * (local_span - 1.20), y_base - 0.10, z_base + 0.15, "C"))
                    atoms_spec.append(("CD2", sign * (local_span - 1.15), y_base + 0.25, z_base - 0.15, "C"))
                elif r_name == "PRO":
                    atoms_spec.append(("CG", x_cb + sign * 0.4, y_base + 0.8, z_base - 0.3, "C"))
                    atoms_spec.append(("CD", x_ca + sign * 0.9, y_base + 0.3, z_base - 0.7, "C"))
                elif r_name == "GLU" and (r_i % 7 == 0):
                    atoms_spec.append(("CG", sign * (local_span - 0.50), y_base + 0.60, z_base, "C"))
                    atoms_spec.append(("CD", sign * (local_span - 1.00), y_base + 1.00, z_base, "C"))
                    atoms_spec.append(("OE1", sign * (local_span - 1.40), y_base + 1.20, z_base, "O"))
                elif r_name == "LYS" and (r_i % 7 == 0):
                    atoms_spec.append(("CG", sign * (local_span - 0.50), y_base + 0.60, z_base, "C"))
                    atoms_spec.append(("CE", sign * (local_span - 1.00), y_base + 1.00, z_base, "C"))
                    atoms_spec.append(("NZ", sign * (local_span - 1.40), y_base + 1.20, z_base, "N"))
                elif r_name == "ASP":
                    atoms_spec.append(("CG", sign * (local_span - 0.60), y_base + 0.30, z_base, "C"))
                    atoms_spec.append(("OD1", sign * (local_span - 1.45), y_base + 0.10, z_base, "O"))
                elif r_name == "ARG":
                    atoms_spec.append(("CZ", sign * (local_span - 0.70), y_base + 0.30, z_base, "C"))
                    atoms_spec.append(("NH1", sign * (local_span - 1.45), y_base + 0.10, z_base, "N"))
                elif r_name == "TYR":
                    atoms_spec.append(("CE1", sign * (local_span - 0.80), y_base, z_base, "C"))
                    atoms_spec.append(("OH", sign * (local_span - 1.45), y_base, z_base, "O"))
                elif r_name in ("SER", "ASN", "GLN", "THR"):
                    og_off = 1.20 if mode_tag == "IDR_DECOY" else 1.45
                    atoms_spec.append(("OG", sign * (local_span - og_off), y_base, z_base, "O"))

                for a_name, ax, ay, az, el in atoms_spec:
                    atom_lines.append(
                        f"ATOM {serial} {el} {a_name} . {r_name} {chain_id} 1 {r_seq} ? "
                        f"{ax:.3f} {ay:.3f} {az:.3f} 1.00 {plddt_val:.2f} {r_seq} {chain_id} 1\n"
                    )
                    serial += 1

        with open(cif_path, "w", encoding="utf-8") as f:
            f.write(f"data_{name}\n#\n_ma_qa_metric_local.type pLDDT\nloop_\n")
            f.write("_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n")
            f.write("_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n")
            f.write("_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n")
            f.write("_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n")
            f.write("_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n")
            f.write("_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n_atom_site.pdbx_PDB_model_num\n")
            f.writelines(atom_lines)
            f.write("#\n")

        total_res = 2 * n_res_chain
        pae_mat: List[List[float]] = []
        for row_i in range(total_res):
            row_vals: List[float] = []
            for col_j in range(total_res):
                if row_i == col_j:
                    row_vals.append(0.8)
                elif (row_i < n_res_chain and col_j < n_res_chain) or (row_i >= n_res_chain and col_j >= n_res_chain):
                    row_vals.append(2.0)
                else:
                    row_vals.append(pae_val)
            pae_mat.append(row_vals)

        token_chains = ["A"] * n_res_chain + ["B"] * n_res_chain
        token_res = list(range(1, n_res_chain + 1)) + list(range(1, n_res_chain + 1))

        with open(json_path, "w", encoding="utf-8") as f:
            json.dump({
                "iptm": iptm_val,
                "ptm": round(iptm_val + 0.04, 2),
                "ranking_score": round(0.8 * iptm_val + 0.2 * (iptm_val + 0.04), 3),
                "token_chain_ids": token_chains,
                "token_res_ids": token_res,
                "pae": pae_mat
            }, f)

def process_structure_file(
    filepath: str,
    out_dir: str,
    temp_k: float = 298.15,
    lambda_pae: float = 8.0,
    ionic_mM: float = 150.0,
    ph_val: float = 7.2,
    cutoff_a: float = 5.5,
    proline_penalty: float = 0.0,
    custom_mutations: Optional[List[str]] = None,
    split_chain: Optional[str] = None
) -> List[Dict]:
    ext = os.path.splitext(filepath)[1].lower()
    if ext in (".cif", ".mmcif"):
        atoms, meta = parse_mmcif_file(filepath)
    else:
        atoms, meta = parse_pdb_file(filepath)

    if not atoms:
        return []

    ordered_residues, chains = build_residues(atoms)
    if split_chain:
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
                        _set_residue_chain(r, c1)
                    for r in d2:
                        _set_residue_chain(r, c2)
                    del chains[sc_cid]
                    chains[c1] = d1
                    chains[c2] = d2
                    kept_ids = {id(r) for cl in chains.values() for r in cl}
                    ordered_residues = [r for r in ordered_residues if id(r) in kept_ids]
        except Exception as _e:
            print(f"[WARN] Could not parse --split-chain '{split_chain}': {_e}")

    if len(chains) < 2:
        return []

    conf_data = parse_companion_json(filepath, ordered_residues)
    mode = detect_confidence_mode(ordered_residues, meta, conf_data)

    sasa_complex = compute_residue_sasa(ordered_residues)
    sasa_free_by_chain: Dict[str, Dict[Tuple[str, int, str], float]] = {
        cid: compute_residue_sasa(c_res) for cid, c_res in chains.items()
    }
    res_to_global_idx = {r.key: idx for idx, r in enumerate(ordered_residues)}

    stem = os.path.splitext(os.path.basename(filepath))[0]
    out_dir = out_dir or "./unippi_results"
    struct_out_dir = os.path.join(out_dir, "converted_structures")
    os.makedirs(struct_out_dir, exist_ok=True)

    chain_remap = export_sanitized_pdb(ordered_residues, os.path.join(struct_out_dir, f"{stem}_sanitized.pdb"))
    export_ligands_sdf(ordered_residues, os.path.join(struct_out_dir, f"{stem}_ligands.sdf"))

    chain_ids = list(chains.keys())
    is_binary_complex = (len(chain_ids) == 2)
    results: List[Dict] = []

    for i in range(len(chain_ids)):
        for j in range(i + 1, len(chain_ids)):
            ca, cb = chain_ids[i], chain_ids[j]
            if not is_binary_complex:
                has_prox = any(
                    math.dist(ra.centroid, rb.centroid) <= 16.0
                    for ra in chains[ca] for rb in chains[cb]
                )
                if not has_prox:
                    continue
                pair_sasa_complex = compute_residue_sasa(chains[ca] + chains[cb])
            else:
                pair_sasa_complex = sasa_complex

            res_pair = analyze_interface_pair(
                ca, cb, chains[ca], chains[cb], res_to_global_idx,
                mode, conf_data, pair_sasa_complex,
                sasa_free_by_chain[ca], sasa_free_by_chain[cb],
                chain_remap=chain_remap,
                temp_k=temp_k, lambda_pae=lambda_pae, ionic_mM=ionic_mM, ph_val=ph_val,
                cutoff_a=cutoff_a, proline_penalty=proline_penalty,
                force_report_if_empty=is_binary_complex,
                custom_mutations=custom_mutations
            )
            if res_pair:
                res_pair["complex_id"] = stem
                res_pair["source_file"] = os.path.basename(filepath)
                results.append(res_pair)

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
                        partner_hits.append((other, rrec.get("lig_burial_pct", 0.0), 0, 0))
                else:
                    _mc_in = rrec.get("metal_coord_bonds", 0)
                    _mc_out = rrec.get("metal_outer_bonds", 0)
                    if o_type in ("PROTEIN", "NUCLEIC", "LIGAND") and (_mc_in >= 1 or _mc_out >= 2):
                        partner_hits.append((other, float(_mc_in) + 0.50 * float(_mc_out), _mc_in, _mc_out))
            if (is_lig_c and len(partner_hits) >= 2) or (not is_lig_c and len(partner_hits) >= 2 and any(h[2] >= 1 for h in partner_hits)):
                partner_hits.sort(key=lambda x: x[1], reverse=True)
                m1, v1, _in1, _out1 = partner_hits[0]
                m2, v2, _in2, _out2 = partner_hits[1]
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
                    cutoff_a=cutoff_a, proline_penalty=proline_penalty,
                    force_report_if_empty=False,
                    custom_mutations=None
                )
                if t_pair:
                    t_pair["complex_id"] = stem
                    t_pair["source_file"] = os.path.basename(filepath)
                    t_pair["ternary_partners"] = comp_id_m
                    t_pair["ternary_burial_split"] = f"{m1}:{v1:.1f}% / {m2}:{v2:.1f}%" if is_lig_c else f"{m1}:{_in1}in+{_out1}out / {m2}:{_in2}in+{_out2}out"
                    t_pair["cooperativity_balance"] = round(min(v1, v2) / max(1e-6, max(v1, v2)), 2)
                    if is_lig_c and t_pair.get("lig_burial_pct", 0.0) >= 50.0 and t_pair["UIS"] >= UIS_DECISION_THRESHOLD:
                        t_pair["pattern_id"] = "Pattern A"
                        t_pair["pattern_desc"] = f"Ternary Molecular Glue / Interface Wedge ({lc} bridges {m1} [{v1:.1f}%] & {m2} [{v2:.1f}%]; composite burial {t_pair['lig_burial_pct']:.1f}%, balance {t_pair['cooperativity_balance']:.2f})"
                    elif not is_lig_c and t_pair.get("metal_coord_bonds", 0) >= 1:
                        if t_pair["UIS"] >= UIS_DECISION_THRESHOLD:
                            t_pair["pattern_id"] = "Pattern A"
                        t_pair["pattern_desc"] = f"Catalytic Metal Bridge ({lc} bridges {m1} [{_in1} inner, {_out1} outer-H2O] & {m2} [{_in2} inner, {_out2} outer-H2O])"
                    results.append(t_pair)
    return results

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

def _hedges_g_py(a: List[float], b: List[float]) -> Optional[float]:
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return None
    m1, m2 = sum(a) / n1, sum(b) / n2
    v1 = sum((x - m1) ** 2 for x in a) / (n1 - 1)
    v2 = sum((x - m2) ** 2 for x in b) / (n2 - 1)
    df = n1 + n2 - 2
    sp = math.sqrt(max(0.0, ((n1 - 1) * v1 + (n2 - 1) * v2) / df))
    if sp < 1e-12:
        return 0.0 if abs(m2 - m1) < 1e-9 else (math.inf if m2 > m1 else -math.inf)
    j_corr = 1.0 - (3.0 / (4.0 * df - 1.0))
    return j_corr * (m2 - m1) / sp

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
    g_dg = _hedges_g_py(dg_a, dg_b)
    g_uis = _hedges_g_py(uis_a, uis_b)
    print(f"\n{'='*96}")
    print(f"[ENSEMBLE DIFF: Condition A (N={len(sub_a)}) vs Condition B (N={len(sub_b)})]")
    _fmt_g = lambda g: "N/A" if g is None else (f"{g:+.2f}" if math.isfinite(g) else ("+inf (0 SD)" if g > 0 else "-inf (0 SD)"))
    print(f"  dUIS = {duis:+.3f} | Hedges' g = {_fmt_g(g_uis)} | Welch p = {f'{p_uis:.4g}' if p_uis is not None else 'N/A'}")
    print(f"  ddG  = {ddg:+.2f} kcal/mol | Kd Fold-Change = {kd_fold:.2g}x | Hedges' g = {_fmt_g(g_dg)} | Welch p = {f'{p_dg:.4g}' if p_dg is not None else 'N/A'}")
    print("  [Statistical Note] Multi-seed ensembles are recommended for formal inferential p-values; single-seed diffusion samples measure sampling consistency.")

    def _norm_res(tok: str) -> str:
        parts = tok.split(":", 1)
        m_num = re.search(r"(\d+)$", parts[-1])
        return f"{parts[0]}:{m_num.group(1)}" if (len(parts) == 2 and m_num) else tok

    def _agg(sub):
        m: Dict[str, Dict] = {}
        for r in sub:
            seen = set()
            for c in r.get("contacts", []):
                if c.get("weight_w", 0.0) < 0.25:
                    continue
                k = f"{_norm_res(c['res_a'])} <-> {_norm_res(c['res_b'])}"
                if k in seen:
                    continue
                seen.add(k)
                rec = m.setdefault(k, {"label": f"{c['res_a']} <-> {c['res_b']}", "dists": []})
                rec["dists"].append(float(c.get("min_dist_A", 0.0)))
        return m
    ma, mb = _agg(sub_a), _agg(sub_b)
    lost, gained = [], []
    diff_rows = ["status,contact_pair,occ_A_pct,occ_B_pct,mean_dist_A_A,mean_dist_B_A"]
    for k in sorted(set(ma.keys()) | set(mb.keys())):
        ra_m, rb_m = ma.get(k), mb.get(k)
        oa = (len(ra_m["dists"]) / len(sub_a)) if ra_m else 0.0
        ob = (len(rb_m["dists"]) / len(sub_b)) if rb_m else 0.0
        lbl = f"{ra_m['label']} [B: {rb_m['label']}]" if (ra_m and rb_m and ra_m["label"] != rb_m["label"]) else (ra_m["label"] if ra_m else rb_m["label"])
        da = f"{sum(ra_m['dists'])/len(ra_m['dists']):.2f}" if ra_m else ""
        db = f"{sum(rb_m['dists'])/len(rb_m['dists']):.2f}" if rb_m else ""
        if oa >= 0.6 and ob < 0.4:
            lost.append(lbl)
            diff_rows.append(f'LOST,"{lbl}",{round(oa*100)},{round(ob*100)},{da},{db}')
        elif oa < 0.4 and ob >= 0.6:
            gained.append(lbl)
            diff_rows.append(f'GAINED,"{lbl}",{round(oa*100)},{round(ob*100)},{da},{db}')
        elif oa >= 0.4 or ob >= 0.4:
            diff_rows.append(f'PERSISTENT,"{lbl}",{round(oa*100)},{round(ob*100)},{da},{db}')
    print(f"  Contact Rewiring: {len(lost)} LOST | {len(gained)} GAINED")
    if lost:
        print("    Top LOST  :", ", ".join(lost[:6]))
    if gained:
        print("    Top GAINED:", ", ".join(gained[:6]))
    diff_csv = os.path.join(out_dir, "diff_contact_rewiring.csv")
    with open(diff_csv, "w", encoding="utf-8") as _df:
        _df.write("\n".join(diff_rows) + "\n")
    print(f"  Saved Diff CSV: {diff_csv}")
    print(f"{'='*96}")

def parse_structure_file(filepath: str) -> Tuple[List[Atom], Dict[str, str]]:
    ext = os.path.splitext(str(filepath))[1].lower()
    if ext in (".cif", ".mmcif"):
        return parse_mmcif_file(str(filepath))
    return parse_pdb_file(str(filepath))

def _collect_struct_files(target: str) -> Tuple[List[str], Optional[str]]:
    import tempfile, zipfile
    tmp_dir = None
    if os.path.isfile(target) and target.lower().endswith(".zip"):
        tmp_dir = tempfile.mkdtemp(prefix="unippi_zip_")
        with zipfile.ZipFile(target, "r") as zf:
            zf.extractall(tmp_dir)
        search_root = tmp_dir
        files = sorted(
            glob.glob(os.path.join(search_root, "**", "*.cif"), recursive=True)
            + glob.glob(os.path.join(search_root, "**", "*.mmcif"), recursive=True)
            + glob.glob(os.path.join(search_root, "**", "*.pdb"), recursive=True)
        )
    elif os.path.isfile(target):
        files = [target]
    else:
        files = sorted(
            glob.glob(os.path.join(target, "**", "*.cif"), recursive=True)
            + glob.glob(os.path.join(target, "**", "*.mmcif"), recursive=True)
            + glob.glob(os.path.join(target, "**", "*.pdb"), recursive=True)
        )
    files = [
        fp for fp in files
        if "_sanitized.pdb" not in fp
        and not any(ig in fp.replace("\\", "/").lower() for ig in ("/msas/", "/templates/", "template_hit"))
    ]
    return files, tmp_dir

def main():
    parser = argparse.ArgumentParser(
        description="UniPPI-OWL CLI v2.2: Universal Multi-Modal Structure-to-Affinity & Interface Auditing Engine"
    )
    parser.add_argument("-i", "--input", type=str, help="Input directory, .zip archive, or single .cif/.pdb file")
    parser.add_argument("-o", "--output", type=str, default="./unippi_results", help="Output directory")
    parser.add_argument("--cutoff", type=float, default=5.5, help="Heavy-atom contact cutoff in A (default: 5.5 A)")
    parser.add_argument("--temp", type=float, default=298.15, help="Temperature in Kelvin (default: 298.15 K)")
    parser.add_argument("--lambda-pae", type=float, default=8.0, help="PAE exponential decay constant (default: 8.0 A)")
    parser.add_argument("--ionic-strength", type=float, default=150.0, help="Ionic strength in mM NaCl (default: 150.0 mM)")
    parser.add_argument("--ph", type=float, default=7.2, help="Solution pH (default: 7.2)")
    parser.add_argument("--proline-penalty", type=float, default=0.0, help="Optional kcal/mol penalty per strained helical Proline (default: 0.0, diagnostic warning only)")
    parser.add_argument("--mutate", type=str, default="", help="Comma-separated custom mutations (e.g. A:LEU5PRO,A:GLU1LYS)")
    parser.add_argument("--split-chain", type=str, default=None, help="Split one chain into 2 domains (e.g. A:1-236,237-406 or A:236)")
    parser.add_argument("--diff-b", type=str, default=None, help="Second structure/.zip/folder (Condition B) for Welch t-test, Hedges' g & Contact Rewiring diff")
    parser.add_argument("--generate-demo", type=str, metavar="DEMO_DIR", help="Generate 6-complex benchmark demo suite")
    args = parser.parse_args()

    if args.generate_demo:
        generate_benchmark_suite(args.generate_demo)
        print(f"[UniPPI-CLI] Generated 6-complex demo suite in: {args.generate_demo}")
        if not args.input:
            args.input = args.generate_demo

    if not args.input:
        parser.print_help()
        sys.exit(0)

    custom_mut_list = [m.strip() for m in args.mutate.split(",") if m.strip()] if args.mutate else None
    os.makedirs(args.output, exist_ok=True)

    struct_files, tmp_a = _collect_struct_files(args.input)
    all_records: List[Dict] = []
    try:
        for fpath in struct_files:
            recs = process_structure_file(
                fpath, args.output,
                temp_k=args.temp, lambda_pae=args.lambda_pae,
                ionic_mM=args.ionic_strength, ph_val=args.ph,
                cutoff_a=args.cutoff, proline_penalty=args.proline_penalty,
                custom_mutations=custom_mut_list, split_chain=args.split_chain
            )
            all_records.extend(recs)
    finally:
        if tmp_a:
            shutil.rmtree(tmp_a, ignore_errors=True)

    all_records.sort(key=lambda x: x["UIS"], reverse=True)

    csv_cols = [
        "complex_id", "chain_pair", "modality", "loci_a", "loci_b", "conf_contacts", "sha256_id", "mode",
        "UIS", "UIS_threshold", "UIS_margin", "threshold_call", "binding_class", "pattern_id",
        "dG_uni_kcal", "dG_prodigy_kcal", "dG_phys_kcal", "Kd_M", "Kd_nM", "ipTM", "LIS",
        "pDockQ_PAE", "pDockQ_coord", "mean_intf_pLDDT", "mean_intf_PAE_A", "BSA_A2", "wBSA_A2",
        "intf_residues", "nis_total", "n_cb8_contacts", "salt_bridges", "h_bonds",
        "hb_base", "hb_2oh", "hb_backbone", "hydrophobic_contacts", "pi_stacking", "cation_pi",
        "halogen_bonds", "metal_coord_bonds", "metal_outer_bonds", "covalent_bonds",
        "disulfides", "steric_clashes", "helical_prolines", "proline_warning", "f_apolar",
        "elec_density_1000A2", "mean_w_intf", "lig_burial_pct", "pocket_quality",
        "rna_specificity_idx", "ternary_partners", "cooperativity_balance",
        "top_hotspot", "top_hotspot_score"
    ]
    csv_path = os.path.join(args.output, "master_interaction_ranking.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(",".join(csv_cols) + "\n")
        for r in all_records:
            f.write(",".join(str(r.get(c, "")) for c in csv_cols) + "\n")

    json_out_path = os.path.join(args.output, "unippi_full_report.json")
    with open(json_out_path, "w", encoding="utf-8") as f:
        json.dump(all_records, f, indent=2)

    print(f"\nModel: {MODEL_VERSION} | Frozen Threshold (UIS): {UIS_DECISION_THRESHOLD:.4f} | Temp: {args.temp:.2f} K | pH: {args.ph:.1f} | NaCl: {args.ionic_strength:.0f} mM")
    print("=" * 128)
    print(f"{'Complex ID':<34} | {'Pair':<16} | {'Modality':<15} | {'UIS':<6} | {'dG (kcal)':<9} | {'Kd (~est)':<11} | {'LIS':<5} | {'Pattern':<9} | {'Class'}")
    print("-" * 128)
    for r in all_records:
        kd_nm = r["Kd_nM"]
        kd_str = ">10mM (NB)" if kd_nm >= 1e7 else (f"{kd_nm/1e6:.2g} mM" if kd_nm >= 1e6 else (f"{kd_nm/1e3:.2g} uM" if kd_nm >= 1e3 else (f"{kd_nm*1e3:.2g} pM" if kd_nm < 0.5 else f"{kd_nm:.2g} nM")))
        print(
            f"{r['complex_id'][:34]:<34} | {r['chain_pair'][:16]:<16} | {r.get('modality','PPI')[:15]:<15} | {r['UIS']:<6.3f} | "
            f"{r['dG_uni_kcal']:<9.2f} | {kd_str:<11} | {r['LIS']:<5.3f} | "
            f"{r['pattern_id']:<9} | {r['binding_class']}"
        )
    print("=" * 128)
    print(f"[UniPPI-CLI] Saved Master CSV: {csv_path}")

    if args.diff_b:
        files_b, tmp_b = _collect_struct_files(args.diff_b)
        recs_b: List[Dict] = []
        try:
            for fpb in files_b:
                recs_b.extend(process_structure_file(
                    fpb, args.output, temp_k=args.temp, lambda_pae=args.lambda_pae,
                    ionic_mM=args.ionic_strength, ph_val=args.ph,
                    cutoff_a=args.cutoff, proline_penalty=args.proline_penalty,
                    split_chain=args.split_chain
                ))
        finally:
            if tmp_b:
                shutil.rmtree(tmp_b, ignore_errors=True)
        run_cli_ensemble_diff(all_records, recs_b, args.output, args.temp)

if __name__ == "__main__":
    main()
