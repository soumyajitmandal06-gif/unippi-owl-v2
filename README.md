# UniPPI-OWL v2.1 (MARVAS) — Universal Multi-Modal Structure-to-Affinity & Interface Auditing Engine

**Live Zero-Server Web App:** [https://soumyajitmandal06-gif.github.io/unippi-owl-v2/](https://soumyajitmandal06-gif.github.io/unippi-owl-v2/)

UniPPI-OWL v2.1-MARVAS is a dual-layer (**100% client-side Browser SPA + pure-Python CLI**) post-prediction biophysical auditing engine for **AlphaFold 3**, **Boltz-1/2**, **Chai-1**, and experimental **PDB/mmCIF** complexes. All calculations run locally with zero external server uploads.

## Supported Biomolecular Modalities (10 Regimes)
- **Protein-Protein (`PPI`)**: Enthalpy-entropy compensated 8-channel solvation + contact-gated PRODIGY hybrid, PAE-weighted `w_ij` IDR damping, Local Interaction Score (`LIS`), and strained helical Proline kink detection.
- **Protein-RNA / DNA & Nucleic Duplexes**: Dedicated 9-channel physics separating phosphate salt bridges (`pKa ~ 1.0`), base-specific vs. ribose `2'-OH` vs. backbone H-bonds, plane-normal nucleobase pi-stacking, `60 deg`-cone cation-pi, and Sequence-Specificity Index (`SSI`).
- **Small Molecules, Drugs & Cofactors**: 3D aromatic ring detector, anionic phosphate/carboxylate & cationic aliphatic amine/guanidinium salt-bridge detectors, directional `>=120 deg` sigma-hole halogen bonds, fractional ligand burial `Q_pocket`, and `LIGAND-LIGAND` support.
- **Catalytic & Structural Metals (`Zn2+`, `Mg2+`, `Ca2+`, `Mn2+`, `Fe`, `Cu`, `Co`, `Ni`)**: Inner-sphere (`<=2.65-2.75 A`) and hydrated outer-sphere (`2.75-4.15 A`) oxygen coordination.
- **Automatic Ternary Molecular-Glue & Metal-Bridge Detection**: Detects composite `A+B-Ligand` wedges and `A+B-Metal` catalytic bridges with burial splits and cooperativity balance.
- **Same-Protein Intra-Chain Domain Slicing**: Custom `--split-chain A:1-236,237-406` or automatic contact-trough monomer splitting.
- **Universal X vs. Y (+ Unlimited Z1..Zn) Studio**: Mix known PDB/UniProt entries with custom lab-synthesized sequences, PDB-numbered mutations, and custom drug SMILES; exports ready-to-run **Boltz-1/2 YAML + 1-click Colab `.sh`**, **AF3 Server JSON**, and **Chai-1 FASTA**.
- **WT vs. Mutant / Apo vs. Holo Ensemble Diff Mode**: Computes `dUIS`, `ddG`, `Kd` fold-change, two-tailed Welch's *t*-test *p*-values, and atomic `LOST`/`GAINED` Contact Rewiring audits.

## Quick CLI Usage
```bash
# Score any AF3 .zip, Boltz-1/Chai-1 .cif + JSON, or experimental PDB/mmCIF
python3 unippi_cli.py -i benchmarks/5ZC9.cif -o results_5zc9

# Custom physiological pH, ionic strength, or intra-chain domain split
python3 unippi_cli.py -i benchmarks/5ZC9.cif --split-chain A:1-236,237-406 --ph 7.2 --ionic-strength 150 -o results_split

# WT vs. Mutant / Apo vs. Holo statistical ensemble diff
python3 unippi_cli.py -i wt_ensemble.zip --diff-b mutant_ensemble.zip -o results_diff
```

## Validated Regression Benchmarks
| Benchmark | Top Interface | Modality | UIS | dG_uni (kcal/mol) | Pattern |
| :--- | :--- | :--- | :---: | :---: | :---: |
| `fold_vip1_bzip61.zip` (AF3 5-model) | `A-B` | `PROTEIN-PROTEIN` | `0.566` | `-10.99` | `Pattern B` |
| `1URN.cif` (U1A-RNA ITC) | `P-A` | `PROTEIN-RNA` | `0.983` | `-13.61` | `Pattern A` |
| `5ZC9.cif` (eIF4A1+RNA+RocA+ANP+Mg) | `A-B` (+ Glue & Mg-Bridge) | `PROTEIN-RNA` | `0.958` | `-12.88` | `Pattern A` |
| `4WNP.cif` (ULK1+3RJ Kinase) | `D-D_3RJ401` | `PROTEIN-LIGAND` | `0.728` | `-8.41` | `Pattern A` |
| `7DIY.cif` (Nsp14-Nsp10+Zn+Mg) | `A-B` (+ 4x Zn2+) | `PROTEIN-PROTEIN` | `0.847` | `-10.62` | `Pattern A` |
| `5ZC9.cif` (`--split-chain A:1-236,237-406`) | `A_D1+A_D2-A_ANP500` | `TERNARY-GLUE` | `0.907` | `-11.80` | `Pattern A` |

## Citation & Author
Architected by **Soumyajit Mandal** (Department of Biochemistry, Central University of Rajasthan, India).
Repository: [https://github.com/soumyajitmandal06-gif/unippi-owl-v2](https://github.com/soumyajitmandal06-gif/unippi-owl-v2) · License: MIT / BSD-3.
