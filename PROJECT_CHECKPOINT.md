# UniPPI-OWL v1.0 — Phase 2 Benchmark & Hardening Checkpoint

## Frozen Engine Physics (`unippi_cli.py` + `index.html`, 52/52 CLI <-> Browser Parity)
- PRODIGY `%NIS_charged` fixed: `0.38100 -> 0.13810` (Vangone & Bonvin 2015).
- Experimental `MODE_C_EXP_CRYSTAL` B-factor damping fixed: `platConf = 1.0`, `rawPlddt` floor `70.0`.
- Nucleic-acid attenuation scalar applied: `0.72` inside `16.5 * tanh((0.72 * E_raw) / 16.5)`, keeping `dG(0) = +2.10 kcal/mol`.

## Completed Phase 2 Benchmarks (Paused During Basics Hardening)
- **Track 1 (PPI — SKEMPI 2.0 WT, N=345):**
  - `dG_phys` (zero-shot): r=+0.303, rho=+0.276, RMSE=3.22, MAE=2.58, Bias=+1.35 (Post-2011 N=91: r=+0.280, RMSE=2.69, Bias=+0.43).
  - `dG_uni` (blend): r=+0.277, rho=+0.258, RMSE=3.10, MAE=2.45, Bias=+1.02 (`PR/PI` N=62: r=+0.569, rho=+0.617).
  - Fixed PRODIGY baseline: r=+0.184, RMSE=3.33 (Post-2011: r=+0.143, RMSE=3.20). Paired p <= 0.0042.
- **Track 2A (Protein-RNA — ProNAB WT, N=276 all-atom):**
  - `dG_RNA` (0.72 M1): r=+0.248, rho=+0.272, RMSE=3.42, MAE=2.86, Bias=-0.57, Recall=77.5%, Mean SSI=0.504 (Raw `-wBSA` r=+0.220).
- **Track 2B (Protein-DNA — ProNAB WT, N=624 all-atom):**
  - `dG_DNA` (0.72 M1): r=+0.230, rho=+0.219, RMSE=3.32, MAE=2.68, Bias=-0.23, Recall=84.0%, Mean SSI=0.326 (Raw `-wBSA` r=+0.163; paired p=0.0142).
- **Track 3B (Ternary Molecular Glues, N=8 + 1):**
  - 8/8 landmark glues detected (`5ZC9`, `1FAP`, `5FQD`, `6H0F`, `6UD7`, `2P1Q`, `1TCO`, `6TD3`); binary control `4WNP` clean.

## Pending Phase 2 Tracks (Resume ONLY After Basics Sanity Suite Passes 100%)
- **Track 3A:** High-resolution Protein-Ligand benchmark (N=285) + 8 multi-chain non-glue controls.
- **Track 4:** AlphaFold 3 Decoy & IDR-Rescue AUROC benchmark (`UIS` & `LIS` vs raw `ipTM`).
