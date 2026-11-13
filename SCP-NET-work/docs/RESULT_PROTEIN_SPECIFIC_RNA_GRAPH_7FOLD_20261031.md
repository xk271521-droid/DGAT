# Protein-specific supervised RNA graph calibration (2026-10-31)

## Purpose

Test whether the strict RNA graph calibration can be improved by making the RNA distance protein-specific. The RNA features used for each protein are selected from source-fit RNA--ADT correlations, then applied to the held-out target using RNA only.

## Protocol

- Frozen checkpoint: the existing `protein_convex_<target>/best_source_validation.pth` for each outer target.
- Outer evaluation: the same seven-target leave-one-section-out protocol.
- For each outer target, RNA--ADT correlations are pooled from the six source sections' fit spots only.
- A spatial 8-neighbor candidate graph is weighted by the protein-specific RNA distance.
- Spatial/RNA blend and per-protein smoothing strength are selected on the six source validation subsets.
- Held-out target ADT is read only for final metric calculation.
- Search setting audited here: top 32 RNA genes per protein, correlation exponent 1.0, blend grid `{0,.25,.5,.75,1}`, smoothing grid `{0,.1,.2,.3,.4,.5}`.

## Results

| mean metric | ProteinConvex base | protein-specific RNA graph |
|---|---:|---:|
| PCC-Pro | 0.586593 | 0.593651 |
| PCC-Spot | 0.629418 | 0.643185 |
| RMSE-Global | 0.285059 | 0.272483 |

Breast_Cancer:

| metric | base | protein-specific RNA graph |
|---|---:|---:|
| PCC-Pro | 0.602035 | 0.609520 |
| PCC-Spot | 0.633267 | 0.635647 |
| RMSE-Global | 0.299574 | 0.291509 |

## Decision

This is a valid strict calibration and improves the frozen base, but it is **not the new best method**. It is below `ProteinConvex + Strict RNA Calibration` on the primary mean PCC-Pro (0.593651 vs 0.594800) and PCC-Spot (0.643185 vs 0.644555), although its mean RMSE is marginally lower (0.272483 vs 0.272659). Keep it as an ablation/negative result; use the strict RNA calibration as the reported best pipeline.

## Artifacts

- Code: `experiments/protein_specific_rna_graph_20261031/protein_specific_rna_graph.py`
- Seven-fold metrics: `experiments/protein_specific_rna_graph_20261031/results_strict_source_validation.csv`
- Per-target predictions: `experiments/protein_specific_rna_graph_20261031/*_predictions.npy`
