# Reproducibility contract

## Frozen registry

All primary cross-section experiments use
`DGAT_final7_breast_ptprc31_muon017_axis0_20261022`:

- 7 sections: A1-LN, D1-LN, A1-TNSL, D1-TNSL, Human Tonsil 10x,
  Human Tonsil AddOns S1, and Breast Cancer.
- 9,764 aligned RNA genes and 31 aligned ADT proteins.
- RNA values follow the frozen registry preprocessing; ADT uses Muon/Seurat
  CLR with `axis=0`.
- PCA and target scaling statistics are fit on source-fit spots only.

## Strict selection rule

Every formal result selects its checkpoint by mean protein-wise PCC on
source-only validation spots. Target ADT is never used to select the epoch,
fit PCA, construct a target RNA neighborhood, or calibrate predictions.
Target metrics recorded at every epoch are diagnostics and are explicitly
labelled oracle-best when reported.

## Local artifact policy

`data/` and `artifacts/` are intentionally ignored by Git. This local
checkout keeps the raw/processed registries only. Large predictions,
attention outputs, per-protein metrics, configs, and split artifacts remain
on the experiment server and are represented in Git by compact summaries and
manifests. No epoch checkpoint is downloaded; if needed, retain only a
source-validation-best checkpoint in `artifacts/checkpoints/`.

## Result inventory

Run the following after syncing local artifacts:

```bash
python scripts/summarize_results.py --results-root artifacts/results \
  --output results/experiment_summary.csv
```

The generated table records the method, source sections, target section,
source-selected epoch, all reported target metrics, and test-oracle metrics.
