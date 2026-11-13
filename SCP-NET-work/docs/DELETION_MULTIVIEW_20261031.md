# MultiViewConsistency cleanup record (2026-10-31)

At the user's request, the active MultiViewConsistency implementation and run artifacts were removed from the workspace.

## Removed

- `SCP-NET-work/experiments/multiview_consistency/`
- `SCP-NET-work/runs/full_lodo_multiview_consistency_20261029/`
- MultiView external-validation JSON and log files under `SCP-NET-work/experiments/pbc_external_validation_20261030/results/`
- `SCP-NET-work/experiments/spatial_smoothing_20261030/analyze_spatial_smoothing.py`, which depended on the removed MultiView run;
- MultiView rows from the active CSV summaries `outputs/retained_models/experiment_summary.csv` and `ST2P/results/experiment_summary.csv`.

## Retained for audit history

- Historical result documents that compare MultiViewConsistency;
- Existing summary snapshots and tables generated before the cleanup;
- Shared preprocessing and helper modules used by other retained models.

MultiViewConsistency is no longer an active retained model track and must not be used in new training, evaluation, fusion, or benchmark claims.
