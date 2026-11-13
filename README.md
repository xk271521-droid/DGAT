# DGAT project code archive

This repository archives the current local research code for spatial RNA-to-ADT
prediction. The current project method is **ProteinConvex + 初始粗糙度尺度门控 +
Strict RNA + Source Affine**. It also contains the retained baseline code,
experiment adapters, historical diagnostic code, and fixed upstream sources.

本次归档日期：2026-11-12。对应私有仓库：
`https://github.com/xk271521-droid/DGAT`。

## Current method and frozen protocol

- Input: 10 datasets, 9384 common RNA genes and 31 common ADT proteins.
- Evaluation: 8 same-tissue targets and 10 cross-tissue targets.
- Checkpoint selection: same-tissue tasks use `alpha=0`; cross-tissue tasks use
  `alpha=15` when initial spatial roughness is below `0.008`, otherwise `alpha=0`.
  The score subtracts the current **absolute** roughness penalty from source
  validation PCC-Pro. The historical directory name `relative_roughness_gated`
  does not mean this archive implements an `R(e)/R(1)` selection rule.
- Strict RNA and Source Affine choices are fitted using source validation data.
  Target ADT is read after selection for final evaluation.
- Code labelled `oracle`, including trajectory analyses, is diagnostic code and
  is not the formal label-free checkpoint-selection method.

## Code map

| Location | Contents |
|---|---|
| `ST2P/src/` | Shared model layers, preprocessing, metrics and retained historical models |
| `ST2P/environment/` | Saved environment specification and history |
| `SCP-NET-work/experiments/relative_roughness_gated_selection/` | Current method, 18-direction launcher and source-only affine calibration |
| `SCP-NET-work/experiments/` | ProteinConvex baseline, cognate anchor, four comparator adapters and diagnostics |
| `SCP-NET-work/runs/` | Retained source snapshots, parameters and compact metrics; large artifacts are external |
| `SCP-NET-work/transfer/` | Retained server launch/adaptation scripts |
| `upstream/` | Fixed original DGAT, DPAS-Graph, NicheTrans and SPINE repositories, with provenance and license files |
| `results/` | Compact CSV/JSON/Markdown result and diagnostic records |
| `figures/` | Retained figure-generation source |
| `tools/export_source_snapshot.py` | Reusable byte-preserving export and verification tool |

The original authors' files remain distinct from our adapters. See
[`upstream/README.md`](upstream/README.md) for fixed commit IDs and the scope of
the public sources. `GIT_SOURCE_MANIFEST.json` records every exported file's
size, SHA-256 and provenance. It deliberately does not hash itself.

## Reproduction entry points

The original scripts and their hyperparameters are archived unchanged. Existing
launchers contain machine-specific Windows or server paths. A fresh clone needs
the frozen data, the appropriate dependencies and local path configuration;
this source upload alone is not a completed fresh-environment reproduction.

1. Restore the frozen data under
   `ST2P/data/processed/DGAT_final10_common9384_muon017_axis0_20261101/`.
   The common panel is recorded in
   `SCP-NET-work/docs/ten_dataset_common_panel_20261101.json`.
2. Consult `ST2P/environment/environment.yml` and the recorded run configuration.
   Comparator dependencies are specified in each original upstream repository.
3. Current training entry:
   `SCP-NET-work/experiments/relative_roughness_gated_selection/scn_pro_relative_roughness_gated.py`.
   The complete parameter and source/target matrix is in
   `run_all_18_background.ps1` in the same directory. That historical launcher
   can remove incomplete output directories, so use a separate output location
   when reproducing it.
4. Current calibration entry:
   `run_source_only_rna_affine_calibration.py` in the same directory. Supply
   `--run-root SCP-NET-work/runs/relative_roughness_gated_runs` explicitly:
   the historical default points at the baseline run root. Use a new
   `--output-root` to avoid overwriting retained results.
5. Consult the recorded parameters and results alongside each source snapshot;
   no checkpoint or prediction matrix is bundled in this source-only archive.

The calibration process dynamically loads `script_snapshot.py` from a run
directory, so the retained run snapshots are intentionally included.

## Results and limitations

The grouped averages are in
[`results/model_grouped_averages_20261106/MODEL_GROUPED_AVERAGES_20261106.csv`](results/model_grouped_averages_20261106/MODEL_GROUPED_AVERAGES_20261106.csv).
For the current method, the 18-target macro averages are PCC-Pro `0.4694787442`,
RMSE-Global `0.2885851786`, and MAE-Global approximately `0.1852496849`.
These are **equal-target macro averages**, including the mean of each target's
RMSE, rather than a new pooled prediction-level RMSE.

Historical documents preserve earlier experiments and may describe earlier
7-section or 9764-gene protocols. The current 10-dataset method is defined by
the entry points above and
[`SCP-NET-work/docs/REPRODUCIBILITY_AUDIT_20261104.md`](SCP-NET-work/docs/REPRODUCIBILITY_AUDIT_20261104.md).
The upload has not rerun training or changed datasets, parameters, or scores.

## Scope of the upload

Raw/preprocessed matrices, model weights, prediction/attention arrays, virtual
environments, downloaded wheels, redundant archives, temporary files, chat dumps
and credentials remain outside this repository. Compact run settings, metrics,
source snapshots, and the complete public tracked files of each fixed upstream
commit are included. When an upstream commit itself tracked bytecode, the
original tracked blob is retained for archive completeness; locally regenerated
bytecode is excluded.

The publishing checkout is `E:/DGAT/git-publication/`. It is independent of the
older nested `ST2P/.git` repository and does not change that repository's remote
or history. To prepare a later snapshot locally:

```powershell
python tools/export_source_snapshot.py --source-root E:/DGAT --output-root E:/DGAT/git-publication
```

The exporter verifies copies, inventories source coverage and scans for common
credential formats; it does not publish anything by itself. Review its report
before committing. Source and model files are copied byte for byte.

Upstream license files and authorship are retained. No new blanket license is
applied to other authors' work.
