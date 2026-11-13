# ST2P: Spatial Transcriptome-to-Protein Prediction

Reproducible implementation and experiment archive for cross-section spatial
RNA-to-protein prediction under the frozen 7-section benchmark.

## Included method for this local task

- **SCN-Pro / SCP-NET** (`src/scn_pro.py`): V10 local niche modeling augmented with an
  RNA-only whole-slice context, bounded niche canonicalization, and a
  low-rank slice-conditioned protein adapter for domain shift.

## Install

```bash
conda env create -f environment/environment.yml
conda activate st2p
```

## Data and artifacts

The raw/processed data and full result archives are not committed to Git.
Expected local layout:

```text
data/raw/DGAT_final7_breast_raw_registry_20261021/
data/processed/DGAT_final7_breast_ptprc31_muon017_axis0_20261022/
```

This local checkout intentionally downloads the two data registries only.
Large prediction/attention archives remain on the experiment server; Git
tracks their compact result summaries, configurations, and audit manifests.

See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) and
[manifests/DATASET_VERSION.md](manifests/DATASET_VERSION.md) for the frozen
protocol and artifact policy.

The formal experiment command and seven-target leave-one-section-out matrix
are maintained in `E:/DGAT/SCP-NET-work/docs/LODO_PROTOCOL_7_SECTION.md`.

## License and data terms

The datasets originate from their respective public providers. Before making
the repository public, retain source URLs, citation requirements, and any
dataset-specific redistribution restrictions in a data card. Do not upload
raw third-party matrices to GitHub without confirming their licenses.
