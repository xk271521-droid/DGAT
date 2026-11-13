#!/usr/bin/env python3
"""Freeze paired raw sections with DGAT's exact Muon-0.1.7 ADT CLR.

This script deliberately performs no train/test split and no panel selection.
It just applies the published DGAT transform to each paired raw-count section
once, so V9 and DGAT can later consume byte-identical transformed inputs.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np


def muon_017_protein_clr_axis0(adt: ad.AnnData) -> ad.AnnData:
    """Replicate ``muon.prot.pp.clr(adata)`` from Muon 0.1.7 exactly.

    DGAT calls this function without an ``axis`` argument.  In that pinned
    version the default is ``axis=0``: each protein is scaled across spots as
    ``log1p(x / exp(mean(log1p(x))))``.  The ADT panel is small, so materialising
    it keeps this implementation transparent and avoids a runtime dependency on
    Muon in the training environment.
    """
    raw = adt.layers["dgat_raw"]
    values = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)
    values = np.asarray(values, dtype=np.float64)
    denominator = np.exp(np.log1p(values).sum(axis=0, keepdims=True) / values.shape[0])
    adt.X = np.log1p(values / denominator).astype(np.float32)
    return adt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw_registry", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--experiment_dir", required=True)
    args = ap.parse_args()

    experiment_dir = Path(args.experiment_dir)
    sys.path.insert(0, str(experiment_dir))
    from triscale_center_protein_v9_dgat_compatible import dgat_preprocess  # noqa: PLC0415

    raw_root = Path(args.raw_registry)
    manifest = json.loads((raw_root / "registry_manifest.json").read_text())
    out_root = Path(args.out_dir)
    if out_root.exists():
        raise FileExistsError(f"Refusing to overwrite prepared registry: {out_root}")
    out_root.mkdir(parents=True)

    records = []
    for record in manifest["datasets"]:
        name = record["dataset"]
        section_dir = raw_root / name
        raw_rna = ad.read_h5ad(section_dir / "rna_raw.h5ad")
        raw_adt = ad.read_h5ad(section_dir / "adt_raw.h5ad")
        # Reuse the previously audited RNA QC/normalisation and ADT panel QC.
        # ``dgat_preprocess`` also harmonizes the documented PTPRC_1/2
        # (CD45RA/CD45RO) vendor identifiers before a cross-section panel is
        # selected; no count values or CLR formula are changed.
        # Its former row-wise ADT transform is immediately replaced below.
        rna, adt = dgat_preprocess(raw_rna, raw_adt)
        adt = muon_017_protein_clr_axis0(adt)
        target_var = np.var(np.asarray(adt.X), axis=0) > 0
        if not np.all(target_var):
            adt = adt[:, target_var].copy()
        destination = out_root / name
        destination.mkdir()
        rna.write_h5ad(destination / "rna_dgat.h5ad", compression="gzip")
        adt.write_h5ad(destination / "adt_dgat_clr.h5ad", compression="gzip")
        audit = {
            "dataset": name,
            "raw_registry": str(section_dir),
            "rna_transform": "DGAT QC; library-size normalize to 10,000; log1p; gene-wise scale; clip 10",
            "adt_transform": "Exact DGAT dependency muon==0.1.7 default CLR: axis=0, log1p(x / exp(mean_spots(log1p(x)))); isotype controls removed",
            "adt_identifier_harmonization": "Vendor gene_ids or feature_id PTPRC_1/PTPRC_2 are used to align the two CD45RA/CD45RO channels across duplicate-display-name and AddOns naming conventions; all other panel names retain their cleaned display labels.",
            "raw_rna_shape": [int(raw_rna.n_obs), int(raw_rna.n_vars)],
            "raw_adt_shape": [int(raw_adt.n_obs), int(raw_adt.n_vars)],
            "post_qc_rna_shape": [int(rna.n_obs), int(rna.n_vars)],
            "post_qc_adt_shape": [int(adt.n_obs), int(adt.n_vars)],
            "rna_file": "rna_dgat.h5ad",
            "adt_file": "adt_dgat_clr.h5ad",
        }
        (destination / "audit.json").write_text(json.dumps(audit, indent=2))
        records.append(audit)
        print(f"prepared {name}: RNA {rna.shape}, ADT {adt.shape}", flush=True)

    common_genes = set(ad.read_h5ad(out_root / records[0]["dataset"] / "rna_dgat.h5ad", backed="r").var_names)
    common_proteins = set(ad.read_h5ad(out_root / records[0]["dataset"] / "adt_dgat_clr.h5ad", backed="r").var_names)
    for record in records[1:]:
        common_genes &= set(ad.read_h5ad(out_root / record["dataset"] / "rna_dgat.h5ad", backed="r").var_names)
        common_proteins &= set(ad.read_h5ad(out_root / record["dataset"] / "adt_dgat_clr.h5ad", backed="r").var_names)
    result = {
        "purpose": "single frozen exact-Muon-0.1.7-DGAT input registry for fair paired-section benchmarking",
        "adt_clr_audit": {"implementation": "muon.prot.pp.clr default semantics", "muon_version": "0.1.7", "axis": 0},
        "source_registry": str(raw_root),
        "datasets": records,
        "all_section_common_gene_count": len(common_genes),
        "all_section_common_protein_count": len(common_proteins),
        "common_panel_is_selected_later": True,
    }
    (out_root / "prepared_manifest.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
