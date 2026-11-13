#!/usr/bin/env python3
"""Build a raw-count seven-section registry with Breast Cancer as target-ready.

The six previously audited LN/Tonsil raw sections are copied only from their
raw-count registry. Breast Cancer is read directly from the original 10x H5;
its Visium array coordinates are aligned by barcode from the supplied spatial
archive. No transformed expression matrices are reused.
"""
from __future__ import annotations

import argparse
import io
import json
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


SECTION_NAMES = (
    "A1_LN", "D1_LN", "A1_TNSL", "D1_TNSL",
    "HumanTonsil_10x", "HumanTonsil_AddOns_S1",
)


def count_audit(matrix) -> dict[str, float | bool]:
    values = matrix.data if hasattr(matrix, "data") else np.asarray(matrix).ravel()
    values = np.asarray(values[: min(100_000, values.size)])
    return {
        "minimum": float(values.min()) if values.size else 0.0,
        "maximum": float(values.max()) if values.size else 0.0,
        "nonnegative": bool(np.all(values >= 0)),
        "integer_fraction_sample": float(np.mean(np.isclose(values, np.round(values)))) if values.size else 1.0,
    }


def write_pair(destination: Path, name: str, rna: ad.AnnData, adt: ad.AnnData, provenance: dict) -> dict:
    if rna.n_obs != adt.n_obs or not rna.obs_names.equals(adt.obs_names):
        raise ValueError(f"{name}: RNA/ADT barcode alignment failed")
    if "spatial" not in rna.obsm or rna.obsm["spatial"].shape != (rna.n_obs, 2):
        raise ValueError(f"{name}: invalid spatial coordinates")
    rna_audit, adt_audit = count_audit(rna.X), count_audit(adt.X)
    if not rna_audit["nonnegative"] or not adt_audit["nonnegative"]:
        raise ValueError(f"{name}: expected non-negative raw counts")
    section = destination / name
    section.mkdir(parents=False)
    rna.write_h5ad(section / "rna_raw.h5ad", compression="gzip")
    adt.write_h5ad(section / "adt_raw.h5ad", compression="gzip")
    audit = {
        "dataset": name,
        "n_spots": int(rna.n_obs),
        "n_rna_features": int(rna.n_vars),
        "n_adt_features": int(adt.n_vars),
        "rna_numeric": rna_audit,
        "adt_numeric": adt_audit,
        **provenance,
    }
    (section / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    return audit


def read_breast_pair(matrix_path: Path, spatial_tar: Path) -> tuple[ad.AnnData, ad.AnnData, dict]:
    multi = sc.read_10x_h5(matrix_path, gex_only=False)
    if "feature_types" not in multi.var:
        raise KeyError("Breast Cancer 10x H5 has no feature_types")
    kind = multi.var["feature_types"].astype(str)
    rna = multi[:, kind.eq("Gene Expression")].copy()
    adt = multi[:, kind.eq("Antibody Capture")].copy()
    with tarfile.open(spatial_tar, "r:gz") as archive:
        member = archive.getmember("spatial/tissue_positions.csv")
        with archive.extractfile(member) as handle:
            positions = pd.read_csv(io.TextIOWrapper(handle, encoding="utf-8"), index_col="barcode")
    if positions.index.has_duplicates:
        raise ValueError("Breast Cancer tissue_positions.csv has duplicate barcodes")
    positions = positions.loc[positions["in_tissue"].astype(int).eq(1)]
    shared = rna.obs_names.intersection(positions.index.astype(str), sort=False)
    if not len(shared):
        raise ValueError("No Breast Cancer 10x barcodes align to in-tissue spatial positions")
    rna, adt = rna[shared].copy(), adt[shared].copy()
    # x=array_col and y=array_row preserve the Visium grid geometry; absolute
    # orientation is immaterial to distance-based neighbourhood construction.
    rna.obsm["spatial"] = positions.loc[shared, ["array_col", "array_row"]].to_numpy(dtype=np.float32)
    return rna, adt, {
        "mode": "raw_10x_h5_plus_vendor_spatial_archive",
        "matrix_source": str(matrix_path),
        "spatial_source": f"{spatial_tar}:spatial/tissue_positions.csv",
        "coordinate_columns": ["array_col", "array_row"],
        "coordinate_scope": "in_tissue==1 only",
        "shared_barcodes": int(len(shared)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing_raw_registry", required=True, type=Path)
    parser.add_argument("--breast_matrix_h5", required=True, type=Path)
    parser.add_argument("--breast_spatial_tar", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite registry: {args.output}")
    args.output.mkdir(parents=True)
    records = []
    for name in SECTION_NAMES:
        source = args.existing_raw_registry / name
        if not (source / "rna_raw.h5ad").exists() or not (source / "adt_raw.h5ad").exists():
            raise FileNotFoundError(f"Missing audited raw source for {name}: {source}")
        rna, adt = ad.read_h5ad(source / "rna_raw.h5ad"), ad.read_h5ad(source / "adt_raw.h5ad")
        records.append(write_pair(args.output, name, rna, adt, {
            "mode": "reused_audited_raw_counts",
            "raw_registry_source": str(source),
        }))
        print(f"registered {name}: RNA={rna.shape}, ADT={adt.shape}", flush=True)
    rna, adt, provenance = read_breast_pair(args.breast_matrix_h5, args.breast_spatial_tar)
    records.append(write_pair(args.output, "Breast_Cancer", rna, adt, provenance))
    print(f"registered Breast_Cancer: RNA={rna.shape}, ADT={adt.shape}", flush=True)
    manifest = {
        "registry": "final seven-section raw-count registry: 2 LN, 4 Tonsil, 1 Breast Cancer",
        "datasets": records,
        "excluded_dataset": "SpaMosaic_HumanTonsil_S1",
    }
    (args.output / "registry_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
