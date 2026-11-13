#!/usr/bin/env python3
import argparse, json, shutil
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import anndata as ad

def write_h5(path, x, barcodes, coords):
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("embeddings", data=x.astype(np.float32), compression="gzip")
        f.create_dataset("barcodes", data=np.asarray(barcodes, dtype="S"))
        f.create_dataset("coords", data=coords.astype(np.float32))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--prepared_root", type=Path, required=True)
    ap.add_argument("--spine_root", type=Path, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--sources", nargs="+", required=True)
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--validation_ratio", type=float, default=0.1)
    args=ap.parse_args()
    names=args.sources+args.targets
    m=json.loads((args.prepared_root/"prepared_manifest.json").read_text())
    all_names=[x["dataset"] for x in m["datasets"]]
    first=args.prepared_root/all_names[0]/"rna_dgat.h5ad"
    gene_sets=[]; protein_sets=[]
    for n in all_names:
        gene_sets.append(set(map(str,ad.read_h5ad(args.prepared_root/n/"rna_dgat.h5ad",backed="r").var_names)))
        protein_sets.append(set(map(str,ad.read_h5ad(args.prepared_root/n/"adt_dgat_clr.h5ad",backed="r").var_names)))
    genes=[g for g in map(str,ad.read_h5ad(first,backed="r").var_names) if all(g in s for s in gene_sets)]
    proteins=[p for p in map(str,ad.read_h5ad(args.prepared_root/all_names[0]/"adt_dgat_clr.h5ad",backed="r").var_names) if all(p in s for s in protein_sets)]
    if len(genes)!=9384 or len(proteins)!=31: raise RuntimeError((len(genes),len(proteins)))
    root=args.spine_root/args.name; shutil.rmtree(root,ignore_errors=True)
    dataset=root/"dataset"/args.name; embed=root/"dataset"/"embed_dataroot"/args.name; split=dataset/"splits"; split.mkdir(parents=True)
    (dataset/"rna_gene_list.json").write_text(json.dumps(genes),encoding="utf-8")
    (dataset/"protein_list.json").write_text(json.dumps(proteins),encoding="utf-8")
    train_ids=[]; val_ids=[]; test_ids=[]
    for n in names:
        rna=ad.read_h5ad(args.prepared_root/n/"rna_dgat.h5ad")[:,genes].copy()
        adt=ad.read_h5ad(args.prepared_root/n/"adt_dgat_clr.h5ad")[:,proteins].copy()
        if not rna.obs_names.equals(adt.obs_names): raise RuntimeError(f"barcode mismatch {n}")
        x=np.asarray(rna.X.toarray() if hasattr(rna.X,"toarray") else rna.X,dtype=np.float32)
        y=np.asarray(adt.X.toarray() if hasattr(adt.X,"toarray") else adt.X,dtype=np.float32)
        coords=np.asarray(rna.obsm["spatial"],dtype=np.float32)
        # Source spots are split internally: SPINE's official test loader is
        # used only on source validation IDs for checkpoint selection.
        if n in args.sources:
            rng=np.random.default_rng(args.seed+args.sources.index(n))
            order=rng.permutation(len(x)); n_val=max(1,int(round(len(x)*args.validation_ratio)))
            idx_map={n+"__fit":order[n_val:],n+"__val":order[:n_val]}
            train_ids.append(n+"__fit"); val_ids.append(n+"__val")
        else:
            idx_map={n:np.arange(len(x))}
            test_ids.append(n)
        for sample,idx in idx_map.items():
            xx=x[idx]; yy=y[idx]; cc=coords[idx]; bc=np.asarray(rna.obs_names).astype(str)[idx]
            (dataset/sample).mkdir(parents=True,exist_ok=True)
            write_h5(embed/sample/"RNA_EMBED"/(sample+".h5"),xx,bc,cc)
            yadata=ad.AnnData(X=yy,obs=pd.DataFrame(index=bc),var=pd.DataFrame(index=proteins))
            yadata.write_h5ad(dataset/sample/"protein.h5ad",compression="gzip")
    (split/"train_0.csv").write_text("sample_id\n"+"\n".join(train_ids)+"\n",encoding="utf-8")
    # Official validation/checkpoint selection sees only source-slice validation.
    (split/"test_0.csv").write_text("sample_id\n"+"\n".join(val_ids)+"\n",encoding="utf-8")
    (root/"test_targets.json").write_text(json.dumps(test_ids,indent=2),encoding="utf-8")
    (root/"manifest.json").write_text(json.dumps({"name":args.name,"source_sections":args.sources,"targets":args.targets,"train_sample_ids":train_ids,"validation_sample_ids":val_ids,"genes":len(genes),"proteins":len(proteins),"source_validation_ratio":args.validation_ratio},indent=2),encoding="utf-8")
    print(root)

if __name__=="__main__": main()
