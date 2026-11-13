#!/usr/bin/env python3
"""Run official SPINE inference on held-out samples after checkpoint selection."""
import argparse, json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import os

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--spine_repo',type=Path,required=True)
    ap.add_argument('--dataset_root',type=Path,required=True)
    ap.add_argument('--embed_root',type=Path,required=True)
    ap.add_argument('--dataset',required=True)
    ap.add_argument('--target_ids',nargs='+',required=True)
    ap.add_argument('--model_dir',type=Path,required=True)
    ap.add_argument('--gpu',default=None,help='physical GPU selected by the queue; use -1 for CPU')
    a=ap.parse_args()
    if a.gpu is not None:
        if str(a.gpu) == '-1':
            os.environ['CUDA_VISIBLE_DEVICES'] = ''
        else:
            os.environ['CUDA_VISIBLE_DEVICES'] = str(a.gpu)
    import torch
    sys.path.insert(0,str(a.spine_repo))
    from spine.data.dataset import RNAToProteinDatasetPath,RNAToProteinDataset
    from spine.model.denoiser import Denoiser
    from spine.flow.interpolant import Interpolant
    from spine.app.flow.test import test
    from types import SimpleNamespace
    cfg=json.loads((a.model_dir/'config.json').read_text())
    cfg['device']=torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    cfg['save_dir']=str(a.model_dir)
    args=SimpleNamespace(**cfg)
    model=Denoiser(args).to(args.device)
    ckpt=torch.load(a.model_dir/'split0'/'checkpoints'/'best.pth',map_location=args.device,weights_only=True)
    model.load_state_dict(ckpt)
    norm_method=None
    if cfg.get('normalize_method') not in (None,'None'):
        from spine.data.normalize_utils import get_normalize_method
        norm_method=get_normalize_method(cfg['normalize_method'])
    feature_normalize=cfg.get('rna_feature_normalize')
    paths=[RNAToProteinDatasetPath(name=n,h5_path=str(a.embed_root/a.dataset/n/'RNA_EMBED'/f'{n}.h5'),protein_h5ad_path=str(a.dataset_root/a.dataset/n/'protein.h5ad'),protein_list_path=str(a.dataset_root/a.dataset/'protein_list.json'),rna_gene_list_path=str(a.dataset_root/a.dataset/'rna_gene_list.json')) for n in a.target_ids]
    loaders=[]
    for p in paths:
        ds=RNAToProteinDataset(p,normalize_method=norm_method,distribution='constant_1.0',sample_times=1,feature_normalize_method=feature_normalize)
        loaders.append(torch.utils.data.DataLoader(ds,batch_size=1,collate_fn=__import__('spine.data.dataset',fromlist=['padding_batcher']).padding_batcher(),num_workers=0))
    interpolant=Interpolant(cfg.get('prior_sampler','gaussian'),device=str(args.device),normalize=cfg.get('prior_sampler','gaussian')!='gaussian')
    summary={}
    summary_root=a.model_dir/'heldout_predictions'; summary_root.mkdir(parents=True,exist_ok=True)
    for p,loader in zip(paths,loaders):
        n=p.name
        one_summary, packed = test(args,interpolant,model,[loader],return_all=True)
        summary[n]=one_summary.get(n, one_summary.get('all'))
        pred=np.asarray(packed['preds_all'],dtype=np.float32)
        truth=np.asarray(packed['targets_all'],dtype=np.float32)
        np.save(summary_root/f'{n}_predictions.npy',pred)
        np.save(summary_root/f'{n}_truth.npy',truth)
        # Analysis-ready files retain spot barcodes and the official protein order.
        with __import__('h5py').File(p.h5_path,'r') as hf:
            barcodes=hf['barcodes'][:].astype(str)
        proteins=json.loads((a.dataset_root/a.dataset/'protein_list.json').read_text())
        ad.AnnData(X=pred,obs=pd.DataFrame(index=barcodes),var=pd.DataFrame(index=proteins)).write_h5ad(summary_root/f'{n}_predictions.h5ad',compression='gzip')
        pd.DataFrame(pred,index=barcodes,columns=proteins).to_csv(summary_root/f'{n}_predictions.csv.gz',compression='gzip')
        (summary_root/f'{n}_metrics.json').write_text(json.dumps(summary[n],indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
    (summary_root/'summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
    print('HELDOUT_PREDICTIONS_SAVED',summary_root,flush=True)

if __name__=='__main__': main()
