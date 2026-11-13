from __future__ import annotations
import argparse, copy, json, random, sys
from pathlib import Path
import anndata as ad
import numpy as np
import scipy.stats
import torch
from torch_geometric.loader import DataLoader

def find_root():
    p=Path(__file__).resolve()
    for q in [p,*p.parents]:
        if (q/'ST2P'/'data'/'processed').exists() and (q/'upstream').exists(): return q
    return p.parents[3]
ROOT=find_root(); UP=ROOT/'upstream'/'DGAT'; sys.path.insert(0,str(UP))
from utils.Graph_utils import MultiGraphDataset, MultiGraphDataset_for_no_protein
from Model.dgat import GATEncoder, Decoder_Protein, Decoder_mRNA

ALPHA,BETA,GAMMA,DELTA,EPSILON=5.0,1.0,1.0,3.0,5.0
HIDDEN=1024; DROPOUT=0.3; EPOCHS=100

def seed_all(s):
    random.seed(s); np.random.seed(s); torch.manual_seed(s)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(s)

def load_panel(path):
    p=json.loads(Path(path).read_text(encoding='utf-8'))
    g=[str(x) for x in p['rna_common_genes']]; pr=[str(x) for x in p['adt_common_proteins']]
    if len(g)!=9384 or len(pr)!=31: raise RuntimeError((len(g),len(pr)))
    return g,pr

def load_pair(root,n,genes,proteins):
    r=ad.read_h5ad(root/n/'rna_dgat.h5ad')[:,genes].copy(); p=ad.read_h5ad(root/n/'adt_dgat_clr.h5ad')[:,proteins].copy()
    common=r.obs_names.intersection(p.obs_names); r=r[common].copy(); p=p[common].copy(); r.uns['name']=n; p.uns['name']=n
    if 'spatial' not in r.obsm: raise RuntimeError(f'{n}: missing spatial')
    return r,p

def split_pair(r,p,seed):
    idx=np.arange(r.n_obs); np.random.default_rng(seed).shuffle(idx); n=max(1,int(round(.1*len(idx))))
    vi=np.sort(idx[:n]); fi=np.sort(idx[n:]); return r[fi].copy(),p[fi].copy(),r[vi].copy(),p[vi].copy()

def batches(ds,shuffle): return DataLoader(ds,batch_size=1,shuffle=shuffle)

def metrics(yh,yt):
    pc=[]
    for j in range(yt.shape[1]):
        if np.std(yt[:,j])>1e-8 and np.std(yh[:,j])>1e-8: pc.append(float(scipy.stats.pearsonr(yh[:,j],yt[:,j]).statistic))
    return {'pcc_pro':float(np.mean(pc)) if pc else float('nan'),'rmse_global':float(np.sqrt(np.mean((yh-yt)**2))),'n_valid_proteins':len(pc)}

def forward(data, em, dm, ep, dp, device):
    mx=data['mRNA'].x.to(device); mei=data[('mRNA','mRNA_knn','mRNA')].edge_index.to(device); mea=data[('mRNA','mRNA_knn','mRNA')].edge_attr.to(device)
    px=data['protein'].x.to(device); pei=data[('protein','protein_knn','protein')].edge_index.to(device); pea=data[('protein','protein_knn','protein')].edge_attr.to(device)
    zm=em(mx,mei,mea); zp=ep(px,pei,pea); return mx,px,zm,zp,dm(zm),dp(zp)

def eval_loader(loader,em,dm,ep,dp,device):
    em.eval();dm.eval();ep.eval();dp.eval(); preds=[]; truths=[]
    with torch.no_grad():
        for d in loader:
            mx,px,zm,zp,_,_=forward(d,em,dm,ep,dp,device); preds.append(dp(zm).cpu().numpy()); truths.append(px.cpu().numpy())
    return np.vstack(preds),np.vstack(truths)

def run_one(name,sources,targets,args,genes,proteins):
    out=Path(args.out_root)/name; out.mkdir(parents=True,exist_ok=True); logp=out/'run.log'
    def log(x): print(x,flush=True); logp.open('a',encoding='utf-8').write(str(x)+'\n')
    seed_all(args.seed); root=Path(args.prepared_root); device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    pairs={n:load_pair(root,n,genes,proteins) for n in sources}; fit_r=[];fit_p=[];val_r=[];val_p=[]
    for i,n in enumerate(sources):
        a,b,c,d=split_pair(*pairs[n],args.seed+i); fit_r.append(a);fit_p.append(b);val_r.append(c);val_p.append(d)
    gdir=out/'graphs'; train=MultiGraphDataset(fit_r,fit_p,save_dir=str(gdir/'fit')); val=MultiGraphDataset(val_r,val_p,save_dir=str(gdir/'val'))
    tl=batches(train,True); vl=batches(val,False)
    em=GATEncoder(len(genes),HIDDEN,dropout=DROPOUT).to(device); dm=Decoder_mRNA(HIDDEN,len(genes),dropout=0).to(device); ep=GATEncoder(len(proteins),HIDDEN,dropout=DROPOUT).to(device); dp=Decoder_Protein(HIDDEN,proteins,dropout=0).to(device)
    oms=[torch.optim.Adam(m.parameters(),lr=lr,weight_decay=wd) for m,lr,wd in [(em,5e-4,2e-5),(dm,1e-4,1e-5),(ep,5e-4,2e-5),(dp,1e-4,1e-5)]]
    mse=torch.nn.MSELoss(); best=-1e9; best_state=None; hist=[]; log(f'[{name}] device={device} genes={len(genes)} proteins={len(proteins)} sources={sources} targets={targets}')
    for epoch in range(1,args.epochs+1):
        em.train();dm.train();ep.train();dp.train(); losses=[]
        for data in tl:
            mx,px,zm,zp,mr,pr=forward(data,em,dm,ep,dp,device); lr=mse(mr,mx).sqrt(); lpr=mse(pr,px).sqrt(); la=mse(zm,zp); lpp=mse(dp(zm),px).sqrt(); lmp=mse(dm(zp),mx).sqrt(); loss=ALPHA*lr+BETA*lpr+GAMMA*la+DELTA*lpp+EPSILON*lmp
            for o in oms: o.zero_grad()
            loss.backward()
            for m in (em,dm,ep,dp): torch.nn.utils.clip_grad_norm_(m.parameters(),1.0)
            for o in oms:o.step()
            losses.append(float(loss.detach().cpu()))
        vp,vt=eval_loader(vl,em,dm,ep,dp,device); vm=metrics(vp,vt); hist.append({'epoch':epoch,'train_loss':float(np.mean(losses)),**{f'val_{k}':v for k,v in vm.items()}}); log(f'[{name}] epoch={epoch}/{args.epochs} loss={np.mean(losses):.5f} val_pcc={vm["pcc_pro"]:.4f}')
        if vm['pcc_pro']>best:
            best=vm['pcc_pro']; best_state={k:copy.deepcopy(v.state_dict()) for k,v in [('em',em),('dm',dm),('ep',ep),('dp',dp)]}; torch.save(best_state,out/'best_state.pt')
        if epoch%10==0:
            for o in oms:
                for pg in o.param_groups: pg['lr']*=.8
    for k,v in [('em',em),('dm',dm),('ep',ep),('dp',dp)]: v.load_state_dict(best_state[k])
    target_metrics={}
    for n in targets:
        r,p=load_pair(root,n,genes,proteins); td=MultiGraphDataset_for_no_protein([r],save_dir=str(gdir/f'target_{n}')); pred=[]
        em.eval();dp.eval()
        with torch.no_grad():
            for d in DataLoader(td,batch_size=1,shuffle=False):
                mx=d['mRNA'].x.to(device); ei=d[('mRNA','mRNA_knn','mRNA')].edge_index.to(device); ea=d[('mRNA','mRNA_knn','mRNA')].edge_attr.to(device); pred.append(dp(em(mx,ei,ea)).cpu().numpy())
        yh=np.vstack(pred).astype(np.float32); yt=np.asarray(p.X.toarray() if hasattr(p.X,'toarray') else p.X,dtype=np.float32); m=metrics(yh,yt); target_metrics[n]=m; np.save(out/f'{n}_predictions.npy',yh); np.save(out/f'{n}_truth.npy',yt); (out/f'{n}_metrics.json').write_text(json.dumps(m,indent=2),encoding='utf-8')
    (out/'metrics_history.json').write_text(json.dumps(hist,indent=2),encoding='utf-8'); summary={'name':name,'sources':sources,'targets':targets,'genes':len(genes),'proteins':proteins,'best_epoch':max(range(1,len(hist)+1),key=lambda i:hist[i-1]['val_pcc_pro']),'best_source_validation_pcc_pro':best,'target_metrics':target_metrics,'official_parameters':{'alpha':ALPHA,'beta':BETA,'gamma':GAMMA,'delta':DELTA,'epsilon':EPSILON,'hidden_dim':HIDDEN,'dropout':DROPOUT,'epochs':args.epochs,'seed':args.seed},'device':str(device)}; (out/'results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8'); log(f'[{name}] COMPLETE targets={targets}')

DIRECTIONS=[('same_LN_A1_to_D1',['A1_LN'],['D1_LN']),('same_LN_D1_to_A1',['D1_LN'],['A1_LN']),('same_TNSL_A1_to_D1',['A1_TNSL'],['D1_TNSL']),('same_TNSL_D1_to_A1',['D1_TNSL'],['A1_TNSL']),('same_HumanTonsil_10x_to_AddOns',['HumanTonsil_10x'],['HumanTonsil_AddOns_S1']),('same_HumanTonsil_AddOns_to_10x',['HumanTonsil_AddOns_S1'],['HumanTonsil_10x']),('same_PBC_6837_to_6835',['PBC_PR_6837'],['PBC_PR_6835_5A']),('same_PBC_6835_to_6837',['PBC_PR_6835_5A'],['PBC_PR_6837']),('cross_holdout_LN',['A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM'],['A1_LN','D1_LN']),('cross_holdout_Tonsil',['A1_LN','D1_LN','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM'],['A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1']),('cross_holdout_Breast_Cancer',['A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','PBC_PR_6837','PBC_PR_6835_5A','GBM'],['Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A']),('cross_holdout_GBM',['A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A'],['GBM'])]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--prepared_root',required=True); ap.add_argument('--panel',required=True); ap.add_argument('--out_root',required=True); ap.add_argument('--epochs',type=int,default=100); ap.add_argument('--seed',type=int,default=0); ap.add_argument('--directions',default='all'); a=ap.parse_args(); genes,proteins=load_panel(a.panel); chosen={n for n,_,_ in DIRECTIONS} if a.directions=='all' else set(a.directions.split(','));
    for n,s,t in DIRECTIONS:
        if n in chosen: run_one(n,s,t,a,genes,proteins)
if __name__=='__main__': main()
