$ErrorActionPreference = 'Stop'
$python = 'C:\Users\xk2714\.conda\envs\st2p\python.exe'
$script = 'E:\DGAT\SCP-NET-work\experiments\nichetrans_10_dataset_20261101\run_nichetrans_10_dataset.py'
$root = 'E:\DGAT\ST2P\data\processed\DGAT_final10_common9384_muon017_axis0_20261101'
$out = 'E:\DGAT\SCP-NET-work\runs\nichetrans_10_dataset_20261101_fixed_panel'
New-Item -ItemType Directory -Force $out | Out-Null
$jobs = @(
 @{Name='same_LN_A1_to_D1'; S=@('A1_LN'); T=@('D1_LN')}, @{Name='same_LN_D1_to_A1'; S=@('D1_LN'); T=@('A1_LN')},
 @{Name='same_TNSL_A1_to_D1'; S=@('A1_TNSL'); T=@('D1_TNSL')}, @{Name='same_TNSL_D1_to_A1'; S=@('D1_TNSL'); T=@('A1_TNSL')},
 @{Name='same_HumanTonsil_10x_to_AddOns'; S=@('HumanTonsil_10x'); T=@('HumanTonsil_AddOns_S1')}, @{Name='same_HumanTonsil_AddOns_to_10x'; S=@('HumanTonsil_AddOns_S1'); T=@('HumanTonsil_10x')},
 @{Name='same_PBC_6837_to_6835'; S=@('PBC_PR_6837'); T=@('PBC_PR_6835_5A')}, @{Name='same_PBC_6835_to_6837'; S=@('PBC_PR_6835_5A'); T=@('PBC_PR_6837')},
 @{Name='cross_holdout_LN'; S=@('A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM'); T=@('A1_LN','D1_LN')},
 @{Name='cross_holdout_Tonsil'; S=@('A1_LN','D1_LN','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM'); T=@('A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1')},
 @{Name='cross_holdout_Breast_Cancer'; S=@('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','GBM'); T=@('Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A')},
 @{Name='cross_holdout_GBM'; S=@('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A'); T=@('GBM')}
)
foreach($j in $jobs){
  $args=@($script,'--prepared_root',$root,'--out_root',$out,'--name',$j.Name,'--sources')+$j.S+@('--targets')+$j.T
  $log=Join-Path $out ($j.Name+'.log')
  & $python @args 2>&1 | Tee-Object -FilePath $log
  if($LASTEXITCODE -ne 0){ throw "failed $($j.Name)" }
}
