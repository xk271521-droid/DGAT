$ErrorActionPreference = 'Stop'
$python = 'C:\Users\xk2714\.conda\envs\st2p\python.exe'
$script = 'E:\DGAT\SCP-NET-work\experiments\nichetrans_10_dataset_20261101\run_nichetrans_10_dataset.py'
$root = 'E:\DGAT\ST2P\data\processed\DGAT_final10_common9384_muon017_axis0_20261101'
$out = 'E:\DGAT\SCP-NET-work\runs\nichetrans_10_dataset_20261102_corrected_breast_holdout'
$name = 'cross_holdout_Breast_Cancer'
$sources = @('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','GBM')
$targets = @('Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A')
$log = Join-Path $out 'queue.log'
New-Item -ItemType Directory -Force -Path $out | Out-Null
& $python $script --prepared_root $root --out_root $out --name $name --sources $sources --targets $targets 2>&1 | Tee-Object -FilePath $log
if($LASTEXITCODE -ne 0){ throw 'NicheTrans corrected Breast Cancer holdout failed' }
