$ErrorActionPreference = 'Stop'
$python = 'C:\Users\xk2714\.conda\envs\st2p\python.exe'
$root = 'E:\DGAT\ST2P\data\processed\DGAT_final10_common9384_muon017_axis0_20261101'
$experiment = 'E:\DGAT\ST2P\src'
$script = 'E:\DGAT\SCP-NET-work\experiments\protein_convex\scn_pro_protein_convex.py'
$runRoot = 'E:\DGAT\SCP-NET-work\runs\protein_convex_10_dataset_20261101'
$logRoot = 'E:\DGAT\SCP-NET-work\experiments\protein_convex_10_dataset_20261101\logs_corrected_cross_tissue'
New-Item -ItemType Directory -Force -Path $runRoot,$logRoot | Out-Null

$jobs = @(
    @{ Name='cross_holdout_LN_A1_LN'; Target='A1_LN'; Sources=@('A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_LN_D1_LN'; Target='D1_LN'; Sources=@('A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_Tonsil_A1_TNSL'; Target='A1_TNSL'; Sources=@('A1_LN','D1_LN','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_Tonsil_D1_TNSL'; Target='D1_TNSL'; Sources=@('A1_LN','D1_LN','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_Tonsil_HumanTonsil_10x'; Target='HumanTonsil_10x'; Sources=@('A1_LN','D1_LN','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_Tonsil_HumanTonsil_AddOns_S1'; Target='HumanTonsil_AddOns_S1'; Sources=@('A1_LN','D1_LN','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_Breast_Cancer_Breast_Cancer'; Target='Breast_Cancer'; Sources=@('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','GBM') },
    @{ Name='cross_holdout_Breast_Cancer_PBC_PR_6837'; Target='PBC_PR_6837'; Sources=@('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','GBM') },
    @{ Name='cross_holdout_Breast_Cancer_PBC_PR_6835_5A'; Target='PBC_PR_6835_5A'; Sources=@('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','GBM') }
)

foreach ($job in $jobs) {
    $out = Join-Path $runRoot $job.Name
    if (Test-Path $out) { throw "Refusing to overwrite existing output: $out" }
    $args = @('-u', $script,
        '--prepared_root', $root, '--out_dir', $out, '--experiment_dir', $experiment,
        '--target', $job.Target, '--source_sections') + $job.Sources + @(
        '--seed','42','--epochs','20','--batch_size','64','--validation_ratio','0.1',
        '--hidden_dim','256','--heads','4','--ffn_dim','256','--dropout','0.3',
        '--protein_embedding_dim','32','--router_rank','32','--lr','0.0005','--weight_decay','0.0001',
        '--aux_weight','0.25','--pcc_weight','0.05','--pcc_min_variance','1e-6',
        '--context_pca_components','32','--context_rank','8','--gene_dropout_rate','0.5','--neighbor_mask_rate','0.25',
        '--local_neighbors','4','--regional_neighbors','8','--program_components','9')
    $log = Join-Path $logRoot ($job.Name + '.log')
    & $python @args 2>&1 | Tee-Object -FilePath $log
    if ($LASTEXITCODE -ne 0) { throw "ProteinConvex failed for $($job.Name) with exit code $LASTEXITCODE" }
}
