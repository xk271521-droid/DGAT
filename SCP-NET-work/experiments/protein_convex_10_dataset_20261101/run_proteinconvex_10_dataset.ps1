$ErrorActionPreference = 'Stop'
$python = 'C:\Users\xk2714\.conda\envs\st2p\python.exe'
$root = 'E:\DGAT\ST2P\data\processed\DGAT_final10_common9384_muon017_axis0_20261101'
$experiment = 'E:\DGAT\ST2P\src'
$script = 'E:\DGAT\SCP-NET-work\experiments\protein_convex\scn_pro_protein_convex.py'
$runRoot = 'E:\DGAT\SCP-NET-work\runs\protein_convex_10_dataset_20261101'
$logRoot = 'E:\DGAT\SCP-NET-work\experiments\protein_convex_10_dataset_20261101\logs'
New-Item -ItemType Directory -Force -Path $runRoot,$logRoot | Out-Null

$jobs = @(
    @{ Name='same_LN_A1_to_D1'; Target='D1_LN'; Sources=@('A1_LN') },
    @{ Name='same_LN_D1_to_A1'; Target='A1_LN'; Sources=@('D1_LN') },
    @{ Name='same_TNSL_A1_to_D1'; Target='D1_TNSL'; Sources=@('A1_TNSL') },
    @{ Name='same_TNSL_D1_to_A1'; Target='A1_TNSL'; Sources=@('D1_TNSL') },
    @{ Name='same_HumanTonsil_10x_to_AddOns'; Target='HumanTonsil_AddOns_S1'; Sources=@('HumanTonsil_10x') },
    @{ Name='same_HumanTonsil_AddOns_to_10x'; Target='HumanTonsil_10x'; Sources=@('HumanTonsil_AddOns_S1') },
    @{ Name='same_PBC_6837_to_6835'; Target='PBC_PR_6835_5A'; Sources=@('PBC_PR_6837') },
    @{ Name='same_PBC_6835_to_6837'; Target='PBC_PR_6837'; Sources=@('PBC_PR_6835_5A') },
    @{ Name='cross_holdout_LN'; Target='A1_LN'; Sources=@('A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_Tonsil'; Target='A1_TNSL'; Sources=@('A1_LN','D1_LN','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_Breast_Cancer'; Target='Breast_Cancer'; Sources=@('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','PBC_PR_6837','PBC_PR_6835_5A','GBM') },
    @{ Name='cross_holdout_GBM'; Target='GBM'; Sources=@('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer','PBC_PR_6837','PBC_PR_6835_5A') }
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
