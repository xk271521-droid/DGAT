$ErrorActionPreference = "Stop"
$python = "C:\Users\xk2714\.conda\envs\st2p\python.exe"
$root = "E:\DGAT\ST2P\data\processed\DGAT_final10_common9384_muon017_axis0_20261101"
$experiment = "E:\DGAT\ST2P\src"
$script = "E:\DGAT\SCP-NET-work\experiments\relative_roughness_gated_selection\scn_pro_relative_roughness_gated.py"
$baseRunRoot = "E:\DGAT\SCP-NET-work\runs\relative_roughness_gated_runs"
$baseLogRoot = "E:\DGAT\SCP-NET-work\experiments\relative_roughness_gated_selection\logs"
New-Item -ItemType Directory -Force -Path $baseRunRoot,$baseLogRoot | Out-Null

$jobs = @(
    @{ Name="cross_holdout_GBM"; Target="GBM"; Sources=@("A1_LN","D1_LN","A1_TNSL","D1_TNSL","HumanTonsil_10x","HumanTonsil_AddOns_S1","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A"); TaskType="cross" },
    @{ Name="same_LN_A1_to_D1"; Target="D1_LN"; Sources=@("A1_LN"); TaskType="same" }
)

foreach ($job in $jobs) {
    $out = Join-Path $baseRunRoot $job.Name
    if (Test-Path $out) { Remove-Item -Recurse -Force $out }
    $log = Join-Path $baseLogRoot ($job.Name + ".log")
    Write-Host "Starting: $($job.Name)..." -ForegroundColor Cyan
    $args = @("-u", $script,
        "--prepared_root", $root, "--out_dir", $out, "--experiment_dir", $experiment,
        "--target", $job.Target, "--source_sections") + $job.Sources + @(
        "--task_type", $job.TaskType, "--growth_threshold", "2.0", "--growth_gamma", "0.010",
        "--seed","42","--epochs","20","--batch_size","64","--validation_ratio","0.1",
        "--hidden_dim","256","--heads","4","--ffn_dim","256","--dropout","0.3",
        "--protein_embedding_dim","32","--router_rank","32","--lr","0.0005","--weight_decay","0.0001",
        "--aux_weight","0.25","--pcc_weight","0.05","--pcc_min_variance","1e-6",
        "--context_pca_components","32","--context_rank","8","--gene_dropout_rate","0.5","--neighbor_mask_rate","0.25",
        "--local_neighbors","4","--regional_neighbors","8","--program_components","9")
    & $python @args 2>&1 | Tee-Object -FilePath $log
    Write-Host "Finished: $($job.Name)" -ForegroundColor Green
}
