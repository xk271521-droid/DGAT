$ErrorActionPreference = "Continue"
$python = "C:\Users\xk2714\.conda\envs\st2p\python.exe"
$root = "E:\DGAT\ST2P\data\processed\DGAT_final10_common9384_muon017_axis0_20261101"
$experiment = "E:\DGAT\ST2P\src"
$script = "E:\DGAT\SCP-NET-work\experiments\relative_roughness_gated_selection\scn_pro_relative_roughness_gated.py"
$baseRunRoot = "E:\DGAT\SCP-NET-work\runs\relative_roughness_gated_runs"
$baseLogRoot = "E:\DGAT\SCP-NET-work\experiments\relative_roughness_gated_selection\logs"
New-Item -ItemType Directory -Force -Path $baseRunRoot,$baseLogRoot | Out-Null

$jobs = @(
    # --- 8 个同组织双向方向 ---
    @{ Name="same_LN_A1_to_D1"; Target="D1_LN"; Sources=@("A1_LN"); TaskType="same" },
    @{ Name="same_LN_D1_to_A1"; Target="A1_LN"; Sources=@("D1_LN"); TaskType="same" },
    @{ Name="same_TNSL_A1_to_D1"; Target="D1_TNSL"; Sources=@("A1_TNSL"); TaskType="same" },
    @{ Name="same_TNSL_D1_to_A1"; Target="A1_TNSL"; Sources=@("D1_TNSL"); TaskType="same" },
    @{ Name="same_HumanTonsil_10x_to_AddOns"; Target="HumanTonsil_AddOns_S1"; Sources=@("HumanTonsil_10x"); TaskType="same" },
    @{ Name="same_HumanTonsil_AddOns_to_10x"; Target="HumanTonsil_10x"; Sources=@("HumanTonsil_AddOns_S1"); TaskType="same" },
    @{ Name="same_PBC_6837_to_6835"; Target="PBC_PR_6835_5A"; Sources=@("PBC_PR_6837"); TaskType="same" },
    @{ Name="same_PBC_6835_to_6837"; Target="PBC_PR_6837"; Sources=@("PBC_PR_6835_5A"); TaskType="same" },

    # --- 4 组严格组织级留出方向（10 个测试目标切片） ---
    @{ Name="cross_holdout_LN_A1_LN"; Target="A1_LN"; Sources=@("A1_TNSL","D1_TNSL","HumanTonsil_10x","HumanTonsil_AddOns_S1","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_LN_D1_LN"; Target="D1_LN"; Sources=@("A1_TNSL","D1_TNSL","HumanTonsil_10x","HumanTonsil_AddOns_S1","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_Tonsil_A1_TNSL"; Target="A1_TNSL"; Sources=@("A1_LN","D1_LN","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_Tonsil_D1_TNSL"; Target="D1_TNSL"; Sources=@("A1_LN","D1_LN","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_Tonsil_HumanTonsil_10x"; Target="HumanTonsil_10x"; Sources=@("A1_LN","D1_LN","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_Tonsil_HumanTonsil_AddOns_S1"; Target="HumanTonsil_AddOns_S1"; Sources=@("A1_LN","D1_LN","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_Breast_Cancer_Breast_Cancer"; Target="Breast_Cancer"; Sources=@("A1_LN","D1_LN","A1_TNSL","D1_TNSL","HumanTonsil_10x","HumanTonsil_AddOns_S1","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_Breast_Cancer_PBC_PR_6837"; Target="PBC_PR_6837"; Sources=@("A1_LN","D1_LN","A1_TNSL","D1_TNSL","HumanTonsil_10x","HumanTonsil_AddOns_S1","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_Breast_Cancer_PBC_PR_6835_5A"; Target="PBC_PR_6835_5A"; Sources=@("A1_LN","D1_LN","A1_TNSL","D1_TNSL","HumanTonsil_10x","HumanTonsil_AddOns_S1","GBM"); TaskType="cross" },
    @{ Name="cross_holdout_GBM"; Target="GBM"; Sources=@("A1_LN","D1_LN","A1_TNSL","D1_TNSL","HumanTonsil_10x","HumanTonsil_AddOns_S1","Breast_Cancer","PBC_PR_6837","PBC_PR_6835_5A"); TaskType="cross" }
)

$total = $jobs.Count
$idx = 0

foreach ($job in $jobs) {
    $idx++
    $out = Join-Path $baseRunRoot $job.Name
    if (Test-Path (Join-Path $out "results.json")) {
        Write-Output "[$idx/$total] $($job.Name) already finished, skipping..."
        continue
    }
    if (Test-Path $out) { Remove-Item -Recurse -Force $out }
    $log = Join-Path $baseLogRoot ($job.Name + ".log")
    Write-Output "[$idx/$total] Starting: $($job.Name)..."
    $args = @("-u", $script,
        "--prepared_root", $root, "--out_dir", $out, "--experiment_dir", $experiment,
        "--target", $job.Target, "--source_sections") + $job.Sources + @(
        "--task_type", $job.TaskType,
        "--seed","42","--epochs","20","--batch_size","64","--validation_ratio","0.1",
        "--hidden_dim","256","--heads","4","--ffn_dim","256","--dropout","0.3",
        "--protein_embedding_dim","32","--router_rank","32","--lr","0.0005","--weight_decay","0.0001",
        "--aux_weight","0.25","--pcc_weight","0.05","--pcc_min_variance","1e-6",
        "--context_pca_components","32","--context_rank","8","--gene_dropout_rate","0.5","--neighbor_mask_rate","0.25",
        "--local_neighbors","4","--regional_neighbors","8","--program_components","9")
    & $python @args 2>&1 | Out-File -FilePath $log -Encoding utf8
    Write-Output "[$idx/$total] Finished: $($job.Name)"
}
Write-Output "=== All 18 targets trained with Domain Scale Adaptive Selection! ==="
