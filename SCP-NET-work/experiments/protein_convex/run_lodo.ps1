$ErrorActionPreference = 'Stop'
$python = 'C:\Users\xk2714\.conda\envs\st2p\python.exe'
$root = 'E:\DGAT\ST2P\data\processed\DGAT_final7_breast_ptprc31_muon017_axis0_20261022'
$experiment = 'E:\DGAT\ST2P\src'
$script = 'E:\DGAT\SCP-NET-work\experiments\protein_convex\scn_pro_protein_convex.py'
$runRoot = 'E:\DGAT\SCP-NET-work\runs'
$targets = @('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer')
$all = @('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer')

foreach ($target in $targets) {
    $out = Join-Path $runRoot ('protein_convex_' + $target)
    $sources = $all | Where-Object { $_ -ne $target }
    $args = @('-u', $script,
        '--prepared_root', $root,
        '--out_dir', $out,
        '--experiment_dir', $experiment,
        '--target', $target,
        '--source_sections') + $sources + @(
        '--seed','42','--epochs','20','--batch_size','64','--validation_ratio','0.1',
        '--hidden_dim','256','--heads','4','--ffn_dim','256','--dropout','0.3',
        '--protein_embedding_dim','32','--router_rank','32','--lr','0.0005','--weight_decay','0.0001',
        '--aux_weight','0.25','--pcc_weight','0.05','--pcc_min_variance','1e-6',
        '--context_pca_components','32','--context_rank','8','--gene_dropout_rate','0.5','--neighbor_mask_rate','0.25',
        '--local_neighbors','4','--regional_neighbors','8','--program_components','9')
    $log = Join-Path $runRoot ('protein_convex_' + $target + '.log')
    & $python @args 2>&1 | Tee-Object -FilePath $log
    if ($LASTEXITCODE -ne 0) { throw "strict LODO failed for $target with exit code $LASTEXITCODE" }
}
