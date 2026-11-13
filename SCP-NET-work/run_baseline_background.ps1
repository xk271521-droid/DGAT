$ErrorActionPreference = 'Stop'
$conda = 'D:\miniconda\condabin\conda.bat'
$prepared = 'E:\DGAT\ST2P\data\processed\DGAT_final7_breast_ptprc31_muon017_axis0_20261022'
$script = 'E:\DGAT\ST2P\src\scn_pro.py'
$experiment = 'E:\DGAT\ST2P\src'
$runroot = 'E:\DGAT\SCP-NET-work\runs'
$targets = @('A1_LN','D1_LN','A1_TNSL','D1_TNSL','HumanTonsil_10x','HumanTonsil_AddOns_S1','Breast_Cancer')

New-Item -ItemType Directory -Force -Path $runroot | Out-Null
foreach ($target in $targets) {
    $sources = @($targets | Where-Object { $_ -ne $target })
    $out = Join-Path $runroot ('baseline_' + $target)
    if (Test-Path (Join-Path $out 'results.json')) {
        Add-Content (Join-Path $runroot 'baseline_background.log') ("SKIP $target already complete $(Get-Date -Format o)")
        continue
    }
    if (Test-Path $out) {
        $out = Join-Path $runroot ('baseline_' + $target + '_retry')
    }
    $log = Join-Path $runroot ('baseline_' + $target + '_background.log')
    Add-Content $log ("START $target $(Get-Date -Format o)")
    & $conda run --no-capture-output -n st2p python $script `
        --prepared_root $prepared --out_dir $out --experiment_dir $experiment `
        --target $target --source_sections $sources --seed 42 --epochs 20 --pcc_weight 0.05 `
        *>> $log
    $code = $LASTEXITCODE
    Add-Content $log ("END $target exit=$code $(Get-Date -Format o)")
    if ($code -ne 0) { exit $code }
}
Add-Content (Join-Path $runroot 'baseline_background.log') ("ALL COMPLETE $(Get-Date -Format o)")
