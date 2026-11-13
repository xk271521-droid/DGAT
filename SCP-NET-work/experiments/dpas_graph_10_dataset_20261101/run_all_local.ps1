$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..')).Path
$python = Join-Path $root '.venv-dpas\Scripts\python.exe'
$script = Join-Path $PSScriptRoot 'run_dpas_graph_10_dataset.py'
$out = Join-Path $root 'SCP-NET-work\runs\dpas_graph_10_dataset_20261101'
New-Item -ItemType Directory -Force -Path $out | Out-Null
& $python $script --out_root $out --epochs 100 --seed 0
