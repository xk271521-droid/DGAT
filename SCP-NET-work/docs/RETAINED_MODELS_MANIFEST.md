# Retained models and result manifest

Updated: 2026-10-31

## Retained improvement code

- `experiments/cognate_anchor/`
- `experiments/protein_convex/`

## Retained corresponding results

- `runs/strict_cognate_anchor_<target>/`
- `runs/protein_convex_<target>/`

The original baseline run directories are also retained for comparison. Temporary candidate branches, pilot/smoke runs, noncanonical ProteinConvex runs, ensembles, audit runs, checkpoints, prediction arrays, logs, and dated documents belonging only to those discarded versions were removed from `SCP-NET-work`.

`ST2P/` and its Git metadata were not touched by this cleanup.

The former MultiViewConsistency implementation and its full LODO run artifacts
were removed on 2026-10-31 at the user's request. Historical result documents
and summary snapshots that mention MultiViewConsistency are retained for audit
history and are not part of the active retained-model scope.
