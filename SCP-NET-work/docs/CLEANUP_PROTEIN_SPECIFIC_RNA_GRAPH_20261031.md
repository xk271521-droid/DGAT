# Cleanup: protein-specific RNA graph experiment (2026-10-31)

The protein-specific supervised RNA graph experiment was removed from the active workspace after comparison with the stricter generic RNA calibration.

- Deleted: `SCP-NET-work/experiments/protein_specific_rna_graph_20261031/`
- Retained: `RESULT_PROTEIN_SPECIFIC_RNA_GRAPH_7FOLD_20261031.md` as an audit record
- Retained as the selected method: `ProteinConvex + Strict RNA Calibration`

The deleted candidate improved the original ProteinConvex baseline but was lower than Strict RNA Calibration on mean PCC-Pro and PCC-Spot.
