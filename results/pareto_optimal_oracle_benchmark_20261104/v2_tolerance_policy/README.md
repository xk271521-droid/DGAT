# Oracle Selection Policy V2
Date: 2026-11-04
Strict dominance: PCC >= max competitor PCC and RMSE <= min competitor RMSE.
PCC-priority: PCC margin >= 0.02 and RMSE gap <= 0.01.
Calibration required: no epoch meets either rule.
This is a development-stage oracle label derived from full target trajectories.
Cross counts: {'strict_dominance': 6, 'pcc_priority_small_rmse_regret': 1, 'calibration_required': 3}
Same counts: {'strict_dominance': 7, 'pcc_priority_small_rmse_regret': 0, 'calibration_required': 1}
