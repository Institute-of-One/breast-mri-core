# Q3 exploratory analysis lock — 2026-09-18

Keep the pre-Q2 88/29 patient split. Full Q2 targets have been exposed; this is
exploratory internal evaluation, never a fresh confirmatory holdout.
Primary target: T1 real/sub retention. T0 secondary, separate model with same
specification; no T0 target as a T1 predictor. Undefined targets are recorded.

Input: precontrast and early images only, within the provided unilateral crop.
Do not read late arrays or the analysis mask for feature extraction. Native
Philips precision corrections follow the already fixed acquisition policy.
Crop localization itself used full-examination knowledge. Thus this is a
crop-conditioned explanatory study, not a validated early-only deployment path.
It cannot close the independent-localization gate in the earlier protocol.

Fixed early support: S0 > 10% of positive-S0 p95; early enhancement >70%.
Features: PE median, PE IQR, PE p90 on early support; fraction of background-
eligible voxels meeting early threshold; median PE divided by early-minus-pre
acquisition interval in minutes; that interval; HR/HER2 categorical group.
The rate is an average enhancement rate, not ultrafast initial kinetics.
No PE threshold tuning, outcome-based voxel selection or late timing predictor.
Empty support is missing; imputation learned in development only.

Comparators: development-mean constant; subtype-only ridge; early-image ridge;
early-image-plus-subtype ridge. Ridge alpha=10 fixed, numeric standardization
and median imputation fitted in development; categorical one-hot unknown-safe.
No hyperparameter or feature search. Do not clip predictions to [0,1]; record
out-of-range predictions. Freeze fitted parameters and test predictions before
reading test targets in the evaluation function.

Primary performance: MAE on retention (report percentage points), paired MAE
difference against constant and subtype baseline. Secondary RMSE and R2.
Patient bootstrap 10,000 seed 20260918 for paired MAE differences. Report all
four comparators and both visits without selecting a successful one after testing.
No success guarantee from Q2. Require independent data before clinical claims.

Q1 diagnostics run separately on seven tolerance-failing exams. Compare native
versus published intensities, float precision, rounded/exact and inclusive/strict
SER thresholds. Record support-count differences separately. Candidate matches
do not establish VOI/connectivity/rounding as the causal mechanism. Do not change
primary volumes or tolerance definitions based on these diagnostics.
