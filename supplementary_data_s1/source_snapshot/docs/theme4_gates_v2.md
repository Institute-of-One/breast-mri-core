# Theme 4 full-cohort gates, v2 — 2026-09-17

User-authorized 117-patient / 234-examination expansion. This amendment precedes
full-cohort Q2 aggregation, but follows the exposed nine-patient pilot. Preserve
v1 and its hashes. No manuscript drafting before the complete cohort flow is known.

Q1 exact agreement is unchanged. Additional engineering tolerance requires BOTH
absolute volume difference <=0.01 cc AND relative difference <=0.1% of stored
volume, separately for FTVpe and FTVser. A visit passes only if both endpoints
pass. Stored zero requires reconstructed zero. These are explicit engineering
limits, not clinically validated equivalence margins and not selected by searching
for 17/18. Report continuous differences even when a tolerance passes. The pilot
has already been viewed; the new margin is prospective only for added cases.

Focus source investigation on the 92-SER-voxel / 344-PE-voxel discrepancy.
Retain the other nine as small unexplained differences; rounding, threshold and
connectivity are candidate mechanisms, not proven causes. New substantial errors
in the expanded cohort require reporting and technical review, not automatic exclusion.
17/18 pilot tolerance passes, if observed, cannot certify an entire independent
FTV implementation. Full Q2 remains explicitly public-reference-mask-conditioned.

Q2 primary data use actual native intensities for Philips when exact physical-grid
mapping is possible; existing native-precision pilot corrections remain unchanged.
Other manufacturers use published DCE crops with source and scaling recorded.
No guessed interpolation or undocumented bit reinterpretation. If native mapping
fails, retain an unresolved examination; do not silently substitute lower-precision
data in the primary result. All 234 planned examinations appear in the flow report.
This is a paired S2-ablation study conditional on the supplied ROI/mask; it does not
reconstruct the full VOI/background/connectivity pipeline from early images alone.

Keep the v1 >=10% absolute relative difference endpoint. Additionally report signed
overestimation >=10%, underestimation >=10%, zero references and discordant zeros;
do not change endpoint because pilot direction is attractive. Patient-bootstrap
intervals and visit-stratified proportions; account for exclusions with worst-case
bounds. Large contribution shows S1 is an inadequate substitute for this metric,
not indispensability of late MRI for all clinical purposes.

Q3 target is retention R=FTVreal/FTVS1sub, with excluded fraction 1-R separately
named. T1 is the primary patient-level target; T0 secondary. Zero denominator is
undefined, not set to zero. No forced clipping to [0,1]. Freeze a manufacturer-
stratified 75/25 patient split now using SHA256 seed theme4-q3-20260917-v2;
all nine exposed pilot patients belong to development. Keep all visits together.
The frozen split precedes full Q2 aggregation, but Q3 remains exploratory: once
full Q2 is viewed, its test targets are no longer fully unobserved. No model tuning
or feature selection on test targets. Early-only inputs and independent validation
requirements from v1 remain. Do not start fitting Q3 in this expansion.

Synthesis and selective-pCR analyses remain frozen. Residual sign analysis, if
performed later, uses frozen predictions and is a mechanism diagnostic, not retraining.
