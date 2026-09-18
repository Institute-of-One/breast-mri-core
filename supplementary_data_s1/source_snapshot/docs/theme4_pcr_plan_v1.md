# Paired pCR discrimination analysis — 2026-09-18

User-authorized post-hoc endpoint pivot after inspecting Q2 and Q3. Local lock
before computing these AUCs; not prospective registration. Original 117 training
patients and labels are previously used data, not an independent validation set.

Compare measured-S2 and S1-substituted reference-mask-conditioned SER volumes.
Use native-precision v2 outputs, without removing the seven Q1 discrepancies.
Primary contrast: absolute change T1 minus T0, AUC difference sub minus real.
Also report T0 and T1 separately and relative change as secondary. Scores are
negative volume at each time and negative signed change (greater reduction means
higher pCR score). Fix direction now; never flip an AUC below 0.5 after inspection.
Relative change score=(T0-T1)/T0; require both methods' T0 volumes >0 and report
paired exclusions. Zero follow-up volumes remain valid. No fitting or covariate
adjustment; these are single-biomarker discrimination AUCs, not calibrated models.

pCR positive uses the existing training CSV pcr==pCR; require exact patient-ID
join, unique IDs, Train split and expected binary categories. No imputation of
missing labels. Subtype strata follow all four source hrher4g categories unchanged.

10,000 paired patient bootstrap replicates, stratified by pCR within each analysis
stratum, seed 20260918. Identical resampled patients for both scores. Percentile
95% intervals for each AUC and sub-real difference. Ties receive half credit.
If either outcome class is absent, AUC/interval unavailable, not 0.5. Report class
counts, flag fewer than 10 of either class, and interpret small strata cautiously.
Intervals are pointwise/exploratory, not multiplicity-adjusted. No result-driven
choice of subtype, time or change definition. No noninferiority margin invented
after results: a CI including zero does not establish preservation/equivalence.

Q1 follow-up stops per user instruction. Mask version/resampling remain possible
explanations, not confirmed causes. Q3 stays frozen and is reported as no clear
improvement by the tested simple features/model, not universal unpredictability.
Q2/BA are secondary measurement comparisons; BA does not validate a clinical
abbreviated protocol. Public localization/masks may incorporate full-DCE knowledge.
Thus this is retrospective late-phase ablation conditional on those masks, not
proof that an autonomous pre-plus-one-early-phase protocol retains clinical utility.

Literature: Kuhl 2014 (PMID24958821) is a screening reader study, not a general
claim of early-phase sufficiency for treatment-response prediction. Li et al.
(PMC5214452) examined PE/SER thresholds and pCR by subtype in ACRIN6657/I-SPY1.
Paired removal/longitudinal framing needs comparison against its full methods
before claiming novelty; threshold removal can itself be a threshold-sweep limit.
