# Theme 4: FTV reproducibility and the contribution of late-phase SER

2026-09-17. User-authorized pivot. This is a timestamped local analysis plan,
not a public preregistration. Existing Theme 1 and Theme 3 results remain frozen.

## Interpretation amendments

- Theme 1: the 73-patient comparison is inconclusive about superiority; it does
  not establish equivalence or refute all selective-information policies. Its CI
  concerns selection versus random selection, not all-DCE versus no-additional-DCE.
  The latter's observed 0.008 improvement cannot establish power for the former.
  Neither 'inherently unidentifiable' nor a required sample size follows from this CI alone.
- Theme 3: the current ridge model saturated the threshold in six validation
  examinations. Squared-error regression may smooth kinetics, but shrinkage
  toward early signal and saturation are not mathematical necessities for all
  voxel regressions. For S2=S1 and S1!=S0, SER=1 exactly; this is the relevant
  deterministic comparator. Denominator reduction only applies when S2>S1>S0.
- Eight of 18 strict reproductions passed. This is a local implementation/data
  discrepancy, not proof of an erroneous public reference or clinical unreliability.
  Q1 promotes its investigation to a primary research question, prospectively
  for expanded analysis; the already observed pilot is explicitly exploratory.

## Cohort and sequencing

Initial full target: the original 117 training patients, T0/T1, 234 examinations.
Inventory every examination, including missing/unprocessable data. No automatic
addition of the consumed 73-patient pCR holdout. The nine engineering patients
and both visits are exposed pilot data, not an independent validation set.
Stage 1 uses the existing 18 complete image/mask pairs and inventories all 234
cached reference headers. Expansion needs a defensible common pipeline and
explicit records of acquisition, exclusions and unresolved conventions.

## Q1: source-to-implementation discrepancy audit

Distinguish FTVpe and FTVser; retain exact stored PE/SER settings and units.
PE=70% and SER=0.9 are standardized sensitivity settings, never assumed source
settings for all patients. Report signed/absolute voxel and cc differences,
relative differences with zero-reference flags, and exact-match status separately.
One to seven voxels near a boundary is not equivalent to a large volume discrepancy.

For each exam record PE/SER precision and strict/inclusive threshold effects;
VOI/omit boundary mapping; background intensity threshold and coordinate domain;
connectivity ordering/neighbour convention; original/rescaled intensity precision;
voxel geometry/resampling; unexpected mask bits; and unresolved residuals.
Use documented settings and one-factor perturbations, not the combination that
best matches the stored volume. Count matches alone cannot identify causation.
Evidence levels: demonstrated by source/voxel correspondence, candidate, unresolved.
Interacting causes must not be forced into an additive mutually exclusive allocation.
Never reinterpret undocumented bit 16 as background solely to improve agreement.

## Q2: paired late-phase ablation

Full-pipeline estimand: process measured (S0,S1,S2) and (S0,S1,S1) identically,
recomputing any late-dependent steps, with acquisition/VOI rules recorded.
Until this is possible, a fixed published mask == 0 analysis is labelled
'reference-mask-conditioned SER-volume', not fully reproduced FTV or an early-only
clinical protocol. Audit whether even that mask/VOI contains late-phase dependence.

Primary threshold: |Vsub-Vreal|/Vreal >= 0.10, fixed now and never optimized.
Report proportions separately at T0 and T1 (patient is sampling unit), plus the
prespecified patient-level either-visit proportion. If Vreal=0 and Vsub>0 classify
as discordant zero-volume, report separately and include as exceeding threshold;
if both zero classify as unchanged. Never replace denominators with epsilon.
Also report (Vsub-Vreal)/Vsub and voxelwise XOR fraction relative to eligible
support; different denominators answer different questions and must be named.

Report signed/absolute cc differences; signed change D=V_T1-V_T0 and difference
Dsub-Dreal; percentage-change difference only when both baselines are positive.
Bland-Altman plots and bias/limits are separate for T0, T1 and longitudinal change.
Wilson 95% intervals for proportions; patient bootstrap (10,000, seed 20260917)
for bias, limits and longitudinal differences. Pilot intervals describe selected
nine-patient data only, not population generalizability. Expanded analysis must
retain all exclusions and report heteroscedasticity and extreme/zero volumes.

## Q3: exploratory prediction, deferred

After Q2, freeze the patient split, predictors and model before fitting Q3.
Because Q2 has exposed targets, this is exploratory internal validation, not
confirmatory preregistration. All nine exposed pilots go to development; all
visits from one patient stay together. New independent patients are required for
confirmatory claims. Do not optimize the 10% definition using Q2's distribution.
Use only early-accessible predictors: PE distribution, time-normalized initial
enhancement (not an ultrafast kinetic slope unless temporal sampling supports it),
and subtype known at the intended decision time. Full-DCE-derived masks/VOIs are
not admissible hidden inputs. No synthesis training or selective-pCR reanalysis.

## Reporting and publication decision

Commit to report either direction and retain null results. Analysis completion
requires Q1 evidence/residual accounting and Q2 patient-level estimates with CIs.
This is NOT automatic proof of novelty, clinical relevance, or paper acceptance.
Small SER contribution supports metric-specific ablation feasibility, not the
safety of omitting a late phase. Large contribution alone is not reproducibility
failure; it may be expected sensitivity to removing real biological information.
Publication strength requires a distinct contribution beyond prior work, enough
precision, defensible measurement, and preferably another collection/implementation.

## Primary literature checked

- Jafri et al., JMRI 2014, doi:10.1002/jmri.24351:
  https://pubmed.ncbi.nlm.nih.gov/24347097/ . Robustness concerns prognostic value
  across thresholds, not invariance of every FTV value.
- Newitt et al., Aegis validation 2014:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3998689/ . Implementation differences
  coexist with similar prognostic performance; these reports are not contradictory.
- Longitudinal variation and pCR prediction, 2023:
  https://pubs.rsna.org/doi/10.1148/rycan.220126 . Further full-text comparison needed.

These sources motivate separating prognostic robustness, numerical reproducibility,
and late-phase information contribution rather than claiming an unresolved conflict.
