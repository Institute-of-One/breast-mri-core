# Prespecified validation plan, ACRIN-6698 independent subset (v1)

Written 2026-11-11, **before any validation volume, score or AUC exists**, in
answer to Referee 1 point 2 of BPEX-107031: "validate the main findings in an
independent cohort using a prespecified analysis plan."

This document is frozen by sha256 before the acquisition script is run. The
recorded hash, the patient list hash and the seed are the evidence that the
plan preceded the result. Nothing in it may be changed after the first
validation AUC is computed; if it is changed, the change and its date are
appended here and the whole analysis is rerun and reported as a second run.

## 1. What is being validated, and what is not

The published analysis (117 BMMR2 training patients) reports that replacing the
late image with the early image inflates the SER-filtered volume substantially
while changing pCR discrimination little, because the substitution reorders few
pCR/non-pCR pairs and the reorderings largely cancel. The emphasis on relative
volume reduction was chosen after the AUC results were inspected. That is the
weakness Referee 1 names, and it is the thing this plan tests: **does the
finding hold on patients whose results were never seen, under an endpoint fixed
in advance?**

**This is not external validation and will not be described as one.** The
validation cohort comes from the same trial (ACRIN-6698 / I-SPY 2), the same
collection, the same acquisition protocol and the same public mask provenance.
It is an independent, non-overlapping subset of the same collection. It cannot
answer transfer to another institution, scanner fleet or protocol.

## 2. Cohort

From `data/raw/acrin6698/Full-Collection-Ancillary-Patient-Information-file.xlsx`,
sheet `ACRIN_6698_Patient_Cohorts_Clin` (385 patients), a patient enters the
validation cohort if and only if **all** of:

1. `T0` = 1 and `T1` = 1;
2. `BMMR2_TRAIN` ≠ 1 (disjoint from the 117 already analysed);
3. `pcr` is exactly `pCR` or `Non-pCR` (not blank, not `NA`);
4. `QC_PASS_FAIL_T0` = `PASS` **and** `QC_PASS_FAIL_T1` = `PASS`;
5. `analysis cohort` = 1;
6. `eligigible` = 1 (the file's spelling).

Criteria 4 to 6 are not a choice made for this validation: all 117 patients of
the published cohort satisfy them, so the rule reproduces the discovery
cohort's own definition rather than inventing one.

**Resulting cohort: 110 patients, 35 pCR and 75 non-pCR, 2625 pCR/non-pCR
pairs** (the published cohort has 117, 36/81, 2916 pairs). 74 of the 110 are
the BMMR2 testing cases whose labels were withheld during the challenge and are
released in this ancillary file; the other 36 are outside both BMMR2 sets.

Composition, recorded now so that any later drift is visible: HR+/HER2− 44,
triple-negative 38, HR+/HER2+ 18, HR−/HER2+ 10; SBR grade III 82, II 26, I 2;
age mean 47.8, SD 9.4, range 25–72.

The patient identifier list is frozen as `docs/theme4_validation_cohort_v1.json`
and hashed in the lock file. No patient may be added or removed afterwards for
any reason connected to its images, volumes, scores or outcome.

### 2.1 Loss after freezing

A patient that cannot be processed is **excluded and counted, never replaced**.
The permitted reasons, fixed here, are exactly:

* TCIA carries no `Analysis Mask` series, or no `original DCE` series, for a
  required visit;
* the acquisition fails after three attempts, or the archive returns a file
  that does not open;
* a geometry, frame-of-reference or patient/visit correspondence check from the
  published pipeline fails;
* the eligible support is empty at a required visit.

Every exclusion is listed with its patient identifier and reason in the result
file. The analysis is reported on the patients that complete, with the attempted
and completed counts stated in the same sentence as the result.

## 3. Procedure

The pipeline is the published one, unchanged: `scripts/theme4_full_cohort.py`
for acquisition and extraction, `scripts/theme3_ftv.py` for the volume
computation, `scripts/theme3_restore_philips.py` for Philips native-grid
recovery. Thresholds, the stored SER value and strict inequality, the
nonfinite-denominator and capping rules, and the engineering tolerance are
those already in the code. **No parameter is re-tuned for this cohort.** The
only new inputs are the patient list and the bootstrap seed.

**Outcome labels are joined last.** Acquisition, extraction and both volume
arms complete for the whole cohort before any pCR label is read into the
analysis. Volumes do not need labels, so this costs nothing and removes the
possibility of a choice made while looking at an outcome.

## 4. Endpoints, scores and estimand

Identical to the published analysis.

* **Primary endpoint: relative volume reduction**, (V_T0 − V_T1) / V_T0,
  computed separately in each arm. It is the primary endpoint *because* it is
  the emphasis the published analysis chose after seeing its results; fixing it
  here in advance is the point of this exercise.
* Secondary: T0 volume and T1 volume.
* Also computed and reported, not interpreted: absolute change V_T0 − V_T1.

Scores: negative volume at T0 and T1, positive reduction for the longitudinal
measures. **Directions are fixed here and will not be reversed**, whatever AUC
results, exactly as in the published analysis.

Estimand: AUC_sub − AUC_real on the same patients. Ties take half credit.

## 5. Interval estimation

10,000 paired, outcome-stratified patient bootstrap resamples; pCR and non-pCR
patients resampled independently within the cohort; the same sampled
identifiers applied to both volume definitions; percentile intervals for each
AUC and for the paired difference.

**Seed: 20261111.** Declared here, used once.

Intervals are pointwise. No multiplicity adjustment is made and none is
claimed; the primary endpoint is single and named above.

## 6. Reproduction criteria, fixed in advance

Published (discovery) values, from `data/derived/theme4_mechanism_v1/rank_mechanism.json`
and the archived AUC results:

| Endpoint | AUC real | AUC sub | Δ | 95% CI | pairs reordered | Spearman |
|---|---:|---:|---:|---|---:|---:|
| Relative reduction (primary) | 0.6461 | 0.6529 | +0.0069 | −0.043 to 0.056 | 350 / 2916 = 12.00% | 0.9211 |
| T0 volume | 0.6591 | 0.6722 | +0.0130 | −0.026 to 0.051 | 246 / 2916 = 8.44% | 0.9459 |
| T1 volume | 0.7126 | 0.7030 | −0.0096 | −0.068 to 0.043 | 356 / 2916 = 12.21% | 0.9139 |

Each criterion below is judged pass or fail **separately** and reported
separately. They are not combined into one verdict, and a failure of any of
them is reported as a failure, not re-analysed.

* **C1 (primary).** The validation Δ for relative reduction lies within the
  discovery 95% CI, −0.043 to 0.056.
* **C2.** |Δ| ≤ 0.05 for all three endpoints. This does not use the discovery
  intervals, so it cannot be satisfied merely by the discovery being imprecise.
* **C3.** The fraction of pCR/non-pCR pairs whose order changes is within
  ±5 percentage points of the discovery fraction for the same endpoint.
* **C4.** Substitution increases volume by at least 10% in at least 90% of
  patients at both visits (discovery: 99.1% at T0, 98.3% at T1).
* **C5.** The Spearman correlation between the measured-late and substituted
  scores is at least 0.85 at all three endpoints.

**What would show the published claim does not hold here.** Failure of C1 or of
C2 means the central claim — that this substitution changes discrimination
little — does not reproduce on patients whose results were not seen. Failure of
C4 would mean the volume inflation itself does not reproduce, which would make
the cohorts not comparable. Failure of C3 or C5 with C1 and C2 passing would
mean the AUC agreement is reproduced but the mechanism offered for it is not.

Each of these will be written in the manuscript in those words if it happens.
No endpoint, threshold, cohort rule or seed will be changed in response to the
result.

## 7. Secondary analyses, prespecified

* Agreement with the stored volumes under the same fixed tolerance
  (|Δ| ≤ 0.01 cc **and** relative ≤ 0.1%), reported as counts.
* Sensitivity analysis repeating the primary estimate with patients outside
  that tolerance excluded, as the discovery cohort's revision will also report.
* The pair decomposition: pairs changed, gained, lost, net concordance, and the
  Spearman correlation, for each endpoint.
* Volume statistics: median volumes, retained fraction quartiles, the standard
  deviations of the logarithms and their ratio, and the correlation between
  log substituted volume and log retained fraction.
* Stratification by site and by scanner manufacturer, reported descriptively.
* Subtype estimates: computed and reported **as exploratory**, with the number
  of pairs beside every one. No conclusion about any subtype will be drawn from
  this cohort either.

## 8. Reporting

Whatever the result, the manuscript will report: the attempted and completed
counts with every exclusion and its reason; all three Δ with intervals; each of
C1 to C5 as pass or fail; the pair decomposition; and the volume statistics.
The discovery and validation cohorts are reported side by side, and the
validation is named as an independent subset of the same collection, never as
external validation.

## 8.1 Amendment 1 — 2026-11-11, before any image was acquired

Made after the metadata census and **before any volume, score or AUC exists**.
The superseded freeze is recorded here so the change is auditable:
plan sha256 `a0f76fbadf34ec5179e05d23f542c29f6d09e67e60915edf53217ab25b38fce5`,
cohort sha256 `037eaaa1e82fde95b732de96c8f90ef06a40f1b4a035ff9cc609120b843df55e`.

**(a) Cohort losses from the archive.** The census of all 110 frozen patients
against TCIA (`scripts/theme4_validation_inventory.py`, 2026-11-11, metadata
only, no outcome read) found that two patients do not carry both of the two
series the pipeline requires at both visits:

* `ACRIN-6698-547405` — one Analysis Mask and one original DCE in the whole
  record, both under the T2 study; nothing at T0 or T1.
* `ACRIN-6698-786753` — Analysis Mask and original DCE at T0 and T2, none at T1.

Both fall under the first reason listed in section 2.1 and are excluded and
counted, not replaced. **The analysed validation cohort is 108 patients.**
These exclusions were decided from the archive's holdings alone; no pCR label
was read to make them, and none has been read at the time of writing.

**(b) Manufacturer strings decide a processing path, so they are normalised.**
The frozen pipeline runs the Philips native-grid intensity recovery only when
the record's manufacturer is exactly `Philips Medical Systems`. The validation
cohort contains 13 patients with that string and **3 with `Philips Healthcare`**,
a value that does not occur anywhere in the published cohort. Left alone, three
Philips patients would silently take the published-crop path while the other
thirteen took the native-grid path, with a free-text DICOM field, not the
vendor, deciding the processing.

The record's manufacturer is therefore normalised by this fixed map before it
reaches the pipeline, and by no other:

```
"Philips Healthcare"     -> "Philips Medical Systems"
"GE MEDICAL SYSTEMS"     -> unchanged
"SIEMENS"                -> unchanged
"Philips Medical Systems"-> unchanged
```

Any manufacturer string outside this map stops the run rather than defaulting
to a path. The map is applied to the record only; the DICOM field is never
rewritten, and both the raw and the normalised strings are reported.

Scanner mix of the 108, from the census, recorded now: GE 71, Siemens 21,
Philips 16 (13 + 3). The published cohort is GE 93, Siemens 13, Philips 11, so
the validation carries proportionally more non-GE examinations than the cohort
it validates. That is reported as a difference between the cohorts, not as a
designed stratification.

## 9. Freeze

```
python scripts/freeze_validation_plan.py           # writes docs/theme4_validation_lock_v1.json
python scripts/freeze_validation_plan.py --check   # refuses if the plan or the list changed
```

The lock records the sha256 of this document, of the frozen patient list, of
the acquisition and volume scripts, and the seed. The acquisition script is not
run until `--check` passes.

## 8.2 Amendment 2 — 2026-11-11, after acquisition and **before any outcome label was read**

The superseded freeze is recorded here so the change is auditable:
plan sha256 `791ce938de6bf70de888441d517e92abf79e718ed32d003d5b921ff51f391efe`,
cohort sha256 `037eaaa1e82fde95b732de96c8f90ef06a40f1b4a035ff9cc609120b843df55e`. The cohort list itself is **unchanged** by
this amendment — only the plan text grew — and the superseded lock is kept
beside the new one as `docs/theme4_validation_lock_v1.superseded-2.json`.

Acquisition completed: 216 examinations attempted, 216 reports, **196 complete**.
Twenty examinations in eleven patients failed, every one of them on the same
check inside the published pipeline:

    assert len(candidates)==1, 'Native series selection not unique'

`theme4_full_cohort.process` identifies the Philips native series as the unique
non-VOLSER MR series at that visit whose image count equals the published
crop's. In this cohort that heuristic finds **zero** candidates, not two: the
native acquisition carries more slices than the uni-lateral crop, so the counts
differ by 2 to 25 per cent. It is a correspondence check from the published
pipeline in the sense of section 2.1, third bullet, and the reading is stated
here rather than assumed.

### The eleven, and the vendor consequence

* `ACRIN-6698-139973` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-169918` (Philips Medical Systems) — T0
* `ACRIN-6698-222135` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-264305` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-343898` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-375954` (Philips Medical Systems) — T0
* `ACRIN-6698-422013` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-689076` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-689452` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-805136` (Philips Medical Systems) — T0, T1
* `ACRIN-6698-973837` (Philips Medical Systems) — T0, T1

All eleven are Philips, which is **11 of the 16 Philips patients**; no GE
(0 of 71) or Siemens (0 of 21) patient was affected. Nine lose both visits and
two lose one; the design is paired, so all eleven are excluded. **97 of the 110
frozen patients complete both visits**, and the analysed cohort is 71 GE,
21 Siemens and 5 Philips.

### Why no replacement rule was written

A hand-written series-selection rule could recover most of them, and it is
**deliberately not written**:

* the native series names are site-specific and share no pattern — `ACRIN DYN`,
  `T1 FS DYN SENSE 5`, `DYNAMIC T1 FS`, `POST GAD DYNAMIC`,
  `Dynamics inject on 2`, `dyn eTHRIVE SENSE`, `VISTA`,
  `T1 FS POST SENSE 5` — so any rule would be a per-site heuristic chosen by
  the author;
* in `ACRIN-6698-973837` the dynamic series is split into separate pre- and
  post-contrast series, so no single native series exists to select;
* section 3 requires the published pipeline **unchanged**. Substituting the
  author's own selection rule for the published one would retire the only
  claim this exercise rests on.

The loss is therefore accepted, counted and reported, and the vendor imbalance
it creates is reported in the same sentence as the result.

### Blinding

This decision was taken from the TCIA series listings and the pipeline's own
assertion alone. **No pCR label had been read when it was made**; the
acquisition summary records `"scope": "volumes only; no outcome label has been
read"`, and the analysis script reads the labels for the first time after this
amendment was written. No endpoint, threshold, cohort rule or seed changed.
