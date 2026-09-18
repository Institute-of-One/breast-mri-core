# breast-mri-core

Analysis code and derived scores for a study of what happens to breast MRI functional
tumor volume, and to its discrimination of pathological complete response (pCR) after
neoadjuvant therapy, when the late post-contrast phase is removed.

The study uses the public ACRIN-6698 / I-SPY 2 collection of The Cancer Imaging Archive:
the 117 participants of its BMMR2 training subset, each imaged before treatment (T0) and
about three weeks after it began (T1). Within the fixed public analysis masks, volumes
filtered by signal enhancement ratio (SER) are computed twice -- once from the acquired
precontrast, early and late images, and once with the late image replaced by the early
one -- and the two versions are compared for pCR discrimination with paired,
outcome-stratified patient bootstrap intervals.

This repository accompanies a manuscript under consideration.

## Contents

`supplementary_data_s1/` is, byte for byte, the manuscript's Supplementary Data S1. Its
`sha256.json` lists the checksum of every file in it; the paths there were written on
Windows and use backslashes as separators.

| Path | What it is |
|---|---|
| `paired_scores.csv` | The 117 public participant IDs, pCR and subtype labels, and the fixed scores used in the reported analysis |
| `reproduce_auc.py` | Recomputes the 20 paired comparisons from `paired_scores.csv` and asserts agreement with the frozen results |
| `auc_results.json`, `verification.json` | The frozen results and the record of their verification |
| `source_snapshot/` | The image-processing and analysis scripts used for the study, with their local import dependencies, the analysis plans and the provenance records |

## Reproducing the reported AUCs and intervals

```
python -m venv .venv
.venv/bin/pip install -r requirements.lock.txt      # .venv\Scripts\pip on Windows
cd supplementary_data_s1
../.venv/bin/python reproduce_auc.py
```

It recomputes all 20 comparisons (seed 20260918, 10,000 outcome-stratified patient
bootstrap resamples, percentile intervals) and prints

```
All 20 AUC comparisons and confidence intervals reproduced within 1e-12.
```

This was checked in a freshly created environment on 2026-09-18 with Python 3.14.3 and
the pinned versions in `requirements.lock.txt`, and CI repeats it, together with the
checksum check, on every push. `requirements-report.lock.txt` adds the plotting stack
used to draw the figures.

## What is not here

No MRI data. The image-processing scripts in `source_snapshot/` read the original TCIA
files and local manifests; they are an archival snapshot rather than a one-command
distribution, and rebuilding the scores from a fresh DICOM download has not been
independently tested. Running them is not needed to reproduce the AUC table.

## Data source and licences

Newitt DC et al. ACRIN 6698/I-SPY2 Breast DWI. The Cancer Imaging Archive, 2021.
https://doi.org/10.7937/tcia.kk02-6d95, CC BY 4.0. The participant IDs and scores in
`paired_scores.csv` are derived from it.

The code is released under the MIT licence (`LICENSE`).

## Citation

See `CITATION.cff`.
