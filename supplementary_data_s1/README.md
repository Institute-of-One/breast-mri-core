# Supplementary Data S1

The paired_scores.csv file contains the 117 public TCIA/BMMR2 participant IDs,
binary pCR labels, subtype labels and fixed scores used in the reported analysis.
T0 and T1 score columns are NEGATIVE volume in mL so that higher scores favor pCR.
Change is T0 minus T1; relative is (T0 minus T1)/T0, without multiplication by 100.

Run `python reproduce_auc.py` with numpy, scipy and scikit-learn installed.
It recomputes the 20 paired comparisons using the original seed 20260918,
10,000 outcome-stratified patient-bootstrap resamples, and percentile intervals.
It asserts agreement with the frozen results to absolute tolerance 1e-12.
It does not train a model or overwrite the frozen results.

The original results.json retains the original primary label all/change.
The manuscript transparently describes the later post-hoc emphasis on T0,
relative reduction and T1. This historical label has not been edited.

The source_snapshot directory contains analysis scripts and their local import
dependencies, retained plans and provenance. It is an archival code snapshot,
not a one-command image-download distribution. Full image processing requires
the original TCIA files and the local manifests described in the scripts.
Reconstruction from a clean DICOM download has not been independently tested.
Do not run acquisition or fitting scripts merely to reproduce the AUC table.

Source data: Newitt DC et al. ACRIN 6698/I-SPY2 Breast DWI. TCIA, 2021.
https://doi.org/10.7937/tcia.kk02-6d95 (CC BY 4.0).
The supplied scores are derived data; no original MRI arrays are included.
Code authorship and reuse terms should be finalized by the submitting author
before a separate public software release. No public repository is claimed.
