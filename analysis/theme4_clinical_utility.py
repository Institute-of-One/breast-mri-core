"""Operating points, reclassification, calibration and net benefit.

    python scripts/theme4_clinical_utility.py

Referee 1 point 8 asks for sensitivity at a fixed specificity and for the
change in patient classification. Referee 2 point 4 asks additionally for
calibration or a decision curve. An AUC is a ranking summary, and two scores
with the same AUC can still send different patients across a decision
boundary, so the request is right.

Three things are computed, all on the published 117-patient cohort, all post
hoc, and all applying the same procedure to both arms so that the comparison
between them stays fair:

1. **Sensitivity at fixed specificity** (0.80 and 0.90), with the threshold
   taken on each arm's own score, and paired patient bootstrap intervals using
   the seed the published analysis used.
2. **Reclassification** at those operating points: how many pCR and how many
   non-pCR participants the substitution moves across the boundary, in each
   direction. Counts, not a net reclassification index, because an index hides
   the direction and the denominators.
3. **Calibration and net benefit**, which need predicted probabilities the
   primary analysis deliberately does not produce. Rather than decline, a
   single-predictor logistic model supplies them out of sample, through the
   **stratified ten-fold** scheme of ``theme4_crossval``. Leave-one-out was
   used first and is wrong here: holding out a participant moves the fitted
   intercept in the direction of their own outcome, which on an independent
   predictor drives the AUC to 0.31 and flattens the calibration slope. The
   module's docstring and ``scripts/test_crossval.py`` carry the measurements.
   This is an added model, labelled as one: the primary comparison does not use
   it, and its apparent performance is not evidence for the primary endpoint.

Writes ``data/derived/theme4_revision_v1/clinical_utility.json``. The published
artefacts are not modified.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

from theme4_crossval import (SEED as CV_SEED, FOLDS, calibration,
                             concordance as auc, cross_validated_probabilities,
                             net_benefit)

ROOT = Path(__file__).resolve().parent.parent
SCORES = ROOT / "data/derived/theme4_pcr_v1/paired_scores.csv"
OUT = ROOT / "data/derived/theme4_revision_v1"

# The seed the published analysis used, so the resampling draws are the same
# ones its intervals were built from.
SEED = 20260918
BOOTSTRAP = 10000
SPECIFICITIES = (0.80, 0.90)
ENDPOINTS = {"T0": "real_T0", "T1": "real_T1", "relative": "real_relative"}
ARMS = {"real": "real_", "sub": "sub_"}

sys.stdout.reconfigure(encoding="utf-8", errors="replace")



def threshold_at_specificity(score: np.ndarray, label: np.ndarray,
                             specificity: float) -> float:
    """The lowest cut-off whose specificity is at least the target.

    Higher score means more likely pCR in every endpoint here, so a
    participant is called positive when the score is at or above the cut-off.
    """
    negatives = np.sort(score[label == 0])
    k = int(np.ceil(specificity * negatives.size))
    k = min(max(k, 1), negatives.size)
    # k negatives must fall below the cut-off.
    return float(np.nextafter(negatives[k - 1], np.inf))


def sens_spec(score: np.ndarray, label: np.ndarray, cut: float) -> tuple[float, float]:
    called = score >= cut
    tp = int((called & (label == 1)).sum())
    fp = int((called & (label == 0)).sum())
    npos = int((label == 1).sum())
    nneg = int((label == 0).sum())
    return tp / npos, 1 - fp / nneg


def paired_bootstrap(real: np.ndarray, sub: np.ndarray, label: np.ndarray,
                     specificity: float) -> dict:
    rng = np.random.default_rng(SEED)
    pos = np.flatnonzero(label == 1)
    neg = np.flatnonzero(label == 0)
    diffs, reals, subs = [], [], []
    for _ in range(BOOTSTRAP):
        take = np.concatenate([rng.choice(pos, pos.size, replace=True),
                               rng.choice(neg, neg.size, replace=True)])
        y = label[take]
        if y.sum() == 0 or y.sum() == y.size:
            continue
        r, s = real[take], sub[take]
        sr = sens_spec(r, y, threshold_at_specificity(r, y, specificity))[0]
        ss = sens_spec(s, y, threshold_at_specificity(s, y, specificity))[0]
        reals.append(sr)
        subs.append(ss)
        diffs.append(ss - sr)
    q = lambda a: [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]
    return {"sensitivity_real_ci": q(reals), "sensitivity_sub_ci": q(subs),
            "difference_ci": q(diffs), "replicates": len(diffs)}


def reclassification(real: np.ndarray, sub: np.ndarray, label: np.ndarray,
                     specificity: float) -> dict:
    cr = threshold_at_specificity(real, label, specificity)
    cs = threshold_at_specificity(sub, label, specificity)
    called_r, called_s = real >= cr, sub >= cs
    out = {}
    for name, mask in (("pCR", label == 1), ("non_pCR", label == 0)):
        out[name] = {
            "n": int(mask.sum()),
            "positive_to_negative": int((mask & called_r & ~called_s).sum()),
            "negative_to_positive": int((mask & ~called_r & called_s).sum()),
            "unchanged": int((mask & (called_r == called_s)).sum()),
        }
    out["any_change"] = int((called_r != called_s).sum())
    return out





def main() -> int:
    rows = list(csv.DictReader(SCORES.open(encoding="utf-8-sig")))
    label = np.array([int(r["pcr"]) for r in rows])
    result = {"n": len(rows), "n_pcr": int(label.sum()),
              "seed": SEED, "bootstrap": BOOTSTRAP,
              "cross_validation": {"scheme": "stratified %d-fold" % FOLDS,
                                   "seed": CV_SEED},
              "note": "post hoc; the decision curve and calibration use a "
                      "cross-validated single-predictor logistic model that "
                      "the primary analysis does not use",
              "endpoints": {}}
    thresholds = np.arange(0.05, 0.75, 0.01)

    for name in ENDPOINTS:
        real = np.array([float(r["real_" + name]) for r in rows])
        sub = np.array([float(r["sub_" + name]) for r in rows])
        entry = {"auc_real": auc(real, label), "auc_sub": auc(sub, label),
                 "operating_points": {}, "calibration": {}, "decision_curve": {}}
        for sp in SPECIFICITIES:
            cr = threshold_at_specificity(real, label, sp)
            cs = threshold_at_specificity(sub, label, sp)
            sr, spr = sens_spec(real, label, cr)
            ss, sps = sens_spec(sub, label, cs)
            entry["operating_points"]["specificity_%.2f" % sp] = {
                "real": {"sensitivity": sr, "specificity": spr, "threshold": cr},
                "sub": {"sensitivity": ss, "specificity": sps, "threshold": cs},
                "sensitivity_difference": ss - sr,
                "bootstrap": paired_bootstrap(real, sub, label, sp),
                "reclassification": reclassification(real, sub, label, sp),
            }
            print("%-9s spec %.2f | sens real %.3f sub %.3f (diff %+.3f) | "
                  "reclassified %d of %d"
                  % (name, sp, sr, ss, ss - sr,
                     entry["operating_points"]["specificity_%.2f" % sp]
                     ["reclassification"]["any_change"], len(rows)))
        for arm, prefix in ARMS.items():
            score = np.array([float(r[prefix + name]) for r in rows])
            prob = cross_validated_probabilities(score, label)
            entry["calibration"][arm] = calibration(prob, label)
            entry["decision_curve"][arm] = net_benefit(prob, label, thresholds)
        c = entry["calibration"]
        print("          calibration slope real %.3f sub %.3f (1.0 is ideal)"
              % (c["real"]["slope"], c["sub"]["slope"]))
        result["endpoints"][name] = entry

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "clinical_utility.json").write_text(json.dumps(result, indent=2),
                                               encoding="utf-8")
    print()
    print("wrote", (OUT / "clinical_utility.json").relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
