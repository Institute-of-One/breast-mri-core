"""Does the volume score discriminate within a subtype, or only between subtypes?

    python scripts/theme4_subtype_adjusted.py

Referee 1 point 6 of BPEX-107031 asks for the analysis adjusted for receptor
subtype. The request is well founded: the pCR rate runs from 16% in
HR+/HER2- to 80% in HR-/HER2+, so any score that happens to track subtype
inherits discrimination it did not earn from the tumour volume.

Three analyses, in increasing strength of claim.

**1. Does the score carry subtype information at all?** A Kruskal-Wallis of
each score across the four subtypes. If it does not, adjustment cannot change
much, and that answers the point directly rather than by assertion.

**2. A within-subtype pooled concordance.** The AUC is the share of
pCR / non-pCR pairs the score orders correctly. Restricting the enumeration to
pairs drawn from the *same* subtype and pooling the counts across subtypes
gives a concordance that no subtype information can inflate. It costs pairs:
744 of the 2916 pairs survive, so the interval is wider, and the smallest cell
(HR-/HER2+ has 2 non-pCR participants) is stated rather than smoothed over.
The bootstrap resamples within subtype-by-outcome cells so each replicate keeps
the design, and both arms are resampled together as the published analysis
does.

**3. A logistic model with subtype terms.** The score, standardised, with and
without three subtype indicators; the coefficient per standard deviation, a
likelihood-ratio test of the score given the subtype terms, and the
out-of-sample AUC of the *combined* model. The log-likelihood is computed
directly from the fitted probabilities, so the test does not depend on a
library's own reporting. The penalty is set very low to approximate an
unpenalised fit; this is stated because it is an approximation, not a
maximum-likelihood fit.

No cross-validated AUC is reported for a single predictor. A one-predictor
logistic model is a monotone reparametrisation of the score, so its AUC is the
score's own AUC by construction; a resampling scheme applied to it reports
leakage, not optimism. ``theme4_crossval`` carries the measurement that found
this -- a leave-one-out AUC of 0.31 on a predictor built independently of the
outcome -- and the stratified folds that remove it.

Writes ``data/derived/theme4_revision_v1/subtype_adjusted.json``. Nothing
published is modified.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2, kruskal

from theme4_crossval import (FOLDS, SEED as CV_SEED, concordance as pair_concordance,
                             cross_validated_probabilities, log_likelihood)

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCORES = ROOT / "data/derived/theme4_pcr_v1/paired_scores.csv"
OUT = ROOT / "data/derived/theme4_revision_v1"

SEED = 20260918          # the seed the published analysis used
BOOTSTRAP = 10000
ENDPOINTS = ("T0", "T1", "change", "relative")
ARMS = ("real", "sub")


def concordance(score: np.ndarray, label: np.ndarray) -> tuple[float, int]:
    """Pairwise concordance and the number of pairs it rests on."""
    pos, neg = score[label == 1], score[label == 0]
    if not pos.size or not neg.size:
        return float("nan"), 0
    d = pos[:, None] - neg[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size), int(d.size)


def pooled_within(score: np.ndarray, label: np.ndarray,
                  group: np.ndarray) -> tuple[float, int]:
    """Concordance over same-subtype pairs only, pooled across subtypes."""
    hits = pairs = 0.0
    for g in np.unique(group):
        m = group == g
        pos, neg = score[m & (label == 1)], score[m & (label == 0)]
        if not pos.size or not neg.size:
            continue
        d = pos[:, None] - neg[None, :]
        hits += (d > 0).sum() + 0.5 * (d == 0).sum()
        pairs += d.size
    return (hits / pairs if pairs else float("nan")), int(pairs)


def cell_bootstrap(real: np.ndarray, sub: np.ndarray, label: np.ndarray,
                   group: np.ndarray) -> dict:
    """Resample within subtype-by-outcome cells; both arms take the same draw."""
    rng = np.random.default_rng(SEED)
    cells = [np.flatnonzero((group == g) & (label == y))
             for g in np.unique(group) for y in (0, 1)]
    cells = [c for c in cells if c.size]
    reals, subs, diffs = [], [], []
    for _ in range(BOOTSTRAP):
        take = np.concatenate([rng.choice(c, c.size, replace=True) for c in cells])
        r, _ = pooled_within(real[take], label[take], group[take])
        s, _ = pooled_within(sub[take], label[take], group[take])
        reals.append(r)
        subs.append(s)
        diffs.append(s - r)
    q = lambda a: [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]
    return {"real_ci": q(reals), "sub_ci": q(subs), "difference_ci": q(diffs),
             "replicates": len(diffs),
             "smallest_cell": int(min(c.size for c in cells))}




def main() -> int:
    rows = list(csv.DictReader(SCORES.open(encoding="utf-8-sig")))
    label = np.array([1 if r["pcr"] in ("1", "pCR") else 0 for r in rows])
    subtype = np.array([r["subtype"] for r in rows])
    names = sorted(set(subtype))
    # Three indicators for four subtypes; the largest is the reference level.
    reference = max(names, key=lambda g: (subtype == g).sum())
    others = [g for g in names if g != reference]
    design = np.column_stack([(subtype == g).astype(float) for g in others])

    result = {
        "scope": "post hoc; Referee 1 point 6",
        "n": len(rows), "n_pcr": int(label.sum()),
        "seed": SEED, "bootstrap": BOOTSTRAP,
        "reference_subtype": reference,
        "subtypes": {g: {"n": int((subtype == g).sum()),
                         "pcr": int(label[subtype == g].sum()),
                         "pcr_rate": float(label[subtype == g].mean())}
                     for g in names},
        "logistic_note": "a very low penalty approximates an unpenalised fit; "
                         "the log-likelihood is computed from the fitted "
                         "probabilities, not read from the library",
        "cross_validation": {"scheme": "stratified %d-fold" % FOLDS,
                             "seed": CV_SEED,
                             "applied_to": "models that combine predictors only"},
        "endpoints": {},
    }

    for endpoint in ENDPOINTS:
        entry = {"score_vs_subtype": {}, "arms": {}}
        for arm in ARMS:
            score = np.array([float(r["%s_%s" % (arm, endpoint)]) for r in rows])
            groups = [score[subtype == g] for g in names]
            k = kruskal(*groups)
            entry["score_vs_subtype"][arm] = {
                "kruskal_statistic": float(k.statistic),
                "p_value": float(k.pvalue),
                "median_by_subtype": {g: float(np.median(score[subtype == g]))
                                      for g in names},
            }

        real = np.array([float(r["real_" + endpoint]) for r in rows])
        sub = np.array([float(r["sub_" + endpoint]) for r in rows])

        unadj_r, all_pairs = concordance(real, label)
        unadj_s, _ = concordance(sub, label)
        adj_r, kept = pooled_within(real, label, subtype)
        adj_s, _ = pooled_within(sub, label, subtype)
        entry["concordance"] = {
            "unadjusted": {"real": unadj_r, "sub": unadj_s,
                           "difference_sub_minus_real": unadj_s - unadj_r,
                           "pairs": all_pairs},
            "within_subtype": {"real": adj_r, "sub": adj_s,
                               "difference_sub_minus_real": adj_s - adj_r,
                               "pairs": kept,
                               "pairs_discarded": all_pairs - kept},
            "within_subtype_bootstrap": cell_bootstrap(real, sub, label, subtype),
            "per_subtype": {g: {"real": concordance(real[subtype == g],
                                                    label[subtype == g])[0],
                                "sub": concordance(sub[subtype == g],
                                                   label[subtype == g])[0],
                                "pairs": concordance(real[subtype == g],
                                                     label[subtype == g])[1]}
                            for g in names},
        }

        for arm, score in (("real", real), ("sub", sub)):
            z = ((score - score.mean()) / score.std(ddof=1)).reshape(-1, 1)
            m_sub, ll_sub = log_likelihood(design, label)
            m_one, ll_one = log_likelihood(z, label)
            m_adj, ll_adj = log_likelihood(np.hstack([z, design]), label)
            statistic = 2 * (ll_adj - ll_sub)
            entry["arms"][arm] = {
                "coefficient_per_sd_unadjusted": float(m_one.coef_[0][0]),
                "odds_ratio_per_sd_unadjusted": float(np.exp(m_one.coef_[0][0])),
                "coefficient_per_sd_adjusted": float(m_adj.coef_[0][0]),
                "odds_ratio_per_sd_adjusted": float(np.exp(m_adj.coef_[0][0])),
                "subtype_coefficients": {g: float(c) for g, c
                                         in zip(others, m_adj.coef_[0][1:])},
                "likelihood_ratio_vs_subtype_only": {
                    "statistic": float(statistic), "df": 1,
                    "p_value": float(chi2.sf(statistic, 1)),
                },
                "auc_score_only": concordance(score, label)[0],
                "auc_score_only_note": "the score's own concordance; a "
                                       "one-predictor logistic model is a "
                                       "monotone transform of it, so there is "
                                       "nothing to cross-validate",
                "cross_validated_auc_subtype_only": pair_concordance(
                    cross_validated_probabilities(design, label), label),
                "cross_validated_auc_score_and_subtype": pair_concordance(
                    cross_validated_probabilities(np.hstack([z, design]), label),
                    label),
            }
        result["endpoints"][endpoint] = entry

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "subtype_adjusted.json").write_text(json.dumps(result, indent=2),
                                               encoding="utf-8")

    print("pCR by subtype (reference level for the model: %s)" % reference)
    for g in names:
        s = result["subtypes"][g]
        print("  %-22s n=%3d  pCR %2d (%.0f%%)" % (g, s["n"], s["pcr"],
                                                   100 * s["pcr_rate"]))
    print()
    for endpoint in ENDPOINTS:
        e = result["endpoints"][endpoint]
        c = e["concordance"]
        print("%s" % endpoint)
        print("  score vs subtype (Kruskal-Wallis): real p=%.4f  sub p=%.4f"
              % (e["score_vs_subtype"]["real"]["p_value"],
                 e["score_vs_subtype"]["sub"]["p_value"]))
        print("  unadjusted      real %.3f  sub %.3f  diff %+.4f   (%d pairs)"
              % (c["unadjusted"]["real"], c["unadjusted"]["sub"],
                 c["unadjusted"]["difference_sub_minus_real"],
                 c["unadjusted"]["pairs"]))
        b = c["within_subtype_bootstrap"]
        print("  within subtype  real %.3f  sub %.3f  diff %+.4f [%+.4f,%+.4f]"
              "   (%d pairs, %d discarded, smallest cell %d)"
              % (c["within_subtype"]["real"], c["within_subtype"]["sub"],
                 c["within_subtype"]["difference_sub_minus_real"],
                 *b["difference_ci"], c["within_subtype"]["pairs"],
                 c["within_subtype"]["pairs_discarded"], b["smallest_cell"]))
        for arm in ARMS:
            a = e["arms"][arm]
            print("  %-4s OR/SD %.3f -> %.3f adjusted | LRT p=%.4f | AUC score "
                  "%.3f | cross-validated: subtype %.3f, score+subtype %.3f"
                  % (arm, a["odds_ratio_per_sd_unadjusted"],
                     a["odds_ratio_per_sd_adjusted"],
                     a["likelihood_ratio_vs_subtype_only"]["p_value"],
                     a["auc_score_only"],
                     a["cross_validated_auc_subtype_only"],
                     a["cross_validated_auc_score_and_subtype"]))
        print()
    print("wrote", (OUT / "subtype_adjusted.json").relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
