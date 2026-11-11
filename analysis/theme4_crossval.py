"""Honest out-of-sample probabilities for the revision's added models.

Both added analyses of the BPEX-107031 revision need predicted probabilities
the primary analysis does not produce: the calibration and decision curve of
Referee 2 point 4, and the combined score-plus-subtype model of Referee 1
point 6. One implementation serves both.

**Leave-one-out is wrong here, and measurably so.** When participant *i* is
held out, the fit sees a sample whose pCR prevalence is lower by one if *i* is
a pCR participant and higher by one if *i* is not. Every pCR participant is
therefore scored by a model with a slightly lower intercept than every non-pCR
participant, which orders the outcomes downwards on its own, with no help from
the predictor. On a predictor generated independently of the outcome at this
sample size the leave-one-out AUC averages **0.31**, not 0.50, and reaches
0.00; dropping each fit's intercept only lifts it to 0.43, because the
coefficient leaks too. A weak predictor's leave-one-out AUC is mostly this
artefact.

**Stratified k-fold removes it**, because stratification holds the outcome
prevalence in every training fold. On the same independent predictor it
averages 0.475 for one predictor and 0.467 for four -- unbiased to within the
ordinary pessimism of training on less data, which errs in the safe direction.
``scripts/test_crossval.py`` keeps both facts as tests.

**For a single monotone predictor there is nothing to cross-validate.** A
one-predictor logistic model is a monotone reparametrisation of the score, so
its AUC *is* the score's AUC by construction, and any difference a resampling
scheme reports is leakage rather than optimism. Cross-validated discrimination
is therefore reported only for models that combine predictors.
"""

from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

SEED = 20260918          # the seed the published analysis used
FOLDS = 10
PENALTY = 1e6            # a very low penalty, approximating an unpenalised fit


def model() -> LogisticRegression:
    return LogisticRegression(C=PENALTY, max_iter=5000)


def cross_validated_probabilities(x: np.ndarray, y: np.ndarray,
                                  seed: int = SEED,
                                  folds: int = FOLDS) -> np.ndarray:
    """Each participant's probability from a fit that did not see them.

    The folds are stratified on the outcome, so every training fold carries
    the cohort's pCR prevalence and no participant is scored by a model whose
    intercept encodes their own outcome.
    """
    x = np.asarray(x, float)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    y = np.asarray(y, int)
    out = np.full(y.size, np.nan)
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    for train, test in splitter.split(x, y):
        fitted = model().fit(x[train], y[train])
        out[test] = fitted.predict_proba(x[test])[:, 1]
    assert not np.isnan(out).any(), "a participant was never in a test fold"
    return out


def concordance(score: np.ndarray, label: np.ndarray) -> float:
    """Pairwise concordance with half credit for ties: the AUC."""
    pos, neg = score[label == 1], score[label == 0]
    if not pos.size or not neg.size:
        return float("nan")
    d = pos[:, None] - neg[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size)


def log_likelihood(x: np.ndarray, y: np.ndarray) -> tuple[LogisticRegression, float]:
    """A near-unpenalised fit and its log-likelihood, computed from the fit."""
    x = np.asarray(x, float)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    fitted = model().fit(x, y)
    p = np.clip(fitted.predict_proba(x)[:, 1], 1e-12, 1 - 1e-12)
    return fitted, float((y * np.log(p) + (1 - y) * np.log(1 - p)).sum())


def calibration(prob: np.ndarray, label: np.ndarray) -> dict:
    """Slope and intercept of the outcome on the predicted log odds."""
    eps = 1e-12
    p = np.clip(prob, eps, 1 - eps)
    logit = np.log(p / (1 - p))
    fitted = model().fit(logit.reshape(-1, 1), label)
    return {"slope": float(fitted.coef_[0][0]),
            "intercept": float(fitted.intercept_[0]),
            "mean_predicted": float(prob.mean()),
            "observed_rate": float(label.mean()),
            "note": "slope 1 and intercept 0 is perfect; a slope below 1 means "
                    "the predictions are too extreme"}


def net_benefit(prob: np.ndarray, label: np.ndarray,
                thresholds: np.ndarray) -> list[dict]:
    """Decision curve: net benefit of treating at each threshold probability."""
    n = label.size
    prevalence = float(label.mean())
    curve = []
    for pt in thresholds:
        called = prob >= pt
        tp = int((called & (label == 1)).sum())
        fp = int((called & (label == 0)).sum())
        w = pt / (1 - pt)
        curve.append({"threshold": float(pt),
                      "net_benefit": float(tp / n - fp / n * w),
                      "treat_all": float(prevalence - (1 - prevalence) * w),
                      "called_positive": int(called.sum())})
    return curve
