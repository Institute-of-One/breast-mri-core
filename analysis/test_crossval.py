"""Leave-one-out leaks the outcome; stratified k-fold does not.

The calibration slope and the combined-model AUC of the BPEX-107031 revision
were first computed leave-one-out, and the ``change`` endpoint returned an AUC
of 0.126 for a score whose own AUC is 0.460. These tests are the
counterexamples that found it and the check that the replacement is sound, so
that the method statement in the manuscript is a demonstrated claim.

    python -m unittest scripts.test_crossval
"""

import unittest

import numpy as np
from sklearn.linear_model import LogisticRegression

from theme4_crossval import (FOLDS, SEED, concordance,
                             cross_validated_probabilities)

# The published cohort's shape: 117 participants, 36 of them pCR.
LABEL = np.array([1] * 36 + [0] * 81)
REPLICATES = 20


def leave_one_out(x, y, keep_intercept=True):
    """The estimator that was wrong, kept so the test can show it is wrong."""
    out = np.empty(y.size)
    for i in range(y.size):
        keep = np.ones(y.size, bool)
        keep[i] = False
        fitted = LogisticRegression(C=1e6, max_iter=5000).fit(x[keep], y[keep])
        out[i] = (fitted.predict_proba(x[i:i + 1])[0, 1] if keep_intercept
                  else float(x[i] @ fitted.coef_[0]))
    return out


def noise(rng, columns=1):
    """A predictor generated without reference to the outcome."""
    return rng.normal(size=(LABEL.size, columns))


class LeaveOneOutLeaksTests(unittest.TestCase):
    def test_leave_one_out_is_biased_far_below_one_half(self):
        """An independent predictor should score 0.50 and scores about 0.31."""
        rng = np.random.default_rng(1)
        got = [concordance(leave_one_out(noise(rng), LABEL), LABEL)
               for _ in range(REPLICATES)]
        self.assertLess(np.mean(got), 0.40)
        self.assertLess(min(got), 0.10)

    def test_the_in_sample_concordance_of_the_same_draws_is_one_half(self):
        """So the bias is the estimator and not the generated predictor."""
        rng = np.random.default_rng(1)
        got = [concordance(noise(rng)[:, 0], LABEL) for _ in range(REPLICATES)]
        self.assertAlmostEqual(float(np.mean(got)), 0.5, delta=0.03)

    def test_dropping_the_intercept_does_not_repair_it(self):
        """The coefficient leaks as well, so the slope-only form is still low."""
        rng = np.random.default_rng(1)
        got = [concordance(leave_one_out(noise(rng), LABEL, keep_intercept=False),
                           LABEL) for _ in range(REPLICATES)]
        self.assertLess(np.mean(got), 0.47)


class StratifiedKFoldTests(unittest.TestCase):
    def test_one_predictor_is_unbiased(self):
        rng = np.random.default_rng(1)
        got = [concordance(cross_validated_probabilities(noise(rng), LABEL), LABEL)
               for _ in range(REPLICATES)]
        self.assertAlmostEqual(float(np.mean(got)), 0.5, delta=0.05)

    def test_four_predictors_are_unbiased(self):
        rng = np.random.default_rng(1)
        got = [concordance(cross_validated_probabilities(noise(rng, 4), LABEL),
                           LABEL) for _ in range(REPLICATES)]
        self.assertAlmostEqual(float(np.mean(got)), 0.5, delta=0.05)

    def test_every_training_fold_keeps_the_prevalence(self):
        """This is the property that removes the intercept leak."""
        from sklearn.model_selection import StratifiedKFold
        x = np.arange(LABEL.size, dtype=float).reshape(-1, 1)
        splitter = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=SEED)
        for train, _ in splitter.split(x, LABEL):
            self.assertAlmostEqual(LABEL[train].mean(), LABEL.mean(), delta=0.01)

    def test_a_real_signal_is_recovered_and_slightly_pessimistic(self):
        """Cross-validation should err downwards, not upwards."""
        rng = np.random.default_rng(1)
        for _ in range(5):
            x = (LABEL * 0.8 + rng.normal(size=LABEL.size)).reshape(-1, 1)
            inside = concordance(x[:, 0], LABEL)
            outside = concordance(cross_validated_probabilities(x, LABEL), LABEL)
            self.assertGreater(outside, 0.55)
            self.assertLessEqual(outside, inside + 0.02)

    def test_every_participant_is_scored_exactly_once(self):
        rng = np.random.default_rng(2)
        prob = cross_validated_probabilities(noise(rng), LABEL)
        self.assertEqual(prob.size, LABEL.size)
        self.assertFalse(np.isnan(prob).any())

    def test_the_result_is_reproducible(self):
        rng = np.random.default_rng(3)
        x = noise(rng)
        a = cross_validated_probabilities(x, LABEL)
        b = cross_validated_probabilities(x, LABEL)
        np.testing.assert_array_equal(a, b)


if __name__ == "__main__":
    unittest.main()
