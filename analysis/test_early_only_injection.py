"""Break the experiment on purpose; the tests must go red.

``scripts/test_early_only.py`` passing proves nothing on its own: a test that
would also pass against broken code is not a check. This file damages the
module in four ways that would each silently flatter the answer, and fails if
the corresponding test stays green.

    python -m unittest scripts.test_early_only_injection
"""

import unittest

import numpy as np

import theme4_early_only as mod
import test_early_only as suite


def run(name: str) -> bool:
    """True when that test passes."""
    loader = unittest.TestLoader()
    tests = [t for t in loader.loadTestsFromModule(suite)._tests
             for t in t._tests if t.id().endswith(name)]
    assert len(tests) == 1, (name, [t.id() for t in tests])
    result = unittest.TestResult()
    tests[0].run(result)
    return result.wasSuccessful()


class InjectedDefectsAreCaught(unittest.TestCase):

    def setUp(self):
        self.saved = {k: getattr(mod, k) for k in
                      ("percent_enhancement", "stored_threshold",
                       "background_mask", "bounding_box")}

    def tearDown(self):
        for k, v in self.saved.items():
            setattr(mod, k, v)

    def test_green_before_any_injection(self):
        for name in ("test_voxels_without_precontrast_signal_count_as_failing",
                     "test_the_wrong_threshold_changes_the_answer",
                     "test_a_support_that_is_not_early_only_is_caught",
                     "test_it_recovers_the_blob"):
            self.assertTrue(run(name), name)

    def test_treating_undefined_enhancement_as_passing_is_caught(self):
        """A voxel with no precontrast signal has no percent enhancement. If it
        were counted as passing, containment would rise for free."""
        def generous(pre, early):
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.divide(100.0 * (early - pre), pre,
                                 out=np.full(pre.shape, np.inf), where=pre > 0)
        mod.percent_enhancement = generous
        self.assertFalse(
            run("test_voxels_without_precontrast_signal_count_as_failing"))

    def test_ignoring_the_stored_threshold_is_caught(self):
        """The defect this one guards against is real: six participants are not
        analysed at 70 %, and an earlier run reported a containment of 0.178 on
        one of them for exactly this reason."""
        mod.stored_threshold = lambda values: 70.0
        self.assertFalse(run("test_the_wrong_threshold_changes_the_answer"))

    def test_a_containment_measure_that_passes_everything_is_caught(self):
        """The planted defect that matters most: a measurement that returns
        1.0 whatever support it is given."""
        mod.percent_enhancement = lambda pre, early: np.full(pre.shape, np.inf)
        self.assertFalse(run("test_a_support_that_is_not_early_only_is_caught"))
        self.assertFalse(run("test_a_shifted_support_is_caught"))

    def test_dropping_the_background_rule_is_caught(self):
        """Low-signal voxels have enormous ratios. Without the published
        background rule the rebuilt support swallows the ones sitting inside
        the lesion's own bounding box, so it stops reproducing the lesion."""
        mod.background_mask = lambda pre, box: np.ones(pre.shape, bool)
        self.assertFalse(run("test_it_recovers_the_blob"))


if __name__ == "__main__":
    unittest.main()
