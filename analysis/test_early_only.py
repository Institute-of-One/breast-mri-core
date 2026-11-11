"""The late-dependence experiment, exercised against planted defects.

A containment fraction of 0.999 is worthless if the measurement would return
0.999 for a support that is *not* an early-only selection, and the guard that
forbids the late image is worthless if nothing has ever tripped it. Each test
here plants the defect the corresponding check exists to catch, and fails if
the check stays quiet.

    python -m unittest scripts.test_early_only
"""

import unittest

import numpy as np

from theme4_early_only import (F1_CONTAINMENT_FLOOR, LateImageTouched, _Guard,
                               automatic_support, bounding_box, containment,
                               dilate, overlap, pe_in_box, percent_enhancement,
                               stored_threshold)

SHAPE = (12, 40, 40)


def phantom(threshold=70.0, seed=0):
    """A precontrast volume and an early volume with one enhancing lesion.

    Inside the lesion the enhancement is above ``threshold`` but below the next
    threshold the archive uses, so reading the wrong threshold changes the
    answer. The lesion is an ellipsoid, not a cube, and the corners of its
    bounding box hold **low-signal voxels whose ratio is enormous** -- which is
    why the specification has a background rule at all. Without that rule a
    rebuilt support swallows them.

    The returned support is exactly the lesion, so a correct containment
    measurement must return 1.0 on it.
    """
    rng = np.random.default_rng(seed)
    pre = rng.uniform(80, 120, SHAPE)
    early = pre * 1.05                                  # 5 %: well below
    z, y, x = np.ogrid[:SHAPE[0], :SHAPE[1], :SHAPE[2]]
    blob = (((z - 5.5) / 2.0) ** 2 + ((y - 16.5) / 5.0) ** 2
            + ((x - 16.5) / 5.0) ** 2) <= 1.0
    # Above this examination's threshold but below the next one up.
    early[blob] = pre[blob] * (1 + 1.2 * threshold / 100)
    box = np.zeros(SHAPE, bool)
    box[4:8, 12:22, 12:22] = True
    faint = box & ~blob
    pre[faint] = 3.0                                    # no real signal
    early[faint] = 30.0                                 # ratio of 900 %
    return pre, early, blob


def values(threshold=70.0, count=None, ser_threshold=None):
    t = ("%g" % threshold)
    s = t if ser_threshold is None else ("%g" % ser_threshold)
    return {"pe": {"label": "FTVpe(PEearly>%s%%)" % t,
                   "lower": 0.0, "upper": 3.0,
                   "count": 0 if count is None else count, "volume_cc": 1.0},
            "ser": {"label": "FTVser(PEearly>%s%%, SER>0.9)" % s,
                    "lower": 0.9, "upper": 3.0, "count": 0, "volume_cc": 1.0}}


class ContainmentDiscriminates(unittest.TestCase):
    """The measurement must not return 1.0 for everything put in front of it."""

    def test_a_truly_early_only_support_is_fully_contained(self):
        pre, early, blob = phantom()
        got = containment(pre, early, blob, values(70.0, blob.sum()))
        self.assertEqual(got["recall_spec"], 1.0)
        self.assertEqual(got["failing_voxels_spec"], 0)
        self.assertTrue(got["support_equals_stored_pe_count"])

    def test_a_support_that_is_not_early_only_is_caught(self):
        """The planted defect: a support holding voxels the early image never
        selected. If containment still reported ~1, every number in part A
        would be an artefact of the measurement rather than a property of the
        archive."""
        pre, early, blob = phantom()
        planted = blob.copy()
        planted[4:8, 28:38, 28:38] = True               # non-enhancing tissue
        got = containment(pre, early, planted, values(70.0))
        self.assertLess(got["recall_spec"], 0.55)
        self.assertLess(got["recall_spec"], F1_CONTAINMENT_FLOOR)
        self.assertGreater(got["failing_voxels_spec"], 0)

    def test_a_shifted_support_is_caught(self):
        pre, early, blob = phantom()
        got = containment(pre, early, np.roll(blob, 14, axis=2),
                          values(70.0))
        self.assertLess(got["recall_spec"], F1_CONTAINMENT_FLOOR)

    def test_the_wrong_threshold_changes_the_answer(self):
        """Six participants carry 40, 50 or 60 %. Reading the stored label is
        not decoration: at 70 % a 40 % examination reports a failure that is
        not there."""
        pre, early, blob = phantom(threshold=40.0)
        right = containment(pre, early, blob, values(40.0))
        wrong = containment(pre, early, blob, values(70.0))
        self.assertEqual(right["recall_spec"], 1.0)
        self.assertEqual(wrong["recall_spec"], 0.0)

    def test_voxels_without_precontrast_signal_count_as_failing(self):
        pre, early, blob = phantom()
        pre = pre.copy()
        pre[blob] = np.where(np.arange(int(blob.sum())).reshape(-1) < 10,
                             0.0, pre[blob])
        got = containment(pre, early, blob, values(70.0))
        self.assertEqual(got["support_voxels_without_precontrast_signal"], 10)
        self.assertEqual(got["failing_voxels_spec"], 10)

    def test_an_empty_supplied_support_is_undefined_not_perfect(self):
        """Two examinations have a stored FTVpe of 0 cc. Reporting 1.0 there
        would add two free successes to the containment distribution, and
        reporting 0.0 would add two free failures (plan §10, amendment 1)."""
        pre, early, _ = phantom()
        got = containment(pre, early, np.zeros(SHAPE, bool), values(70.0))
        self.assertEqual(got["status"], "empty supplied support")
        self.assertIsNone(got["recall_spec"])
        self.assertIsNone(got["recall_floored"])
        self.assertEqual(got["failing_voxels_spec"], 0)

    def test_a_non_empty_support_is_marked_measured(self):
        pre, early, blob = phantom()
        self.assertEqual(containment(pre, early, blob,
                                     values(70.0))["status"], "measured")

    def test_percent_enhancement_never_passes_an_undefined_voxel(self):
        pe = percent_enhancement(np.array([0.0, 1.0]), np.array([5.0, 5.0]))
        self.assertEqual(pe[0], -np.inf)
        self.assertEqual(pe[1], 400.0)


class StoredThresholdIsRead(unittest.TestCase):

    def test_both_labels_must_name_the_same_threshold(self):
        with self.assertRaises(AssertionError):
            stored_threshold(values(70.0, ser_threshold=50.0))

    def test_each_declared_threshold_parses(self):
        for t in (40.0, 50.0, 60.0, 70.0):
            self.assertEqual(stored_threshold(values(t)), t)

    def test_an_unrecognised_label_is_refused(self):
        bad = values(70.0)
        bad["pe"]["label"] = "FTVpe(peak enhancement)"
        with self.assertRaises(AssertionError):
            stored_threshold(bad)


class TheGuardStopsTheLateImage(unittest.TestCase):

    def test_reading_late_before_the_unlock_raises(self):
        guard = _Guard(pre=1, early=2, late=3)
        self.assertEqual(guard["pre"], 1)
        with self.assertRaises(LateImageTouched):
            guard["late"]

    def test_a_support_builder_that_reads_late_is_caught(self):
        """The planted defect: a rule that peeks at the late image. This is the
        one thing the experiment must not do, so it must not be possible to do
        it quietly."""
        pre, early, blob = phantom()
        guard = _Guard(pre=pre, early=early, late=early * 1.1)

        def leaky(g):
            return pe_in_box(g["pre"], g["late"], blob, 70.0, 0)

        with self.assertRaises(LateImageTouched):
            leaky(guard)
        guard.unlock()
        self.assertIsNotNone(leaky(guard))      # only after the SER step opens

    def test_an_unknown_key_is_not_silently_allowed(self):
        with self.assertRaises(KeyError):
            _Guard(pre=1)["enhanced"]


class PeInBoxUsesOnlyWhatItClaims(unittest.TestCase):

    def test_it_recovers_the_blob(self):
        pre, early, blob = phantom()
        support = pe_in_box(pre, early, blob, 70.0, 0)
        dice, _ = overlap(support, blob)
        self.assertGreater(dice, 0.99)

    def test_a_margin_can_only_grow_the_box(self):
        pre, early, blob = phantom()
        got = [pe_in_box(pre, early, blob, 70.0, m) for m in (0, 2, 4, 8)]
        for small, big in zip(got, got[1:]):
            self.assertTrue(np.all(big | small == big))
            self.assertGreaterEqual(np.count_nonzero(big),
                                    np.count_nonzero(small))

    def test_the_box_is_the_bounding_box_of_the_supplied_support(self):
        _, _, blob = phantom()
        box = bounding_box(blob)
        self.assertEqual([(s.start, s.stop) for s in box],
                         [(4, 8), (12, 22), (12, 22)])

    def test_a_margin_is_clipped_to_the_array(self):
        box = dilate(bounding_box(np.ones(SHAPE, bool)), 8, SHAPE)
        self.assertEqual([(s.start, s.stop) for s in box],
                         [(0, 12), (0, 40), (0, 40)])

    def test_an_empty_support_is_reported_and_not_invented(self):
        pre, early, blob = phantom()
        self.assertIsNone(pe_in_box(pre, early, blob, 1000.0, 0))


class TheAutomaticRedrawFindsTheBrighterSide(unittest.TestCase):
    """Part C's negative result, as a test rather than as an anecdote."""

    def test_it_lands_on_the_wrong_blob_when_that_one_is_larger(self):
        pre, early, blob = phantom()
        decoy = np.zeros(SHAPE, bool)
        decoy[3:9, 10:30, 26:38] = True                 # larger, other side
        early = early.copy()
        early[decoy] = pre[decoy] * 3.0
        support = automatic_support(pre, early, 70.0)
        self.assertIsNotNone(support)
        # It returns a support, and that support is disjoint from the lesion.
        self.assertEqual(overlap(support, blob)[0], 0.0)
        centre = np.argwhere(support).mean(axis=0)
        self.assertLess(abs(centre[2] - np.argwhere(decoy).mean(axis=0)[2]),
                        abs(centre[2] - np.argwhere(blob).mean(axis=0)[2]))


if __name__ == "__main__":
    unittest.main()
