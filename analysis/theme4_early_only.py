"""Is the analysis support late-dependent? Three measurements, one load per examination.

    python scripts/theme4_early_only.py --smoke --limit 4   # execution test
    python scripts/theme4_early_only.py --freeze            # pin plan + script
    python scripts/theme4_early_only.py                     # run it

Referee 1 point 3 and Referee 2 point 2 of BPEX-107031 ask whether the paired
comparison is biased because both arms are measured inside an analysis support
that the archive built with the full dynamic series, late image included. The
rule is written down, before any number from part B or part C existed, in
``docs/theme4_early_only_support_plan_v1.md``.

Part A -- **containment**. Does the archive's support contain a single voxel
that the late image could have selected? Each voxel of the supplied support is
tested against percent enhancement computed from the precontrast and early
images alone, at the threshold the archive itself stored for that examination.
There is no constant to choose here: the threshold is read per examination from
the archive's own label, which is not always 70 %. **This measurement was made
before the plan was written** and is recomputed here on the full cohort; the
plan says so.

Part B -- **the paired comparison inside a late-independent support**. A support
rebuilt as percent enhancement above the stored threshold, inside the supplied
support's bounding box dilated by a margin, with the published background rule.
Every step reads only ``pre`` and ``early``; the margin is not chosen but run
over the declared grid.

Part C -- **the fully automatic redraw**, reported as the negative result it is.
The box is found from orthogonal projections with no reader input at all. It
lands on the contralateral breast when that side enhances more.

``_Guard`` wraps the loaded arrays and raises on any access to ``late`` while a
support is being built, so parts B and C cannot quietly use it. The late image
is unlocked once, for the SER measurement, which is the step that is supposed
to use it.

Writes ``data/derived/theme4_revision_v1/early_only_support.json``. Nothing
published is modified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PLAN = ROOT / "docs/theme4_early_only_support_plan_v1.md"
LOCK = ROOT / "docs/theme4_early_only_lock_v1.json"
SOURCE = ROOT / "data/derived/theme4_full_v2"
SCORES = ROOT / "data/derived/theme4_pcr_v1/paired_scores.csv"
OUT = ROOT / "data/derived/theme4_revision_v1"

# Published specification constants. None of these is a choice made here, and
# the enhancement threshold is not even a constant: it is read per examination
# from the archive's own stored label.
BACKGROUND_FRACTION = 0.60
PROJECTION_PERCENTILE = 95
BACKGROUND_PERCENTILE = 95
# The specification gives no margin, so no margin is selected. Declared grid:
MARGINS = (0, 2, 4, 8)
VISITS = ("T0", "T1")
ENDPOINTS = ("T0", "T1", "change", "relative")
# Falsification thresholds for part B, declared in the plan before any result.
F1_CONTAINMENT_FLOOR = 0.95
F2_EMPTY_SHARE = 0.10
F3_MEDIAN_DICE = 0.30
PE_LABEL = re.compile(r"^FTVpe\(PEearly>([0-9]+(?:\.[0-9]+)?)%\)$")
SER_LABEL = re.compile(r"^FTVser\(PEearly>([0-9]+(?:\.[0-9]+)?)%, SER>")


class LateImageTouched(RuntimeError):
    pass


class _Guard:
    """An array view that refuses the late image until it is unlocked."""

    def __init__(self, **arrays):
        self._arrays = arrays
        self._unlocked = False

    def __getitem__(self, key):
        if key == "late" and not self._unlocked:
            raise LateImageTouched(
                "the late image was read while a late-independent support was "
                "being built; that is the one thing this experiment must not do")
        return self._arrays[key]

    def unlock(self):
        self._unlocked = True


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_freeze() -> dict:
    if not LOCK.exists():
        raise SystemExit("no lock: run with --freeze once, before any result "
                         "from part B exists, and keep the lock")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    for name, path in (("plan_sha256", PLAN), ("script_sha256", Path(__file__))):
        now = digest(path)
        if lock[name] != now:
            raise SystemExit(
                "refusing to run: %s has changed since the freeze\n"
                "  frozen %s\n  now    %s\n"
                "A rule edited after a result is not a prespecified rule. If the "
                "change is intended, record it as a named amendment in the plan "
                "and re-freeze deliberately." % (path.name, lock[name], now))
    return lock


def freeze() -> int:
    if LOCK.exists():
        raise SystemExit("a lock already exists; refusing to overwrite it")
    lock = {"plan": str(PLAN.relative_to(ROOT)).replace("\\", "/"),
            "plan_sha256": digest(PLAN),
            "script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
            "script_sha256": digest(Path(__file__)),
            "margins": list(MARGINS),
            "specification_constants": {
                "enhancement_threshold": "read per examination from the "
                                         "archive's stored FTVpe label",
                "background_fraction": BACKGROUND_FRACTION,
                "projection_percentile": PROJECTION_PERCENTILE,
                "background_percentile": BACKGROUND_PERCENTILE},
            "falsification": {"F1_containment_floor": F1_CONTAINMENT_FLOOR,
                              "F2_empty_share": F2_EMPTY_SHARE,
                              "F3_median_dice": F3_MEDIAN_DICE},
            "already_observed_before_the_freeze":
                "part A, the containment fraction, was measured in an "
                "exploratory run before this plan was written; the plan records "
                "the observed value. Parts B and C had produced no comparison "
                "and no overlap statistic at the freeze.",
            "scope": "the rule for testing whether the analysis support is "
                     "late-dependent, pinned before any result from part B"}
    LOCK.write_text(json.dumps(lock, indent=2), encoding="utf-8")
    print("froze the plan and the script:")
    print("  plan   %s" % lock["plan_sha256"])
    print("  script %s" % lock["script_sha256"])
    return 0


# --------------------------------------------------------------------------
# the archive's own threshold


def stored_threshold(values: dict) -> float:
    """The enhancement threshold this examination was analysed at, in percent.

    Not a constant. Six participants of the published cohort carry 40, 50 or
    60 % rather than the nominal 70 %, and testing them at 70 % would report a
    containment failure that is an error of mine, not a property of the
    archive.
    """
    pe = PE_LABEL.match(values["pe"]["label"])
    ser = SER_LABEL.match(values["ser"]["label"])
    assert pe and ser, (values["pe"]["label"], values["ser"]["label"])
    assert pe.group(1) == ser.group(1), (
        "the stored SER label names a different enhancement threshold from the "
        "stored PE label: %s vs %s" % (values["pe"]["label"],
                                       values["ser"]["label"]))
    return float(pe.group(1))


def percent_enhancement(pre: np.ndarray, early: np.ndarray) -> np.ndarray:
    """``100 (early - pre) / pre`` where ``pre > 0``; ``-inf`` elsewhere.

    Undefined is never treated as zero and never as passing: a voxel with no
    precontrast signal is counted as failing the test, which is the reading
    that cannot flatter the result.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(100.0 * (early - pre), pre,
                         out=np.full(pre.shape, -np.inf, dtype=float),
                         where=pre > 0)


def containment(pre, early, eligible, values) -> dict:
    """Part A. What share of the supplied support is early-only selectable?

    Two examinations of the published cohort have an **empty** supplied
    support: the archive found no enhancing voxel at the follow-up visit and
    stored an FTVpe of 0 cc. The containment fraction is 0/0 there. It is
    reported as undefined rather than as 1.0 or 0.0, and those examinations are
    excluded from the distribution (plan §10, amendment 1).
    """
    threshold = stored_threshold(values)
    n = int(np.count_nonzero(eligible))
    if not n:
        return {"status": "empty supplied support",
                "stored_threshold_percent": threshold,
                "stored_label": values["pe"]["label"],
                "support_voxels": 0,
                "stored_pe_count": int(values["pe"]["count"]),
                "support_equals_stored_pe_count":
                    int(values["pe"]["count"]) == 0,
                "recall_spec": None, "recall_floored": None,
                "failing_voxels_spec": 0,
                "support_voxels_without_precontrast_signal": 0}
    pe = percent_enhancement(pre, early)
    inside = pe[eligible]
    # The specification's reading: PE above the threshold.
    passes_spec = int(np.count_nonzero(inside > threshold))
    # The archive's PE maps are stored floored to whole percent, so a support
    # built on them satisfies the stricter test too. Reported, not relied on.
    floored = np.floor(np.maximum(inside, 0.0))
    passes_floored = int(np.count_nonzero(np.isfinite(inside)
                                          & (floored > threshold)))
    return {"status": "measured",
            "stored_threshold_percent": threshold,
            "stored_label": values["pe"]["label"],
            "support_voxels": n,
            "stored_pe_count": int(values["pe"]["count"]),
            "support_equals_stored_pe_count":
                int(values["pe"]["count"]) == n,
            "recall_spec": passes_spec / n,
            "recall_floored": passes_floored / n,
            "failing_voxels_spec": n - passes_spec,
            "support_voxels_without_precontrast_signal":
                int(np.count_nonzero(~np.isfinite(inside)))}


# --------------------------------------------------------------------------
# the supports


def background_mask(pre: np.ndarray, box: tuple) -> np.ndarray:
    """The published background rule, applied inside the given box."""
    cut = BACKGROUND_FRACTION * np.percentile(pre[box], BACKGROUND_PERCENTILE)
    return pre > cut


def dilate(box: tuple, margin: int, shape: tuple) -> tuple:
    if not margin:
        return box
    return tuple(slice(max(0, s.start - margin), min(n, s.stop + margin))
                 for s, n in zip(box, shape))


def bounding_box(mask: np.ndarray) -> tuple:
    where = [np.flatnonzero(mask.any(axis=tuple(a for a in range(3) if a != k)))
             for k in range(3)]
    assert all(w.size for w in where), "cannot bound an empty mask"
    return tuple(slice(int(w.min()), int(w.max()) + 1) for w in where)


def pe_in_box(pre, early, eligible, threshold, margin):
    """Part B's support. Reads ``pre``, ``early`` and the supplied box only.

    The box is the one thing taken from the archive, and it is the one thing
    the specification says a person supplies: the VOI. By part A the supplied
    support is itself an early-only selection, so its bounding box carries no
    information the late image could have put there.
    """
    box = dilate(bounding_box(eligible), margin, pre.shape)
    inside = np.zeros(pre.shape, bool)
    inside[box] = True
    enhancing = percent_enhancement(pre, early) > threshold
    support = inside & enhancing & background_mask(pre, box)
    return support if support.any() else None


def box_from_projections(subtraction: np.ndarray):
    """Part C's box. No reader input of any kind.

    Each projection bounds the two axes it keeps; the three sets of bounds are
    intersected, so each axis is bounded twice. The component kept on each
    projection is the largest one, which is why this lands on the contralateral
    breast when that side enhances more. That failure is the result, not a bug
    to be engineered away.
    """
    lower = [None, None, None]
    upper = [None, None, None]
    for axis in range(3):
        kept = [a for a in range(3) if a != axis]
        flat = subtraction.max(axis=axis)
        seed = flat > np.percentile(flat, PROJECTION_PERCENTILE)
        if not seed.any():
            return None
        labels, count = ndimage.label(seed, structure=np.array(
            [[0, 1, 0], [1, 1, 1], [0, 1, 0]]))          # 4-connectivity
        sizes = ndimage.sum(seed, labels, range(1, count + 1))
        where = np.argwhere(labels == int(np.argmax(sizes)) + 1)
        for i, a in enumerate(kept):
            lo, hi = int(where[:, i].min()), int(where[:, i].max()) + 1
            lower[a] = lo if lower[a] is None else max(lower[a], lo)
            upper[a] = hi if upper[a] is None else min(upper[a], hi)
    if any(u is None or u <= l for l, u in zip(lower, upper)):
        return None
    return tuple(slice(l, u) for l, u in zip(lower, upper))


def automatic_support(pre, early, threshold):
    box = box_from_projections(early - pre)
    if box is None:
        return None
    inside = np.zeros(pre.shape, bool)
    inside[box] = True
    support = (inside & (percent_enhancement(pre, early) > threshold)
               & background_mask(pre, box))
    return support if support.any() else None


def overlap(a: np.ndarray, b: np.ndarray):
    inter = float(np.count_nonzero(a & b))
    sa, sb = float(np.count_nonzero(a)), float(np.count_nonzero(b))
    union = sa + sb - inter
    return (2 * inter / (sa + sb) if sa + sb else float("nan"),
            inter / union if union else float("nan"))


# --------------------------------------------------------------------------


def analyse(pid: str, visit: str) -> dict:
    """One examination, one load, all three parts."""
    import theme3_ftv as ftv
    from theme3_ftv import load_case, volume

    # The published cohort's extractions live under theme4_full_v2, and
    # theme3_ftv binds OUT at import time to the pilot directory. The
    # discovery run patches it the same way in theme4_full_cohort.main.
    ftv.OUT = SOURCE
    pre, early, late, eligible, voxel_cc, values, _v, _u = load_case(pid, visit)
    guard = _Guard(pre=pre, early=early, late=late)
    threshold = stored_threshold(values)
    lower, upper = values["ser"]["lower"], values["ser"]["upper"]

    record = {"patient_id": pid, "visit": visit, "voxel_cc": voxel_cc,
              "supplied_support_voxels": int(np.count_nonzero(eligible)),
              "supplied_support_cc": float(np.count_nonzero(eligible) * voxel_cc),
              "containment": containment(guard["pre"], guard["early"],
                                         eligible, values)}

    supplied_empty = not np.count_nonzero(eligible)
    supports = {}
    if supplied_empty:
        # No box exists, so no support can be rebuilt. The archive measured
        # 0 cc here in both of its own components, and the published analysis
        # carries this examination with a volume of zero; parts B and C do the
        # same, so the examination stays in the paired comparison with the
        # value the published analysis gives it (plan §10, amendment 1).
        supports = {m: None for m in MARGINS}
        supports["automatic"] = None
    else:
        for margin in MARGINS:
            supports[margin] = pe_in_box(guard["pre"], guard["early"],
                                         eligible, threshold, margin)
        supports["automatic"] = automatic_support(guard["pre"],
                                                  guard["early"], threshold)
    guard.unlock()                                  # the SER step may use late

    def arms(support):
        if support is None and supplied_empty:
            return {"status": "complete", "supplied_support_was_empty": True,
                    "support_voxels": 0, "support_cc": 0.0,
                    "dice_vs_supplied": None, "jaccard_vs_supplied": None,
                    "contains_supplied_fraction": None,
                    "real_cc": 0.0, "sub_cc": 0.0, "retained": None}
        if support is None:
            return {"status": "empty support"}
        real = volume(pre, early, late, support, voxel_cc, lower, upper)
        sub = volume(pre, early, early, support, voxel_cc, lower, upper)
        # The support is built from PE above the threshold, so early exceeds
        # precontrast everywhere on it and the substituted ratio is 1, which is
        # above every stored lower bound. Asserted, not assumed.
        assert sub["count"] == int(np.count_nonzero(support)), (
            "%s %s: a support voxel has no early enhancement" % (pid, visit))
        dice, jaccard = overlap(support, eligible)
        return {"status": "complete",
                "support_voxels": int(np.count_nonzero(support)),
                "support_cc": float(np.count_nonzero(support) * voxel_cc),
                "dice_vs_supplied": dice, "jaccard_vs_supplied": jaccard,
                "contains_supplied_fraction":
                    float(np.count_nonzero(support & eligible)
                          / np.count_nonzero(eligible)),
                "real_cc": real["volume_cc"], "sub_cc": sub["volume_cc"],
                "retained": (real["volume_cc"] / sub["volume_cc"]
                             if sub["volume_cc"] else None)}

    record["by_margin"] = {str(m): arms(supports[m]) for m in MARGINS}
    record["automatic"] = arms(supports["automatic"])
    return record


def describe(values) -> dict:
    a = np.asarray([v for v in values if v is not None], dtype=float)
    a = a[np.isfinite(a)]
    if not a.size:
        return {"n": 0}
    return {"n": int(a.size), "min": float(a.min()),
            "q1": float(np.percentile(a, 25)), "median": float(np.median(a)),
            "q3": float(np.percentile(a, 75)), "max": float(a.max()),
            "mean": float(a.mean())}


def paired_comparison(rows, records) -> dict:
    """The published estimator on the published four endpoints."""
    from theme4_pcr_auc import compare

    by_key = {(r["patient_id"], r["visit"]): r for r in records}
    usable = [r for r in rows
              if all(by_key.get((r["patient_id"], v), {}).get("status")
                     == "complete" for v in VISITS)]
    if not usable:
        return {"status": "no participant has both visits"}
    out = {"participants": len(usable),
           "participant_ids": [r["patient_id"] for r in usable]}
    arms = {}
    for arm in ("real_cc", "sub_cc"):
        v0 = np.array([by_key[(r["patient_id"], "T0")][arm] for r in usable])
        v1 = np.array([by_key[(r["patient_id"], "T1")][arm] for r in usable])
        # The published sign convention, from theme4_pcr_auc.main: the visit
        # scores are negative volumes, change is T0 minus T1, and the relative
        # reduction is undefined -- not zero -- when the T0 volume is zero.
        arms[arm] = {"T0": -v0, "T1": -v1, "change": v0 - v1,
                     "relative": np.where(v0 > 0, (v0 - v1)
                                          / np.where(v0 > 0, v0, 1), np.nan)}
    y = np.array([int(r["pcr"] in ("1", "pCR")) for r in usable])
    for endpoint in ENDPOINTS:
        real, sub = arms["real_cc"][endpoint], arms["sub_cc"][endpoint]
        keep = np.isfinite(real) & np.isfinite(sub)
        res = compare(y[keep], real[keep], sub[keep])
        res["excluded"] = [r["patient_id"]
                           for r, k in zip(usable, keep) if not k]
        out[endpoint] = res
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--freeze", action="store_true",
                    help="pin the plan and this script, once")
    ap.add_argument("--limit", type=int, default=0,
                    help="first N participants only (smoke test; the reported "
                         "run uses all of them)")
    ap.add_argument("--smoke", action="store_true",
                    help="execution test before the freeze: runs without a "
                         "lock and writes nothing. Debugging the code is not "
                         "tuning the rule, but the boundary is made explicit "
                         "rather than left to trust.")
    args = ap.parse_args()
    if args.freeze:
        return freeze()

    import csv

    lock = ({"plan_sha256": "smoke run: not frozen",
             "script_sha256": "smoke run: not frozen"} if args.smoke
            else require_freeze())
    rows = list(csv.DictReader(SCORES.open(encoding="utf-8-sig")))
    assert len(rows) == 117, (
        "expected the published 117 participants, got %d" % len(rows))
    if args.limit:
        rows = rows[:args.limit]

    records = []
    for i, r in enumerate(rows, 1):
        for visit in VISITS:
            records.append(analyse(r["patient_id"], visit))
        c = [x["containment"] for x in records[-2:]]
        show = lambda v: "undefined" if v is None else "%.5f" % v
        print("%3d/%d %s  containment %s %s  threshold %g %g"
              % (i, len(rows), r["patient_id"], show(c[0]["recall_spec"]),
                 show(c[1]["recall_spec"]), c[0]["stored_threshold_percent"],
                 c[1]["stored_threshold_percent"]), flush=True)
    assert len(records) == 2 * len(rows), "an examination was skipped"

    # ---- part A ----------------------------------------------------------
    cont = [r["containment"] for r in records]
    measured = [c for c in cont if c["status"] == "measured"]
    total = sum(c["support_voxels"] for c in cont)
    failing = sum(c["failing_voxels_spec"] for c in cont)
    part_a = {
        "examinations": len(cont),
        "examinations_measured": len(measured),
        "examinations_with_an_empty_supplied_support":
            [[c2["patient_id"], c2["visit"]] for c2, c in zip(records, cont)
             if c["status"] == "empty supplied support"],
        "thresholds_used": {str(t): sum(c["stored_threshold_percent"] == t
                                        for c in cont)
                            for t in sorted({c["stored_threshold_percent"]
                                             for c in cont})},
        "recall_spec": describe([c["recall_spec"] for c in cont]),
        "recall_floored": describe([c["recall_floored"] for c in cont]),
        "support_voxels_total": total,
        "failing_voxels_total": failing,
        "failing_voxel_share": failing / total,
        "examinations_below_0_99": sum(c["recall_spec"] < 0.99
                                       for c in measured),
        "worst_examinations": sorted(
            ([c["recall_spec"], r["patient_id"], r["visit"]]
             for r, c in zip(records, cont) if c["status"] == "measured"))[:5],
        "examinations_with_support_equal_to_stored_pe_count":
            sum(c["support_equals_stored_pe_count"] for c in cont),
        "support_voxels_without_precontrast_signal":
            sum(c["support_voxels_without_precontrast_signal"] for c in cont),
        "statement": "every voxel of the supplied support is tested against "
                     "percent enhancement computed from the precontrast and "
                     "early images only, at the threshold the archive stored "
                     "for that examination",
    }

    # ---- part B ----------------------------------------------------------
    part_b = {}
    for margin in MARGINS:
        per_exam = [dict(r["by_margin"][str(margin)],
                         patient_id=r["patient_id"], visit=r["visit"])
                    for r in records]
        done = [m for m in per_exam if m["status"] == "complete"]
        empty = {v: sum(m["status"] != "complete" for m in per_exam
                        if m["visit"] == v) for v in VISITS}
        part_b[str(margin)] = {
            "examinations": len(per_exam), "complete": len(done),
            "empty_support_by_visit": empty,
            "carried_as_zero_because_the_supplied_support_was_empty":
                sum(m.get("supplied_support_was_empty", False)
                    for m in per_exam),
            "dice_vs_supplied": describe([m["dice_vs_supplied"] for m in done]),
            "jaccard_vs_supplied": describe([m["jaccard_vs_supplied"]
                                             for m in done]),
            "contains_supplied_fraction":
                describe([m["contains_supplied_fraction"] for m in done]),
            "support_cc": describe([m["support_cc"] for m in done]),
            "supplied_support_cc": describe([r["supplied_support_cc"]
                                             for r in records]),
            "retained_fraction": describe([m["retained"] for m in done]),
            "comparison": paired_comparison(rows, per_exam),
            "per_examination": done,
        }
        e = part_b[str(margin)]
        print("\nmargin %d voxels | complete %d/%d | empty T0 %d T1 %d | "
              "Dice %.3f [%.3f-%.3f] | support %.1f cc vs supplied %.1f cc"
              % (margin, len(done), len(per_exam), empty["T0"], empty["T1"],
                 e["dice_vs_supplied"]["median"], e["dice_vs_supplied"]["q1"],
                 e["dice_vs_supplied"]["q3"], e["support_cc"]["median"],
                 e["supplied_support_cc"]["median"]), flush=True)
        for endpoint in ENDPOINTS:
            c = e["comparison"].get(endpoint)
            if c and c.get("status") == "estimated":
                print("    %-9s real %.3f sub %.3f  diff %+.4f [%+.4f,%+.4f]"
                      % (endpoint, c["real_auc"], c["sub_auc"],
                         c["difference_sub_minus_real"], *c["difference_ci"]))

    # ---- part C ----------------------------------------------------------
    auto = [dict(r["automatic"], patient_id=r["patient_id"], visit=r["visit"])
            for r in records]
    done = [m for m in auto if m["status"] == "complete"]
    part_c = {
        "examinations": len(auto), "complete": len(done),
        "dice_vs_supplied": describe([m["dice_vs_supplied"] for m in done]),
        "disjoint_from_supplied": sum(m["dice_vs_supplied"] == 0
                                      for m in done),
        "disjoint_examinations": [[m["patient_id"], m["visit"]] for m in done
                                  if m["dice_vs_supplied"] == 0],
        "comparison": paired_comparison(rows, auto),
        "statement": "the box is found from orthogonal projections with no "
                     "reader input; the largest bright component is kept, "
                     "which is the contralateral breast whenever that side "
                     "enhances more",
        "per_examination": done,
    }
    print("\nautomatic redraw | complete %d/%d | Dice %.3f [%.3f-%.3f] | "
          "disjoint from the supplied support in %d examinations"
          % (len(done), len(auto), part_c["dice_vs_supplied"]["median"],
             part_c["dice_vs_supplied"]["q1"], part_c["dice_vs_supplied"]["q3"],
             part_c["disjoint_from_supplied"]), flush=True)

    # ---- falsification, judged in the plan's own words -------------------
    verdicts = {}
    verdicts["F1_the_supplied_support_is_an_early_only_selection"] = {
        "pass": part_a["recall_spec"]["min"] >= F1_CONTAINMENT_FLOOR,
        "threshold": F1_CONTAINMENT_FLOOR,
        "worst_examination_recall": part_a["recall_spec"]["min"],
        "statement": "at its worst examination the supplied support is %.3f "
                     "early-only selectable" % part_a["recall_spec"]["min"],
    }
    flips = [{"margin": m, "endpoint": e,
              "difference": part_b[str(m)]["comparison"][e]
                  ["difference_sub_minus_real"],
              "ci": part_b[str(m)]["comparison"][e]["difference_ci"]}
             for m in MARGINS for e in ENDPOINTS
             if part_b[str(m)]["comparison"].get(e, {}).get("status")
             == "estimated"
             and (part_b[str(m)]["comparison"][e]["difference_ci"][0] > 0
                  or part_b[str(m)]["comparison"][e]["difference_ci"][1] < 0)]
    verdicts["F2_an_interval_excluding_zero_in_the_rebuilt_support"] = {
        "pass": not flips,
        "statement": "at no margin does the AUC difference have an interval "
                     "that excludes zero" if not flips else
                     "an interval excludes zero; the manuscript's conclusion "
                     "would depend on the supplied support",
        "occurrences": flips,
    }
    worst_empty = max(max(part_b[str(m)]["empty_support_by_visit"].values())
                      / len(rows) for m in MARGINS)
    verdicts["F3_empty_support_share"] = {
        "pass": worst_empty <= F2_EMPTY_SHARE,
        "threshold": F2_EMPTY_SHARE, "worst_share": worst_empty,
        "statement": "the rule produces a support for all but %.1f%% of "
                     "examinations at the worst margin" % (100 * worst_empty),
    }
    best_dice = max(part_b[str(m)]["dice_vs_supplied"]["median"]
                    for m in MARGINS)
    verdicts["F4_median_dice"] = {
        "pass": best_dice >= F3_MEDIAN_DICE,
        "threshold": F3_MEDIAN_DICE, "best_median_dice": best_dice,
        "statement": "the rebuilt and supplied supports overlap at a median "
                     "Dice of %.3f at the best margin" % best_dice,
    }

    result = {
        "scope": "post hoc; Referee 1 point 3 and Referee 2 point 2",
        "plan_sha256": lock["plan_sha256"],
        "script_sha256": lock["script_sha256"],
        "participants": len(rows), "margins": list(MARGINS),
        "limits": "an algorithmic reconstruction of the documented rule, not a "
                  "reader study; for non-Philips participants the arrays are "
                  "already the archive's published DCE crop, so parts B and C "
                  "test the support within that crop and not the crop itself; "
                  "OMIT regions and the specification's neighbour operation "
                  "are not reconstructed",
        "containment": part_a,
        "pe_in_box_by_margin": part_b,
        "automatic_redraw": part_c,
        "falsification": verdicts,
    }

    print()
    for name, v in verdicts.items():
        print("%-56s %s  %s" % (name, "PASS" if v["pass"] else "FAIL",
                                v["statement"]))
    if args.smoke:
        print("\nsmoke run: nothing written")
        return 0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "early_only_support.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    print("\nwrote", (OUT / "early_only_support.json").relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
