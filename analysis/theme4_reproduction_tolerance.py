"""Why seven examinations miss the reproduction tolerance, and what happens without them.

    python scripts/theme4_reproduction_tolerance.py

Referee 1 point 4 of BPEX-107031 asks how large the seven reproduction
mismatches are, what causes them, and whether the result survives their
removal. Three things are computed.

**1. The criterion and which leg of it fails.** ``theme4_full_cohort.tolerance``
calls an examination reproduced when, for *both* the PE and the SER component,

    |calculated - stored| <= 0.01 cc   **and**   |calculated - stored| / stored <= 0.001

The two legs are conjunctive and neither is scale-free: 0.01 cc is about 23
voxels at this resolution, so for a 78 cc volume the absolute leg is a
single-voxel criterion, and for a 0.36 cc volume the relative leg is one. Each
failure is therefore attributed to the leg that failed, and the voxel-count
difference is reported beside the volume difference.

**2. Whether the eligible set reproduced.** The SER volume is counted over the
voxels of the fixed reference mask whose early percent enhancement exceeds the
archive's threshold. When that eligible count reproduces exactly and only the
SER count differs, the difference is voxels sitting on the SER cut-off. When
the eligible count itself differs, the enhancement map differs, which is a
larger claim and is reported as one.

**3. The sensitivity analysis.** The endpoint comparison does not read the
stored volume: both arms are computed from the images, and the stored value is
a check on the pipeline, not an input. A reproduction mismatch therefore cannot
propagate into the AUC arithmetic, and what the sensitivity analysis asks is
the other question -- whether the six participants whose images my pipeline
reproduces least faithfully are carrying the result. The published comparison
is rerun on the remaining 111 by calling ``theme4_pcr_auc.compare`` itself, so
the estimator, the seed and the 10 000 replicates are the published ones.

Writes ``data/derived/theme4_revision_v1/reproduction_tolerance.json``. Nothing
published is modified.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SOURCE = ROOT / "data/derived/theme4_full_v2"
SCORES = ROOT / "data/derived/theme4_pcr_v1/paired_scores.csv"
OUT = ROOT / "data/derived/theme4_revision_v1"

ABSOLUTE_CC = 0.01
RELATIVE = 0.001
COMPONENTS = ("pe", "ser")
ENDPOINTS = ("T0", "T1", "change", "relative")


def classify(report: dict) -> dict:
    """Attribute one examination's failure to a leg and to a stage."""
    legs, detail = set(), {}
    for key in COMPONENTS:
        stored, calc = report["stored"][key], report["calculated"][key]
        delta = calc["volume_cc"] - stored["volume_cc"]
        rel = delta / stored["volume_cc"] if stored["volume_cc"] else float("nan")
        fails_abs = abs(delta) > ABSOLUTE_CC
        fails_rel = abs(rel) > RELATIVE
        if fails_abs:
            legs.add("absolute")
        if fails_rel:
            legs.add("relative")
        detail[key] = {
            "stored_count": stored["count"], "calculated_count": calc["count"],
            "count_difference": calc["count"] - stored["count"],
            "stored_cc": stored["volume_cc"], "calculated_cc": calc["volume_cc"],
            "difference_cc": delta, "relative_difference": rel,
            "fails_absolute_leg": fails_abs, "fails_relative_leg": fails_rel,
        }
    eligible = report["calculated"]["ser"]["eligible_voxels"]
    stage = ("the eligible set differs: the early-enhancement decision changed"
             if detail["pe"]["count_difference"] != 0 else
             "the eligible set reproduces exactly; the difference is voxels on "
             "the SER cut-off")
    voxel_cc = (detail["pe"]["stored_cc"] / detail["pe"]["stored_count"]
                if detail["pe"]["stored_count"] else float("nan"))
    return {
        "patient_id": report["patient_id"], "visit": report["visit"],
        "precision": report["precision"],
        "failed_legs": sorted(legs), "stage": stage,
        "eligible_voxels": eligible,
        "voxel_volume_cc": voxel_cc,
        "absolute_leg_in_voxels": ABSOLUTE_CC / voxel_cc if voxel_cc else float("nan"),
        "zero_denominator_in_mask": report["calculated"]["ser"]["zero_denominator_in_mask"],
        "negative_denominator_in_mask": report["calculated"]["ser"]["negative_denominator_in_mask"],
        "components": detail,
    }


def census() -> tuple[list[dict], dict]:
    """Every examination's worst relative difference, and the failures in full."""
    worst, failures, zero_reference = [], [], []
    for path in sorted(SOURCE.glob("*_q2.json")):
        report = json.loads(path.read_text(encoding="utf-8"))
        if report["status"] != "complete":
            continue
        if not report.get("within_tolerance"):
            failures.append(classify(report))
        ratios = [abs(report["calculated"][k]["volume_cc"]
                      - report["stored"][k]["volume_cc"]) / report["stored"][k]["volume_cc"]
                  for k in COMPONENTS if report["stored"][k]["volume_cc"]]
        if not ratios:
            # The archive stores zero: a complete response on imaging. The
            # criterion for a zero reference is an exact zero, and both of
            # these recompute to exactly zero, so no ratio exists to report.
            zero_reference.append({"patient_id": report["patient_id"],
                                   "visit": report["visit"],
                                   "calculated_cc": [report["calculated"][k]["volume_cc"]
                                                     for k in COMPONENTS]})
            continue
        worst.append(max(ratios))
    a = np.array(worst)
    spread = {
        "examinations_with_a_non_zero_reference": int(a.size),
        "examinations_with_a_zero_reference": zero_reference,
        "bit_identical": int((a == 0).sum()),
        "median": float(np.median(a)),
        "p90": float(np.percentile(a, 90)),
        "p99": float(np.percentile(a, 99)),
        "maximum": float(a.max()),
        "exceeding_0.1_percent": int((a > 1e-3).sum()),
        "exceeding_1_percent": int((a > 1e-2).sum()),
    }
    return failures, spread


def main() -> int:
    from theme4_pcr_auc import compare

    failures, spread = census()
    excluded = sorted({f["patient_id"] for f in failures})

    rows = list(csv.DictReader(SCORES.open(encoding="utf-8-sig")))
    keep = [r for r in rows if r["patient_id"] not in excluded]
    assert len(rows) - len(keep) == len(excluded), "an excluded id is not in the scores"

    published = {r["metric"]: r for r in json.loads(
        (ROOT / "data/derived/theme4_pcr_v1/results.json").read_text(encoding="utf-8")
    )["results"] if r["group"] == "all"}

    comparisons = {}
    for name, subset in (("all", rows), ("reproduced", keep)):
        label = np.array([1 if r["pcr"] in ("1", "pCR") else 0 for r in subset])
        entry = {}
        for endpoint in ENDPOINTS:
            real = np.array([float(r["real_" + endpoint]) for r in subset])
            sub = np.array([float(r["sub_" + endpoint]) for r in subset])
            entry[endpoint] = compare(label, real, sub)
            if name == "all":
                # The whole-cohort arm must return the published numbers, or the
                # subset arm beside it is not the published estimator either.
                p = published[endpoint]
                for field in ("real_auc", "sub_auc", "difference_sub_minus_real"):
                    assert entry[endpoint][field] == p[field], (endpoint, field)
                assert entry[endpoint]["difference_ci"] == p["difference_ci"], endpoint
        comparisons[name] = entry

    result = {
        "scope": "post hoc; Referee 1 point 4",
        "criterion": {
            "absolute_cc": ABSOLUTE_CC, "relative": RELATIVE,
            "form": "both components must satisfy both legs",
            "note": "neither leg is scale-free, so a failure is reported against "
                    "the leg that failed",
        },
        "note": "the endpoint comparison never reads the stored volume: both arms "
                "are computed from the images, so a reproduction mismatch cannot "
                "enter the AUC arithmetic. The sensitivity analysis asks instead "
                "whether these six participants carry the result.",
        "spread_over_all_examinations": spread,
        "failures": failures,
        "excluded_patients": excluded,
        "comparisons": comparisons,
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "reproduction_tolerance.json").write_text(json.dumps(result, indent=2),
                                                     encoding="utf-8")

    print("reproduction difference over the %d examinations with a non-zero "
          "reference (worst component):"
          % spread["examinations_with_a_non_zero_reference"])
    print("  bit-identical %d | median %.2e | p90 %.2e | p99 %.2e | max %.2e"
          % (spread["bit_identical"], spread["median"], spread["p90"],
             spread["p99"], spread["maximum"]))
    print("  above 0.1%%: %d examinations | above 1%%: %d"
          % (spread["exceeding_0.1_percent"], spread["exceeding_1_percent"]))
    print("  %d further examinations store a zero reference and recompute to "
          "exactly zero" % len(spread["examinations_with_a_zero_reference"]))
    print()
    print("the %d examinations outside tolerance, in %d participants:"
          % (len(failures), len(excluded)))
    for f in failures:
        s = f["components"]["ser"]
        print("  %s %s  %-26s ser %+d voxels, %+.5f cc (%+.4f%%)  legs: %s"
              % (f["patient_id"][-6:], f["visit"], f["precision"],
                 s["count_difference"], s["difference_cc"],
                 100 * s["relative_difference"], ", ".join(f["failed_legs"])))
        print("      %s   (0.01 cc = %.0f voxels here)"
              % (f["stage"], f["absolute_leg_in_voxels"]))
    print()
    for endpoint in ENDPOINTS:
        a, b = comparisons["all"][endpoint], comparisons["reproduced"][endpoint]
        print("%-9s n=%3d diff %+.4f [%+.4f,%+.4f]   ->  n=%3d diff %+.4f [%+.4f,%+.4f]"
              % (endpoint, a["n"], a["difference_sub_minus_real"], *a["difference_ci"],
                 b["n"], b["difference_sub_minus_real"], *b["difference_ci"]))
    print()
    print("wrote", (OUT / "reproduction_tolerance.json").relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
