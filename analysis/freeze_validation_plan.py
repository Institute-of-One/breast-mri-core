"""Freeze the prespecified validation plan before the validation is run.

    python scripts/freeze_validation_plan.py           # write the lock
    python scripts/freeze_validation_plan.py --check    # refuse if anything moved

Referee 1 asked for an independent cohort analysed under a prespecified plan.
A plan is only prespecified if it can be shown to have existed first. This
derives the cohort from the rule in the plan rather than from a hand-written
list -- so the rule and the list cannot disagree -- and records the sha256 of
the plan, the derived list, the seed and every script that will produce a
number, in one lock file.

The acquisition script refuses to start unless ``--check`` passes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parent.parent
XLSX = ROOT / "data/raw/acrin6698/Full-Collection-Ancillary-Patient-Information-file.xlsx"
SHEET = "ACRIN_6698_Patient_Cohorts_Clin"
PLAN = ROOT / "docs/theme4_validation_plan_v1.md"
COHORT = ROOT / "docs/theme4_validation_cohort_v1.json"
LOCK = ROOT / "docs/theme4_validation_lock_v1.json"

# Fixed in the plan. Changing either here is changing the plan.
SEED = 20261111
BOOTSTRAP = 10000

# Every script that will turn an image into a number in this validation.
PIPELINE = [
    "scripts/theme4_full_cohort.py",
    "scripts/theme3_ftv.py",
    "scripts/theme3_restore_philips.py",
    "scripts/theme3_extract_development.py",
]

# Section 6 of the plan, so that the criteria are in the lock and not only in
# prose that could be softened after the result is known.
DISCOVERY = {
    "relative": {"auc_real": 0.6461, "auc_sub": 0.6529, "difference": 0.0069,
                 "ci": [-0.043, 0.056], "fraction_pairs_changed": 0.1200,
                 "spearman": 0.9211},
    "T0": {"auc_real": 0.6591, "auc_sub": 0.6722, "difference": 0.0130,
           "ci": [-0.026, 0.051], "fraction_pairs_changed": 0.0844,
           "spearman": 0.9459},
    "T1": {"auc_real": 0.7126, "auc_sub": 0.7030, "difference": -0.0096,
           "ci": [-0.068, 0.043], "fraction_pairs_changed": 0.1221,
           "spearman": 0.9139},
}
# Plan amendment 1(b). The frozen pipeline selects the Philips native-grid
# recovery by an exact manufacturer string, and the validation cohort contains
# a Philips spelling the published cohort does not. The map is pinned here so
# that it cannot be widened after a result is seen; anything outside it stops
# the run rather than defaulting to a processing path.
MANUFACTURER_MAP = {
    "Philips Healthcare": "Philips Medical Systems",
    "Philips Medical Systems": "Philips Medical Systems",
    "GE MEDICAL SYSTEMS": "GE MEDICAL SYSTEMS",
    "SIEMENS": "SIEMENS",
}

CRITERIA = {
    "C1_primary_difference_within_discovery_ci": {"endpoint": "relative",
                                                  "ci": [-0.043, 0.056]},
    "C2_difference_magnitude": {"endpoints": ["relative", "T0", "T1"], "max_abs": 0.05},
    "C3_pairs_changed_band": {"endpoints": ["relative", "T0", "T1"],
                              "tolerance_percentage_points": 5.0},
    "C4_volume_inflation": {"min_fraction_ge_10pct": 0.90, "visits": ["T0", "T1"]},
    "C5_spearman": {"endpoints": ["relative", "T0", "T1"], "min": 0.85},
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def derive_cohort() -> dict:
    """The rule of plan section 2, applied to the ancillary file."""
    import openpyxl

    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    rows = list(wb[SHEET].iter_rows(values_only=True))
    head = [str(h) for h in rows[0]]
    records = [dict(zip(head, r)) for r in rows[1:]]

    def flag(value) -> bool:
        return str(value).strip() in ("1", "1.0", "True")

    selected, counts = [], {"total": len(records)}
    for r in records:
        if not (flag(r["T0"]) and flag(r["T1"])):
            continue
        counts["both_visits"] = counts.get("both_visits", 0) + 1
        if flag(r["BMMR2_TRAIN"]):
            continue
        counts["outside_bmmr2_training"] = counts.get("outside_bmmr2_training", 0) + 1
        if str(r["pcr"]).strip() not in ("pCR", "Non-pCR"):
            continue
        counts["outcome_known"] = counts.get("outcome_known", 0) + 1
        if str(r["QC_PASS_FAIL_T0"]).strip() != "PASS":
            continue
        if str(r["QC_PASS_FAIL_T1"]).strip() != "PASS":
            continue
        counts["qc_pass_both"] = counts.get("qc_pass_both", 0) + 1
        if not flag(r["analysis cohort"]) or not flag(r["eligigible"]):
            continue
        counts["analysis_cohort_and_eligible"] = counts.get(
            "analysis_cohort_and_eligible", 0) + 1
        selected.append(r)

    patients = sorted(str(r["TCIA PATIENT ID"]) for r in selected)
    assert len(patients) == len(set(patients)), "duplicate patient identifier"
    pcr = sorted(str(r["TCIA PATIENT ID"]) for r in selected
                 if str(r["pcr"]).strip() == "pCR")
    return {
        "rule": "plan section 2: T0 and T1 present; not BMMR2_TRAIN; pcr in "
                "{pCR, Non-pCR}; QC_PASS_FAIL_T0 and _T1 both PASS; analysis "
                "cohort = 1; eligigible = 1",
        "source_file": XLSX.name,
        "source_sheet": SHEET,
        "funnel": counts,
        "n": len(patients),
        "n_pcr": len(pcr),
        "n_non_pcr": len(patients) - len(pcr),
        "pairs": len(pcr) * (len(patients) - len(pcr)),
        "patients": patients,
    }


def build_lock(cohort: dict) -> dict:
    return {
        "plan": str(PLAN.relative_to(ROOT)).replace("\\", "/"),
        "plan_sha256": sha256_file(PLAN),
        "cohort_file": str(COHORT.relative_to(ROOT)).replace("\\", "/"),
        "cohort_sha256": sha256_file(COHORT),
        "n": cohort["n"], "n_pcr": cohort["n_pcr"],
        "n_non_pcr": cohort["n_non_pcr"], "pairs": cohort["pairs"],
        "seed": SEED, "bootstrap_samples": BOOTSTRAP,
        "pipeline_sha256": {p: sha256_file(ROOT / p) for p in PIPELINE},
        "ancillary_file_sha256": sha256_file(XLSX),
        "discovery": DISCOVERY,
        "criteria": CRITERIA,
        "manufacturer_map": MANUFACTURER_MAP,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="verify the lock instead of writing it")
    args = ap.parse_args()

    if not PLAN.exists():
        print("FAIL the plan is missing: %s" % PLAN)
        return 1

    cohort = derive_cohort()
    if args.check:
        if not LOCK.exists() or not COHORT.exists():
            print("FAIL no lock: run without --check before acquiring any image")
            return 1
        on_disk = json.loads(COHORT.read_text(encoding="utf-8"))
        if on_disk["patients"] != cohort["patients"]:
            print("FAIL the rule no longer yields the frozen patient list "
                  "(%d frozen, %d now)" % (len(on_disk["patients"]), cohort["n"]))
            return 1
        lock = json.loads(LOCK.read_text(encoding="utf-8"))
        want = build_lock(cohort)
        bad = [k for k in ("plan_sha256", "cohort_sha256", "pipeline_sha256",
                           "ancillary_file_sha256", "seed", "bootstrap_samples",
                           "discovery", "criteria", "manufacturer_map")
               if lock.get(k) != want[k]]
        if bad:
            print("FAIL the frozen plan has moved: %s" % ", ".join(bad))
            for k in bad:
                if isinstance(want[k], dict):
                    for kk in want[k]:
                        if lock.get(k, {}).get(kk) != want[k][kk]:
                            print("   %s/%s changed" % (k, kk))
            return 1
        print("OK   plan, cohort (%d patients, %d pairs), seed %d and %d pipeline "
              "scripts are unchanged since the freeze"
              % (lock["n"], lock["pairs"], lock["seed"], len(lock["pipeline_sha256"])))
        return 0

    if LOCK.exists():
        print("FAIL a lock already exists. A prespecified plan is frozen once; "
              "to change it, append the change and its date to the plan and say "
              "so in the manuscript.")
        return 1
    COHORT.write_text(json.dumps(cohort, indent=2) + "\n", encoding="utf-8")
    LOCK.write_text(json.dumps(build_lock(cohort), indent=2) + "\n", encoding="utf-8")
    print("froze %d patients (%d pCR, %d non-pCR, %d pairs), seed %d"
          % (cohort["n"], cohort["n_pcr"], cohort["n_non_pcr"], cohort["pairs"], SEED))
    print("plan   sha256 %s" % sha256_file(PLAN))
    print("cohort sha256 %s" % sha256_file(COHORT))
    print("wrote %s" % LOCK.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
