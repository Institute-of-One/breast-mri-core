"""Does the substitution behave differently by site, vendor or late-phase index?

    python scripts/theme4_strata.py

Referee 1 point 5 of BPEX-107031 asks for the results stratified by acquisition
site, by scanner manufacturer and by which temporal index the archive nominates
as the late phase. Two analyses are reported, because the obvious one is the
weaker one.

**The per-stratum AUC difference is reported, and is close to uninformative.**
The smallest vendor stratum holds 11 participants and the smallest site holds
2, so a stratum AUC is estimated from a handful of pCR participants. Every
stratum carries its own participant counts and the published estimator's
small-class warning, and no stratum difference is interpreted as a finding. No
threshold for stratum dependence was prespecified, and none is claimed.

**The magnitude of the substitution is analysed instead, and that analysis is
properly powered**, because it needs no outcome: the retained fraction
V_real / V_sub is one number per examination, so all 234 examinations
contribute, not just the 36 pCR participants. If the substitution's effect were
concentrated in one vendor or in one phase index, it would show here. The
pairwise reversal rate of section 2.6 is reported per stratum on the same
footing.

Sites are reported individually down to 10 participants and pooled below that;
the rule is fixed here rather than chosen from the counts, and the pooled
stratum is named with its members so the reader can see what went into it.

Writes ``data/derived/theme4_revision_v1/strata.json``. Nothing published is
modified.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import kruskal

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SOURCE = ROOT / "data/derived/theme4_full_v2"
SCORES = ROOT / "data/derived/theme4_pcr_v1/paired_scores.csv"
OUT = ROOT / "data/derived/theme4_revision_v1"

ENDPOINTS = ("T0", "T1", "change", "relative")
VISITS = ("T0", "T1")
SITE_FLOOR = 10          # fixed before looking at the site counts
MANUFACTURER = {"GE MEDICAL SYSTEMS": "GE", "SIEMENS": "Siemens",
                "Philips Medical Systems": "Philips"}


def metadata(pid: str) -> dict:
    """Site, vendor and the archive's own phase indices, per participant."""
    d = json.loads((SOURCE / ("%s_extraction.json" % pid)).read_text(encoding="utf-8"))
    sites = {d["visits"][v]["site"] for v in VISITS}
    indices = {tuple(d["visits"][v]["indices"]) for v in VISITS}
    phases = {d["visits"][v]["phases"] for v in VISITS}
    assert len(sites) == 1, "%s is imaged at more than one site" % pid
    late = sorted({i[-1] for i in indices})
    return {
        "site": sites.pop(),
        "manufacturer": MANUFACTURER[d["manufacturer"]],
        # The late phase is the last of the three indices the archive supplies.
        "late_index": ("index %d" % late[0] if len(late) == 1
                       else "differs between visits (%s)"
                            % ", ".join(str(i) for i in late)),
        "phases": ("%d phases" % sorted(phases)[0] if len(phases) == 1
                   else "differs between visits"),
    }


def exam_measurements() -> list[dict]:
    """One row per examination: the retained fraction, with no outcome attached.

    Two examinations have no retained fraction because both arms measure
    exactly zero -- a complete response on imaging, where the ratio is 0/0 and
    not a small number. They are carried with ``retained=None`` and counted
    where they are dropped, never silently skipped.
    """
    rows = []
    for path in sorted(SOURCE.glob("*_q2.json")):
        r = json.loads(path.read_text(encoding="utf-8"))
        if r["status"] != "complete":
            continue
        retained = r["retention_real_over_sub"]
        assert retained is not None or (r["real_cc"] == 0 and r["sub_cc"] == 0), \
            "%s %s has no retained fraction but a non-zero volume" % (
                r["patient_id"], r["visit"])
        rows.append({"patient_id": r["patient_id"], "visit": r["visit"],
                     "real_cc": r["real_cc"], "sub_cc": r["sub_cc"],
                     "retained": retained})
    assert len(rows) == 234, "expected 234 complete examinations, found %d" % len(rows)
    return rows


def describe(values: np.ndarray) -> dict:
    return {"n": int(values.size), "median": float(np.median(values)),
            "q1": float(np.percentile(values, 25)),
            "q3": float(np.percentile(values, 75)),
            "min": float(values.min()), "max": float(values.max())}


def site_label(site: str, counts: dict) -> str:
    return site if counts[site] >= SITE_FLOOR else "pooled sites below %d" % SITE_FLOOR


def main() -> int:
    from theme4_pcr_auc import compare
    from theme4_rank_preservation import pair_table

    scores = list(csv.DictReader(SCORES.open(encoding="utf-8-sig")))
    meta = {r["patient_id"]: metadata(r["patient_id"]) for r in scores}

    site_counts = {}
    for m in meta.values():
        site_counts[m["site"]] = site_counts.get(m["site"], 0) + 1
    pooled = sorted(s for s, n in site_counts.items() if n < SITE_FLOOR)

    def key_of(pid: str, factor: str) -> str:
        return (site_label(meta[pid]["site"], site_counts) if factor == "site"
                else meta[pid][factor])

    factors = ("site", "manufacturer", "late_index", "phases")
    result = {
        "scope": "post hoc exploratory stratification; Referee 1 point 5",
        "site_rule": "reported individually at %d participants or more, pooled "
                     "below that; the floor was fixed before the counts were "
                     "read" % SITE_FLOOR,
        "pooled_sites": pooled,
        "caution": "a stratum AUC rests on that stratum's pCR participants "
                   "alone, so the smallest strata are not informative about the "
                   "substitution; no stratum difference is claimed, and no "
                   "threshold for one was prespecified. The retained-fraction "
                   "analysis beside it uses every examination and no outcome.",
        "sites_in_cohort": sorted(site_counts),
        "discrimination": {}, "magnitude": {},
    }

    # --- Discrimination, by stratum. The weak analysis, reported as weak. ---
    for factor in factors:
        groups = {}
        for r in scores:
            groups.setdefault(key_of(r["patient_id"], factor), []).append(r)
        entry = {}
        for name, subset in sorted(groups.items()):
            label = np.array([1 if r["pcr"] in ("1", "pCR") else 0 for r in subset])
            per_endpoint = {}
            for endpoint in ENDPOINTS:
                real = np.array([float(r["real_" + endpoint]) for r in subset])
                sub = np.array([float(r["sub_" + endpoint]) for r in subset])
                per_endpoint[endpoint] = compare(label, real, sub)
            entry[name] = per_endpoint
        result["discrimination"][factor] = entry

    # --- Magnitude, by stratum. Every examination, no outcome. ---
    exams = exam_measurements()
    for factor in factors:
        entry = {}
        for visit in VISITS:
            groups = {}
            for e in exams:
                if e["visit"] != visit or e["patient_id"] not in meta:
                    continue
                groups.setdefault(key_of(e["patient_id"], factor), []).append(e)
            per_stratum, samples, counted = {}, [], 0
            for name, subset in sorted(groups.items()):
                defined = [e for e in subset if e["retained"] is not None]
                counted += len(subset)
                retained = np.array([e["retained"] for e in defined])
                real = np.array([e["real_cc"] for e in subset])
                sub = np.array([e["sub_cc"] for e in subset])
                keep = real > 0
                pairs = (pair_table(np.log(sub[keep]), np.log(real[keep] / sub[keep]))
                         if keep.sum() > 1 else None)
                per_stratum[name] = {
                    "retained_fraction": describe(retained),
                    "examinations_without_a_retained_fraction":
                        len(subset) - len(defined),
                    "examinations_with_a_zero_measured_volume": int((~keep).sum()),
                    "pairs": pairs["pairs"] if pairs else 0,
                    "reversed": pairs["reversed"] if pairs else 0,
                    "fraction_reversed": pairs["fraction_reversed"] if pairs else None,
                }
                samples.append(retained)
            assert counted == len(scores), \
                "%s/%s covered %d of %d participants" % (factor, visit, counted,
                                                         len(scores))
            test = (kruskal(*samples) if len(samples) > 1
                    and all(s.size > 0 for s in samples) else None)
            entry[visit] = {
                "strata": per_stratum,
                "kruskal_wallis": ({"statistic": float(test.statistic),
                                    "p_value": float(test.pvalue),
                                    "groups": len(samples),
                                    "note": "across strata, on the retained "
                                            "fraction; descriptive, not a "
                                            "prespecified test"}
                                   if test else None),
            }
        result["magnitude"][factor] = entry

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "strata.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    for factor in factors:
        print("=== %s" % factor)
        for name in sorted(result["discrimination"][factor]):
            d = result["discrimination"][factor][name]["T1"]
            m0 = result["magnitude"][factor]["T0"]["strata"][name]
            m1 = result["magnitude"][factor]["T1"]["strata"][name]
            warn = " (small class)" if d.get("small_class_warning") else ""
            diff = ("%+.4f [%+.4f,%+.4f]" % (d["difference_sub_minus_real"],
                                             *d["difference_ci"])
                    if d["status"] == "estimated" else d["status"])
            print("  %-28s n=%3d pCR=%2d | T1 AUC diff %s%s"
                  % (name, d["n"], d.get("pcr", 0), diff, warn))
            print("      retained fraction  T0 %.3f [%.3f-%.3f]   T1 %.3f [%.3f-%.3f]"
                  % (m0["retained_fraction"]["median"], m0["retained_fraction"]["q1"],
                     m0["retained_fraction"]["q3"], m1["retained_fraction"]["median"],
                     m1["retained_fraction"]["q1"], m1["retained_fraction"]["q3"]))
            print("      reversed pairs     T0 %s   T1 %s"
                  % ("%.2f%%" % (100 * m0["fraction_reversed"])
                     if m0["fraction_reversed"] is not None else "n/a",
                     "%.2f%%" % (100 * m1["fraction_reversed"])
                     if m1["fraction_reversed"] is not None else "n/a"))
        for visit in VISITS:
            t = result["magnitude"][factor][visit]["kruskal_wallis"]
            if t:
                print("  %s retained fraction across strata: H=%.3f p=%.4f"
                      % (visit, t["statistic"], t["p_value"]))
        print()

    print("wrote", (OUT / "strata.json").relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
