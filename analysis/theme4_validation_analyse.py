"""The prespecified independent validation, run once.

    python scripts/theme4_validation_analyse.py --root D:/DevGit/breast-mri-validation

Everything this script does was fixed in
``docs/theme4_validation_plan_v1.md`` before any volume existed, and the plan,
the cohort list, the four pipeline scripts and the seed are hash-pinned in
``docs/theme4_validation_lock_v1.json``. The script refuses to run if any of
them moved.

Order matters and is enforced:

1. verify the freeze;
2. read the volumes, which carry no outcome;
3. **then** join the pCR labels, for the first time in the whole exercise;
4. estimate, with seed 20261111 and 10 000 replicates;
5. judge C1 to C5, each separately, each in the plan's own words.

**The estimator is the published one, and that is proved rather than claimed.**
``theme4_pcr_auc.compare`` hard-codes the discovery seed, so the bootstrap is
reimplemented here with the seed as an argument -- and
``check_reproduces_published`` runs the reimplementation on the discovery
cohort with the discovery seed and asserts it returns the published AUCs and
intervals exactly. If it does not, this script stops before touching the
validation data.

Writes ``<root>/data/derived/theme4_validation_v1/validation_results.json``.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

COHORT = PROJECT / "docs/theme4_validation_cohort_v1.json"
LOCK = PROJECT / "docs/theme4_validation_lock_v1.json"
ANCILLARY = PROJECT / ("data/raw/acrin6698/"
                       "Full-Collection-Ancillary-Patient-Information-file.xlsx")
SHEET = "ACRIN_6698_Patient_Cohorts_Clin"
DISCOVERY = PROJECT / "data/derived/theme4_pcr_v1"
VISITS = ("T0", "T1")
ENDPOINTS = ("relative", "T0", "T1", "change")
PRIMARY = "relative"

# Section 6 of the plan. Discovery values and the thresholds, fixed in advance.
DISCOVERY_VALUES = {
    "relative": {"real": 0.6461, "sub": 0.6529, "delta": 0.0069,
                 "ci": [-0.043, 0.056], "reordered": 0.1200, "spearman": 0.9211},
    "T0": {"real": 0.6591, "sub": 0.6722, "delta": 0.0130,
           "ci": [-0.026, 0.051], "reordered": 0.0844, "spearman": 0.9459},
    "T1": {"real": 0.7126, "sub": 0.7030, "delta": -0.0096,
           "ci": [-0.068, 0.043], "reordered": 0.1221, "spearman": 0.9139},
}
C2_LIMIT = 0.05
C3_POINTS = 0.05
C4_SHARE = 0.90
C4_INFLATION = 0.10
C5_SPEARMAN = 0.85


def require_freeze() -> dict:
    p = subprocess.run([sys.executable, "-B",
                        str(PROJECT / "scripts/freeze_validation_plan.py"), "--check"],
                       cwd=PROJECT, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    sys.stdout.write(p.stdout or "")
    if p.returncode != 0:
        sys.stdout.write(p.stderr or "")
        raise SystemExit("refusing to analyse: the frozen plan does not verify")
    return json.loads(LOCK.read_text(encoding="utf-8"))


def auc_boot(score: np.ndarray, pos: np.ndarray, neg: np.ndarray) -> np.ndarray:
    """Vectorised bootstrap AUC. Transcribed from theme4_pcr_auc.auc_boot."""
    values = np.concatenate([score[pos], score[neg]], axis=1)
    ranks = rankdata(values, axis=1, method="average")
    npos, nneg = pos.shape[1], neg.shape[1]
    return (ranks[:, :npos].sum(axis=1) - npos * (npos + 1) / 2) / (npos * nneg)


def compare(y: np.ndarray, real: np.ndarray, sub: np.ndarray, seed: int,
            replicates: int = 10000) -> dict:
    """theme4_pcr_auc.compare, with the seed as an argument. Same checks."""
    positives = np.flatnonzero(y == 1)
    negatives = np.flatnonzero(y == 0)
    base = {"n": int(len(y)), "pcr": int(len(positives)),
            "non_pcr": int(len(negatives)),
            "small_class_warning": bool(min(len(positives), len(negatives)) < 10)}
    if not len(positives) or not len(negatives):
        return dict(base, status="not_estimable")
    ar, ass = roc_auc_score(y, real), roc_auc_score(y, sub)
    for s, a in ((real, ar), (sub, ass)):
        d = s[positives, None] - s[negatives]
        assert np.isclose(np.mean((d > 0) + .5 * (d == 0)), a)
    rng = np.random.default_rng(seed)
    pos = rng.choice(positives, (replicates, len(positives)), replace=True)
    neg = rng.choice(negatives, (replicates, len(negatives)), replace=True)
    br, bs = auc_boot(real, pos, neg), auc_boot(sub, pos, neg)
    for i in range(3):
        ix = np.r_[pos[i], neg[i]]
        assert np.isclose(br[i], roc_auc_score(y[ix], real[ix]))
        assert np.isclose(bs[i], roc_auc_score(y[ix], sub[ix]))
    return dict(base, status="estimated", real_auc=float(ar), sub_auc=float(ass),
                real_ci=np.quantile(br, [.025, .975]).tolist(),
                sub_ci=np.quantile(bs, [.025, .975]).tolist(),
                difference_sub_minus_real=float(ass - ar),
                difference_ci=np.quantile(bs - br, [.025, .975]).tolist())


def check_reproduces_published() -> None:
    """The reimplementation must return the published numbers, or stop."""
    published = {r["metric"]: r for r in json.loads(
        (DISCOVERY / "results.json").read_text(encoding="utf-8"))["results"]
        if r["group"] == "all"}
    rows = list((DISCOVERY / "paired_scores.csv").open(encoding="utf-8-sig"))
    rows = list(csv.DictReader(rows))
    y = np.array([1 if r["pcr"] in ("1", "pCR") else 0 for r in rows])
    for metric in ("T0", "T1", "relative", "change"):
        got = compare(y, np.array([float(r["real_" + metric]) for r in rows]),
                      np.array([float(r["sub_" + metric]) for r in rows]),
                      seed=20260918)
        want = published[metric]
        for field in ("real_auc", "sub_auc", "difference_sub_minus_real"):
            assert got[field] == want[field], (metric, field, got[field], want[field])
        assert got["difference_ci"] == want["difference_ci"], metric
    print("estimator: reproduces the published discovery results exactly")


def pair_decomposition(y: np.ndarray, real: np.ndarray,
                       sub: np.ndarray) -> dict:
    """Pairs changed, gained and lost, as section 7 of the plan requires."""
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    cr = np.sign(real[pos, None] - real[neg]) * .5 + .5
    cs = np.sign(sub[pos, None] - sub[neg]) * .5 + .5
    changed = cr != cs
    return {"pairs": int(cr.size), "changed": int(changed.sum()),
            "fraction_changed": float(changed.sum() / cr.size),
            "gained": int(((cs > cr)).sum()), "lost": int(((cs < cr)).sum()),
            "net_concordance": float((cs - cr).sum()),
            "tied_real": int((cr == .5).sum()), "tied_sub": int((cs == .5).sum())}


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def read_volumes(vroot: Path) -> tuple[list[dict], dict]:
    """Volumes only. No outcome is read in this function."""
    out = vroot / "data/derived/theme4_validation_v1"
    reports = {}
    for path in sorted(out.glob("*_q2.json")):
        r = json.loads(path.read_text(encoding="utf-8"))
        if r["status"] == "complete":
            reports[(r["patient_id"], r["visit"])] = r
    plan = json.loads((out / "validation_plan.json").read_text(encoding="utf-8"))
    rows, excluded = [], []
    for record in plan["records"]:
        pid = record["patient_id"]
        if not all((pid, v) in reports for v in VISITS):
            missing = [v for v in VISITS if (pid, v) not in reports]
            excluded.append({"patient_id": pid, "visits_missing": missing,
                             "manufacturer": record["manufacturer"]})
            continue
        a, b = reports[(pid, "T0")], reports[(pid, "T1")]
        row = {"patient_id": pid, "manufacturer": record["manufacturer"],
               "site": record["visits"]["T1"]["site"],
               "within_tolerance": bool(a.get("within_tolerance")
                                        and b.get("within_tolerance"))}
        for arm in ("real", "sub"):
            t0, t1 = a[arm + "_cc"], b[arm + "_cc"]
            row[arm + "_T0"] = -t0
            row[arm + "_T1"] = -t1
            row[arm + "_change"] = t0 - t1
            row[arm + "_relative"] = (t0 - t1) / t0 if t0 else 0.0
            row[arm + "_V0"], row[arm + "_V1"] = t0, t1
        rows.append(row)
    return rows, {"planned": plan["planned"], "prepared": plan["prepared"],
                  "complete_both_visits": len(rows), "excluded": excluded}


def join_labels(rows: list[dict]) -> dict:
    """The first and only place an outcome enters this exercise.

    The sheet is read exactly as ``freeze_validation_plan.derive_cohort``
    reads it -- same workbook, same sheet, same header handling, same
    identifier column -- so the labels join to the cohort the rule selected
    and not to some other reading of the same file.
    """
    import openpyxl

    book = openpyxl.load_workbook(ANCILLARY, read_only=True, data_only=True)
    sheet = list(book[SHEET].iter_rows(values_only=True))
    book.close()
    header = [str(h) for h in sheet[0]]
    table = {str(r["TCIA PATIENT ID"]): r
             for r in (dict(zip(header, values)) for values in sheet[1:])
             if r.get("TCIA PATIENT ID") is not None}

    labelled, missing = 0, []
    for row in rows:
        entry = table.get(row["patient_id"])
        if entry is None or str(entry.get("pcr")).strip() not in ("pCR", "Non-pCR"):
            missing.append(row["patient_id"])
            continue
        row["pcr"] = 1 if str(entry.get("pcr")).strip() == "pCR" else 0
        row["subtype"] = str(entry.get("hrher4g"))
        row["sbr_grade"] = str(entry.get("SBRgrade"))
        labelled += 1
    assert not missing, ("the frozen cohort rule required a pCR label, so every "
                         "patient must carry one: %s" % missing)
    return {"labelled": labelled,
            "identifier_column": "TCIA PATIENT ID",
            "source": ANCILLARY.name, "sheet": SHEET}


def judge(results: dict, decomposition: dict, volumes: dict,
          spearmans: dict) -> dict:
    """C1 to C5, separately, in the plan's own words."""
    verdicts = {}

    d = results[PRIMARY]["difference_sub_minus_real"]
    lo, hi = DISCOVERY_VALUES[PRIMARY]["ci"]
    verdicts["C1"] = {
        "statement": "the validation delta for relative reduction lies within "
                     "the discovery 95% CI, -0.043 to 0.056",
        "value": d, "interval": [lo, hi], "pass": bool(lo <= d <= hi)}

    worst = max(abs(results[e]["difference_sub_minus_real"])
                for e in DISCOVERY_VALUES)
    verdicts["C2"] = {
        "statement": "|delta| <= 0.05 for all three endpoints",
        "largest_absolute_difference": worst, "limit": C2_LIMIT,
        "per_endpoint": {e: results[e]["difference_sub_minus_real"]
                         for e in DISCOVERY_VALUES},
        "pass": bool(worst <= C2_LIMIT)}

    gaps = {e: decomposition[e]["fraction_changed"] - DISCOVERY_VALUES[e]["reordered"]
            for e in DISCOVERY_VALUES}
    verdicts["C3"] = {
        "statement": "the fraction of pCR/non-pCR pairs whose order changes is "
                     "within 5 percentage points of the discovery fraction for "
                     "the same endpoint",
        "per_endpoint": {e: {"validation": decomposition[e]["fraction_changed"],
                             "discovery": DISCOVERY_VALUES[e]["reordered"],
                             "difference": gaps[e]} for e in DISCOVERY_VALUES},
        "pass": bool(all(abs(g) <= C3_POINTS for g in gaps.values()))}

    verdicts["C4"] = {
        "statement": "substitution increases volume by at least 10% in at least "
                     "90% of patients at both visits",
        "per_visit": volumes, "required_share": C4_SHARE,
        "pass": bool(all(v["share_inflated_at_least_10_percent"] >= C4_SHARE
                         for v in volumes.values()))}

    verdicts["C5"] = {
        "statement": "the Spearman correlation between the measured-late and "
                     "substituted scores is at least 0.85 at all three endpoints",
        "per_endpoint": spearmans, "floor": C5_SPEARMAN,
        "pass": bool(all(spearmans[e] >= C5_SPEARMAN for e in DISCOVERY_VALUES))}
    return verdicts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path)
    args = ap.parse_args()

    lock = require_freeze()
    check_reproduces_published()
    seed = int(lock["seed"])
    replicates = int(lock["bootstrap_samples"])

    vroot = args.root.resolve()
    rows, flow = read_volumes(vroot)
    print("volumes read: %d of %d prepared patients complete both visits"
          % (flow["complete_both_visits"], flow["prepared"]))
    assert all("pcr" not in r for r in rows), "an outcome leaked in before step 3"

    joined = join_labels(rows)
    y = np.array([r["pcr"] for r in rows])
    print("outcome labels joined: %d patients, %d pCR, %d non-pCR"
          % (joined["labelled"], int(y.sum()), int((1 - y).sum())))
    print()

    results, decomposition, spearmans = {}, {}, {}
    for endpoint in ENDPOINTS:
        real = np.array([r["real_" + endpoint] for r in rows])
        sub = np.array([r["sub_" + endpoint] for r in rows])
        results[endpoint] = compare(y, real, sub, seed=seed, replicates=replicates)
        decomposition[endpoint] = pair_decomposition(y, real, sub)
        spearmans[endpoint] = spearman(real, sub)

    volumes = {}
    for visit in VISITS:
        real = np.array([r["real_V" + visit[-1]] for r in rows])
        sub = np.array([r["sub_V" + visit[-1]] for r in rows])
        inflated = np.where(real > 0, (sub - real) / np.where(real > 0, real, 1),
                            np.where(sub > 0, np.inf, 0.0))
        retained = np.divide(real, sub, out=np.full(real.shape, np.nan),
                             where=sub > 0)
        volumes[visit] = {
            "median_real_cc": float(np.median(real)),
            "median_sub_cc": float(np.median(sub)),
            "share_inflated_at_least_10_percent":
                float(np.mean(inflated >= C4_INFLATION)),
            "patients_inflated": int((inflated >= C4_INFLATION).sum()),
            "n": int(real.size),
            "retained_fraction_quartiles":
                np.nanpercentile(retained, [25, 50, 75]).tolist(),
        }

    verdicts = judge(results, decomposition, volumes, spearmans)

    tolerant = [r for r in rows if r["within_tolerance"]]
    sensitivity = {}
    if len(tolerant) < len(rows):
        yt = np.array([r["pcr"] for r in tolerant])
        for endpoint in ENDPOINTS:
            sensitivity[endpoint] = compare(
                yt, np.array([r["real_" + endpoint] for r in tolerant]),
                np.array([r["sub_" + endpoint] for r in tolerant]),
                seed=seed, replicates=replicates)

    strata = {}
    for factor in ("manufacturer", "site"):
        groups = {}
        for r in rows:
            groups.setdefault(r[factor], []).append(r)
        strata[factor] = {}
        for name, subset in sorted(groups.items()):
            ys = np.array([r["pcr"] for r in subset])
            strata[factor][name] = {
                "n": len(subset), "pcr": int(ys.sum()),
                PRIMARY: compare(ys,
                                 np.array([r["real_" + PRIMARY] for r in subset]),
                                 np.array([r["sub_" + PRIMARY] for r in subset]),
                                 seed=seed, replicates=replicates)}

    out = vroot / "data/derived/theme4_validation_v1"
    payload = {
        "scope": "prespecified independent validation of BPEX-107031, run once",
        "plan_sha256": lock["plan_sha256"], "cohort_sha256": lock["cohort_sha256"],
        "seed": seed, "bootstrap": replicates,
        "flow": flow, "n": len(rows), "n_pcr": int(y.sum()),
        "results": results, "pair_decomposition": decomposition,
        "spearman": spearmans, "volumes": volumes,
        "criteria": verdicts,
        "sensitivity_excluding_out_of_tolerance": sensitivity,
        "strata": strata,
        "discovery_values": DISCOVERY_VALUES,
    }
    (out / "validation_results.json").write_text(json.dumps(payload, indent=2),
                                                 encoding="utf-8")

    print("%-9s %-28s %-28s" % ("endpoint", "validation (n=%d)" % len(rows),
                                "discovery (n=117)"))
    for endpoint in ENDPOINTS:
        r = results[endpoint]
        d = DISCOVERY_VALUES.get(endpoint)
        print("%-9s real %.3f sub %.3f  %+.4f [%+.4f,%+.4f]   %s"
              % (endpoint, r["real_auc"], r["sub_auc"],
                 r["difference_sub_minus_real"], *r["difference_ci"],
                 ("%+.4f [%+.3f,%+.3f]" % (d["delta"], *d["ci"])) if d else ""))
    print()
    for name in ("C1", "C2", "C3", "C4", "C5"):
        v = verdicts[name]
        print("%s  %-4s  %s" % (name, "PASS" if v["pass"] else "FAIL",
                                v["statement"]))
    print()
    print("wrote", out / "validation_results.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
