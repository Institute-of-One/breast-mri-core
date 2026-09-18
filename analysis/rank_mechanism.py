"""Why a large volume change leaves pCR discrimination nearly unchanged: the pair decomposition.

An AUC is the fraction of (pCR, non-pCR) pairs that a score orders correctly, ties counting
one half (Hanley and McNeil 1982). Replacing the late image with the early one therefore
changes the AUC only through pairs whose ordering it changes, and the change in AUC is
exactly the mean change in pair concordance. This script counts those pairs from the frozen
paired scores of Supplementary Data S1, and reports alongside them how much the volumes
themselves moved and how well their ranking was preserved.

Descriptive and post hoc: nothing is fitted, resampled, thresholded or tuned, and the
frozen AUCs are recomputed first and required to match `auc_results.json` to 1e-12, so
the decomposition is of the reported numbers and not of a variant.

    python analysis/rank_mechanism.py            # writes analysis/rank_mechanism.json
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
S1 = ROOT / "supplementary_data_s1"
OUT = Path(__file__).resolve().with_suffix(".json")
METRICS = {"T0": ("real_T0", "sub_T0"), "T1": ("real_T1", "sub_T1"),
           "relative": ("real_relative", "sub_relative"), "change": ("real_change", "sub_change")}


def concordance(pos: np.ndarray, neg: np.ndarray) -> np.ndarray:
    """Pair matrix: 1 where the pCR score is higher, 0.5 on ties, 0 otherwise."""
    diff = pos[:, None] - neg[None, :]
    return (diff > 0) + 0.5 * (diff == 0)


def main() -> int:
    with open(S1 / "paired_scores.csv", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    pcr = np.array([int(r["pcr"]) for r in rows]) == 1
    frozen = {r["metric"]: r for r in json.loads((S1 / "auc_results.json").read_text(encoding="utf-8"))["results"]
              if r["group"] == "all"}

    out: dict = {"n": len(rows), "pcr": int(pcr.sum()), "non_pcr": int((~pcr).sum()),
                 "pairs": int(pcr.sum() * (~pcr).sum()), "metrics": {}}
    for metric, (real_col, sub_col) in METRICS.items():
        real = np.array([float(r[real_col]) for r in rows])
        sub = np.array([float(r[sub_col]) for r in rows])
        c_real = concordance(real[pcr], real[~pcr])
        c_sub = concordance(sub[pcr], sub[~pcr])
        auc_real, auc_sub = float(c_real.mean()), float(c_sub.mean())

        # The decomposition is only meaningful if it is of the reported AUCs.
        assert abs(auc_real - frozen[metric]["real_auc"]) < 1e-12, (metric, auc_real)
        assert abs(auc_sub - frozen[metric]["sub_auc"]) < 1e-12, (metric, auc_sub)

        delta = c_sub - c_real
        gained, lost = float(delta[delta > 0].sum()), float(-delta[delta < 0].sum())
        assert abs((gained - lost) / delta.size - (auc_sub - auc_real)) < 1e-12
        rho = spearmanr(real, sub).statistic
        out["metrics"][metric] = {
            "auc_real": auc_real, "auc_sub": auc_sub, "auc_difference": auc_sub - auc_real,
            "pairs_changed": int((delta != 0).sum()),
            "fraction_pairs_changed": float((delta != 0).mean()),
            "concordance_gained": gained, "concordance_lost": lost,
            "net_concordance": gained - lost,
            "pairs_gaining": int((delta > 0).sum()), "pairs_losing": int((delta < 0).sum()),
            "tied_pairs_real": int((c_real == 0.5).sum()), "tied_pairs_sub": int((c_sub == 0.5).sum()),
            "spearman_real_vs_sub": float(rho),
        }

    # How far the volumes themselves moved. Scores are negative volumes in mL.
    for visit in ("T0", "T1"):
        v_real = -np.array([float(r[f"real_{visit}"]) for r in rows])
        v_sub = -np.array([float(r[f"sub_{visit}"]) for r in rows])
        ok = v_sub > 0
        retained = v_real[ok] / v_sub[ok]
        q1, med, q3 = np.percentile(retained, [25, 50, 75])
        out[f"volume_{visit}"] = {
            "median_real_ml": float(np.median(v_real)), "median_sub_ml": float(np.median(v_sub)),
            "retained_fraction_median": float(med), "retained_fraction_q1": float(q1),
            "retained_fraction_q3": float(q3), "n_with_positive_sub_volume": int(ok.sum()),
            "n_real_zero": int((v_real == 0).sum()),
        }
        # Outcome-free account of rank preservation: log Vreal = log Vsub + log(retained fraction).
        # Where the retained fraction varies much less than the volume, ordering survives.
        pos = (v_real > 0) & (v_sub > 0)
        log_sub, log_rf = np.log(v_sub[pos]), np.log(v_real[pos] / v_sub[pos])
        out[f"volume_{visit}"].update({
            "n_log": int(pos.sum()),
            "sd_log_sub_volume": float(np.std(log_sub, ddof=1)),
            "sd_log_retained_fraction": float(np.std(log_rf, ddof=1)),
            "ratio_sd_log": float(np.std(log_rf, ddof=1) / np.std(log_sub, ddof=1)),
            "corr_log_sub_vs_log_retained": float(np.corrcoef(log_sub, log_rf)[0, 1]),
        })

    # The same accounting within each subtype, where few pairs carry each AUC difference.
    subtype = np.array([r["subtype"] for r in rows])
    frozen_all = json.loads((S1 / "auc_results.json").read_text(encoding="utf-8"))["results"]
    out["subtypes"] = {}
    for group in sorted(set(subtype)):
        g = subtype == group
        entry = {"n": int(g.sum()), "pcr": int((pcr & g).sum()), "pairs": int((pcr & g).sum() * (~pcr & g).sum())}
        for metric, (real_col, sub_col) in METRICS.items():
            real = np.array([float(r[real_col]) for r in rows])[g]
            sub = np.array([float(r[sub_col]) for r in rows])[g]
            pg = pcr[g]
            delta = concordance(sub[pg], sub[~pg]) - concordance(real[pg], real[~pg])
            ref = next(r for r in frozen_all if r["group"] == group and r["metric"] == metric)
            assert abs(delta.mean() - ref["difference_sub_minus_real"]) < 1e-12, (group, metric)
            entry[metric] = {"auc_difference": float(delta.mean()), "net_concordance": float(delta.sum()),
                             "pairs_changed": int((delta != 0).sum()),
                             "auc_real": ref["real_auc"], "auc_sub": ref["sub_auc"],
                             "difference_ci_low": ref["difference_ci"][0], "difference_ci_high": ref["difference_ci"][1]}
        out["subtypes"][group] = entry

    OUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    for metric, m in out["metrics"].items():
        print(f"{metric:>8}: AUC {m['auc_real']:.3f} -> {m['auc_sub']:.3f}  "
              f"pairs changed {m['pairs_changed']}/{out['pairs']} ({m['fraction_pairs_changed']:.1%})  "
              f"gained {m['concordance_gained']:.1f} lost {m['concordance_lost']:.1f}  "
              f"Spearman {m['spearman_real_vs_sub']:.3f}")
    for visit in ("T0", "T1"):
        v = out[f"volume_{visit}"]
        print(f"{visit}: median volume {v['median_real_ml']:.2f} -> {v['median_sub_ml']:.2f} mL; retained fraction "
              f"median {v['retained_fraction_median']:.3f} (IQR {v['retained_fraction_q1']:.3f}-{v['retained_fraction_q3']:.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
