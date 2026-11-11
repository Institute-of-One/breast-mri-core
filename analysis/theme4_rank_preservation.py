"""What actually governs whether the SER criterion preserves the volume ordering.

    python scripts/theme4_rank_preservation.py

Referee 1 point 7 of BPEX-107031 objects that section 2.6 offers a smaller
spread of the retained fraction as a necessary condition for rank preservation,
and that ordering also depends on the correlation between the components and on
the spacing between individual patients. **The referee is right.** This computes
the condition that does hold, and measures it.

Write X for log V_sub and Y for log(V_real / V_sub), so that
log V_real = X + Y. For a pair of participants with differences dX and dY, the
substituted and measured-late volumes order that pair the same way exactly when

    dX * (dX + dY) > 0,

that is, when dY does not both oppose dX and exceed it in size. Two consequences
follow immediately, and both contradict a condition stated on the spread of Y
alone:

* **Not necessary.** If Y = aX with a > 0, then sd(Y) = a sd(X) can be as large
  as one likes while log V_real = (1 + a) X is strictly increasing in X, so no
  pair reverses.
* **Not sufficient.** If Y is independent of X with a small spread, pairs whose
  volumes lie close together still reverse; the number of them depends on the
  spacing of X, not on sd(Y) alone.

So the quantity to report is the share of pairs that actually reverse, which is
computed here exactly over all pairs, with no resampling and no outcome.

Reads the archived paired scores. Writes
``data/derived/theme4_revision_v1/rank_preservation.json``; nothing published is
modified.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SCORES = ROOT / "data/derived/theme4_pcr_v1/paired_scores.csv"
OUT = ROOT / "data/derived/theme4_revision_v1"

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def volumes(rows: list[dict], visit: str) -> tuple[np.ndarray, np.ndarray]:
    """The scores are negative volumes at T0 and T1; return the volumes."""
    real = np.array([-float(r["real_" + visit]) for r in rows])
    sub = np.array([-float(r["sub_" + visit]) for r in rows])
    return real, sub


def pair_table(x: np.ndarray, y: np.ndarray) -> dict:
    """Exact enumeration over all unordered pairs. No sampling."""
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    upper = np.triu(np.ones_like(dx, dtype=bool), k=1)
    dx, dy = dx[upper], dy[upper]
    total = dx.size
    product = dx * (dx + dy)
    reversed_ = product < 0
    tied = product == 0
    opposed = np.sign(dy) == -np.sign(dx)
    bigger = np.abs(dy) > np.abs(dx)
    return {
        "pairs": int(total),
        "reversed": int(reversed_.sum()),
        "fraction_reversed": float(reversed_.sum() / total),
        "tied_after": int(tied.sum()),
        "pairs_with_dY_opposed_to_dX": int(opposed.sum()),
        "fraction_dY_opposed": float(opposed.sum() / total),
        "pairs_with_abs_dY_exceeding_abs_dX": int(bigger.sum()),
        "fraction_abs_dY_exceeds": float(bigger.sum() / total),
        # The spread-only reading would predict reversal whenever |dY| > |dX|.
        # Reversal also needs the signs to oppose, which is the correlation the
        # referee names; the gap between these two numbers is that requirement.
        "predicted_by_magnitude_alone_but_not_reversed":
            int((bigger & ~reversed_).sum()),
    }


def kendall_tau_b(a: np.ndarray, b: np.ndarray) -> float:
    da = a[:, None] - a[None, :]
    db = b[:, None] - b[None, :]
    upper = np.triu(np.ones_like(da, dtype=bool), k=1)
    da, db = da[upper], db[upper]
    concordant = int(((da * db) > 0).sum())
    discordant = int(((da * db) < 0).sum())
    ties_a = int((da == 0).sum())
    ties_b = int((db == 0).sum())
    n = da.size
    denominator = np.sqrt((n - ties_a) * (n - ties_b))
    return float((concordant - discordant) / denominator) if denominator else float("nan")


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    def rank(v):
        order = np.argsort(v, kind="mergesort")
        r = np.empty(len(v), float)
        r[order] = np.arange(1, len(v) + 1)
        # average ranks for ties
        _, inverse, counts = np.unique(v, return_inverse=True, return_counts=True)
        sums = np.zeros(len(counts))
        np.add.at(sums, inverse, r)
        return (sums / counts)[inverse]
    ra, rb = rank(a), rank(b)
    return float(np.corrcoef(ra, rb)[0, 1])


def main() -> int:
    rows = list(csv.DictReader(SCORES.open(encoding="utf-8-sig")))
    result = {
        "source": str(SCORES.relative_to(ROOT)).replace("\\", "/"),
        "n_participants": len(rows),
        "condition": "a pair keeps its order iff dX*(dX+dY) > 0, with X = log "
                     "V_sub and Y = log(V_real/V_sub); the spread of Y alone is "
                     "neither necessary nor sufficient",
        "visits": {},
    }
    for visit in ("T0", "T1"):
        real, sub = volumes(rows, visit)
        keep = real > 0
        n_zero = int((~keep).sum())
        x = np.log(sub[keep])
        y = np.log(real[keep] / sub[keep])
        table = pair_table(x, y)
        entry = {
            "n_log": int(keep.sum()),
            "n_excluded_zero_measured_volume": n_zero,
            "sd_log_sub_volume": float(x.std(ddof=1)),
            "sd_log_retained_fraction": float(y.std(ddof=1)),
            "ratio_sd_log": float(y.std(ddof=1) / x.std(ddof=1)),
            "corr_log_sub_vs_log_retained": float(np.corrcoef(x, y)[0, 1]),
            "spearman_real_vs_sub": spearman(real[keep], sub[keep]),
            "kendall_tau_b_real_vs_sub": kendall_tau_b(real[keep], sub[keep]),
        }
        entry.update(table)
        result["visits"][visit] = entry
        print("%s  n=%d  sd(X)=%.4f sd(Y)=%.4f ratio=%.4f corr=%+.4f"
              % (visit, entry["n_log"], entry["sd_log_sub_volume"],
                 entry["sd_log_retained_fraction"], entry["ratio_sd_log"],
                 entry["corr_log_sub_vs_log_retained"]))
        print("      pairs %d | reversed %d (%.2f%%) | Kendall tau_b %.4f | Spearman %.4f"
              % (entry["pairs"], entry["reversed"], 100 * entry["fraction_reversed"],
                 entry["kendall_tau_b_real_vs_sub"], entry["spearman_real_vs_sub"]))
        print("      |dY|>|dX| in %d pairs (%.2f%%), of which %d do NOT reverse "
              "because the signs agree"
              % (entry["pairs_with_abs_dY_exceeding_abs_dX"],
                 100 * entry["fraction_abs_dY_exceeds"],
                 entry["predicted_by_magnitude_alone_but_not_reversed"]))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rank_preservation.json").write_text(json.dumps(result, indent=2),
                                                encoding="utf-8")
    print()
    print("wrote", (OUT / "rank_preservation.json").relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
