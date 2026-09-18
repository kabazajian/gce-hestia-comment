"""Stage 5 summarizer (local, light): per-variant deltas vs baseline.

For each variant and each statistic: median and max |delta| over the 36
views, plus whether any view's gnfw2-vs-combo_1 Pearson preference flips.
The verdict per axis is a counting statement: which Stage 4 conclusions
move, and by how much, under each systematic.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

PATH = Path(sys.argv[1] if len(sys.argv) > 1
            else "/private/nfs/Data2/hestia_muru/stage5_systematics.json")
STATS = ("corr_gnfw2", "corr_combo_1", "corr_f98", "q30", "c4_30", "conc_3_10")
PAIRS = (("base", ("k16", "k64", "smax50", "rsun8p0", "half_a", "half_b")),
         ("g31_base", ("g31_center_ahf",)))


def main():
    d = json.load(open(PATH))
    for base_key, variants in PAIRS:
        base = d[base_key]
        pref0 = [r["corr_gnfw2"] > r["corr_combo_1"] for r in base]
        print(f"== vs {base_key} ({sum(pref0)}/36 views prefer gnfw2) ==")
        hdr = f"  {'variant':<14s}" + "".join(f" {s:>12s}" for s in STATS) + "   flips"
        print(hdr + "\n  " + " " * 14 + "  (median |d| over 36 views; max in parens)")
        for v in variants:
            rows = d[v]
            cells = []
            for s in STATS:
                dd = np.abs([r[s] - b[s] for r, b in zip(rows, base)])
                cells.append(f"{np.median(dd):.4f}({dd.max():.3f})")
            flips = sum((r["corr_gnfw2"] > r["corr_combo_1"]) != p0
                        for r, p0 in zip(rows, pref0))
            print(f"  {v:<14s}" + "".join(f" {c:>12s}" for c in cells)
                  + f"   {flips}/36")
        print()


if __name__ == "__main__":
    main()
