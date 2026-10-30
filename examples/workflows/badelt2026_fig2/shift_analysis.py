
#!/usr/bin/env python3
"""
fourplot.py — 2x2 correlation plots split by label (tw/fw) and column pair.

Input format (whitespace-separated, one row per line):
    <label> <a> <b> <c> [hue_string]

The optional 5th column can be any string and is used as hue (one colour per unique value).

Usage:
    python fourplot.py data.txt
    python fourplot.py data.txt -o plot.png
    cat data.txt | python fourplot.py
"""

import argparse
import sys
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
def load_data(source):
    rows = defaultdict(list)
    src = source if hasattr(source, "read") else open(source)
    for lineno, line in enumerate(src, 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            print(f"Warning: skipping line {lineno} (expected ≥4 fields): {line!r}",
                  file=sys.stderr)
            continue
        label = parts[0]
        try:
            a, b, c = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            print(f"Warning: skipping line {lineno} (non-numeric values): {line!r}",
                  file=sys.stderr)
            continue
        hue = parts[4] if len(parts) >= 5 else None
        rows[label].append((a, b, c, hue))
    return rows

def make_palette(hue_values):
    """Assign a colour to each unique hue string, None/other → grey."""
    unique = sorted(v for v in hue_values if v is not None and v != "other")
    colors = ["steelblue", "tomato", "seagreen", "darkorange",
              "mediumpurple", "saddlebrown", "deeppink", "teal"]
    palette = {v: colors[i % len(colors)] for i, v in enumerate(unique)}
    palette["add_stack"] = "steelblue"
    palette["del_stack"] = "darkorange"
    palette["neutral_2s"] = "deeppink"
    palette["neutral_1s"] = "teal"
    palette["neutral_trans"] = "deeppink"
    palette["neutral_cis"] = "teal" # cheaper?
    palette["other"] = "lightgrey"
    palette[None] = "lightgrey"
    return palette

def scatter_panel(ax, rows, col_idx, xlabel, ylabel, title, palette):
    """col_idx: 1 → b, 2 → c (index into the (a,b,c,hue) tuple)."""
    if not rows:
        ax.set_title(f"{title}\n(no data)")
        return

    has_hue = any(r[3] is not None for r in rows)

    def plot_group(subset, hue_key):
        if not subset:
            return
        x = np.array([r[0]/100 for r in subset])
        y = np.array([r[col_idx]/100 for r in subset])
        ax.scatter(x, y, alpha=0.6, edgecolors="white", linewidths=0.4,
                   s=40, color=palette[hue_key], zorder=3,
                   label=str(hue_key) if hue_key is not None else "?")

    if has_hue:
        # group by hue value; plot None first, then sorted hue values (on top)
        groups = defaultdict(list)
        for r in rows:
            groups[r[3]].append(r)
        plot_group(groups.get(None, []), None)
        plot_group(groups.get("other", []), "other")
        for hue_key in sorted(k for k in groups if k is not None and k != "other"):
            plot_group(groups[hue_key], hue_key)
    else:
        plot_group(rows, None)

    x_all = np.array([r[0]       for r in rows])
    y_all = np.array([r[col_idx] for r in rows])
    lo = min(x_all.min(), y_all.min())
    hi = max(x_all.max(), y_all.max())
    lo = min(-10, -10)
    hi = min(15, 15)
    pad = (hi - lo) * 0.05 or 1.0
    pad = (hi - lo) * 0.05 or 1.0
    lims = (lo - pad, hi + pad)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.plot(lims, lims, color="grey", linewidth=1.0, linestyle="--", zorder=1, label="y = x")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(f"{title}  (n={len(x_all)})")
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8)
def main():
    parser = argparse.ArgumentParser(
        description="2x2 correlation plots split by label and column pair.",
        epilog="Reads from stdin if no file is given.",
    )
    parser.add_argument("input", nargs="?", default=None)
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("--col1-label", default="a", help="Name for 1st numeric column")
    parser.add_argument("--col2-label", default="b", help="Name for 2nd numeric column")
    parser.add_argument("--col3-label", default="c", help="Name for 3rd numeric column")
    args = parser.parse_args()

    source = args.input if args.input is not None else sys.stdin
    data = load_data(source)

    tw = data.get("tw", [])
    fw = data.get("fw", [])

    fig, axes = plt.subplots(2, 1, figsize=(4.5, 7.5))
    fig.suptitle(args.col3_label, fontsize=13, y=0.99)

    a, b, c = args.col1_label, args.col2_label, args.col3_label

    tw_palette = make_palette([r[3] for r in tw])
    fw_palette = make_palette([r[3] for r in fw])

    scatter_panel(axes[0], tw, 1, a, b, f"three-way shift evaluation", tw_palette)
    #scatter_panel(axes[0, 1], tw, 2, a, c, f"three-way shift evaluation", tw_palette)
    scatter_panel(axes[1], fw, 1, a, b, f"four-way shift evaluation", fw_palette)
    #scatter_panel(axes[1, 1], fw, 2, a, c, f"four-way shift evaluation", fw_palette)

    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"Saved: {args.output}")
    else:
        plt.show()
if __name__ == "__main__":
    main()
