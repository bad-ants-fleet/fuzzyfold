#!/usr/bin/env python3
"""Evaluate a structure like ff-eval, and draw it as an arc diagram.

Input (file or stdin), as for ff-eval:

    >optional header
    SEQUENCE        (strands separated by '+' or '&')
    STRUCTURE       (dot-bracket; ()[]{}<> and Aa..Zz for pseudoknots)

Each stem of a pseudoknot gets its own color, as the energy evaluation sees
the stems (see `fuzzyfold.stems`): strained stems are solid, relaxed stems
dashed. Other pairs are grey. Strand breaks are dotted vertical lines.

Requires the fuzzyfold Python module (crates/ff_python) and matplotlib.

    python ff_eval_arcs.py input.fa -v -o xr.svg
"""

import sys
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.lines import Line2D

import fuzzyfold as ff


def read_input(handle):
    """Header (or None), sequence and structure of the first record."""
    lines = [l.strip() for l in handle if l.strip()]
    header = lines.pop(0)[1:].strip() if lines and lines[0].startswith(">") else None
    if len(lines) < 2:
        sys.exit("Expected a sequence and a structure.")
    return header, lines[0], lines[1].split()[0]


def arc_diagram(sequence, structure, title, path):
    n = len(structure)
    pairs = ff.pairs(structure)
    stems = ff.stems(structure)
    colors = plt.get_cmap("tab10").colors

    # The stem of each pair, if it belongs to a pseudoknot.
    style = {}
    legend = []
    for b, (stem, strained) in enumerate(stems):
        color = colors[b % len(colors)]
        line = "-" if strained else "--"
        for pair in stem:
            style[pair] = (color, line)
        (i, j), (p, q) = stem[0], stem[-1]
        label = f"({i}, {j})..({p}, {q}), {'strained' if strained else 'relaxed'}"
        legend.append(Line2D([0], [0], color=color, linestyle=line, label=label))

    span = max((j - i for i, j in pairs), default=1)
    fig, ax = plt.subplots(figsize=(max(6, 0.15 * n), max(2, 0.08 * span + 1)))
    for i, j in pairs:
        color, line = style.get((i, j), ("0.6", "-"))
        ax.add_patch(Arc(((i + j) / 2, 0), j - i, j - i, theta1=0, theta2=180,
                         color=color, linestyle=line, linewidth=1.5))
    for k, c in enumerate(sequence):
        if c in "+&":
            ax.axvline(k, color="0.3", linestyle=":", linewidth=1)
        elif n <= 200:
            ax.text(k, -0.6, c, ha="center", va="top", fontsize=7, family="monospace")

    ax.set_xlim(-1, n)
    ax.set_ylim(-2, span / 2 + 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=9)
    if legend:
        ax.legend(handles=legend, fontsize=7, loc="upper left", bbox_to_anchor=(1, 1), frameon=False)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("input", nargs="?", default="-", help="input file, or - for stdin")
    parser.add_argument("--params", default="rna_turner_2004",
                        choices=["rna_turner_2004", "rna_extended", "dna"])
    parser.add_argument("--celsius", type=float, default=37.0)
    parser.add_argument("-v", "--verbose", action="store_true", help="list all loops")
    parser.add_argument("-o", "--output", help="arc diagram (.svg, .png, .pdf); "
                        "default: <header>.svg")
    args = parser.parse_args()

    with (sys.stdin if args.input == "-" else open(args.input)) as handle:
        header, sequence, structure = read_input(handle)
    if len(sequence) != len(structure):
        sys.exit("Sequence and structure differ in length.")

    model = ff.Fuzzyfold(params=args.params, celsius=args.celsius)
    try:
        terms = model.energy_terms(sequence, structure)
    except ValueError as e:
        sys.exit(str(e))
    energy = sum(e for level, _, e in terms if level == 0)

    if header:
        print(f">{header}")
    if args.verbose:
        for level, term, e in terms:
            indent = "  " * level + ("  " if level else "")
            print(f"{indent}{term + ':' if level == 0 else term:<{41 - len(indent)}} {e:>6.2f}")
    print(f"{sequence}\n{structure} {energy:>6.2f}")

    path = args.output or f"{(header or 'structure').split()[0]}.svg"
    title = f"{header + ': ' if header else ''}{energy:.2f} kcal/mol"
    arc_diagram(sequence, structure, title, path)
    print(f"Arc diagram: {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
