"""Plots of macrostate occupancies over time (needs matplotlib).

    python -m fuzzyfold.plotting run.tln -o run -f pdf png

The x-axis is linear up to t_sep and logarithmic after it, or logarithmic
only without t_sep. A Timecourse or a .tln file gives t_sep and the
standard errors (shown as bands) itself; for a .nxy file, give t_sep and
num_sims. With t_ext and seq_length, a top axis shows the sequence length
of a co-transcriptional simulation.
"""

import sys
import argparse
from dataclasses import dataclass

import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as ticker
    import matplotlib.transforms as mtransforms
    from matplotlib.patches import ConnectionPatch
except ImportError as e:
    raise ImportError("fuzzyfold.plotting needs matplotlib (pip install matplotlib)") from e

from .timecourse import Timecourse

# Series with these names (case-insensitive) are light grey, behind the
# other curves, and last in the legend (unless given in labels).
GREY_NAMES = {"unassigned"}
GREY_COLOR = "0.75"
LABEL_COLOR = "black"
# Minimum gap between neighboring tick labels (points).
LABEL_PAD_PT = 4
# Minimum spacing of the sequence length minor ticks (points).
EXT_MIN_TICK_PT = 1.5
# Opacity of the standard error bands.
BAND_ALPHA = 0.2


@dataclass
class Style:
    """Figure width (inches) and font sizes (points).

    Line widths scale with fig_w, font sizes do not.
    """
    fig_w: float = 7.0
    aspect: float = 3.0 / 7.0
    title: float = 12
    labels: float = 8
    legend: float = 6.5
    ticks: float = 7
    split_label: float = 8

    def lw(self, base):
        return base * self.fig_w / 7.0


# Tick labels

def _fmt_sci(x):
    """'0', '$10^{n}$' or '$k{\\cdot}10^{n}$' (k with up to 3 digits)."""
    if x == 0:
        return "0"
    sign = "-" if x < 0 else ""
    x = abs(x)
    exp = int(np.floor(np.log10(x)))
    coeff = float(f"{x / 10**exp:.3g}")
    if coeff >= 10:
        coeff /= 10
        exp += 1
    if coeff == 1:
        return r"$%s10^{%d}$" % (sign, exp)
    return r"$%s%s{\cdot}10^{%d}$" % (sign, f"{coeff:g}", exp)


def _fmt_val(x):
    """Plain number for moderate magnitudes, scientific otherwise."""
    if x == 0:
        return "0"
    if 1e-3 <= abs(x) < 1e4:
        return f"{x:.4g}"
    return _fmt_sci(x)


def _fmt_log(x):
    """Powers of ten as 10^n, everything else as _fmt_val."""
    if x > 0:
        e = np.log10(x)
        if abs(e - round(e)) < 1e-9:
            return r"$10^{%d}$" % round(e)
    return _fmt_val(x)


def _log_candidates(lo, hi):
    """Major ticks inside (lo, hi): decades, or with 2 and 5 times decades
    if fewer than two decades fall inside."""
    e0 = int(np.floor(np.log10(lo))) - 1
    e1 = int(np.ceil(np.log10(hi))) + 1

    def inside(v):
        return lo * (1 + 1e-9) < v < hi * (1 - 1e-9)

    decades = [10.0**e for e in range(e0, e1 + 1) if inside(10.0**e)]
    if len(decades) >= 2:
        return decades
    extra = sorted(k * 10.0**e for e in range(e0, e1 + 1) for k in (2, 5)
                   if inside(k * 10.0**e))
    return decades + extra


def _log_tick_entries(ax, lo, hi, n):
    """Tick entries of a log axis: automatic, or n forced ticks evenly
    spaced in log space (rounded to 2 digits)."""
    if n is None:
        return [(ax, v, _fmt_log(v), False) for v in _log_candidates(lo, hi)]
    values = [float(f"{v:.2g}") for v in np.geomspace(lo, hi, int(n) + 2)[1:-1]]
    return [(ax, v, _fmt_log(v), True) for v in values if lo < v < hi]


def _log_minor(ax, enabled):
    """Unlabelled minor ticks at 2..9 times 10^n, or none."""
    if enabled:
        ax.xaxis.set_minor_locator(
            ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    else:
        ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())


# A tick entry is (axes, x, label, forced). Forced labels are always shown,
# others in order if they do not overlap a shown label. Overlap is measured
# in display space, so it works across the two panels of a split axis.

def _declutter(fig, entries, fontsize):
    renderer = fig.canvas.get_renderer()
    pad = LABEL_PAD_PT * fig.dpi / 72.0
    spans = []
    for ax, x, label, _ in entries:
        if not label:
            spans.append(None)
            continue
        cx = ax.transData.transform((x, 0.0))[0]
        t = fig.text(0, 0, label, fontsize=fontsize)
        w = t.get_window_extent(renderer=renderer).width
        t.remove()
        spans.append((cx - w / 2 - pad / 2, cx + w / 2 + pad / 2))

    keep = [False] * len(entries)
    taken = []
    order = ([i for i, e in enumerate(entries) if e[3]] +
             [i for i, e in enumerate(entries) if not e[3]])
    for i in order:
        if spans[i] is None:
            continue
        a, b = spans[i]
        if entries[i][3] or all(b <= c or a >= d for c, d in taken):
            keep[i] = True
            taken.append((a, b))
    return keep


def _apply_ticks(fig, entries, fontsize):
    """Fixed major ticks per axes; dropped labels are empty."""
    keep = _declutter(fig, entries, fontsize)
    by_ax = {}
    for (ax, x, label, _), k in zip(entries, keep):
        ticks = by_ax.setdefault(ax, {})
        if x not in ticks or (k and label):
            ticks[x] = label if k else ""
    for ax, ticks in by_ax.items():
        xs = sorted(ticks)
        ax.set_xticks(xs)
        ax.set_xticklabels([ticks[x] for x in xs], fontsize=fontsize)


def _offset(fig, base, dy_pt):
    """base shifted vertically by dy_pt points."""
    return base + mtransforms.ScaledTranslation(0, dy_pt / 72.0, fig.dpi_scale_trans)


# Sequence length axis

def _thin_minor(ax, lengths, times):
    """Every step, or every 2nd, 5th, ... if the ticks would be too dense."""
    if len(times) < 2:
        return times
    px = ax.transData.transform(np.column_stack([times, np.zeros_like(times)]))[:, 0]
    limit = EXT_MIN_TICK_PT * ax.figure.dpi / 72.0
    for step in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000):
        sel = lengths % step == 0
        if sel.sum() < 2 or np.min(np.diff(px[sel])) >= limit:
            return times[sel]
    return times[lengths % 1000 == 0]


def _add_ext_axes(fig, panels, ext, style):
    """A top axis with the sequence length on each panel that needs one.

    panels: (ax, lo, hi, lo_inclusive, spines_to_hide)
    ext:    (lengths, times, every)
    """
    lengths, times, every = ext
    first, last = int(lengths[0]), int(lengths[-1])
    entries, twins = [], []
    for ax, lo, hi, lo_incl, hide in panels:
        tol = 1e-9 * max(abs(lo), abs(hi), 1e-300)
        m = (times >= lo - tol) if lo_incl else (times > lo + tol)
        m &= times <= hi + tol
        if ax.get_xscale() == "log":
            m &= times > 0
        if not m.any():
            continue
        tw = ax.twiny()
        tw.set_xscale(ax.get_xscale())
        tw.set_xlim(lo, hi)
        for sp in hide:
            tw.spines[sp].set_visible(False)
        tw.tick_params(axis="x", which="both", labelsize=style.ticks)
        tw.xaxis.set_minor_locator(ticker.FixedLocator(_thin_minor(tw, lengths[m], times[m])))
        tw.xaxis.set_minor_formatter(ticker.NullFormatter())
        for length, t in zip(lengths[m], times[m]):
            length = int(length)
            if length in (first, last) or length % every == 0:
                entries.append((tw, float(t), str(length), length in (first, last)))
        twins.append(tw)
    if entries:
        _apply_ticks(fig, entries, style.ticks)
    return twins


# Input

def read_nxy(path):
    """Times, names and occupancies of a .nxy file (header row: time and
    the names)."""
    names, rows = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if names is None:
                names = line.split()[1:]
            else:
                rows.append([float(x) for x in line.split()])
    if names is None or not rows:
        raise ValueError(f"{path}: no data")
    data = np.array(rows)
    return data[:, 0], names, data[:, 1:]


@dataclass
class _Data:
    times: np.ndarray
    names: list
    occupancy: np.ndarray
    error: np.ndarray = None
    t_sep: float = None
    t_ext: float = None
    seq_length: int = None
    start_length: int = None


def _load(data, t_sep, num_sims, t_ext, seq_length, start_length):
    """The data to plot; arguments that are given override the data."""
    if not isinstance(data, Timecourse):
        path = str(data)
        if path.endswith(".nxy"):
            times, names, occupancy = read_nxy(path)
            error = None
            if num_sims:
                error = np.sqrt(occupancy * (1 - occupancy) / num_sims)
            return _Data(times, names, occupancy, error, t_sep, t_ext, seq_length,
                         start_length or 1)
        data = Timecourse.read(path)
    return _Data(data.times, data.names, data.occupancy(), data.standard_error(),
                 t_sep if t_sep is not None else data.t_sep,
                 t_ext if t_ext is not None else data.t_ext,
                 seq_length if seq_length is not None else data.seq_length,
                 start_length or data.start_length or 1)


# Plot

def plot_timecourse(data, output=None, *,
                    t_sep=None, num_sims=None, bands=True,
                    t_ext=None, seq_length=None, start_length=None,
                    title="", plim=1e-2, labels=None, labels_strict=False,
                    split_pos=0.5, lin_ticks=None, log_ticks=None,
                    xlim=None, ylim=None, rescale=1.0,
                    xlabel="time [s]", ylabel="occupancy",
                    ext_label="sequence length [nt]", ext_every=10,
                    legend_loc=None, style=None):
    """Plot macrostate occupancies over time; returns the figure.

    Args:
        data:          a Timecourse, or the path of a .tln or .nxy file.
        output:        file name(s) to save to (format from the extension).
        t_sep:         split the x-axis here: linear before, log after;
                       default from the data (not for .nxy), else log only.
        num_sims:      the number of simulations of a .nxy file, for bands.
        bands:         show the standard errors as bands.
        t_ext:         time per extension; with seq_length, a top axis shows
                       the sequence length (start_length at t=0, one more
                       every t_ext). Default from the data.
        title:         plot title.
        plim:          minimum peak occupancy to plot a series (not for
                       series in labels).
        labels:        series first in the legend, in this order, always
                       plotted. Other series follow in the order they reach
                       plim, grey series last.
        labels_strict: only series in labels get a color and a legend
                       entry, the others are thin and grey.
        split_pos:     fraction of the width for the linear panel.
        lin_ticks:     number of evenly spaced, always labelled ticks on the
                       linear axis (default automatic).
        log_ticks:     number of always labelled ticks on the log axis,
                       evenly spaced in log space, without minor ticks
                       (default decades).
        xlim, ylim:    axis limits (xlim in rescaled units).
        rescale:       factor for all times (also t_sep and t_ext).
        ext_every:     label every n-th sequence length (and the first and
                       last).
        legend_loc:    matplotlib legend location; default 'upper right'
                       (split) or 'best' (log only).
        style:         a Style (figure width, font sizes).
    """
    style = style or Style()
    d = _load(data, t_sep, num_sims, t_ext, seq_length, start_length)
    time = d.times * rescale
    t_sep = None if d.t_sep is None else d.t_sep * rescale
    error = d.error if bands else None
    series = _series(time, d.names, d.occupancy, error, plim, labels, labels_strict)

    ext = None
    if d.t_ext is not None and d.seq_length is not None:
        if d.t_ext <= 0:
            raise ValueError("t_ext must be > 0")
        lengths = np.arange(d.start_length, int(d.seq_length) + 1)
        ext = (lengths, (lengths - d.start_length) * d.t_ext * rescale, max(1, int(ext_every)))

    fig_h = style.fig_w * style.aspect
    ylim = ylim if ylim is not None else (-0.02, 1.02)
    common = dict(time=time, series=series, title=title, ylim=ylim, xlim=xlim,
                  style=style, fig_h=fig_h, log_ticks=log_ticks, xlabel=xlabel,
                  ylabel=ylabel, ext=ext, ext_label=ext_label, legend_loc=legend_loc)
    if t_sep is None:
        fig = _plot_log(**common)
    else:
        fig = _plot_split(t_sep=t_sep, split_pos=split_pos, lin_ticks=lin_ticks, **common)

    for path in [output] if isinstance(output, str) else (output or []):
        fig.savefig(path, bbox_inches="tight")
    return fig


@dataclass
class _Series:
    name: str
    values: np.ndarray
    error: np.ndarray
    color: str
    lw: float
    alpha: float
    label: str
    zorder: float


def _series(time, names, occupancy, error, plim, labels, labels_strict):
    """The series to plot, in legend order, with their colors."""
    columns = {name: i for i, name in enumerate(names)}
    labels = list(labels or [])
    missing = [n for n in labels if n not in columns]
    if missing:
        raise ValueError(f"labels not in the data: {' '.join(missing)} "
                         "(a .tln file has only the occupied macrostates)")

    def first_seen(name):
        vals = occupancy[:, columns[name]]
        hit = np.nonzero(vals >= plim if plim > 0 else vals > 0)[0]
        return time[hit[0]] if hit.size else np.inf

    rest = [n for n in names
            if n not in labels and not (plim and occupancy[:, columns[n]].max() < plim)]
    rest.sort(key=first_seen)
    rest = ([n for n in rest if n.lower() not in GREY_NAMES] +
            [n for n in rest if n.lower() in GREY_NAMES])

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    out, color_idx = [], 0
    for name in labels + rest:
        vals = occupancy[:, columns[name]]
        err = None if error is None else error[:, columns[name]]
        primary = not labels or not labels_strict or name in labels
        label = name if primary else "_nolegend_"
        if name.lower() in GREY_NAMES:
            out.append(_Series(name, vals, err, GREY_COLOR, 1.2, 1.0, label, 1.9))
        elif primary:
            out.append(_Series(name, vals, err, colors[color_idx % len(colors)],
                               2, 1.0, label, 2))
            color_idx += 1
        else:
            out.append(_Series(name, vals, None, "0.65", 0.8, 0.8, label, 1.8))
    return out


def _draw(ax, time, s, style, sl=slice(None), label=True):
    ax.plot(time[sl], s.values[sl], "-", lw=style.lw(s.lw), color=s.color,
            alpha=s.alpha, label=s.label if label else None, zorder=s.zorder)
    if s.error is not None:
        ax.fill_between(time[sl], (s.values - s.error)[sl], (s.values + s.error)[sl],
                        color=s.color, alpha=BAND_ALPHA, lw=0, zorder=s.zorder - 0.5)


def _legend(ax, series, style, ncols, loc):
    if not any(not s.label.startswith("_") for s in series):
        return
    leg = ax.legend(fontsize=style.legend, ncols=ncols, loc=loc,
                    facecolor="white", framealpha=0.8, edgecolor="0.8",
                    handlelength=1.2, handletextpad=0.4, borderpad=0.3,
                    labelspacing=0.2, columnspacing=0.8, borderaxespad=0.3)
    for line in leg.get_lines():
        line.set_linewidth(style.lw(1.2))


def _plot_log(time, series, title, ylim, xlim, style, fig_h, xlabel, ylabel,
              ext, ext_label, log_ticks, legend_loc):
    """One panel, log scale."""
    positive = time[time > 0]
    if xlim is not None:
        t_min, t_max = xlim
    else:
        if positive.size == 0:
            raise ValueError("no positive times to show on a log axis")
        t_min, t_max = positive[0], time[-1]
    if t_min <= 0:
        raise ValueError(f"a log axis needs xmin > 0 (got {t_min})")
    if t_max <= t_min:
        raise ValueError(f"xmax ({t_max}) must be > xmin ({t_min})")

    fig, ax = plt.subplots(figsize=(style.fig_w, fig_h))
    ax.set_xscale("log")
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(*ylim)
    _log_minor(ax, enabled=log_ticks is None)
    ax.tick_params(labelsize=style.ticks)
    ax.grid(axis="both", which="major", alpha=0.5, color="gray",
            linestyle="--", linewidth=style.lw(0.5))
    for s in series:
        _draw(ax, time, s, style)

    entries = [(ax, t_min, _fmt_log(t_min), True), (ax, t_max, _fmt_log(t_max), True)]
    entries += _log_tick_entries(ax, t_min, t_max, log_ticks)
    _apply_ticks(fig, entries, style.ticks)

    if ext is not None:
        for tw in _add_ext_axes(fig, [(ax, t_min, t_max, True, ())], ext, style):
            tw.set_xlabel(ext_label, fontsize=style.labels, color=LABEL_COLOR, labelpad=2)

    ax.set_ylabel(ylabel, fontsize=style.labels, color=LABEL_COLOR)
    ax.set_xlabel(xlabel, fontsize=style.labels, color=LABEL_COLOR)
    if title:
        ax.set_title(title, fontsize=style.title)
    _legend(ax, series, style, ncols=1, loc=legend_loc or "best")
    return fig


def _plot_split(time, series, title, ylim, xlim, t_sep, split_pos, lin_ticks,
                style, fig_h, xlabel, ylabel, ext, ext_label, log_ticks, legend_loc):
    """Linear panel (left) up to t_sep, log panel (right) after it."""
    t_min = xlim[0] if xlim is not None else time[0]
    t_max = xlim[1] if xlim is not None else time[-1]
    if not t_min < t_sep < t_max:
        raise ValueError(f"t_sep ({t_sep}) must be between {t_min} and {t_max}")
    if t_sep <= 0:
        raise ValueError(f"t_sep ({t_sep}) must be > 0 for the log panel")
    split_pos = float(np.clip(split_pos, 0.05, 0.95))

    fig = plt.figure(figsize=(style.fig_w, fig_h))
    gs = gridspec.GridSpec(1, 2, width_ratios=[split_pos, 1.0 - split_pos], wspace=0.0)
    ax_lin = fig.add_subplot(gs[0])
    ax_log = fig.add_subplot(gs[1])
    for ax in (ax_lin, ax_log):
        ax.set_ylim(*ylim)
        ax.grid(axis="both", which="major", alpha=0.5, color="gray",
                linestyle="--", linewidth=style.lw(0.5))
        ax.axvline(x=t_sep, color="black", lw=style.lw(2), zorder=5)
        ax.tick_params(labelsize=style.ticks)
    ax_lin.set_xlim(t_min, t_sep)
    ax_lin.spines["right"].set_visible(False)
    ax_log.set_xscale("log")
    ax_log.set_xlim(t_sep, t_max)
    _log_minor(ax_log, enabled=log_ticks is None)
    ax_log.spines["left"].set_visible(False)
    ax_log.yaxis.set_tick_params(left=False, right=False)
    ax_log.set_yticklabels([])

    # Always labelled: the beginning, t_sep (once, on the linear panel), the end.
    entries = [(ax_lin, t_min, _fmt_val(t_min), True),
               (ax_lin, t_sep, _fmt_val(t_sep), True),
               (ax_log, t_sep, "", False),
               (ax_log, t_max, _fmt_log(t_max), True)]
    if lin_ticks is not None:
        for v in np.linspace(t_min, t_sep, int(lin_ticks) + 1)[1:-1]:
            entries.append((ax_lin, float(v), _fmt_val(v), True))
    else:
        loc = ticker.MaxNLocator(nbins=max(1, int(style.fig_w * split_pos / 0.7)),
                                 steps=[1, 2, 2.5, 5, 10])
        span = t_sep - t_min
        for v in loc.tick_values(t_min, t_sep):
            if t_min + 1e-9 * span < v < t_sep - 1e-9 * span:
                entries.append((ax_lin, float(v), _fmt_val(v), False))
    entries += _log_tick_entries(ax_log, t_sep, t_max, log_ticks)
    _apply_ticks(fig, entries, style.ticks)

    # Each panel gets the points up to and one past its range.
    k = max(1, min(int(np.searchsorted(time, t_sep)), len(time) - 1))
    for s in series:
        _draw(ax_lin, time, s, style, slice(None, k + 1), label=False)
        _draw(ax_log, time, s, style, slice(k - 1, None))

    twins = []
    if ext is not None:
        twins = _add_ext_axes(fig, [(ax_lin, t_min, t_sep, True, ("right",)),
                                    (ax_log, t_sep, t_max, False, ("left",))], ext, style)

    # Above: sequence length ticks, its label, the title. Below: time
    # ticks, lin/log arrows, the time label (in points).
    top_ticks = 7 + 1.3 * style.ticks
    bottom_ticks = 7 + 1.5 * style.ticks
    y_arrow = -(bottom_ticks + 3)
    y_xlabel = y_arrow - (2 + 1.3 * style.split_label + 3)
    p0, p1 = ax_lin.get_position(), ax_log.get_position()
    x_center = 0.5 * (p0.x0 + p1.x1)
    top = mtransforms.blended_transform_factory(fig.transFigure, ax_lin.transAxes)

    ax_lin.set_ylabel(ylabel, fontsize=style.labels, color=LABEL_COLOR)
    for ax, label in ((ax_lin, "lin"), (ax_log, "log")):
        tf = _offset(fig, ax.transAxes, y_arrow)
        fig.add_artist(ConnectionPatch(
            xyA=(0.0, 0.0), coordsA=tf, xyB=(1.0, 0.0), coordsB=tf,
            arrowstyle="->", color=LABEL_COLOR, lw=style.lw(0.8), clip_on=False))
        ax.annotate(label, xy=(0.0, 0.0), xycoords=tf, xytext=(2, -2),
                    textcoords="offset points", ha="left", va="top",
                    fontsize=style.split_label, color=LABEL_COLOR, annotation_clip=False)
    fig.text(x_center, 0.0, xlabel, ha="center", va="top", fontsize=style.labels,
             color=LABEL_COLOR, transform=_offset(fig, top, y_xlabel))

    y_next = 4
    if twins:
        fig.text(x_center, 1.0, ext_label, ha="center", va="bottom",
                 fontsize=style.labels, color=LABEL_COLOR,
                 transform=_offset(fig, top, top_ticks + 1))
        y_next = top_ticks + 1 + 1.3 * style.labels + 4
    if title:
        fig.text(x_center, 1.0, title, ha="center", va="bottom",
                 fontsize=style.title, transform=_offset(fig, top, y_next))
    _legend(ax_log, series, style, ncols=2, loc=legend_loc or "upper right")
    return fig


# Command line

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="python -m fuzzyfold.plotting",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot macrostate occupancies over time from .tln or .nxy files.")
    p.add_argument("inputs", nargs="+", metavar="FILE",
                   help="a .nxy file, or .tln files (merged)")
    p.add_argument("-o", "--output", default=None, metavar="BASE",
                   help="output file name without extension (default: the first input)")
    p.add_argument("-f", "--formats", nargs="+", default=["pdf"], metavar="FMT",
                   help="output formats: pdf, svg, png, ..., show")
    p.add_argument("--merged", default=None, metavar="FILE",
                   help="also write the merged .tln inputs to this .tln or .nxy file")
    p.add_argument("--t-sep", type=float, default=None, metavar="T",
                   help="split the x-axis here, in input units (default for .tln: "
                        "the end of the linear output times; for .nxy: log only)")
    p.add_argument("--num-sims", type=int, default=None, metavar="N",
                   help="number of simulations of a .nxy file, for error bands")
    p.add_argument("--no-bands", action="store_true", help="no standard error bands")
    p.add_argument("--t-ext", type=float, default=None, metavar="T",
                   help="time per extension, in input units; with --seq-length, a top "
                        "axis shows the sequence length")
    p.add_argument("--seq-length", type=int, default=None, metavar="N",
                   help="full sequence length")
    p.add_argument("--start-length", type=int, default=None, metavar="N",
                   help="sequence length at time 0 (default 1)")
    p.add_argument("--ext-every", type=int, default=10, metavar="N",
                   help="label every N-th sequence length")
    p.add_argument("--ext-label", default="sequence length [nt]", metavar="STR")
    p.add_argument("--title", default="", metavar="STR")
    p.add_argument("--plim", type=float, default=1e-2,
                   help="minimum peak occupancy to plot a series")
    p.add_argument("--labels", nargs="+", default=None, metavar="NAME",
                   help="series first in the legend, in this order, always plotted")
    p.add_argument("--labels-strict", action="store_true",
                   help="only --labels get a color and a legend entry")
    p.add_argument("--split-pos", type=float, default=0.5, metavar="FRAC",
                   help="fraction of the width for the linear panel")
    p.add_argument("--lin-ticks", type=int, default=None, metavar="N",
                   help="number of evenly spaced ticks on the linear axis")
    p.add_argument("--log-ticks", type=int, default=None, metavar="N",
                   help="number of ticks on the log axis, evenly spaced in log space")
    p.add_argument("--rescale", type=float, default=1.0, metavar="FACTOR",
                   help="multiply all times (also --t-sep and --t-ext) by this factor")
    p.add_argument("--xlabel", default="time [s]", metavar="STR")
    p.add_argument("--ylabel", default="occupancy", metavar="STR")
    p.add_argument("--xlim", nargs=2, type=float, default=None, metavar=("XMIN", "XMAX"),
                   help="x-axis limits, in rescaled units")
    p.add_argument("--ylim", nargs=2, type=float, default=None, metavar=("YMIN", "YMAX"))
    p.add_argument("--legend-loc", default=None, metavar="LOC",
                   help="matplotlib legend location (default: 'upper right' split, "
                        "'best' log only)")
    defaults = Style()
    p.add_argument("--fig-w", type=float, default=defaults.fig_w, metavar="INCHES",
                   help="figure width; line widths scale with it, font sizes do not")
    p.add_argument("--title-size", type=float, default=defaults.title, metavar="PT")
    p.add_argument("--label-size", type=float, default=defaults.labels, metavar="PT")
    p.add_argument("--legend-size", type=float, default=defaults.legend, metavar="PT")
    p.add_argument("--tick-size", type=float, default=defaults.ticks, metavar="PT")
    args = p.parse_args(argv)

    if args.inputs[0].endswith(".nxy"):
        if len(args.inputs) > 1:
            p.error(".nxy files cannot be merged, use .tln files")
        if args.merged:
            p.error("--merged needs .tln inputs")
        data = args.inputs[0]
    else:
        data = Timecourse.read(args.inputs[0]).merge(
            *(Timecourse.read(f) for f in args.inputs[1:]))
        if args.merged:
            data.save(args.merged)
            print(f"Wrote: {args.merged}")

    style = Style(fig_w=args.fig_w, title=args.title_size, labels=args.label_size,
                  legend=args.legend_size, ticks=args.tick_size)
    fig = plot_timecourse(
        data, t_sep=args.t_sep, num_sims=args.num_sims, bands=not args.no_bands,
        t_ext=args.t_ext, seq_length=args.seq_length, start_length=args.start_length,
        title=args.title, plim=args.plim, labels=args.labels,
        labels_strict=args.labels_strict, split_pos=args.split_pos,
        lin_ticks=args.lin_ticks, log_ticks=args.log_ticks,
        xlim=tuple(args.xlim) if args.xlim else None,
        ylim=tuple(args.ylim) if args.ylim else None, rescale=args.rescale,
        xlabel=args.xlabel, ylabel=args.ylabel, ext_label=args.ext_label,
        ext_every=args.ext_every, legend_loc=args.legend_loc, style=style)

    base = args.output or args.inputs[0].rsplit(".", 1)[0]
    for fmt in args.formats:
        if fmt != "show":
            fig.savefig(f"{base}.{fmt}", bbox_inches="tight")
            print(f"Wrote: {base}.{fmt}")
    if "show" in args.formats:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
