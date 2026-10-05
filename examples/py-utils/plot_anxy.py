#!/usr/bin/env python3
#
# plot_anxy.py
#
# Plot annotated nxy files: whitespace-separated with a header row
# giving column names (first column = time, rest = state occupancies).
#
# Features
#   * Pure log x-axis, or a split linear/log x-axis via --t-split/--split-pos.
#   * Log axes get decade major ticks (10^n labels) and 2..9 minor ticks,
#     or exactly N log-spaced ticks via --log-ticks (cf. --lin-ticks).
#   * The beginning, the lin/log split and the end of the time axis are
#     always labelled; other tick labels are dropped only if they would
#     overlap one of those (or each other).
#   * Optional secondary (top) axis for the growing sequence length via
#     --t-ext/--seq-length: length 1 at t=0, +1 every t_ext, a tick per
#     step, labels every --ext-every steps plus the first and full length.
#   * --fig-w sets the figure width; font sizes stay fixed so text remains
#     readable at any size, while line widths scale proportionally.
#
# Example input:
#
#      time    Unassigned    LM1    LM2    LM3
#   0.000e0    0.000e0       0.0    0.0    1.0
#   1.111e-8   ...
#

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.transforms as mtransforms
from matplotlib.patches import ConnectionPatch

# Default figure size in inches (width, height).
DEFAULT_FIG_W = 7.0
DEFAULT_ASPECT = 3.0 / 7.0   # height = width * aspect

# Font sizes (points) — fixed, not scaled with figure width.
TITLE = 12
AXISLABELS = 8      # axis titles: time, sequence length, occupancy
LEGEND = 6.5        # legend entries
SPLIT_LABEL = 8     # 'lin' / 'log' arrow labels
TICKS = 7          # tick labels: time and sequence length

# Colour of the axis titles and the lin/log markers.
LABEL_COLOR = 'black'

# Series with these names (case-insensitive) are drawn light grey instead of
# taking a colour from the cycle, sit behind the other curves, and are listed
# last in the legend (unless --labels places them explicitly).
GREY_NAMES = {'unassigned'}
GREY_COLOR = '0.75'

# Minimum horizontal gap between two neighbouring tick labels (points).
LABEL_PAD_PT = 4
# Minimum spacing between drawn sequence-length minor ticks (points).
EXT_MIN_TICK_PT = 1.5


# ─────────────────────────────────────────────────────────────────────────────
# Small helpers
# ─────────────────────────────────────────────────────────────────────────────

def _s(fig_w, base, ref=DEFAULT_FIG_W):
    """Scale a line/geometry value proportionally to fig_w.
    Font sizes are NOT scaled — they stay readable at any figure size.
    """
    return base * (fig_w / ref)


def _fmt_sci(x):
    """Scientific label: '0', '$10^{n}$' or '$k{\\cdot}10^{n}$' (k up to 3 sig. digits)."""
    if x == 0:
        return '0'
    sign = '-' if x < 0 else ''
    x = abs(x)
    exp = int(np.floor(np.log10(x)))
    coeff = float(f'{x / 10**exp:.3g}')
    if coeff >= 10:            # rounding pushed us to the next decade
        coeff /= 10
        exp += 1
    if coeff == 1:
        return r'$%s10^{%d}$' % (sign, exp)
    return r'$%s%s{\cdot}10^{%d}$' % (sign, f'{coeff:g}', exp)


def _fmt_val(x):
    """Plain number for moderate magnitudes, scientific otherwise."""
    if x == 0:
        return '0'
    if 1e-3 <= abs(x) < 1e4:
        return f'{x:.4g}'
    return _fmt_sci(x)


def _fmt_log(x):
    """Log-axis label: powers of ten as 10^n, everything else via _fmt_val."""
    if x > 0:
        e = np.log10(x)
        if abs(e - round(e)) < 1e-9:
            return r'$10^{%d}$' % round(e)
    return _fmt_val(x)


def _log_candidates(lo, hi):
    """Candidate major ticks strictly inside (lo, hi) for a log axis.

    Decades first (they take label priority).  If fewer than two decades fall
    inside the range, 2·10^n and 5·10^n are added so short ranges still get
    readable ticks.
    """
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


def _log_even(lo, hi, n):
    """n ticks evenly spaced in log space strictly inside (lo, hi).

    Values are rounded to 2 significant digits so tick and label agree.
    """
    vals = np.geomspace(lo, hi, int(n) + 2)[1:-1]
    return [float(f'{v:.2g}') for v in vals]


def _log_tick_entries(ax, lo, hi, n):
    """Tick entries for a log axis: automatic (optional) or n forced ticks."""
    if n is None:
        return [(ax, v, _fmt_log(v), False) for v in _log_candidates(lo, hi)]
    return [(ax, v, _fmt_log(v), True) for v in _log_even(lo, hi, n)
            if lo < v < hi]


def _setup_log_minor(ax, enabled=True):
    """Unlabelled minor ticks at 2..9 × 10^n on a log x-axis, or none."""
    if enabled:
        ax.xaxis.set_minor_locator(
            ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    else:
        ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())


def _offset(fig, base, dy_pt):
    """`base` transform shifted vertically by dy_pt points."""
    return base + mtransforms.ScaledTranslation(0, dy_pt / 72.0,
                                                fig.dpi_scale_trans)


def _legend(ax, fig_w, ncols, loc):
    """Compact legend: small font, short handles, tight spacing."""
    leg = ax.legend(fontsize=LEGEND, ncols=ncols, loc=loc,
                    facecolor='white', framealpha=0.8, edgecolor='0.8',
                    handlelength=1.2, handletextpad=0.4, borderpad=0.3,
                    labelspacing=0.2, columnspacing=0.8, borderaxespad=0.3)
    for line in leg.get_lines():
        line.set_linewidth(_s(fig_w, 1.2))
    return leg


def _has_legend_entries(series_list):
    return any(not s[5].startswith('_') for s in series_list)


# ─────────────────────────────────────────────────────────────────────────────
# Tick-label decluttering
# ─────────────────────────────────────────────────────────────────────────────
#
# A tick "entry" is (axes, x, label, forced).  Forced labels are always shown;
# optional labels are shown in order of appearance if they do not overlap a
# label that is already shown.  Overlap is measured in display space, so this
# works across the two panels of a split axis.

def _declutter(fig, entries, fontsize=None, pad_pt=LABEL_PAD_PT):
    fontsize = TICKS if fontsize is None else fontsize
    renderer = fig.canvas.get_renderer()
    pad = pad_pt * fig.dpi / 72.0
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


def _apply_ticks(entries, keep, fontsize=None):
    """Set fixed major ticks per axes; dropped labels become empty strings."""
    fontsize = TICKS if fontsize is None else fontsize
    by_ax = {}
    for (ax, x, label, _), k in zip(entries, keep):
        by_ax.setdefault(ax, {})
        # Same position twice: keep whichever has a visible label.
        if x not in by_ax[ax] or (k and label):
            by_ax[ax][x] = label if k else ''
    for ax, ticks in by_ax.items():
        xs = sorted(ticks)
        ax.set_xticks(xs)
        ax.set_xticklabels([ticks[x] for x in xs], fontsize=fontsize)


# ─────────────────────────────────────────────────────────────────────────────
# Sequence-length (extension) axis
# ─────────────────────────────────────────────────────────────────────────────

def _thin_minor(ax, L, T, min_pt=EXT_MIN_TICK_PT):
    """Drop per-step minor ticks if they would be packed tighter than min_pt."""
    if len(T) < 2:
        return T
    px = ax.transData.transform(np.column_stack([T, np.zeros_like(T)]))[:, 0]
    limit = min_pt * ax.figure.dpi / 72.0
    for step in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000):
        sel = (L % step == 0)
        if sel.sum() < 2 or np.min(np.diff(px[sel])) >= limit:
            return T[sel]
    return T[L % 1000 == 0]


def _add_ext_axes(fig, panels, ext):
    """Add a top axis showing sequence length on each panel that needs one.

    panels: list of (ax, lo, hi, lo_inclusive, spines_to_hide)
    ext:    (lengths, times, every)
    Returns the list of twin axes created.
    """
    lengths, times, every = ext
    n_full = int(lengths[-1])
    entries, twins = [], []

    for ax, lo, hi, lo_incl, hide in panels:
        tol = 1e-9 * max(abs(lo), abs(hi), 1e-300)
        if lo_incl:
            m = (times >= lo - tol) & (times <= hi + tol)
        else:
            m = (times > lo + tol) & (times <= hi + tol)
        if ax.get_xscale() == 'log':
            m &= times > 0
        if not m.any():
            continue

        tw = ax.twiny()
        tw.set_xscale(ax.get_xscale())
        tw.set_xlim(lo, hi)
        for sp in hide:
            tw.spines[sp].set_visible(False)
        tw.tick_params(axis='x', which='both', labelsize=TICKS)

        L, T = lengths[m], times[m]
        tw.xaxis.set_minor_locator(ticker.FixedLocator(_thin_minor(tw, L, T)))
        tw.xaxis.set_minor_formatter(ticker.NullFormatter())

        for l, t in zip(L, T):
            l = int(l)
            if l == 1 or l == n_full or l % every == 0:
                entries.append((tw, float(t), str(l), l in (1, n_full)))
        twins.append(tw)

    if entries:
        _apply_ticks(entries, _declutter(fig, entries))
    return twins


# ─────────────────────────────────────────────────────────────────────────────
# Input
# ─────────────────────────────────────────────────────────────────────────────

def parse_anxy(stream):
    """Parse an annotated nxy file.

    Returns:
        headers (list[str]): column names (first entry is 'time')
        data    (np.ndarray): shape (n_timepoints, n_columns)
    """
    headers = None
    rows = []
    for line in stream:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if headers is None:
            headers = line.split()
            continue
        rows.append(list(map(float, line.split())))
    if headers is None or not rows:
        raise ValueError("Empty or header-only input.")
    return headers, np.array(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Main plotting entry point
# ─────────────────────────────────────────────────────────────────────────────

def plot_anxy(stream, basename, formats,
              title='',
              plim=1e-2,
              labels=None,
              labels_strict=False,
              t_split=None,
              split_pos=0.5,
              lin_ticks=None,
              log_ticks=None,
              xlim=None,
              ylim=None,
              rescale=1.0,
              xlabel='time [s]',
              ylabel='occupancy',
              fig_w=DEFAULT_FIG_W,
              t_ext=None,
              seq_length=None,
              ext_label='sequence length [nt]',
              ext_every=10,
              tick_size=None,
              title_size=None,
              legend_size=None,
              legend_loc=None,
              label_size=None):
    """Plot an annotated nxy file, optionally with a split linear/log x-axis.

    Args:
        stream:        Readable stream (stdin or file).
        basename:      Output file base name (no extension).
        formats:       List of extensions, e.g. ['pdf', 'png'].
                       Pass [] or ['show'] to display interactively.
        title:         Plot title string.
        plim:          Minimum peak occupancy to plot a series.  Series listed
                       in labels are always plotted regardless of plim.
        labels:        If given, these series come first in the legend, in
                       this order, and are always included (overrides plim).
                       Remaining series that pass plim follow in order of
                       first appearance; 'Unassigned' goes last.
        labels_strict: If True (requires labels), only the listed series get
                       colour and legend entries.  Everything else is plotted
                       thin gray with no label.
        t_split:       If given, split the x-axis: linear in [t_min, t_split],
                       log in [t_split, t_max].  If None, use a pure log scale.
        split_pos:     Fraction of figure width given to the linear panel (0–1).
        lin_ticks:     If given (int), place exactly this many evenly-spaced,
                       always-labelled ticks on the linear axis, starting at
                       t_min (t_split excluded).  If None, automatic.
        log_ticks:     If given (int), place exactly this many always-labelled
                       ticks, evenly spaced in log space, strictly between the
                       start and the end of the log axis (the ends are always
                       labelled anyway).  Minor ticks are switched off.
                       If None, decades with 2..9 minor ticks are used.
        xlim:          (xmin, xmax) override for the combined x range.
                       Applied after rescaling.
        ylim:          (ymin, ymax) override.
        rescale:       Multiply every time value by this factor before plotting.
        xlabel:        X-axis label.
        ylabel:        Y-axis label.
        fig_w:         Figure width in inches.
        t_ext:         Time per sequence extension step, in the units of the
                       input file (rescale is applied to it as well).  If set,
                       a top axis shows the sequence length: 1 at t=0, L at
                       (L-1)*t_ext.
        seq_length:    Full sequence length (required with t_ext).
        ext_label:     Label of the sequence-length axis.
        ext_every:     Label every n-th length (plus 1 and seq_length).
        tick_size:     Font size (pt) of all tick labels (time, sequence
                       length, occupancy).  Default: TICKS.
        title_size:    Font size (pt) of the title.  Default: TITLE.
        legend_size:   Font size (pt) of the legend.  Default: LEGEND.
        label_size:    Font size (pt) of the axis titles.  Default: AXISLABELS.
        legend_loc:    Matplotlib legend location, e.g. 'upper left',
                       'center right', 'best'.  Default: 'upper right' on
                       split plots, 'best' on pure log plots.
    """
    global TICKS, TITLE, LEGEND, AXISLABELS
    if label_size is not None:
        AXISLABELS = label_size
    if legend_size is not None:
        LEGEND = legend_size
    if tick_size is not None:
        TICKS = tick_size
    if title_size is not None:
        TITLE = title_size

    headers, data = parse_anxy(stream)
    time = data[:, 0] * rescale

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    label_set = set(labels) if labels else set()

    all_series = {name: data[:, i] for i, name in enumerate(headers[1:], start=1)}

    if labels:
        missing = [l for l in labels if l not in all_series]
        if missing:
            raise ValueError(f"--labels not found in header: {' '.join(missing)}")

    # Order (= legend order and colour order):
    # 1. --labels first, exactly in the given order, always included.
    # 2. Remaining series that pass plim, in order of first appearance: the
    #    first time point at which their occupancy reaches plim (> 0 if plim
    #    is 0).  Ties keep file order.
    # 3. Grey series (GREY_NAMES, e.g. 'Unassigned') last, unless listed in 1.
    threshold = plim if plim else 0.0

    def first_seen(name):
        vals = all_series[name]
        hit = np.nonzero(vals >= threshold if threshold > 0 else vals > 0)[0]
        return time[hit[0]] if hit.size else np.inf

    explicit = list(labels) if labels else []
    rest = [n for n in headers[1:]
            if n not in explicit and not (plim and all_series[n].max() < plim)]
    rest.sort(key=first_seen)                        # stable: ties → file order
    rest = ([n for n in rest if n.lower() not in GREY_NAMES] +
            [n for n in rest if n.lower() in GREY_NAMES])
    ordered_names = explicit + rest

    # Primary series get colours in order; secondaries are thin gray.
    # Grey-named series are light grey, drawn behind the others, and do not
    # use up a colour.
    series_list = []  # (name, vals, color, lw_s, alpha, legend_name, zorder)
    color_idx = 0
    for name in ordered_names:
        vals = all_series[name]
        is_primary = (not label_set) or (name in label_set) or (not labels_strict)
        if name.lower() in GREY_NAMES:
            series_list.append((name, vals, GREY_COLOR, 1.2, 1.0,
                                name if is_primary else '_nolegend_', 1.9))
        elif is_primary:
            color = colors[color_idx % len(colors)]
            color_idx += 1
            series_list.append((name, vals, color, 2, 1.0, name, 2))
        else:
            series_list.append((name, vals, '0.65', 0.8, 0.8, '_nolegend_', 1.8))

    ymin, ymax = ylim if ylim is not None else (-0.02, 1.02)
    fig_h = fig_w * DEFAULT_ASPECT

    ext = None
    if t_ext is not None:
        if seq_length is None or seq_length < 1:
            raise ValueError("t_ext requires seq_length >= 1")
        if t_ext <= 0:
            raise ValueError("t_ext must be > 0")
        lengths = np.arange(1, int(seq_length) + 1)
        ext = (lengths, (lengths - 1) * t_ext * rescale, max(1, int(ext_every)))

    common = dict(time=time, series_list=series_list, basename=basename,
                  formats=formats, title=title, ymin=ymin, ymax=ymax,
                  xlim=xlim, fig_w=fig_w, fig_h=fig_h, log_ticks=log_ticks,
                  xlabel=xlabel, ylabel=ylabel, ext=ext, ext_label=ext_label,
                  legend_loc=legend_loc)
    if t_split is None:
        _plot_log(**common)
    else:
        _plot_split(t_split=t_split, split_pos=split_pos,
                    lin_ticks=lin_ticks, **common)


def _plot_log(time, series_list, basename, formats, title, ymin, ymax, xlim,
              fig_w, fig_h, xlabel, ylabel, ext, ext_label, log_ticks=None,
              legend_loc=None):
    """Single-panel log-scale plot."""
    positive = time[time > 0]
    if xlim is not None:
        t_min, t_max = xlim
    else:
        if positive.size == 0:
            raise ValueError("No positive time values to show on a log axis.")
        t_min, t_max = positive[0], time[-1]
    if t_min <= 0:
        raise ValueError(f"Log axis needs xmin > 0 (got {t_min}).")
    if t_max <= t_min:
        raise ValueError(f"xmax ({t_max}) must be > xmin ({t_min}).")

    lw_grid = _s(fig_w, 0.5)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xscale('log')
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(ymin, ymax)
    _setup_log_minor(ax, enabled=log_ticks is None)
    ax.tick_params(labelsize=TICKS)
    ax.grid(axis='both', which='major', alpha=0.5,
            color='gray', linestyle='--', linewidth=lw_grid)

    for name, vals, color, lw_s, alpha, legend_name, z in series_list:
        ax.plot(time, vals, '-', lw=_s(fig_w, lw_s), color=color,
                alpha=alpha, label=legend_name, zorder=z)

    # Time ticks: start and end always labelled, decades in between.
    entries = [(ax, t_min, _fmt_log(t_min), True),
               (ax, t_max, _fmt_log(t_max), True)]
    entries += _log_tick_entries(ax, t_min, t_max, log_ticks)
    _apply_ticks(entries, _declutter(fig, entries))

    if ext is not None:
        for tw in _add_ext_axes(fig, [(ax, t_min, t_max, True, ())], ext):
            tw.set_xlabel(ext_label, fontsize=AXISLABELS, color=LABEL_COLOR, labelpad=2)

    ax.set_ylabel(ylabel, fontsize=AXISLABELS, color=LABEL_COLOR)
    ax.set_xlabel(xlabel, fontsize=AXISLABELS, color=LABEL_COLOR)
    if title:
        ax.set_title(title, fontsize=TITLE)
    if _has_legend_entries(series_list):
        _legend(ax, fig_w, ncols=1, loc=legend_loc or 'best')

    _save(fig, basename, formats)


def _plot_split(time, series_list, basename, formats, title, ymin, ymax, xlim,
                t_split, split_pos, lin_ticks, fig_w, fig_h, xlabel, ylabel,
                ext, ext_label, log_ticks=None, legend_loc=None):
    """Split linear (left) / log (right) plot."""
    t_min = xlim[0] if xlim is not None else time[0]
    t_max = xlim[1] if xlim is not None else time[-1]

    if t_split <= t_min:
        raise ValueError(f"t_split ({t_split}) must be > t_min ({t_min})")
    if t_split >= t_max:
        raise ValueError(f"t_split ({t_split}) must be < t_max ({t_max})")
    if t_split <= 0:
        raise ValueError(f"t_split ({t_split}) must be > 0 for the log panel")

    split_pos = float(np.clip(split_pos, 0.05, 0.95))

    lw_grid = _s(fig_w, 0.5)
    lw_split = _s(fig_w, 2)
    lw_arrow = _s(fig_w, 0.8)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, 2, width_ratios=[split_pos, 1.0 - split_pos],
                           wspace=0.0)
    ax_lin = fig.add_subplot(gs[0])
    ax_log = fig.add_subplot(gs[1])

    for ax in (ax_lin, ax_log):
        ax.set_ylim(ymin, ymax)
        ax.grid(axis='both', which='major', alpha=0.5,
                color='gray', linestyle='--', linewidth=lw_grid)
        ax.axvline(x=t_split, color='black', lw=lw_split, zorder=5)
        ax.tick_params(labelsize=TICKS)

    ax_lin.set_xlim(t_min, t_split)
    ax_lin.spines['right'].set_visible(False)

    ax_log.set_xscale('log')
    ax_log.set_xlim(t_split, t_max)
    _setup_log_minor(ax_log, enabled=log_ticks is None)
    ax_log.spines['left'].set_visible(False)
    ax_log.yaxis.set_tick_params(left=False, right=False)
    ax_log.set_yticklabels([])

    # ── Time ticks ───────────────────────────────────────────────────────
    # Forced: beginning, split (labelled once, on the linear panel), end.
    entries = [(ax_lin, t_min, _fmt_val(t_min), True),
               (ax_lin, t_split, _fmt_val(t_split), True),
               (ax_log, t_split, '', False),        # tick only, no 2nd label
               (ax_log, t_max, _fmt_log(t_max), True)]

    if lin_ticks is not None:
        for v in np.linspace(t_min, t_split, int(lin_ticks) + 1)[1:-1]:
            entries.append((ax_lin, float(v), _fmt_val(v), True))
    else:
        width_in = fig_w * split_pos
        loc = ticker.MaxNLocator(nbins=max(1, int(width_in / 0.7)),
                                 steps=[1, 2, 2.5, 5, 10])
        span = t_split - t_min
        for v in loc.tick_values(t_min, t_split):
            if t_min + 1e-9 * span < v < t_split - 1e-9 * span:
                entries.append((ax_lin, float(v), _fmt_val(v), False))

    entries += _log_tick_entries(ax_log, t_split, t_max, log_ticks)
    _apply_ticks(entries, _declutter(fig, entries))

    # ── Plot series ──────────────────────────────────────────────────────
    split_idx = int(np.searchsorted(time, t_split))
    split_idx = max(1, min(split_idx, len(time) - 1))

    for name, vals, color, lw_s, alpha, legend_name, z in series_list:
        ax_lin.plot(time[:split_idx + 1], vals[:split_idx + 1],
                    '-', lw=_s(fig_w, lw_s), color=color, alpha=alpha, zorder=z)
        ax_log.plot(time[split_idx - 1:], vals[split_idx - 1:],
                    '-', lw=_s(fig_w, lw_s), color=color, alpha=alpha,
                    label=legend_name, zorder=z)

    # ── Sequence-length axis (top) ───────────────────────────────────────
    twins = []
    if ext is not None:
        twins = _add_ext_axes(fig, [
            (ax_lin, t_min, t_split, True, ('right',)),
            (ax_log, t_split, t_max, False, ('left',)),
        ], ext)

    # ── Layout of everything above / below the axes (in points) ─────────
    # Above:  [ext tick labels] → [ext axis label, right at the ticks] → title
    # Below:  time tick labels → lin/log arrows with labels → time axis label
    top_ticks = 7 + 1.3 * TICKS               # tick + pad + label height
    bottom_ticks = 7 + 1.5 * TICKS            # taller: 10^n superscripts
    y_arrow = -(bottom_ticks + 3)
    y_xlabel = y_arrow - (2 + 1.3 * SPLIT_LABEL + 3)

    p0, p1 = ax_lin.get_position(), ax_log.get_position()
    x_center = 0.5 * (p0.x0 + p1.x1)
    top = mtransforms.blended_transform_factory(fig.transFigure, ax_lin.transAxes)

    ax_lin.set_ylabel(ylabel, fontsize=AXISLABELS, color=LABEL_COLOR)

    for ax, label in ((ax_lin, 'lin'), (ax_log, 'log')):
        tf = _offset(fig, ax.transAxes, y_arrow)
        fig.add_artist(ConnectionPatch(
            xyA=(0.0, 0.0), coordsA=tf, xyB=(1.0, 0.0), coordsB=tf,
            arrowstyle='->', color=LABEL_COLOR, lw=lw_arrow, clip_on=False))
        ax.annotate(label, xy=(0.0, 0.0), xycoords=tf,
                    xytext=(2, -2), textcoords='offset points',
                    ha='left', va='top', fontsize=SPLIT_LABEL,
                    color=LABEL_COLOR, annotation_clip=False)

    fig.text(x_center, 0.0, xlabel, ha='center', va='top', fontsize=AXISLABELS,
             color=LABEL_COLOR, transform=_offset(fig, top, y_xlabel))

    y_next = 4
    if twins:
        fig.text(x_center, 1.0, ext_label, ha='center', va='bottom',
                 fontsize=AXISLABELS, color=LABEL_COLOR,
                 transform=_offset(fig, top, top_ticks + 1))
        y_next = top_ticks + 1 + 1.3 * AXISLABELS + 4

    if title:
        fig.text(x_center, 1.0, title, ha='center', va='bottom',
                 fontsize=TITLE, transform=_offset(fig, top, y_next))

    if _has_legend_entries(series_list):
        _legend(ax_log, fig_w, ncols=2, loc=legend_loc or 'upper right')

    _save(fig, basename, formats)


def _save(fig, basename, formats):
    """Save fig to each format, or show interactively."""
    written = []
    show = False
    for fmt in formats:
        if fmt == 'show':
            show = True
            continue
        pfile = f'{basename}.{fmt}'
        fig.savefig(pfile, bbox_inches='tight')
        written.append(pfile)
    if show or not formats:
        plt.show()
    plt.close(fig)
    for f in written:
        print(f'Wrote: {f}')


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='plot_anxy: plot annotated nxy occupancy files.')

    parser.add_argument(
        'infile', nargs='?', default='-',
        help='Input file (annotated nxy). Use - or omit for stdin.')
    parser.add_argument(
        '-o', '--output', default='plot_anxy', metavar='PATH',
        help='Output file base name (without extension).')
    parser.add_argument(
        '--title', default='', metavar='STR',
        help='Plot title, shown above everything else. Quote it if it '
             'contains spaces; matplotlib mathtext ($...$) works.')
    parser.add_argument(
        '--title-size', type=float, default=TITLE, metavar='PT',
        help='Title font size in points.')
    parser.add_argument(
        '--label-size', type=float, default=AXISLABELS, metavar='PT',
        help='Font size in points of the axis titles (time, sequence length, '
             'occupancy).')
    parser.add_argument(
        '--legend-size', type=float, default=LEGEND, metavar='PT',
        help='Legend font size in points.')
    parser.add_argument(
        '--legend-loc', default=None, metavar='LOC',
        help="Legend position, e.g. 'upper left', 'center right', 'best'. "
             "Default: 'upper right' (split) / 'best' (log).")
    parser.add_argument(
        '--tick-size', type=float, default=TICKS, metavar='PT',
        help='Font size in points of the tick labels (time, sequence length '
             'and occupancy).')
    parser.add_argument(
        '-f', '--formats', nargs='+', default=['pdf'], metavar='FMT',
        help='Output formats: pdf, svg, png, eps, show.')
    parser.add_argument(
        '--plim', type=float, default=1e-2,
        help='Minimum peak occupancy to plot a series.')
    parser.add_argument(
        '--labels', nargs='+', default=None, metavar='COL',
        help='Series to highlight: listed first in the legend in this order '
             'with full colour, and always included regardless of --plim. '
             'Without it, series are ordered by when they first reach --plim, '
             'with Unassigned last.')
    parser.add_argument(
        '--labels-strict', action='store_true',
        help='With --labels: everything outside --labels is drawn thin gray '
             'with no legend entry.')
    parser.add_argument(
        '--t-split', type=float, default=None, metavar='T',
        help='Split x-axis here (rescaled units): linear left, log right. '
             'If omitted, a single log-scale axis is used.')
    parser.add_argument(
        '--split-pos', type=float, default=0.5, metavar='FRAC',
        help='Fraction of figure width for the linear panel (0–1).')
    parser.add_argument(
        '--lin-ticks', type=int, default=None, metavar='N',
        help='Number of evenly spaced, always-labelled ticks on the linear '
             'axis, starting at t_min (t_split excluded). Automatic if omitted.')
    parser.add_argument(
        '--log-ticks', type=int, default=None, metavar='N',
        help='Number of always-labelled ticks on the log axis, evenly spaced '
             'in log space between its start and end (both ends are always '
             'labelled). Values are rounded to 2 significant digits and '
             'minor ticks are switched off. If omitted, powers of ten '
             'with 2..9 minor ticks are used.')
    parser.add_argument(
        '--rescale', type=float, default=1.0, metavar='FACTOR',
        help='Multiply all time values (and --t-ext) by this factor.')
    parser.add_argument(
        '--xlabel', default='time [s]', metavar='STR',
        help='X-axis label.')
    parser.add_argument(
        '--ylabel', default='occupancy', metavar='STR',
        help='Y-axis label.')
    parser.add_argument(
        '--xlim', nargs=2, type=float, default=None, metavar=('XMIN', 'XMAX'),
        help='X-axis limits (combined range, rescaled units).')
    parser.add_argument(
        '--ylim', nargs=2, type=float, default=None, metavar=('YMIN', 'YMAX'),
        help='Y-axis limits.')
    parser.add_argument(
        '--fig-w', type=float, default=DEFAULT_FIG_W, metavar='INCHES',
        help='Figure width in inches. Height is set automatically. '
             'Line widths scale with this value; font sizes stay fixed.')
    parser.add_argument(
        '--t-ext', type=float, default=None, metavar='T',
        help='Time per extension step, in input-file units. Adds a top axis '
             'with the sequence length (1 at t=0, +1 every T). '
             'Requires --seq-length.')
    parser.add_argument(
        '--seq-length', type=int, default=None, metavar='N',
        help='Full sequence length for the --t-ext axis.')
    parser.add_argument(
        '--ext-every', type=int, default=10, metavar='N',
        help='Label every N-th length on the --t-ext axis '
             '(1 and the full length are always labelled).')
    parser.add_argument(
        '--ext-label', default='sequence length [nt]', metavar='STR',
        help='Label of the --t-ext axis.')

    args = parser.parse_args()
    if args.t_ext is not None and args.seq_length is None:
        parser.error('--t-ext requires --seq-length')

    stream = sys.stdin if args.infile == '-' else open(args.infile)
    try:
        plot_anxy(
            stream,
            basename=args.output,
            formats=args.formats,
            title=args.title,
            plim=args.plim,
            labels=args.labels,
            labels_strict=args.labels_strict,
            t_split=args.t_split,
            split_pos=args.split_pos,
            lin_ticks=args.lin_ticks,
            log_ticks=args.log_ticks,
            xlim=tuple(args.xlim) if args.xlim else None,
            ylim=tuple(args.ylim) if args.ylim else None,
            rescale=args.rescale,
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            fig_w=args.fig_w,
            t_ext=args.t_ext,
            seq_length=args.seq_length,
            ext_label=args.ext_label,
            ext_every=args.ext_every,
            tick_size=args.tick_size,
            title_size=args.title_size,
            legend_size=args.legend_size,
            legend_loc=args.legend_loc,
            label_size=args.label_size,
        )
    finally:
        if args.infile != '-':
            stream.close()


if __name__ == '__main__':
    main()