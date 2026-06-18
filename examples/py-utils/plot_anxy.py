#!/usr/bin/env python3
#
# plot_anxy.py
#
# Plot annotated nxy files: whitespace-separated with a header row
# giving column names (first column = time, rest = state occupancies).
#
# Supports a split linear/log x-axis via --t-split and --split-pos.
# --fig-w sets the figure width; font sizes stay fixed so text remains
# readable at any size, while line widths scale proportionally.
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
from matplotlib.patches import ConnectionPatch

# Default figure size in inches (width, height).
# All style values are expressed as multiples of width so they scale cleanly.
DEFAULT_FIG_W = 7.0
DEFAULT_ASPECT = 3.0 / 7.0   # height = width * aspect

# Font Sizes
TITLE = 12
AXISLABELS = 10
LEGEND = 8
TICKS = 9

def _s(fig_w, base, ref=DEFAULT_FIG_W):
    """Scale a line/geometry value proportionally to fig_w.
    Font sizes are NOT scaled — they stay readable at any figure size.
    """
    return base * (fig_w / ref)

def _fmt_val(x):
    """Per-tick scientific notation label without a shared multiplier.

    Produces '0', '$10^{n}$', or '$k{\\cdot}10^{n}$'.
    """
    if x == 0:
        return '0'
    if 1e-4 <= abs(x) < 1e4:
        return f'{x:g}'   # plain float: 0.025, 0.1, 1000
    exp = int(np.floor(np.log10(abs(x))))
    coeff = x / 10**exp
    if abs(coeff - 1.0) < 1e-9:
        return r'$10^{%d}$' % exp
    if abs(coeff - round(coeff)) < 1e-9:
        return r'$%d{\cdot}10^{%d}$' % (int(round(coeff)), exp)
    return r'$%.1f{\cdot}10^{%d}$' % (coeff, exp)


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


def plot_anxy(stream, basename, formats,
              title='',
              plim=1e-2,
              labels=None,
              labels_strict=False,
              t_split=None,
              split_pos=0.5,
              lin_ticks=None,
              xlim=None,
              ylim=None,
              rescale=1.0,
              xlabel='time [s]',
              ylabel='occupancy',
              fig_w=DEFAULT_FIG_W):
    """Plot an annotated nxy file, optionally with a split linear/log x-axis.

    Args:
        stream:        Readable stream (stdin or file).
        basename:      Output file base name (no extension).
        formats:       List of extensions, e.g. ['pdf', 'png'].
                       Pass [] or ['show'] to display interactively.
        title:         Plot title string.
        plim:          Minimum peak occupancy to plot a series.  Series listed
                       in labels are always plotted regardless of plim.
        labels:        If given, these series are plotted first in this order
                       and always included (overrides plim).  Remaining series
                       that pass plim follow after.
        labels_strict: If True (requires labels), only the listed series get
                       colour and legend entries.  Everything else is plotted
                       thin gray with no label.
        t_split:       If given, split the x-axis: linear in [t_min, t_split],
                       log in [t_split, t_max].  If None, use a pure log scale.
        split_pos:     Fraction of figure width given to the linear panel (0–1).
        lin_ticks:     If given (int), place exactly this many evenly-spaced ticks
                       on the linear axis (t_split excluded).  If None, automatic.
        xlim:          (xmin, xmax) override for the combined x range.
                       Applied after rescaling.
        ylim:          (ymin, ymax) override.
        rescale:       Multiply every time value by this factor before plotting.
        xlabel:        X-axis label (default: 'time [s]').
        ylabel:        Y-axis label (default: 'occupancy').
        fig_w:         Figure width in inches.
    """
    headers, data = parse_anxy(stream)
    time = data[:, 0] * rescale

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    label_set = set(labels) if labels else set()

    # Build name→values for every column in the file.
    all_series = {}
    for i, name in enumerate(headers[1:], start=1):
        all_series[name] = data[:, i]

    # Decide which series to include and in what order.
    # 1. Labels first (in the given order), always included.
    # 2. Remaining series that pass plim, in file order.
    ordered_names = list(labels) if labels else []
    for name in headers[1:]:
        if name in ordered_names:
            continue
        vals = all_series[name]
        if plim and vals.max() < plim:
            continue
        ordered_names.append(name)

    # Assign styles. Primary (in label_set, or all if no labels given) get
    # colours in order; secondaries get thin gray and no legend label.
    series_list = []  # (name, vals, color, lw_s, alpha, legend_name)
    color_idx = 0
    for name in ordered_names:
        vals = all_series[name]
        is_primary = (not label_set) or (name in label_set) or (not labels_strict)
        if is_primary:
            color = colors[color_idx % len(colors)]
            color_idx += 1
            lw_s, alpha = 2, 1.0
            legend_name = name
        else:
            # secondary: thin gray, no legend entry
            color = '0.65'
            lw_s, alpha = 0.8, 0.8
            legend_name = '_nolegend_'
        series_list.append((name, vals, color, lw_s, alpha, legend_name))

    ymin, ymax = (-0.02, 1.02)
    if ylim is not None:
        ymin, ymax = ylim

    fig_h = fig_w * DEFAULT_ASPECT

    if t_split is None:
        _plot_log(time, series_list, basename, formats,
                  title, ymin, ymax, xlim, fig_w, fig_h, xlabel, ylabel)
    else:
        _plot_split(time, series_list, basename, formats,
                    title, ymin, ymax, xlim, t_split, split_pos, lin_ticks,
                    fig_w, fig_h, xlabel, ylabel)


def _plot_log(time, series_list, basename, formats,
              title, ymin, ymax, xlim, fig_w, fig_h, xlabel, ylabel):
    """Single-panel log-scale plot."""
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xscale('log')
    ax.set_ylim([ymin, ymax])
    if xlim is not None:
        ax.set_xlim(xlim)

    lw = _s(fig_w, 2)
    for name, vals, color, lw_s, alpha, legend_name in series_list:
        ax.plot(time, vals, '-', lw=_s(fig_w, lw_s), color=color,
                alpha=alpha, label=legend_name)

    ax.set_ylabel(ylabel, fontsize=AXISLABELS)
    ax.set_xlabel(xlabel, fontsize=AXISLABELS)
    ax.tick_params(labelsize=TICKS)
    if title:
        ax.set_title(title, fontsize=TITLE)
    ax.legend(facecolor='white', framealpha=1, fontsize=LEGEND)

    _save(fig, basename, formats)


def _plot_split(time, series_list, basename, formats,
                title, ymin, ymax, xlim, t_split, split_pos, lin_ticks,
                fig_w, fig_h, xlabel, ylabel):
    """Split linear (left) / log (right) plot."""
    t_min = xlim[0] if xlim is not None else time[0]
    t_max = xlim[1] if xlim is not None else time[-1]

    if t_split <= t_min:
        raise ValueError(f"t_split ({t_split}) must be > t_min ({t_min})")
    if t_split >= t_max:
        raise ValueError(f"t_split ({t_split}) must be < t_max ({t_max})")

    split_pos = float(np.clip(split_pos, 0.05, 0.95))

    # Line/geometry values scale with fig_w; font sizes stay fixed.
    lw_grid   = _s(fig_w, 0.5)
    lw_split  = _s(fig_w, 2)
    lw_arrow  = _s(fig_w, 0.8)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, 2,
                           width_ratios=[split_pos, 1.0 - split_pos],
                           wspace=0.0)
    ax_lin = fig.add_subplot(gs[0])
    ax_log = fig.add_subplot(gs[1])

    for ax in (ax_lin, ax_log):
        ax.set_ylim([ymin, ymax])
        ax.grid(axis='y', which='major', alpha=0.5,
                color='gray', linestyle='--', linewidth=lw_grid)
        ax.grid(axis='x', which='major', alpha=0.5,
                color='gray', linestyle='--', linewidth=lw_grid)
        ax.axvline(x=t_split, color='black', lw=lw_split, zorder=5)
        ax.tick_params(labelsize=TICKS)

    ax_lin.set_xlim([t_min, t_split])
    ax_lin.spines['right'].set_visible(False)

    ax_log.set_xscale('log')
    ax_log.set_xlim([t_split, t_max])
    ax_log.spines['left'].set_visible(False)
    ax_log.yaxis.set_tick_params(left=False, right=False)
    ax_log.set_yticklabels([])

    # ── Linear axis ticks ────────────────────────────────────────────────
    if lin_ticks is not None:
        ticks = np.linspace(t_min, t_split, lin_ticks + 2)[:-1]
    elif split_pos < 0.20:
        ticks = [t_min] if t_min == 0.0 else []
    else:
        auto_max = 2 if split_pos < 0.35 else 4
        loc = ticker.MaxNLocator(nbins=auto_max, prune=None)
        loc.set_axis(ax_lin.xaxis)
        clearance = (t_split - t_min) * 0.20
        ticks = [v for v in loc.tick_values(t_min, t_split)
                 if t_min - 1e-12 <= v <= t_split - clearance]

    ax_lin.set_xticks(ticks)
    ax_lin.set_xticklabels([_fmt_val(v) for v in ticks], fontsize=TICKS)

    # ── Log axis: hide the leftmost tick label to avoid crowding the split ──
    fig.canvas.draw()
    log_labels = ax_log.xaxis.get_majorticklabels()
    if log_labels:
        log_labels[0].set_visible(False)

    # ── Plot series ──────────────────────────────────────────────────────
    split_idx = int(np.searchsorted(time, t_split))
    split_idx = max(1, min(split_idx, len(time) - 1))

    for name, vals, color, lw_s, alpha, legend_name in series_list:
        ax_lin.plot(time[:split_idx + 1], vals[:split_idx + 1],
                    '-', lw=_s(fig_w, lw_s), color=color, alpha=alpha)
        ax_log.plot(time[split_idx - 1:], vals[split_idx - 1:],
                    '-', lw=_s(fig_w, lw_s), color=color, alpha=alpha,
                    label=legend_name)

    # ── Labels ───────────────────────────────────────────────────────────
    ax_lin.set_ylabel(ylabel, fontsize=AXISLABELS)
    fig.text(0.5, -0.02, xlabel, ha='center', va='top', fontsize=AXISLABELS)

    # Two left-to-right arrows above the plot, labelled at their starts.
    for ax, label in ((ax_lin, 'lin'), (ax_log, 'log')):
        fig.add_artist(ConnectionPatch(
            xyA=(0.0, 1.0), coordsA='axes fraction', axesA=ax,
            xyB=(1.0, 1.0), coordsB='axes fraction', axesB=ax,
            arrowstyle='->', color='dimgray', lw=lw_arrow, clip_on=False,
        ))
        ax.annotate(label, xy=(0.0, 1.0), xycoords='axes fraction',
                    xytext=(2, 4), textcoords='offset points',
                    ha='left', va='bottom', fontsize=LEGEND,
                    color='dimgray', annotation_clip=False)

    ax_log.legend(facecolor='white', framealpha=0.8, ncols=2,
                  loc='upper right', fontsize=LEGEND)
    if title:
        fig.suptitle(title, fontsize=TITLE, y=1.0, va='bottom')

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
        help='Plot title. If omitted, no title is shown.')
    parser.add_argument(
        '-f', '--formats', nargs='+', default=['pdf'], metavar='FMT',
        help='Output formats: pdf, svg, png, eps, show.')
    parser.add_argument(
        '--plim', type=float, default=1e-2,
        help='Minimum peak occupancy to plot a series.')
    parser.add_argument(
        '--labels', nargs='+', default=None, metavar='COL',
        help='Series to highlight: plotted first in this order with full '
             'colour, and always included regardless of --plim.')
    parser.add_argument(
        '--labels-strict', action='store_true',
        help='With --labels: everything outside --labels is drawn thin gray '
             'with no legend entry.')
    parser.add_argument(
        '--t-split', type=float, default=None, metavar='T',
        help='Split x-axis here: linear left, log right. '
             'If omitted, a single log-scale axis is used.')
    parser.add_argument(
        '--split-pos', type=float, default=0.5, metavar='FRAC',
        help='Fraction of figure width for the linear panel (0–1).')
    parser.add_argument(
        '--lin-ticks', type=int, default=None, metavar='N',
        help='Number of ticks on the linear axis (t_split excluded). '
             'If omitted, chosen automatically.')
    parser.add_argument(
        '--rescale', type=float, default=1.0, metavar='FACTOR',
        help='Multiply all time values by this factor before plotting.')
    parser.add_argument(
        '--xlabel', default='time [s]', metavar='STR',
        help='X-axis label.')
    parser.add_argument(
        '--ylabel', default='occupancy', metavar='STR',
        help='Y-axis label.')
    parser.add_argument(
        '--xlim', nargs=2, type=float, default=None, metavar=('XMIN', 'XMAX'),
        help='X-axis limits (combined range).')
    parser.add_argument(
        '--ylim', nargs=2, type=float, default=None, metavar=('YMIN', 'YMAX'),
        help='Y-axis limits.')
    parser.add_argument(
        '--fig-w', type=float, default=DEFAULT_FIG_W, metavar='INCHES',
        help='Figure width in inches. Height is set automatically. '
             'All fonts and line widths scale with this value.')

    args = parser.parse_args()

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
            xlim=tuple(args.xlim) if args.xlim else None,
            ylim=tuple(args.ylim) if args.ylim else None,
            rescale=args.rescale,
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            fig_w=args.fig_w,
        )
    finally:
        if args.infile != '-':
            stream.close()


if __name__ == '__main__':
    main()
