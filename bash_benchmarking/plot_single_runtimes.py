#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

tag = "FF_silent_kawasaki"

# Load VKF benchmark files
files = {
    10:     f"simulate_benchmark_{tag}_t10.csv",
    100:    f"simulate_benchmark_{tag}_t100.csv",
    1000:   f"simulate_benchmark_{tag}_t1000.csv",
    10000:  f"simulate_benchmark_{tag}_t10000.csv",
    100000: f"simulate_benchmark_{tag}_t100000.csv",
    1000000:f"simulate_benchmark_{tag}_t1000000.csv",
}

dfs = []
for t_end, path in files.items():
    df = pd.read_csv(path)
    df["t_end"] = t_end
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

# Compute stats per (seq_len, t_end)
summary = (
    df.groupby(["seq_len", "t_end"])["elapsed_seconds"]
    .agg(
        mean="mean",
        median="median",
        q1=lambda x: x.quantile(0.25),
        q3=lambda x: x.quantile(0.75),
        min="min",
        max="max",
    )
    .reset_index()
)

plt.figure(figsize=(8, 5))
color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

for i, (seq_len, g) in enumerate(summary.groupby("seq_len")):
    if seq_len == 10:
        continue;
    color = color_cycle[i % len(color_cycle)]

    # Mean (solid)
    plt.plot(g["t_end"], g["mean"], color=color, linestyle="-", marker="o", label=f"L={seq_len} mean")

    # Median (dashed)
    plt.plot(g["t_end"], g["median"], color=color, linestyle="--", marker=".", label=f"L={seq_len} median")

    # IQR
    plt.vlines(g["t_end"], g["q1"], g["q3"], color=color, linewidth=3)

    # Min–max
    plt.vlines(g["t_end"], g["min"], g["max"], color=color, linewidth=1, linestyles=":")

plt.xscale("log")
plt.yscale("log")
plt.ylim(5e-3, 2e3)
plt.xlim(5e0, 2e6)
plt.xlabel("Simulation time (--t-end)")
plt.ylabel("Wall-clock runtime [s]")
plt.title(f"Runtime {tag} (100 repeats)")
plt.legend(ncol=2, fontsize=8)
plt.savefig(f"FF_{tag}.pdf")
