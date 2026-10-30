"""Macrostate occupancies over time, from many stochastic simulations."""

import json
import numpy as np

UNASSIGNED = "Unassigned"


class Timecourse:
    """The macrostate counts of many simulations at common output times.

    Attributes:
        times:        output times, shape (n,).
        counts:       {(length, name): counts at each time}, one entry per
                      macrostate; co-transcriptional macrostates have one
                      entry per sequence length.
        counter:      number of simulations at each time, shape (n,).
        t_sep:        end of the linear output times (inferred from the
                      times if not given), or None.
        t_ext:        time per extension (co-transcriptional), or None.
        seq_length:   full sequence length, or None.
        start_length: length of the start structure, or None.

    Simulator.simulate_timecourse returns a Timecourse with all attributes
    set.
    """

    def __init__(self, times, counts, counter, t_sep=None, t_ext=None,
                 seq_length=None, start_length=None):
        self.times = np.asarray(times, dtype=float)
        self.counts = {(int(length), str(name)): np.asarray(c, dtype=np.int64)
                       for (length, name), c in counts.items()}
        self.counter = np.asarray(counter, dtype=np.int64)
        self.t_sep = t_sep if t_sep is not None else infer_t_sep(self.times)
        self.t_ext = t_ext
        self.seq_length = seq_length
        self.start_length = start_length
        n = len(self.times)
        if self.counter.shape != (n,) or any(c.shape != (n,) for c in self.counts.values()):
            raise ValueError("times, counts and counter must have the same length")

    @property
    def names(self):
        """Macrostate names, Unassigned first, then in the order of counts."""
        names = list(dict.fromkeys(name for _, name in self.counts))
        if UNASSIGNED in names:
            names.remove(UNASSIGNED)
            names.insert(0, UNASSIGNED)
        return names

    @property
    def num_sims(self):
        """The number of simulations."""
        return int(self.counter.max()) if len(self.counter) else 0

    def occupancy(self):
        """Occupancy of each macrostate (summed over lengths), shape (n, len(names))."""
        return self._name_counts() / np.maximum(self.counter, 1)[:, None]

    def standard_error(self):
        """Standard error of the occupancies, sqrt(p (1 - p) / counter)."""
        p = self.occupancy()
        return np.sqrt(p * (1 - p) / np.maximum(self.counter, 1)[:, None])

    def _name_counts(self):
        names = self.names
        out = np.zeros((len(self.times), len(names)), dtype=np.int64)
        for (_, name), c in self.counts.items():
            out[:, names.index(name)] += c
        return out

    def merge(self, *others):
        """A new Timecourse with the counts of this one and others."""
        merged = Timecourse(self.times, self.counts, self.counter, self.t_sep,
                            self.t_ext, self.seq_length, self.start_length)
        for other in others:
            if len(other.times) != len(merged.times) or not np.allclose(
                    other.times, merged.times, rtol=1e-9, atol=0):
                raise ValueError("cannot merge time courses with different output times")
            for attr in ("t_ext", "seq_length", "start_length"):
                a, b = getattr(merged, attr), getattr(other, attr)
                if a is not None and b is not None and a != b:
                    raise ValueError(f"cannot merge time courses with different {attr}")
                if a is None:
                    setattr(merged, attr, b)
            merged.counter = merged.counter + other.counter
            for key, c in other.counts.items():
                merged.counts[key] = merged.counts.get(key, 0) + c
        return merged

    @classmethod
    def read(cls, path, **attributes):
        """Read a .tln file (as written by ff-timecourse).

        Further attributes (t_sep, t_ext, seq_length, start_length) can be
        given.
        """
        path = str(path)
        if not path.endswith(".tln"):
            raise ValueError(f"{path}: expected a .tln file (.nxy files have no counts)")
        with open(path) as f:
            points = json.load(f)["points"]
        times = [p["time"] for p in points]
        counter = [p["counter"] for p in points]
        counts = {}
        for i, p in enumerate(points):
            for length, name, c in p["ensemble"]:
                counts.setdefault((length, name), np.zeros(len(points), dtype=np.int64))[i] += c
        counts = dict(sorted(counts.items(), key=lambda kv: (kv[0][1], kv[0][0])))
        return cls(times, counts, counter, **attributes)

    def save(self, path):
        """Write a .tln or .nxy file (the formats of ff-timecourse)."""
        path = str(path)
        if path.endswith(".tln"):
            points = [{
                "time": float(t),
                "ensemble": [[length, name, int(c[i])]
                             for (length, name), c in self.counts.items() if c[i] > 0],
                "counter": int(self.counter[i]),
            } for i, t in enumerate(self.times)]
            text = json.dumps({"points": points}, indent=2)
        elif path.endswith(".nxy"):
            occupancy = self.occupancy()
            lines = [f"{'time':>13}" + "".join(f" {n:>13}" for n in self.names)]
            for t, row in zip(self.times, occupancy):
                lines.append(_exp(t) + "".join(" " + _exp(x) for x in row))
            text = "\n".join(lines) + "\n"
        else:
            raise ValueError(f"{path}: expected a .tln or .nxy file")
        with open(path, "w") as f:
            f.write(text)

    def __repr__(self):
        return (f"Timecourse({len(self.times)} times, {len(self.names)} macrostates, "
                f"{self.num_sims} simulations)")


def infer_t_sep(times):
    """The end of the evenly spaced times from 0, or None if all are."""
    times = np.asarray(times, dtype=float)
    if len(times) < 3 or times[0] != 0:
        return None
    step = np.diff(times)
    even = np.abs(step - step[0]) <= 1e-6 * step[0]
    if even.all():
        return None
    return float(times[np.argmin(even)])


def _exp(x):
    """x as {:13.6e} in Rust: no '+' and no leading zeros in the exponent."""
    mantissa, exponent = f"{x:.6e}".split("e")
    return f"{mantissa}e{int(exponent)}".rjust(13)
