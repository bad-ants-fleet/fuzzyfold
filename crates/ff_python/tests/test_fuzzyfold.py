"""Tests of the Python bindings: energies, structures and kinetics.

    cd crates/ff_python && maturin develop --release && pytest tests
"""

import pytest

import fuzzyfold as ff

HAIRPIN_SEQ = "GGGGAAACCCC"
HAIRPIN = "((((...))))"
OPEN = "..........."


# Energies

def test_viennarna_default_energy():
    # Turner 2004 at 37 C, as RNAeval -d2.
    assert ff.ViennaRNA().energy_of_structure(HAIRPIN_SEQ, HAIRPIN) == pytest.approx(-4.50)


def test_pseudoknot_energy_and_terms():
    model = ff.Fuzzyfold()
    seq, pk = "GCGAUUUCUGACCGCUUUUUUGUCAG", "[[[....(((((]]]......)))))"
    energy = model.energy_of_structure(seq, pk)
    assert energy == pytest.approx(-3.57)
    terms = model.energy_terms(seq, pk)
    assert sum(e for level, _, e in terms if level == 0) == pytest.approx(energy)


def test_strand_order_does_not_matter():
    model = ff.Fuzzyfold()
    a = model.energy_of_structure("GGGGAAA+UUUCCCC", "((((...+...))))")
    b = model.energy_of_structure("UUUCCCC+GGGGAAA", "...((((+))))...")
    assert a == pytest.approx(b)


def test_invalid_structure_raises():
    with pytest.raises(ValueError):
        ff.Fuzzyfold().energy_of_structure(HAIRPIN_SEQ, "((((...)))")


# Structures

def test_structure_functions():
    s = "((..[[..))..]]"
    assert ff.pairs(s) == [(0, 9), (1, 8), (4, 13), (5, 12)]
    assert ff.normalize("[(.)]") == "((.))"
    assert ff.density(s) == 2


# Kinetics

def test_trajectory():
    sim = ff.Simulator(k0=1e5)
    states = list(sim.simulate_trajectory(HAIRPIN_SEQ, start=OPEN, t_end=1e-3))
    assert states[0][0] == OPEN
    # (structure, energy, time, waiting time, flux); the times add up.
    for (_, _, t, tinc, _), (_, _, t_next, _, _) in zip(states, states[1:]):
        assert t + tinc == pytest.approx(t_next)
    assert states[-1][2] + states[-1][3] == pytest.approx(1e-3)


def test_cotranscriptional_trajectory_reaches_full_length():
    sim = ff.Simulator(k0=1e5)
    states = list(sim.simulate_trajectory(HAIRPIN_SEQ, t_ext=1e-4, t_end=1e-3))
    assert len(states[0][0]) < len(HAIRPIN_SEQ)
    assert len(states[-1][0]) == len(HAIRPIN_SEQ)


def test_ensemble_counts_all_simulations():
    times = ff.lin_log_times(t_sep=1e-5, t_end=1e-2, t_lin=5, t_log=5)
    ensemble = ff.Simulator().simulate_ensemble(HAIRPIN_SEQ, start=OPEN, times=times, num_sims=20)
    assert [t for t, _ in ensemble] == pytest.approx(times)
    assert all(sum(counts.values()) == 20 for _, counts in ensemble)


def test_timecourse_occupancies(tmp_path):
    ms = tmp_path / "hairpin.ms"
    ms.write_text(f">hairpin\n{HAIRPIN_SEQ}\n{HAIRPIN}\n.(((...))).\n")
    tc = ff.Simulator().simulate_timecourse(HAIRPIN_SEQ, start=OPEN, num_sims=10, macrostates=[ms])
    assert tc.names[:2] == ["Unassigned", "hairpin"]
    assert tc.occupancy().sum(axis=1) == pytest.approx(1.0)
    # Written and read again.
    tc.save(tmp_path / "run.tln")
    again = ff.Timecourse.read(tmp_path / "run.tln")
    assert again.occupancy() == pytest.approx(tc.occupancy())


def test_a_to_i_edit():
    # The A at 3 opposite the C at 14 in a 1x1 interior loop. With a rate
    # far above all others, the edit is the first move.
    seq, start = "GCGACGCAAAAGCGCCGC", "(((.(((....))).)))"
    sim = ff.Simulator(k0=1.0, k_atoi=1e12)
    trajectory = sim.simulate_trajectory(seq, start=start, t_end=1e-6)
    states = list(trajectory)
    assert states[1][0] == "(((((((....)))))))"
    assert trajectory.sequence == "GCGICGCAAAAGCGCCGC"


def test_explorer_includes_start():
    neighbors = ff.Explorer().explore(HAIRPIN_SEQ, HAIRPIN, delta=3.0)
    structures = [s for s, _ in neighbors]
    assert HAIRPIN in structures
    assert all(e <= neighbors[structures.index(HAIRPIN)][1] + 3.0 + 1e-9 for _, e in neighbors)
