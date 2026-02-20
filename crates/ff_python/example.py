import fuzzyfold as ff

seq = "UGCCUAGAGAGUCAGGUGAU"
db1 = ".((((.((...))))))..."

# Testing energy evaluation.
emodel = ff.ViennaRNA()
print(emodel.energy_of_structure(seq, db1))

# Let's get a standard simulator (25C, why not?).
ssa = ff.Simulator(celsius = 25, k0 = 1)


# Simulate from start stucture:
print(f"{seq}")
for ss, en, gtime, _tinc, _wtime in ssa.simulate(seq, db1, t_end = 4000):
    print(f"{ss} {en/100:6.2f} {gtime:12.8e}")
print(f"{seq}")

# Co-transcriptional mode:
print(f"{seq}")
for ss, en, gtime, _tinc, _wtime in ssa.simulate(seq, None, t_ext = 4000, t_end = 4000):
    print(f"{ss} {en/100:6.2f} {gtime:12.8e}")
print(f"{seq}")

# Co-transcriptional mode with shift moves:
ssa = ff.Simulator(k0 = 1, k3ws = 1, k4ws = 1)
print(f"{seq}")
for ss, en, gtime, _tinc, _wtime in ssa.simulate(seq, None, t_ext = 4000, t_end = 4000):
    print(f"{ss} {en/100:6.2f} {gtime:12.8e}")
print(f"{seq}")


