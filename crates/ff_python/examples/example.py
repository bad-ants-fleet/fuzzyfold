import fuzzyfold as ff

# Default, ViennaRNA with Turner 2004 parameters.
# Convenient for comparisons to ViennaRNA's thermodynamic folding.

seq = "UGCCUAGAGAGUCAGGUGAU"
db1 = ".((((.((...))))))..."

# Testing energy evaluation.
emodel = ff.ViennaRNA()
print("{}\n{} {:>6.2f}".format(seq, db1, emodel.energy_of_structure(seq, db1)))

emodel = ff.ViennaRNA(celsius = 25)
print("{}\n{} {:>6.2f}".format(seq, db1, emodel.energy_of_structure(seq, db1)))

# Using an expanded parameterset with special hairpin energies and pseudouridine.
seq = "UGCCUAGAGAGPCAGGPGAP"
db1 = ".((((.((...))))))..."

emodel = ff.ViennaRNA(params = "rna_extended", celsius = 25)
print("{}\n{} {:>6.2f}".format(seq, db1, emodel.energy_of_structure(seq, db1)))


# Simulate from start stucture:
seq = "UGCCUAGAGAGUCAGGUGAU"
db1 = ".((((.((...))))))..."

ssa = ff.Simulator(k0 = 1)
print(f"{seq}")
for ss, en, gtime, _tinc, _wtime in ssa.simulate(seq, db1, t_end = 40):
    print(f"{ss} {en/100:6.2f} {gtime:12.8e}")
print(f"{seq}")

# Co-transcriptional mode:
print(f"{seq}")
for ss, en, gtime, _tinc, _wtime in ssa.simulate(seq, None, t_ext = 40, t_end = 40):
    print(f"{ss} {en/100:6.2f} {gtime:12.8e}")
print(f"{seq}")

# Co-transcriptional mode modfications and shift moves:
seq = "UGCCUAGAGAGPCAGGPGAU"
db1 = ".((((.((...))))))..."
ssa = ff.Simulator(params = "rna_extended", k0 = 1, k3ws = 1, k4ws = 1)
print(f"{seq}")
for ss, en, gtime, _tinc, _wtime in ssa.simulate(seq, None, t_ext = 40, t_end = 40):
    print(f"{ss} {en/100:6.2f} {gtime:12.8e}")
print(f"{seq}")

