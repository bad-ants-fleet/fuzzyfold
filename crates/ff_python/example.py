import fuzzyfold as ff

seq = "UGCCUAGAGAGUCAGGUGAU"
db1 = ".((((.((...))))))..."

emodel = ff.ViennaRNA()

print(emodel.energy_of_structure(seq, db1))

ssa = ff.Simulator(k0 = 1)

for line in ssa.simulate(seq, None, t_ext = 4000, t_end = 4000, silent = False):
    print(line)

