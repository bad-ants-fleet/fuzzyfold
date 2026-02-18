import fuzzyfold as ff

seq = "UGCCUAGAGAGUCAGGUGAU"
db1 = ".((((.((...))))))."

#emodel = ff.ViennaRNA()
#print(emodel.energy_of_structure(seq, db1))

ssa = ff.Simulator(k0 = 1)

print(f"{seq}")
for i, (ss, en, gt, ti, wt) in enumerate(ssa.simulate(seq, None, t_ext = 4, t_end = 4)):
    print(f"{ss} {en/100:6.2f} {gt:12.8e} {ti:12.8e}")
print(f"{seq}")

# for line in ssa.simulate(seq, None, t_ext = 4000, t_end = 4000):
#     print(line)

