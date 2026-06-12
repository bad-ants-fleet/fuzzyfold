#!/usr/bin/bash

# Reproduce results -- fuzzyfold

# ## Generate macrostate files from local minima. 
# cat dld1_lm1.na | ff-explore --rna turner2004 --delta 3.00 --sorted > dld1_lm1.ms
# cat dld1_lm2.na | ff-explore --rna turner2004 --delta 3.00 --sorted > dld1_lm2.ms
# cat dld1_lm3.na | ff-explore --rna turner2004 --delta 3.00 --sorted > dld1_lm3.ms
# cat dld1_lm4.na | ff-explore --rna turner2004 --delta 3.00 --sorted > dld1_lm4.ms
# cat dld1_lm5.na | ff-explore --rna turner2004 --delta 3.00 --sorted > dld1_lm5.ms
#  
# ## Generate data from stochastic simulations (run repeatedly or adjust number of simulations).
# cat dld1_lm2.na | ff-timecourse --rna turner2004 --macrostates dld1*.ms --k0 1e5 --t-ext 1e-5 --t-end 100 -n 100 --output dld1_lm2_t100
# cat dld1_lm4.na | ff-timecourse --rna turner2004 --macrostates dld1*.ms --k0 1e5 --t-ext 1e-5 --t-end 100 -n 100 --output dld1_lm4_t100
# cat dld1_lm5.na | ff-timecourse --rna turner2004 --macrostates dld1*.ms --k0 1e5 --t-ext 1e-5 --t-end 100 -n 100 --output dld1_lm5_t100

## Publication plots
python ../../py-utils/plot_anxy.py dld1_lm2_t100.nxy --t-split 1e-5 --split-pos 0.125 --title "10,000 aggregated, stochastic simulations" \
    -o dld1_lm2_t100 -f pdf --fig-w 5 --plim 0 --labels Unassigned LM1 LM2 --labels-strict
python ../../py-utils/plot_anxy.py dld1_lm4_t100.nxy --t-split 1e-5 --split-pos 0.125 --title "10,000 aggregated, stochastic simulations" \
    -o dld1_lm4_t100 -f pdf --fig-w 5 --plim 0 --labels Unassigned LM1 LM2 LM3 LM4 --labels-strict
python ../../py-utils/plot_anxy.py dld1_lm5_t100.nxy --t-split 1e-5 --split-pos 0.125 --title "10,000 aggregated, stochastic simulations" \
    -o dld1_lm5_t100 -f pdf --fig-w 5 --plim 0 --labels Unassigned LM1 LM2 LM5 --labels-strict

