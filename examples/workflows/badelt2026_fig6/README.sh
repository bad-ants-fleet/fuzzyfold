#!/usr/bin/bash

# Reproduce results -- fuzzyfold

## Generate macrostate files from local minima. 
cat dld2_lm1.na | ff-explore --rna turner2004 --delta 4.00 --sorted > dld2_lm1.ms
cat dld2_lm2.na | ff-explore --rna turner2004 --delta 4.00 --sorted > dld2_lm2.ms
cat dld2_lm3.na | ff-explore --rna turner2004 --delta 4.00 --sorted > dld2_lm3.ms
cat dld2_lm4.na | ff-explore --rna turner2004 --delta 4.00 --sorted > dld2_lm4.ms
cat dld2_lm5.na | ff-explore --rna turner2004 --delta 4.00 --sorted > dld2_lm5.ms
 
## Generate data from stochastic simulations (run repeatedly or adjust number of simulations).
cat dld2_lm4.na | ff-timecourse --rna turner2004 --macrostates dld2*.ms --k0 1e5 --k3ws 0 --k4ws 0 --t-ext 0.1 --t-end 100 -n 100 --output dld2_lm4_t100_k0_1e5_k3ws_0_k4ws_0
cat dld2_lm4.na | ff-timecourse --rna turner2004 --macrostates dld2*.ms --k0 1e5 --k3ws 1e4 --k4ws 0 --t-ext 0.1 --t-end 100 -n 100 --output dld2_lm4_t100_k0_1e5_k3ws_1e4_k4ws_0
cat dld2_lm4.na | ff-timecourse --rna turner2004 --macrostates dld2*.ms --k0 1e5 --k3ws 1e4 --k4ws 1e4 --t-ext 0.1 --t-end 100 -n 100 --output dld2_lm4_t100_k0_1e5_k3ws_1e4_k4ws_1e4

## Publication plots
python ../../py-utils/plot_anxy.py dld2_lm4_t100_k0_1e5_k3ws_0_k4ws_0.nxy --t-split 0.1 --split-pos 0.5 --title "ff-timecourse (open / close)" \
    -o dld2_lm4_t100_k0_1e5_k3ws_0_k4ws_0 -f pdf --fig-w 5 
python ../../py-utils/plot_anxy.py dld2_lm4_t100_k0_1e5_k3ws_1e4_k4ws_0.nxy --t-split 0.1 --split-pos 0.5 --title "ff-timecourse (open / close / 3-way shift)" \
    -o dld2_lm4_t100_k0_1e5_k3ws_1e4_k4ws_0 -f pdf --fig-w 5 
python ../../py-utils/plot_anxy.py dld2_lm4_t100_k0_1e5_k3ws_1e4_k4ws_1e4.nxy --t-split 0.1 --split-pos 0.5 --title "ff-timecourse (open / close / 3-way & 4-way shift)" \
    -o dld2_lm4_t100_k0_1e5_k3ws_1e4_k4ws_1e4 -f pdf --fig-w 5

