#!/usr/bin/bash

# Reproduce results -- ViennaRNA-v2.7.0, barriers-v1.7.0, DrTransformer-v1.0

## Generate rate matrix
#tail -n 1 ../dld1.na | RNAsubopt -s -e 15 | barriers --minh 4 --max 6 --rates > dld1.bar

## Simulate rate matrix (using DrSimulate from DrTransformer)
header="time LM1 LM2 LM3 LM4 Unassigned LM5" 
DrSimulate --transpose --p0 2=1 --t1 0.01 --t8 1e7 --t-log 100 -r rates.out > dld1_lm2_t100_bar.nxy
sed -i "1s/^/${header}\n/" dld1_lm2_t100_bar.nxy
python ../../../py-utils/plot_anxy.py dld1_lm2_t100_bar.nxy --rescale 1e-5 --t-split 1e-5 \
    --split-pos 0.125 --title "coarse-grained, deterministic simulation" -o dld1_lm2_t100_bar --fig-w 5 \
    --labels Unassigned LM1 LM2 --plim 0 --labels-strict

DrSimulate --transpose --p0 4=1 --t1 0.01 --t8 1e7 --t-log 100 -r rates.out > dld1_lm4_t100_bar.nxy
sed -i "1s/^/${header}\n/" dld1_lm4_t100_bar.nxy
python ../../../py-utils/plot_anxy.py dld1_lm4_t100_bar.nxy --rescale 1e-5 --t-split 1e-5 \
    --split-pos 0.125 --title "coarse-grained, deterministic simulation" -o dld1_lm4_t100_bar --fig-w 5 \
    --labels Unassigned LM1 LM2 LM3 LM4 --plim 0 --labels-strict

# Note local minimum 6 in the barrier tree corresponds to LM5 in our analysis.
DrSimulate --transpose --p0 6=1 --t1 0.01 --t8 1e7 --t-log 100 -r rates.out > dld1_lm5_t100_bar.nxy
sed -i "1s/^/${header}\n/" dld1_lm5_t100_bar.nxy
python ../../../py-utils/plot_anxy.py dld1_lm5_t100_bar.nxy --rescale 1e-5 --t-split 1e-5 \
    --split-pos 0.125 --title "coarse-grained, deterministic simulation" -o dld1_lm5_t100_bar --fig-w 5 \
    --labels Unassigned LM1 LM2 LM5 --plim 0 --labels-strict


