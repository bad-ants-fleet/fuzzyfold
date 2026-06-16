# Reprode the statistics on three-way and four-way shift evaluation.

# NOTE: The following analysis depends on a feature called shift-analysis,
# which you need to specify at compile-time.

## Generate macrostate files from local minima. 
cat dld1_lm1.na | cargo run --bin ff-explore --features ff_kinetics/shift_analysis -- --rna turner2004 --delta 12.00 --three-way-shifts --four-way-shifts > dld1_lm1_d12.mss
cat dld2_lm1.na | cargo run --bin ff-explore --features ff_kinetics/shift_analysis -- --rna turner2004 --delta 12.00 --three-way-shifts --four-way-shifts > dld2_lm1_d12.mss

## Filter the relevant informations.
cat dld1_lm1_d12.mss | grep -v "\." > dld1_lm1_d12.se
cat dld2_lm1_d12.mss | grep -v "\." > dld2_lm1_d12.se

## Plot the two panels for comparison.
cat dld1_lm1_d12.se | python ../../py-utils/shift_analysis.py \
    --col1-label "$\Delta E_1 + \Delta E_2$ [kcal/mol]" \
    --col2-label "$\max(\Delta E_1, \Delta E_2)$ [kcal/mol]" \
    --col3-label "RNA design dld1 ($\Delta E = 12$ kcal/mol)" \
    -o dld1_de12.png

cat dld2_lm1_d12.se | python ../../py-utils/shift_analysis.py \
    --col1-label "$\Delta E_1 + \Delta E_2$ [kcal/mol]" \
    --col2-label "$\max(\Delta E_1, \Delta E_2)$ [kcal/mol]" \
    --col3-label "RNA design dld2 ($\Delta E = 12$ kcal/mol)" \
    -o dld2_de12.png

