#!/usr/bin/env bash
set -euo pipefail

# Parameters
N=100    # number of sequences
LENGTHS=(10 50 100 250 500 750 1000 2500 5000)

for L in "${LENGTHS[@]}"; do
    SEQFILE="benchmark_random_sequences_len${L}.fa"
    STRFILE="benchmark_random_structures_len${L}.vrna"

    # Clear or create the output file
    > "$SEQFILE"
    > "$STRFILE"
    
    # Generate and format as FASTA
    for i in $(seq 1 $N); do
        seq=$(ff-randseq --length "$L" --num 1)
        str=$(echo "${seq}" | RNAsubopt -p 1)
        echo ">random_seq_$i" >> "$SEQFILE"
        echo "$seq" >> "$SEQFILE"
        echo ">random_str_$i" >> "$STRFILE"
        echo "$str" >> "$STRFILE"
    done
    echo "✅ Wrote $N sequences of length $L to $SEQFILE"
    echo "✅ Wrote $N sequences of length $L to $STRFILE"
done


