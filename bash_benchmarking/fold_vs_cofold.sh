#!/usr/bin/env bash

set -euo pipefail

# Usage: ./benchmark_folding.sh input1.fa input2.fa ...
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 input_files..."
    echo "Example: $0 rand*.fa"
    exit 1
fi

TEXT="100000"

# CALL="RNAfold"
# TAG="RF"

CALL="ff-trajectory --k0 1 --t-ext ${TEXT} --t-end ${TEXT} --silent"
TAG="FF_silent_metropolis_cotr"

# Output CSV
RESULTS="simulate_benchmark_${TAG}_t${TEXT}.csv"
echo "program,input_file,seq_index,seq_len,elapsed_seconds" > "$RESULTS"

# Iterate over programs and input files
echo "$TAG: $CALL"

BIN="${CALL%% *}" # take everything before first space
if ! command -v "$BIN" &> /dev/null; then
    echo "$BIN not found in PATH!"
fi

for infile in "$@"; do
    if [[ ! -f $infile ]]; then
        echo "WARNING: no $infile found!"
        continue
    fi
    numseq=$(grep -c '^>' "$infile")
    echo " - Running $CALL on $infile ($numseq sequences)..."

    idx=0
    while true; do
        read -r header || break
        read -r seq || break
        read -r struct || break
        idx=$((idx + 1))

        # Program input: seq on first line, structure on second line
        input="${seq}"
        seq_len=${#seq}

        start=$(date +%s.%N)
        echo -e $input | $CALL >/dev/null 
        end=$(date +%s.%N)

        runtime=$(awk -v s="$start" -v e="$end" 'BEGIN {printf "%.9f", e - s}')
        echo "$CALL,$infile,$idx,$seq_len,$runtime" >> "$RESULTS"
    done < "$infile"
done

echo "✅ Benchmark completed. Results in $RESULTS"

