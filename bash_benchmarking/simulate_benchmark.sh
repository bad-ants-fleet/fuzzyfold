#!/usr/bin/env bash

set -euo pipefail

# Usage: ./benchmark_folding.sh input1.fa input2.fa ...
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 input_files..."
    echo "Example: $0 rand*.fa"
    exit 1
fi

TIME="10"

KF="Kinfold --noShift --fpt --met --time ${TIME} --start --logML --cut 9999" 
FF="ff-trajectory --t-end ${TIME}"

# Programs to benchmark (space-separated)
PROGRAMS=("$FF")

# Output CSV
RESULTS="simulate_benchmark_results_FF_t${TIME}.csv"
echo "program,input_file,seq_index,seq_len,elapsed_seconds" > "$RESULTS"

# Iterate over programs and input files
for prog in "${PROGRAMS[@]}"; do
    echo $prog
    prog_bin="${prog%% *}"   # take everything before first space

    if ! command -v "$prog_bin" &> /dev/null; then
        echo "⚠️ Skipping $prog_bin (not found in PATH)"
        continue
    fi

    for infile in "$@"; do
        if [[ ! -f $infile ]]; then
            echo "⚠️ Skipping $infile (file not found)"
            continue
        fi
        numseq=$(grep -c '^>' "$infile")
        echo "⏱  Running $prog on $infile ($numseq sequences)..."

        idx=0
        while true; do
            read -r header || break
            read -r seq || break
            read -r struct || break
            idx=$((idx + 1))

            # Program input: seq on first line, structure on second line
            input="${seq}\n${struct}"
            seq_len=${#seq}

            start=$(date +%s.%N)
            echo -e $input | $prog >/dev/null 
            end=$(date +%s.%N)

            runtime=$(awk -v s="$start" -v e="$end" 'BEGIN {printf "%.9f", e - s}')
            echo "$prog,$infile,$idx,$seq_len,$runtime" >> "$RESULTS"
        done < "$infile"
    done
done

echo "✅ Benchmark completed. Results in $RESULTS"

