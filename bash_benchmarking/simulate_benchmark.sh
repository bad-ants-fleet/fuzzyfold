#!/usr/bin/env bash

set -euo pipefail

# Usage: ./benchmark_folding.sh input1.fa input2.fa ...
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 input_files..."
    echo "Example: $0 rand*.fa"
    exit 1
fi

# Programs to benchmark (space-separated)
PROGRAMS=("Kinfold --fpt --met --time 10 --start" \
          "ff-trajectory --t-end 10")

# Output CSV
RESULTS="simulate_benchmark_results.csv"
echo "program,input_file,num_sequences,elapsed_seconds" > "$RESULTS"

# Iterate over programs and input files
for prog in "${PROGRAMS[@]}"; do
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
        total_time=0
        count=0

        while true; do
            read -r header || break
            read -r seq || break
            read -r struct || break
            count=$((count + 1))

            # Program input: seq on first line, structure on second line
            input="${seq}\n${struct}"
            #echo -e $input

            start=$(date +%s.%N)
            bash -c "$prog <<< \$'${input}'" >/dev/null 2>&1
            end=$(date +%s.%N)
            runtime=$(awk -v s="$start" -v e="$end" 'BEGIN {printf "%.9f", e - s}')

            #runtime=$(/usr/bin/time -f "%e" bash -c "$prog <<< \$'${input}'" 2>&1 > /dev/null)
            total_time=$(awk -v a="$total_time" -v b="$runtime" 'BEGIN {printf "%.6f", a + b}')
            #echo "$runtime"
        done < "$infile"

        echo "$prog,$infile,$count,$total_time" >> "$RESULTS"
    done
done

echo "✅ Benchmark completed. Results in $RESULTS"

