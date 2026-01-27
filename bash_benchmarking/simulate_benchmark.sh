#!/usr/bin/env bash

set -euo pipefail

# Usage: ./benchmark_folding.sh input1.fa input2.fa ...
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 input_files..."
    echo "Example: $0 rand*.fa"
    exit 1
fi

TIME="1000"

CALL="ff-trajectory --k0 1 --ks 1 --t-end ${TIME} "
TAG="FF_metropolis"

# CALL="ff-trajectory --k0 1 --ks 1 --t-end ${TIME} --silent "
# TAG="FF_silent_metropolis"
# 
# CALL="ff-trajectory --k0 1 --ks 1 --t-end ${TIME} --silent --model kawasaki"
# TAG="FF_silent_kawasaki"
# 
# CALL="ff-trajectory --k0 1 --ks 1 --t-end ${TIME} --silent --shift-moves "
# TAG="FF_silent_metropolis_shift"
#   
# CALL="ff-trajectory --k0 1 --ks 1 --t-end ${TIME} --silent --model kawasaki --shift-moves"
# TAG="FF_silent_kawasaki_shift"
 
# CALL="Kinfold --noShift --fpt --met --time ${TIME} --start --logML --cut 99999" 
# TAG="KF_metropolis"
 
# CALL="Kinfold --noShift --fpt --met --time ${TIME} --start --logML --silent" 
# TAG="KF_silent_metropolis"

# Output CSV
RESULTS="simulate_benchmark_${TAG}_t${TIME}.csv"
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
        input="${seq}\n${struct}"
        seq_len=${#seq}

        start=$(date +%s.%N)
        echo -e $input | $CALL >/dev/null 
        end=$(date +%s.%N)

        runtime=$(awk -v s="$start" -v e="$end" 'BEGIN {printf "%.9f", e - s}')
        echo "$CALL,$infile,$idx,$seq_len,$runtime" >> "$RESULTS"
    done < "$infile"
done

echo "✅ Benchmark completed. Results in $RESULTS"

