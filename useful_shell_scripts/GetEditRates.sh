#!/bin/bash

set -euo pipefail

input=""
output=""

while getopts "i:o:" flag; do
    case "$flag" in
        i) input=$OPTARG ;;
        o) output=$OPTARG ;;
        ?) echo "Usage ./GetEditRates.sh -i input_directory -o output_directory" ; exit 1 ;;
    esac
done

if [[ -z "$input" || -z "$output" ]]; then
    echo "Usage ./GetEditRates.sh -i input_directory -o output_directory"
    exit 1
fi



check_file() {
    local path=$1

    if ! [[  -f "$path" || -d "$path" ]]; then
        echo "${path} not found"
        exit 1
    fi
}

check_files() {
    for thing in "$@"; do
        check_file "$thing"
    done
}

log () {
    echo "[$(date '+%H:%M:%S')] $1"
}

inputs=("$input" "$output")

log "Checking inputs..."
check_files "${inputs[@]}"
log "Inputs checked"

log "Grabbing gene name and date from input"
gene_name=$(echo "$input" | cut -d'.' -f2)
sample_date=$(echo "$input" | cut -d'/' -f8 | cut -d'.' -f1)

input_directory="${input}/${gene_name}_X*_R*_D05*readstats.tsv"
output_file="${output}/${gene_name}.editrates.${sample_date}.tsv"

log "Writing edit rate file"
for f in $input_directory; do
    tail -1 "$f" | awk -v OFS='\t' '{ print $1, ($10 + $11) / $3 }' >> "$output_file"
done
log "Done!"

