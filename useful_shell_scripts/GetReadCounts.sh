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
