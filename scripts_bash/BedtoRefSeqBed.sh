#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 -I <input.bed> -O <output.bed>"
    exit 1
}

# Parse command-line arguments
while getopts ":I:O:" opt; do
    case "${opt}" in
        I)
            INPUT_BED="${OPTARG}"
            ;;
        O)
            OUTPUT_BED="${OPTARG}"
            ;;
        *)
            usage
            ;;
    esac
done

# Check if both input and output files are specified
if [[ -z "$INPUT_BED" || -z "$OUTPUT_BED" ]]; then
    usage
fi

# Use sed to replace chromosome names with RefSeq accession numbers
sed -e 's/^chr1/NC_000067.7/' \
    -e 's/^chr2/NC_000068.8/' \
    -e 's/^chr3/NC_000069.7/' \
    -e 's/^chr4/NC_000070.7/' \
    -e 's/^chr5/NC_000071.7/' \
    -e 's/^chr6/NC_000072.7/' \
    -e 's/^chr7/NC_000073.7/' \
    -e 's/^chr8/NC_000074.7/' \
    -e 's/^chr9/NC_000075.7/' \
    -e 's/^chr10/NC_000076.7/' \
    -e 's/^chr11/NC_000077.7/' \
    -e 's/^chr12/NC_000078.7/' \
    -e 's/^chr13/NC_000079.7/' \
    -e 's/^chr14/NC_000080.7/' \
    -e 's/^chr15/NC_000081.7/' \
    -e 's/^chr16/NC_000082.7/' \
    -e 's/^chr17/NC_000083.7/' \
    -e 's/^chr18/NC_000084.7/' \
    -e 's/^chr19/NC_000085.7/' \
    -e 's/^chrX/NC_000086.8/' \
    -e 's/^chrY/NC_000087.8/' \
    "$INPUT_BED" > "$OUTPUT_BED"

# Notify the user
echo "New BED file with RefSeq accession numbers created: $OUTPUT_BED"

