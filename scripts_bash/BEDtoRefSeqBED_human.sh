#!/bin/bash

# BedtoRefSeqBed_human.sh
#
# Description:
#   Converts BED files with UCSC-style chromosome names (chr1, chr2, etc.)
#   to RefSeq accession numbers (NC_000001.11, NC_000002.12, etc.) for GRCh38.
#   This is useful when working with tools that require RefSeq chromosome naming.
#
# Usage:
#   ./BedtoRefSeqBed_human.sh -I input.bed -O output.bed
#
# Parameters:
#   -I  Input BED file with UCSC-style chromosome names
#   -O  Output BED file with RefSeq accession numbers
#
# Example:
#   ./BedtoRefSeqBed_human.sh -I ucsc.bed -O refseq.bed
#
# Author: Anton Zhelonkin
# Date: 2025

# Usage function to display help
usage() {
    echo "Usage: $0 -i <input.bed> -o <output.bed>"
    exit 1
}

# Parse command-line arguments
while getopts ":i:o:" opt; do
    case "${opt}" in
        i)
            INPUT_BED="${OPTARG}"
            ;;
        o)
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

# Check if input file exists
if [ ! -f "$INPUT_BED" ]; then
    echo "Error: Input file '$INPUT_BED' does not exist."
    exit 1
fi

# Use sed to replace chromosome names with RefSeq accession numbers for GRCh38
sed -e 's/^chr1\t/NC_000001.11\t/' \
    -e 's/^chr2\t/NC_000002.12\t/' \
    -e 's/^chr3\t/NC_000003.12\t/' \
    -e 's/^chr4\t/NC_000004.12\t/' \
    -e 's/^chr5\t/NC_000005.10\t/' \
    -e 's/^chr6\t/NC_000006.12\t/' \
    -e 's/^chr7\t/NC_000007.14\t/' \
    -e 's/^chr8\t/NC_000008.11\t/' \
    -e 's/^chr9\t/NC_000009.12\t/' \
    -e 's/^chr10\t/NC_000010.11\t/' \
    -e 's/^chr11\t/NC_000011.10\t/' \
    -e 's/^chr12\t/NC_000012.12\t/' \
    -e 's/^chr13\t/NC_000013.11\t/' \
    -e 's/^chr14\t/NC_000014.9\t/' \
    -e 's/^chr15\t/NC_000015.10\t/' \
    -e 's/^chr16\t/NC_000016.10\t/' \
    -e 's/^chr17\t/NC_000017.11\t/' \
    -e 's/^chr18\t/NC_000018.10\t/' \
    -e 's/^chr19\t/NC_000019.10\t/' \
    -e 's/^chr20\t/NC_000020.11\t/' \
    -e 's/^chr21\t/NC_000021.9\t/' \
    -e 's/^chr22\t/NC_000022.11\t/' \
    -e 's/^chrX\t/NC_000023.11\t/' \
    -e 's/^chrY\t/NC_000024.10\t/' \
    -e 's/^chrM\t/NC_012920.1\t/' \
    "$INPUT_BED" > "$OUTPUT_BED"

# Notify the user
echo "New BED file with RefSeq accession numbers created: $OUTPUT_BED"
