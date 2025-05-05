#!/bin/bash

# GTFtoRefFlat.sh
# 
# Description:
#   Converts GTF annotation files to RefFlat format, which is required by some
#   RNA-seq analysis tools like Picard CollectRnaSeqMetrics. The script handles
#   both regular GTF files and gzipped GTF files (.gz).
#
# Usage:
#   ./GTFtoRefFlat.sh -i input.gtf[.gz] -o output.refFlat
#
# Parameters:
#   -i  Input GTF file (can be gzipped)
#   -o  Output RefFlat file
#
# Example:
#   ./GTFtoRefFlat.sh -i gencode.v38.annotation.gtf.gz -o gencode.v38.refFlat
#
# Author: Anton Zhelonkin
# Date: 2025

# Usage message function
usage() {
    echo "Usage: $0 -i input.gtf[.gz] -o output.refFlat"
    exit 1
}

# Parse command-line arguments
while getopts "i:o:" opt; do
    case ${opt} in
        i )
            GTF_FILE=$OPTARG
            ;;
        o )
            REFLAT_OUTPUT=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Ensure both input and output arguments are provided
if [ -z "$GTF_FILE" ] || [ -z "$REFLAT_OUTPUT" ]; then
    usage
fi

# Check if input file exists
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: Input file '$GTF_FILE' does not exist."
    exit 1
fi

# Check if the GTF file is gzipped
if [[ "$GTF_FILE" == *.gz ]]; then
    echo "Detected gzipped GTF file. Processing..."
    GTF_COMMAND="zcat"
else
    echo "Processing GTF file..."
    GTF_COMMAND="cat"
fi

# Extracting fields from the GTF file and formatting them into RefFlat format
$GTF_COMMAND "$GTF_FILE" | awk '
BEGIN {OFS="\t"}
{
    if ($3 == "exon") {
        # Extract relevant fields from GTF
        chrom = $1;  # Chromosome
        start = $4;  # Exon start position
        end = $5;    # Exon end position
        strand = $7; # Strand (+ or -)

        # Parse gene name and transcript name (adjust based on GTF format)
        geneName = ""; transcriptName = "";
        for (i = 9; i <= NF; i++) {
            if ($i ~ /gene_id/) {
                geneName = $(i+1); gsub(/"|;/, "", geneName);
            }
            if ($i ~ /transcript_id/ && $(i+1) != "\"\"") {
                transcriptName = $(i+1); gsub(/"|;/, "", transcriptName);
            }
        }

        # Skip lines where transcript name is empty
        if (transcriptName == "") {
            next;
        }

        # Store information for each transcript
        key = geneName "\t" transcriptName "\t" chrom "\t" strand;

        # Initialize or update transcript details
        if (!(key in transcripts)) {
            transcripts[key] = geneName "\t" transcriptName "\t" chrom "\t" strand "\t" start "\t" end "\t" start "\t" end "\t1\t" start "," "\t" end ",";
        } else {
            split(transcripts[key], fields, "\t");
            txStart = (start < fields[5]) ? start : fields[5];
            txEnd = (end > fields[6]) ? end : fields[6];
            exonCount = fields[9] + 1;
            exonStarts = fields[10] "" start ",";
            exonEnds = fields[11] "" end ",";
            transcripts[key] = geneName "\t" transcriptName "\t" chrom "\t" strand "\t" txStart "\t" txEnd "\t" fields[7] "\t" fields[8] "\t" exonCount "\t" exonStarts "\t" exonEnds;
        }
    }
}

# Output in RefFlat format
END {
    for (key in transcripts) {
        print transcripts[key];
    }
}' > "$REFLAT_OUTPUT"

# Confirm output file creation
echo "RefFlat file created: $REFLAT_OUTPUT"
