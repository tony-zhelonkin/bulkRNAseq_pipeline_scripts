#!/bin/bash

# GTFtoBED.sh
# 
# Description:
#   Converts GTF annotation files to BED format. The script can handle
#   both regular GTF files and gzipped GTF files (.gz).
#
# Usage:
#   ./GTFtoBED.sh -i input.gtf[.gz] -o output.bed [-f feature_type]
#
# Parameters:
#   -i  Input GTF file (can be gzipped)
#   -o  Output BED file
#   -f  Feature type to extract (default: exon)
#        Common values: exon, CDS, gene, transcript, etc.
#
# Example:
#   ./GTFtoBED.sh -i gencode.v38.annotation.gtf.gz -o gencode.v38.exons.bed
#   ./GTFtoBED.sh -i gencode.v38.annotation.gtf -o gencode.v38.genes.bed -f gene
#
# Author: Tony Zhelonkin
# Date: 2023

# Usage message function
usage() {
    echo "Usage: $0 -i input.gtf[.gz] -o output.bed [-f feature_type]"
    exit 1
}

# Default feature type
FEATURE_TYPE="exon"

# Parse command-line arguments
while getopts "i:o:f:" opt; do
    case ${opt} in
        i )
            GTF_FILE=$OPTARG
            ;;
        o )
            BED_OUTPUT=$OPTARG
            ;;
        f )
            FEATURE_TYPE=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Ensure both input and output arguments are provided
if [ -z "$GTF_FILE" ] || [ -z "$BED_OUTPUT" ]; then
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

echo "Extracting $FEATURE_TYPE features from GTF and converting to BED format..."

# Convert GTF to BED format
$GTF_COMMAND "$GTF_FILE" | awk -v feature="$FEATURE_TYPE" 'BEGIN {OFS="\t"}
{
    if ($3 == feature) {
        chrom = $1;
        start = $4 - 1;  # BED format is 0-based for start position
        end = $5;
        strand = $7;
        
        # Extract gene_id and gene_name
        gene_id = ""; gene_name = "";
        for (i = 9; i <= NF; i++) {
            if ($i ~ /gene_id/) {
                gene_id = $(i+1); gsub(/"|;/, "", gene_id);
            }
            if ($i ~ /gene_name/) {
                gene_name = $(i+1); gsub(/"|;/, "", gene_name);
            }
        }
        
        # If gene_name is not available, use gene_id as name
        name = (gene_name != "") ? gene_name : gene_id;
        
        # Output in BED format: chrom, start, end, name, score, strand
        # Using "." as a placeholder for score (column 5)
        print chrom, start, end, name, ".", strand;
    }
}' > "$BED_OUTPUT"

# Count the number of features extracted
FEATURE_COUNT=$(wc -l < "$BED_OUTPUT")
echo "BED file created: $BED_OUTPUT with $FEATURE_COUNT $FEATURE_TYPE features"
