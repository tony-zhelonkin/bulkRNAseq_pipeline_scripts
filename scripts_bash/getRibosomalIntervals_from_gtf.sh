#!/bin/bash

# getRibosomalIntervals_from_gtf.sh
# 
# Description:
#   Extracts ribosomal RNA intervals from a GTF annotation file and creates a
#   Picard-compatible interval list file. This is useful for RNA-seq analysis
#   to identify and potentially filter out ribosomal RNA reads.
#
# Usage:
#   ./getRibosomalIntervals_from_gtf.sh -i input.gtf -r reference.fna[.gz] -o output.interval_list
#
# Parameters:
#   -i  Input GTF annotation file
#   -r  Reference genome FASTA file (can be gzipped)
#   -o  Output interval list file
#
# Example:
#   ./getRibosomalIntervals_from_gtf.sh -i GRCm39_genomic.gtf -r GRCm39_genomic.fna -o ribosomal.interval_list
#
# Author: Tony Zhelonkin
# Date: 2023

# Usage message function
usage() {
    echo "Usage: $0 -i input.gtf -r reference.fna[.gz] -o output.interval_list"
    exit 1
}

# Parse command-line arguments
while getopts "i:r:o:" opt; do
    case ${opt} in
        i )
            GTF_FILE=$OPTARG
            ;;
        r )
            REFERENCE_GENOME=$OPTARG
            ;;
        o )
            INTERVAL_LIST_OUTPUT=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Ensure all required arguments are provided
if [ -z "$GTF_FILE" ] || [ -z "$REFERENCE_GENOME" ] || [ -z "$INTERVAL_LIST_OUTPUT" ]; then
    usage
fi

# Check if input files exist
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file '$GTF_FILE' does not exist."
    exit 1
fi

if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "Error: Reference genome file '$REFERENCE_GENOME' does not exist."
    exit 1
fi

# Check if the reference genome is gzipped
if [[ "$REFERENCE_GENOME" == *.gz ]]; then
    echo "Detected gzipped reference genome. Uncompressing for processing..."
    UNCOMPRESSED_REF="${REFERENCE_GENOME%.gz}"
    gunzip -c "$REFERENCE_GENOME" > "$UNCOMPRESSED_REF"
    REFERENCE_GENOME="$UNCOMPRESSED_REF"
    CLEANUP_NEEDED=true
else
    CLEANUP_NEEDED=false
fi

# Step 1: Generate the header with chromosome sizes from the .fai file
# Check if .fai index exists, if not, create it
if [ ! -f "$REFERENCE_GENOME.fai" ]; then
    echo "FASTA index (.fai) not found. Generating index..."
    samtools faidx "$REFERENCE_GENOME"
fi

# Step 2: Add header to the output interval list file
echo "@HD    VN:1.5  SO:coordinate" > "$INTERVAL_LIST_OUTPUT"

# Generate @SQ lines from the .fai index and append them to the interval list
awk '{print "@SQ\tSN:"$1"\tLN:"$2}' "$REFERENCE_GENOME.fai" >> "$INTERVAL_LIST_OUTPUT"

# Step 3: Extract ribosomal RNA intervals from the GTF file
# The GTF file is expected to have "rRNA" annotations in the feature field or gene names
echo "Extracting ribosomal RNA intervals from GTF file..."

grep -i "rRNA" "$GTF_FILE" | awk 'BEGIN {OFS="\t"} 
{
    if ($3 == "exon") {
        chr = $1;         # Chromosome
        start = $4;       # Start position
        end = $5;         # End position
        strand = $7;      # Strand (+/-)
        
        # Extract the gene/transcript name (adjust according to the GTF format)
        gene_id = ""; 
        for(i=9; i<=NF; i++) {
            if($i ~ /gene_name/ || $i ~ /gene_id/) {
                gene_id = $(i+1); 
                gsub(/[";]/, "", gene_id);  # Remove quotes and semicolons
                break;
            }
        }
        if(gene_id == "") { gene_id = "rRNA"; }  # Default to rRNA if no gene name found

        # Output in Picard interval list format (chr, start, end, strand, gene_name)
        print chr, start, end, strand, gene_id;
    }
}' >> "$INTERVAL_LIST_OUTPUT"

# Clean up temporary uncompressed reference if needed
if [ "$CLEANUP_NEEDED" = true ]; then
    echo "Cleaning up temporary uncompressed reference..."
    rm "$REFERENCE_GENOME"
    # Also remove the index file if we created it for the temporary file
    if [ -f "$REFERENCE_GENOME.fai" ]; then
        rm "$REFERENCE_GENOME.fai"
    fi
fi

echo "Ribosomal interval list created successfully: $INTERVAL_LIST_OUTPUT"

