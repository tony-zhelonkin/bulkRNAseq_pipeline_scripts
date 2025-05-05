#!/bin/bash

# GTFtoBED12.sh
# 
# Description:
#   Converts GTF annotation files to BED12 format for use with RSeQC tools
#   like infer_experiment.py. The script creates whole-transcript records with
#   exon structure encoded in the block fields, as required by RSeQC.
#
# Usage:
#   ./GTFtoBED12.sh -i input.gtf[.gz] -o output.bed
#
# Parameters:
#   -i  Input GTF file (can be gzipped)
#   -o  Output BED12 file
#
# Example:
#   ./GTFtoBED12.sh -i gencode.v38.annotation.gtf.gz -o gencode.v38.bed12
#
# Author: Anton Zhelonkin
# Date: 2025

# Usage message function
usage() {
    echo "Usage: $0 -i input.gtf[.gz] -o output.bed"
    exit 1
}

# Parse command-line arguments
while getopts "i:o:" opt; do
    case ${opt} in
        i )
            GTF_FILE=$OPTARG
            ;;
        o )
            BED_OUTPUT=$OPTARG
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

echo "Converting GTF to BED12 format for RSeQC..."

# Create a temporary file for intermediate processing
TEMP_FILE=$(mktemp)

# First pass: Extract transcript information and exon/CDS coordinates
$GTF_COMMAND "$GTF_FILE" | awk 'BEGIN {OFS="\t"}
{
    # Extract transcript ID
    transcript_id = "";
    for (i = 9; i <= NF; i++) {
        if ($i ~ /transcript_id/) {
            transcript_id = $(i+1);
            gsub(/"|;/, "", transcript_id);
            break;
        }
    }
    
    if (transcript_id == "") next;
    
    # Extract gene ID and name
    gene_id = ""; gene_name = "";
    for (i = 9; i <= NF; i++) {
        if ($i ~ /gene_id/) {
            gene_id = $(i+1);
            gsub(/"|;/, "", gene_id);
        }
        if ($i ~ /gene_name/) {
            gene_name = $(i+1);
            gsub(/"|;/, "", gene_name);
        }
    }
    
    # Store chromosome and strand information
    chrom = $1;
    strand = $7;
    
    # Store transcript info
    if (!(transcript_id in tx_chrom)) {
        tx_chrom[transcript_id] = chrom;
        tx_strand[transcript_id] = strand;
        tx_gene_id[transcript_id] = gene_id;
        tx_gene_name[transcript_id] = gene_name;
    }
    
    # Process different feature types
    if ($3 == "exon") {
        start = $4;
        end = $5;
        
        # Update transcript boundaries
        if (!(transcript_id in tx_start) || start < tx_start[transcript_id]) {
            tx_start[transcript_id] = start;
        }
        if (!(transcript_id in tx_end) || end > tx_end[transcript_id]) {
            tx_end[transcript_id] = end;
        }
        
        # Store exon coordinates
        if (!(transcript_id in exon_starts)) {
            exon_starts[transcript_id] = start;
            exon_ends[transcript_id] = end;
            exon_count[transcript_id] = 1;
        } else {
            exon_starts[transcript_id] = exon_starts[transcript_id] "," start;
            exon_ends[transcript_id] = exon_ends[transcript_id] "," end;
            exon_count[transcript_id]++;
        }
    } else if ($3 == "CDS") {
        start = $4;
        end = $5;
        
        # Store CDS coordinates for thick start/end
        if (!(transcript_id in cds_start) || start < cds_start[transcript_id]) {
            cds_start[transcript_id] = start;
        }
        if (!(transcript_id in cds_end) || end > cds_end[transcript_id]) {
            cds_end[transcript_id] = end;
        }
    }
}

# Output collected transcript information
END {
    for (tx_id in tx_chrom) {
        # Set thick start/end (CDS region)
        thick_start = (tx_id in cds_start) ? cds_start[tx_id] : tx_start[tx_id];
        thick_end = (tx_id in cds_end) ? cds_end[tx_id] : tx_end[tx_id];
        
        # For non-coding RNAs, set thick start/end equal to transcript start/end
        if (!(tx_id in cds_start)) {
            thick_start = tx_start[tx_id];
            thick_end = tx_start[tx_id];  # Zero-length CDS for non-coding
        }
        
        # Sort exons by start position
        n = split(exon_starts[tx_id], starts, ",");
        split(exon_ends[tx_id], ends, ",");
        
        # Create arrays for sorting
        for (i = 1; i <= n; i++) {
            pos[i] = i;
            start_pos[i] = starts[i];
        }
        
        # Simple bubble sort
        for (i = 1; i <= n; i++) {
            for (j = i + 1; j <= n; j++) {
                if (start_pos[pos[i]] > start_pos[pos[j]]) {
                    temp = pos[i];
                    pos[i] = pos[j];
                    pos[j] = temp;
                }
            }
        }
        
        # Reconstruct sorted exon lists
        sorted_starts = "";
        sorted_ends = "";
        for (i = 1; i <= n; i++) {
            if (i > 1) {
                sorted_starts = sorted_starts ",";
                sorted_ends = sorted_ends ",";
            }
            sorted_starts = sorted_starts starts[pos[i]];
            sorted_ends = sorted_ends ends[pos[i]];
        }
        
        # Use gene_name if available, otherwise gene_id
        name = (tx_gene_name[tx_id] != "") ? tx_gene_name[tx_id] : tx_gene_id[tx_id];
        if (name == "") name = tx_id;
        
        # Output transcript info
        print tx_chrom[tx_id], tx_start[tx_id], tx_end[tx_id], tx_id, 0, tx_strand[tx_id], 
              thick_start, thick_end, "0", exon_count[tx_id], 
              sorted_starts, sorted_ends, name;
    }
}' > "$TEMP_FILE"

# Second pass: Convert to BED12 format (0-based start coordinates)
awk 'BEGIN {OFS="\t"}
{
    # BED is 0-based for start position
    $2 = $2 - 1;
    if ($7 != $2) $7 = $7 - 1;  # Only adjust thick_start if it differs from start
    
    # Convert exon starts and ends to BED12 block format
    split($11, starts, ",");
    split($12, ends, ",");
    
    # Calculate block sizes and relative starts
    block_sizes = "";
    block_starts = "";
    
    for (i = 1; i <= $10; i++) {
        # Calculate block size
        size = ends[i] - starts[i] + 1;
        
        # Calculate relative start position (0-based)
        rel_start = starts[i] - ($2 + 1);
        
        if (i > 1) {
            block_sizes = block_sizes ",";
            block_starts = block_starts ",";
        }
        
        block_sizes = block_sizes size;
        block_starts = block_starts rel_start;
    }
    
    # Output BED12 format
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, block_sizes, block_starts, $13;
}' "$TEMP_FILE" > "$BED_OUTPUT"

# Clean up temporary file
rm "$TEMP_FILE"

# Count the number of transcripts
TRANSCRIPT_COUNT=$(wc -l < "$BED_OUTPUT")
echo "BED12 file created: $BED_OUTPUT with $TRANSCRIPT_COUNT transcripts"
