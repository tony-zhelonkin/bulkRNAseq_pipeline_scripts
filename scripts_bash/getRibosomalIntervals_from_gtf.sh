#!/bin/bash

# Variables
GTF_FILE="/home/mogilenko_lab/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf"  # Path to your GTF annotation file
INTERVAL_LIST_OUTPUT="ribosomal.interval_list"  # Output file for Picard
REFERENCE_GENOME="/home/mogilenko_lab/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.fna"

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

echo "Ribosomal interval list created successfully: $INTERVAL_LIST_OUTPUT"

