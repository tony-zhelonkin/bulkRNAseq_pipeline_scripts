#!/bin/bash

# Define the directory containing your BAM files
BAM_DIR="."  # Adjust this if needed
# Define the directory to store the output for Picard
PICARD_OUTPUT_DIR="./QC/picard_out/dedup"
# Define the directory for the MultiQC report
MULTIQC_DIR="./QC/multiQC/dedup"
# Define the path to the reference genome file (ensure this path is correct)
REFERENCE_GENOME="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.fna"
# Define the path to the refFlat file (ensure this path is correct)
REF_FLAT_FILE="/mouse_genome/ref_genome/GRCm39_refFlat"
# Number of jobs to run in parallel
NUM_JOBS=2

# Create the output directories if they don't exist
mkdir -p "$PICARD_OUTPUT_DIR"
mkdir -p "$MULTIQC_DIR"

# Function to process BAM file: remove duplicates and run RNA-seq QC
process_bam() {
    BAM_FILE="$1"
    BASE_NAME=$(basename "$BAM_FILE" .bam)
    
    DEDUP_BAM="$PICARD_OUTPUT_DIR/${BASE_NAME}_dedup.bam"
    METRICS_FILE="$PICARD_OUTPUT_DIR/${BASE_NAME}_dedup_metrics.txt"
    RNA_METRICS_FILE="$PICARD_OUTPUT_DIR/${BASE_NAME}_rna_metrics.txt"

    echo "Processing $BAM_FILE - Removing Duplicates..."
    
    # Step 1: Remove duplicates with Picard MarkDuplicates
    picard MarkDuplicates \
        I="$BAM_FILE" \
        O="$DEDUP_BAM" \
        M="$METRICS_FILE" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true

    echo "Processing $DEDUP_BAM - Running RNA-seq QC with Picard..."

    # Step 2: Run Picard CollectRnaSeqMetrics on the deduplicated BAM file
    picard CollectRnaSeqMetrics \
        I="$DEDUP_BAM" \
        O="$RNA_METRICS_FILE" \
        R="$REFERENCE_GENOME" \
        REF_FLAT="$REF_FLAT_FILE" \
        STRAND_SPECIFICITY="NONE"
}

export -f process_bam  # Export the function so GNU parallel can access it

# Step 3: Use GNU parallel to process all BAM files in the directory
find "$BAM_DIR" -name "*.bam" | parallel -j "$NUM_JOBS" process_bam {}

# Step 4: Run MultiQC to aggregate all Picard results
echo "Running MultiQC..."
multiqc "$PICARD_OUTPUT_DIR" -o "$MULTIQC_DIR" -n "PostAlignmentPicardStats"

echo "Processing complete. Duplicates removed, QC metrics collected, and reports aggregated."

