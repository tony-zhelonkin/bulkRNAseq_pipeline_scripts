#!/bin/bash

# TMP dir set to docker user accessible location. May be commented if user::sudo
export TMPDIR="/workspaces/Yasmine-retroT/.tmp"

# Number of parallel jobs (for GNU parallel)
JOB_NUMBER=6  # Number of parallel Picard, RSeQC, and Samtools processes to run

# Directories (use absolute paths)
BAM_DIR="."  # Directory containing BAM files
PICARD_OUTPUT_DIR="../post_align_QC/picard"  # Output directory for Picard metrics
RSEQ_OUTPUT_DIR="../post_align_QC/rseqc"  # Output directory for RSeQC metrics
SAMTOOLS_OUTPUT_DIR="../post_align_QC/samtools"  # Output directory for Samtools flagstat
MULTIQC_DIR="../post_align_QC/multiQC"  # Output directory for MultiQC report
LOG_DIR="../post_align_QC/logs"  # Directory to store logs

# Define the reference genome and refFlat file (make sure these paths are correct)
REFERENCE_GENOME="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GCF_000001635.27_GRCm39_genomic.fna"  # Path to reference genome
REF_FLAT_FILE="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GCF_000001635.27_GRCm39_genomic.reflat"  # Path to refFlat file for Picard
BED_ANNOTATION_FILE="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GRCm39_RfSeq.bed"  # BED file for RSeQC infer_experiment.py
RIBOSOMAL_INTERVALS="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/ribosomal.interval_list" # Ribosomal intervals list for Picard estimation


# Create output directories if they don't exist
mkdir -p "$PICARD_OUTPUT_DIR" "$RSEQ_OUTPUT_DIR" "$SAMTOOLS_OUTPUT_DIR" "$MULTIQC_DIR" "$LOG_DIR"

run_picard() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local log_file="$LOG_DIR/${base_name}_picard.log"
    echo "Running Picard CollectRnaSeqMetrics on $bam_file..."
    picard CollectRnaSeqMetrics \
        I="$bam_file" \
        O="$PICARD_OUTPUT_DIR/${base_name}_rna_metrics.txt" \
        R="$REFERENCE_GENOME" \
        REF_FLAT="$REF_FLAT_FILE" \
        STRAND_SPECIFICITY="NONE" \
        RIBOSOMAL_INTERVALS="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/ribosomal.interval_list" \
        > "$log_file" 2>&1
    if [[ $? -ne 0 ]]; then
        echo "Picard CollectRnaSeqMetrics failed for $bam_file. Check log for details: $log_file"
    fi
}

# Function to run Picard MarkDuplicates
run_picard_mark_duplicates() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local output_file="$PICARD_OUTPUT_DIR/${base_name}_marked_duplicates_metrics.txt"
    local log_file="$LOG_DIR/${base_name}_mark_duplicates.log"
    
    # Create temp directory with full permissions
    mkdir -p tmp && chmod 777 tmp
    
    # Set JAVA_TOOL_OPTIONS environment variable for this command
    JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$PWD/tmp" \
    picard MarkDuplicates \
        I="$bam_file" \
        O="$PICARD_OUTPUT_DIR/${base_name}_marked_duplicates.bam" \
        METRICS_FILE="$output_file" \
        CREATE_INDEX=true \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR="$PWD/tmp" \
        > "$log_file" 2>&1
        
    if [[ $? -ne 0 ]]; then
        echo "Picard MarkDuplicates failed for $bam_file. Check log for details: $log_file"
    fi
}

# Function to run RSeQC infer_experiment.py
run_rseqc() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local log_file="$LOG_DIR/${base_name}_rseqc.log"
    local output_dir="$RSEQ_OUTPUT_DIR/$base_name"
    mkdir -p "$output_dir"
    echo "Running RSeQC infer_experiment.py on $bam_file..."
    infer_experiment.py -i "$bam_file" -r "$BED_ANNOTATION_FILE" -s 200000  > "$log_file" 2>&1
    if [[ $? -ne 0 ]]; then
        echo "RSeQC failed for $bam_file. Check log for details: $log_file"
    fi
}

# Function to run Samtools flagstat
run_samtools_flagstat() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local output_file="$SAMTOOLS_OUTPUT_DIR/${base_name}_flagstat.txt"
    local log_file="$LOG_DIR/${base_name}_samtools_flagstat.log"
    echo "Running Samtools flagstat on $bam_file..."
    samtools flagstat "$bam_file" > "$output_file" 2> "$log_file"
    if [[ $? -ne 0 ]]; then
        echo "Samtools flagstat failed for $bam_file. Check log for details: $log_file"
    fi
}

# Export the functions and variables so they can be used by GNU parallel
export -f run_picard
export -f run_picard_mark_duplicates
export -f run_rseqc
export -f run_samtools_flagstat
export PICARD_OUTPUT_DIR RSEQ_OUTPUT_DIR SAMTOOLS_OUTPUT_DIR REFERENCE_GENOME REF_FLAT_FILE BED_ANNOTATION_FILE LOG_DIR

# Find all BAM files and run Picard, RSeQC, and Samtools flagstat in parallel
find "$BAM_DIR" -type f -name "*.bam" | parallel -j "$JOB_NUMBER" run_picard {}
find "$BAM_DIR" -type f -name "*.bam" | parallel -j 2 run_picard_mark_duplicates {}
find "$BAM_DIR" -type f -name "*.bam" | parallel -j "$JOB_NUMBER" run_rseqc {}
find "$BAM_DIR" -type f -name "*.bam" | parallel -j "$JOB_NUMBER" run_samtools_flagstat {}

# Prepare directories for STAR metrics
STAR_METRICS_DIRS=$(find "$BAM_DIR" -type d -name "*_STARpass1")

# Run MultiQC to aggregate the results
echo "Running MultiQC..."
multiqc "$PICARD_OUTPUT_DIR" "$RSEQ_OUTPUT_DIR" "$SAMTOOLS_OUTPUT_DIR" $STAR_METRICS_DIRS -o "$MULTIQC_DIR" -n "PostAlignmentQC"
if [[ $? -ne 0 ]]; then
    echo "MultiQC failed. Check output in: $MULTIQC_DIR"
fi

echo "Picard, RSeQC, Samtools flagstat, STAR metrics, and MultiQC processing complete."


