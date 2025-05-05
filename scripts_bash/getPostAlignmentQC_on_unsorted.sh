#!/bin/bash

# ------------------------------------------------------------------------------
# getPostAlignmentQC.sh
# ------------------------------------------------------------------------------

# This script:
# 1) Creates a coordinate-sorted copy of each .bam file (retaining the original).
# 2) Collects post-alignment QC metrics (Picard, RSeQC, Samtools, etc.).
# 3) Aggregates results via MultiQC.

# ------------------------------------------------------------------------------
# Environment and parallel settings
# ------------------------------------------------------------------------------
# Environment and parallel settings
export TMPDIR="/workspaces/Yasmine-retroT/.tmp"

# Number of parallel jobs
JOB_NUMBER=6

# Directories (use absolute paths)
BAM_DIR="."  # Directory containing BAM files
PICARD_OUTPUT_DIR="../post_align_QC/picard"
RSEQ_OUTPUT_DIR="../post_align_QC/rseqc"
SAMTOOLS_OUTPUT_DIR="../post_align_QC/samtools"
MULTIQC_DIR="../post_align_QC/multiQC"
LOG_DIR="../post_align_QC/logs"

# Reference files
REFERENCE_GENOME="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GCF_000001635.27_GRCm39_genomic.fna"
REF_FLAT_FILE="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GCF_000001635.27_GRCm39_genomic.reflat"
BED_ANNOTATION_FILE="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GRCm39_RfSeq.bed"
#RIBOSOMAL_INTERVALS="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/ribosomal.interval_list" # Using as variable leads to CollectRNAseqMetrics fail

# Create output directories
mkdir -p "$PICARD_OUTPUT_DIR" "$RSEQ_OUTPUT_DIR" "$SAMTOOLS_OUTPUT_DIR" "$MULTIQC_DIR" "$LOG_DIR"

# ------------------------------------------------------------------------------
# Function: sort_bam
# Create a coordinate-sorted copy of the BAM and preserve the original
# ------------------------------------------------------------------------------
sort_bam() {
    local bam_file="$1"
    local sorted_bam="${bam_file%.bam}.sorted.bam"

    if [[ -f "$sorted_bam" ]]; then
        echo "[sort_bam] $sorted_bam already exists, using existing file."
        return 0
    fi

    if samtools view -H "$bam_file" 2>/dev/null | grep -q "@HD.*SO:coordinate"; then
        echo "[sort_bam] $bam_file is already coordinate-sorted."
        cp "$bam_file" "$sorted_bam"
    else
        echo "[sort_bam] Sorting $bam_file => $sorted_bam ..."
        samtools sort -o "$sorted_bam" "$bam_file"
        echo "[sort_bam] Created $sorted_bam. (Original $bam_file retained.)"
    fi
}

# ------------------------------------------------------------------------------
# Function: run_picard
# Run Picard CollectRnaSeqMetrics on each sorted BAM
# ------------------------------------------------------------------------------
run_picard() {
    local bam_file="$1"
    local base_name=$(basename "${bam_file%.sorted.bam}" .bam)
    local sorted_bam="${bam_file}"
    local log_file="$LOG_DIR/${base_name}_picard.log"

    if [[ ! -f "$sorted_bam" ]]; then
        echo "[run_picard] Could not find BAM file: $sorted_bam. Skipping."
        return 1
    fi

    echo "Running Picard CollectRnaSeqMetrics on $sorted_bam..."
    picard CollectRnaSeqMetrics \
        INPUT="$sorted_bam" \
        OUTPUT="$PICARD_OUTPUT_DIR/${base_name}_rna_metrics.txt" \
        REF_FLAT="$REF_FLAT_FILE" \
        STRAND_SPECIFICITY=NONE \
        REFERENCE_SEQUENCE="$REFERENCE_GENOME" \
        RIBOSOMAL_INTERVALS="/workspaces/Yasmine-retroT/0_Data/ref/ref_genome/ribosomal.interval_list" \
        > "$log_file" 2>&1

    if [[ $? -ne 0 ]]; then
        echo "Picard CollectRnaSeqMetrics failed for $sorted_bam. See log: $log_file"
    fi
}


# ------------------------------------------------------------------------------
# Function: run_picard_mark_duplicates
# Run Picard MarkDuplicates on each sorted BAM
# ------------------------------------------------------------------------------
run_picard_mark_duplicates() {
    local bam_file="$1"
    local base_name=$(basename "${bam_file%.sorted.bam}" .bam)
    local sorted_bam="${bam_file}"
    local output_bam="$PICARD_OUTPUT_DIR/${base_name}_marked_duplicates.bam"
    local output_metrics="$PICARD_OUTPUT_DIR/${base_name}_marked_duplicates_metrics.txt"
    local log_file="$LOG_DIR/${base_name}_mark_duplicates.log"

    mkdir -p tmp && chmod 777 tmp

    if [[ ! -f "$sorted_bam" ]]; then
        echo "[run_picard_mark_duplicates] Could not find BAM file: $sorted_bam. Skipping."
        return 1
    fi

    echo "Running Picard MarkDuplicates on $sorted_bam..."
    JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$PWD/tmp" \
    picard MarkDuplicates \
        -I "$sorted_bam" \
        -O "$output_bam" \
        -METRICS_FILE "$output_metrics" \
        -CREATE_INDEX true \
        -REMOVE_DUPLICATES false \
        -VALIDATION_STRINGENCY LENIENT \
        -TMP_DIR "$PWD/tmp" \
        > "$log_file" 2>&1

    if [[ $? -ne 0 ]]; then
        echo "Picard MarkDuplicates failed for $sorted_bam. See log: $log_file"
    fi
}

# ------------------------------------------------------------------------------
# Function: run_rseqc
# Run RSeQC infer_experiment.py on each sorted BAM
# ------------------------------------------------------------------------------
run_rseqc() {
    local bam_file="$1"
    local base_name=$(basename "${bam_file%.sorted.bam}" .bam)
    local sorted_bam="${bam_file}"
    local log_file="$LOG_DIR/${base_name}_rseqc.log"
    local output_dir="$RSEQ_OUTPUT_DIR/$base_name"

    mkdir -p "$output_dir"

    if [[ ! -f "$sorted_bam" ]]; then
        echo "[run_rseqc] Could not find BAM file: $sorted_bam. Skipping."
        return 1
    fi

    echo "Running RSeQC infer_experiment.py on $sorted_bam..."
    infer_experiment.py \
        -i "$sorted_bam" \
        -r "$BED_ANNOTATION_FILE" \
        -s 200000 \
        > "$output_dir/${base_name}_infer_experiment.txt" 2> "$log_file"

    if [[ $? -ne 0 ]]; then
        echo "RSeQC failed for $sorted_bam. Check log: $log_file"
    fi
}

# ------------------------------------------------------------------------------
# Function: run_samtools_flagstat
# Run Samtools flagstat on each sorted BAM
# ------------------------------------------------------------------------------
run_samtools_flagstat() {
    local bam_file="$1"
    local base_name=$(basename "${bam_file%.sorted.bam}" .bam)
    local sorted_bam="${bam_file}"
    local output_file="$SAMTOOLS_OUTPUT_DIR/${base_name}_flagstat.txt"
    local log_file="$LOG_DIR/${base_name}_samtools_flagstat.log"

    if [[ ! -f "$sorted_bam" ]]; then
        echo "[run_samtools_flagstat] Could not find BAM file: $sorted_bam. Skipping."
        return 1
    fi

    echo "Running Samtools flagstat on $sorted_bam..."
    samtools flagstat "$sorted_bam" > "$output_file" 2> "$log_file"

    if [[ $? -ne 0 ]]; then
        echo "Samtools flagstat failed for $sorted_bam. Check log: $log_file"
    fi
}

# Export functions and variables for parallel execution
export -f sort_bam
export -f run_picard
export -f run_picard_mark_duplicates
export -f run_rseqc
export -f run_samtools_flagstat
export PICARD_OUTPUT_DIR RSEQ_OUTPUT_DIR SAMTOOLS_OUTPUT_DIR REFERENCE_GENOME REF_FLAT_FILE BED_ANNOTATION_FILE LOG_DIR

# Process BAM files
find "$BAM_DIR" -type f -name "*.bam" -not -name "*.sorted.bam" | parallel -j "$JOB_NUMBER" sort_bam {}

# Run QC on sorted BAMs
find "$BAM_DIR" -type f -name "*.sorted.bam" | parallel -j "$JOB_NUMBER" run_picard
find "$BAM_DIR" -type f -name "*.sorted.bam" | parallel -j 2 run_picard_mark_duplicates
find "$BAM_DIR" -type f -name "*.sorted.bam" | parallel -j "$JOB_NUMBER" run_rseqc
find "$BAM_DIR" -type f -name "*.sorted.bam" | parallel -j "$JOB_NUMBER" run_samtools_flagstat

# Get STAR metrics directories
STAR_METRICS_DIRS=$(find "$BAM_DIR" -type d -name "*_STARpass1")

# Run MultiQC
echo "Running MultiQC..."
multiqc "$PICARD_OUTPUT_DIR" "$RSEQ_OUTPUT_DIR" "$SAMTOOLS_OUTPUT_DIR" $STAR_METRICS_DIRS -o "$MULTIQC_DIR" -n "PostAlignmentQC" --force

if [[ $? -ne 0 ]]; then
    echo "MultiQC failed. Check output in: $MULTIQC_DIR"
fi

echo "Post-alignment QC processing complete."
