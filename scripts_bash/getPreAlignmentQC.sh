#!/bin/bash

# Number of parallel jobs (for GNU parallel)
JOB_NUMBER=6  # Number of parallel FastQC processes to run

# Directories (use absolute paths)
FASTQ_DIR="."  # Directory where FASTQ files are located
FASTQC_OUTPUT_DIR="/workspaces/Yasmine-retroT/0_Data/Preprocessing/pre_align_QC/fastQC"  # Output directory for FastQC reports
MULTIQC_OUTPUT_DIR="/workspaces/Yasmine-retroT/0_Data/Preprocessing/pre_align_QC/multiQC"  # Output directory for MultiQC report
FASTQ_LOGS="/workspaces/Yasmine-retroT/0_Data/Preprocessing/pre_align_QC/logs"  # Directory for logs

# Create output directories if they don't exist
mkdir -p "$FASTQC_OUTPUT_DIR" "$MULTIQC_OUTPUT_DIR" "$FASTQ_LOGS"

# Find all .fq.gz, .fastq.gz, and .fasta files in the specified FASTQ_DIR
find "$FASTQ_DIR" -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fasta" \) > "$FASTQ_LOGS/fastq_files_list.txt"

# Function to run FastQC and log output
run_fastqc() {
    local file=$1
    local base_name=$(basename "$file")
    local log_file="$FASTQ_LOGS/${base_name}_fastqc.log"

    echo "Running FastQC on $file..."
    fastqc "$file" --outdir "$FASTQC_OUTPUT_DIR" --quiet > "$log_file" 2>&1
}

export -f run_fastqc
export FASTQC_OUTPUT_DIR FASTQ_LOGS

# Run FastQC on each file using GNU parallel with logging
cat "$FASTQ_LOGS/fastq_files_list.txt" | parallel -j "$JOB_NUMBER" run_fastqc {}

# Run MultiQC to aggregate the FastQC reports
multiqc "$FASTQC_OUTPUT_DIR" -o "$MULTIQC_OUTPUT_DIR"

# Optional: Logging completion
echo "FastQC and MultiQC analysis completed." >> "$FASTQ_LOGS/fastqc.log"#
