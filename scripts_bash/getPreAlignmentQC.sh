#!/bin/bash

# getPreAlignmentQC.sh
#
# Description:
#   Performs pre-alignment quality control on FASTQ files using FastQC and MultiQC.
#   This script processes all FASTQ files in the input directory, runs FastQC on each,
#   and then aggregates the results using MultiQC.
#
# Usage:
#   ./getPreAlignmentQC.sh -i input_dir -o output_dir [-j num_jobs]
#
# Parameters:
#   -i  Input directory containing FASTQ files (.fq.gz, .fastq.gz, or .fasta)
#   -o  Output directory for FastQC and MultiQC reports
#   -j  Number of parallel FastQC processes to run (default: 6)
#
# Example:
#   ./getPreAlignmentQC.sh -i /path/to/fastq -o /path/to/qc -j 8
#
# Author: Tony Zhelonkin
# Date: 2023

# Default number of parallel jobs
JOB_NUMBER=6

# Usage message function
usage() {
    echo "Usage: $0 -i input_dir -o output_dir [-j num_jobs]"
    exit 1
}

# Parse command-line arguments
while getopts "i:o:j:" opt; do
    case ${opt} in
        i )
            FASTQ_DIR=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        j )
            JOB_NUMBER=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Ensure required arguments are provided
if [ -z "$FASTQ_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

# Check if input directory exists
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Input directory '$FASTQ_DIR' does not exist."
    exit 1
fi

# Define output subdirectories
FASTQC_OUTPUT_DIR="$OUTPUT_DIR/fastQC"
MULTIQC_OUTPUT_DIR="$OUTPUT_DIR/multiQC"
FASTQ_LOGS="$OUTPUT_DIR/logs"
# Create output directories if they don't exist
mkdir -p "$FASTQC_OUTPUT_DIR" "$MULTIQC_OUTPUT_DIR" "$FASTQ_LOGS"

# Find all .fq.gz, .fastq.gz, and .fasta files in the specified FASTQ_DIR
find "$FASTQ_DIR" -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fasta" \) > "$FASTQ_LOGS/fastq_files_list.txt"

# Check if any FASTQ files were found
if [ ! -s "$FASTQ_LOGS/fastq_files_list.txt" ]; then
    echo "Error: No FASTQ files found in '$FASTQ_DIR'"
    exit 1
fi

# Count the number of files to process
FILE_COUNT=$(wc -l < "$FASTQ_LOGS/fastq_files_list.txt")
echo "Found $FILE_COUNT FASTQ files to process"


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
echo "Starting FastQC analysis with $JOB_NUMBER parallel processes..."
cat "$FASTQ_LOGS/fastq_files_list.txt" | parallel -j "$JOB_NUMBER" run_fastqc {}

# Run MultiQC to aggregate the FastQC reports

# Run MultiQC to aggregate the FastQC reports
echo "Running MultiQC to aggregate FastQC reports..."
multiqc "$FASTQC_OUTPUT_DIR" -o "$MULTIQC_OUTPUT_DIR" > "$FASTQ_LOGS/multiqc.log" 2>&1

# Logging completion
echo "FastQC and MultiQC analysis completed." | tee -a "$FASTQ_LOGS/fastqc.log"
echo "FastQC reports: $FASTQC_OUTPUT_DIR"
echo "MultiQC report: $MULTIQC_OUTPUT_DIR/multiqc_report.html"
