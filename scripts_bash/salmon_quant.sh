#!/bin/bash

# Constants
FASTQ_DIR="."                                # Directory where FASTQ files are located
SALMON_OUTPUT_DIR="./salmon"                 # Base directory for Salmon output
INDEX_DIR="/mouse_genome/ref_genome/transcripts/salmon_index"  # Prespecified Salmon index directory
EXPERIMENT_OUTPUT_DIR="$SALMON_OUTPUT_DIR/experiment_type"      # Directory for experiment strandedness output
LOG_DIR="$SALMON_OUTPUT_DIR/logs"            # Directory for logs

# Number of parallel jobs for GNU parallel
JOB_NUMBER=1

# Create directories if they don't exist
mkdir -p "$EXPERIMENT_OUTPUT_DIR" "$LOG_DIR"

# Path to the transcriptome FASTA file
TRANSCRIPTOME_FASTA="/mouse_genome/ref_genome/transcripts/GCF_000001635.27/rna.fna"

# Build Salmon index if it does not exist
if [ ! -d "$INDEX_DIR" ]; then
    echo "Building Salmon index in $INDEX_DIR..."
    salmon index -t "$TRANSCRIPTOME_FASTA" -i "$INDEX_DIR"
else
    echo "Salmon index already exists in $INDEX_DIR."
fi

# Function to infer strandedness using Salmon
run_salmon() {
    local fastq_file=$1
    local base_name=$(basename "$fastq_file")
    local sample_name="${base_name%.*}"  # Remove file extension for sample name
    local result_dir="$EXPERIMENT_OUTPUT_DIR/${sample_name}_salmon_output"
    local log_file="$LOG_DIR/${sample_name}_salmon.log"

    mkdir -p "$result_dir"

    echo "Running Salmon on $fastq_file..."
    salmon quant -i "$INDEX_DIR" -l A -r "$fastq_file" -o "$result_dir" > "$log_file" 2>&1
}

export -f run_salmon
export INDEX_DIR EXPERIMENT_OUTPUT_DIR LOG_DIR

# Find all FASTQ files and run Salmon sequentially with logging
find "$FASTQ_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" \) | while read -r fastq_file; do
    run_salmon "$fastq_file"
done

# Run MultiQC to gather reports
echo "Running MultiQC to gather Salmon reports..."
multiqc "$SALMON_OUTPUT_DIR"

echo "Salmon quantification and MultiQC report generation completed."

