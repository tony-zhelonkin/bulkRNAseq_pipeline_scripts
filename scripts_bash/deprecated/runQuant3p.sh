#!/bin/bash

# Number of parallel jobs (for GNU parallel)
JOB_NUMBER=8  # Set the desired number of parallel jobs

# Directories (use absolute paths)
BAM_DIR="."  # Directory containing BAM files
GTF_FILE="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf"  # Path to GTF file
OUTPUT_DIR="./quant3p/counts"  # Output directory for quant3p results
LOG_DIR="./quant3p/logs"  # Directory for logs

# Create output and log directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Function to run quant3p and save logs
run_quant3p() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .bam)
    local result_file="$OUTPUT_DIR/${base_name}_quant3p_results.txt"
    local log_file="$LOG_DIR/${base_name}_quant3p.log"

    echo "Running quant3p on $bam_file..."
    quant3p -n "$base_name" -g "$GTF_FILE" "$bam_file" > "$result_file" 2> "$log_file"
}

export -f run_quant3p
export GTF_FILE OUTPUT_DIR LOG_DIR

# Find all BAM files and run quant3p in parallel with logging
find "$BAM_DIR" -type f -name "*.bam" | parallel -j "$JOB_NUMBER" run_quant3p {}

# Optional: Logging completion
echo "quant3p analysis completed." >> "$LOG_DIR/quant3p_analysis.log"
