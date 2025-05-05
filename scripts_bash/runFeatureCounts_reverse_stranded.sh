#!/bin/bash

# Directory containing BAM files
BAM_DIR="."  # Only the current directory

# Annotation file
ANNOTATION_FILE="/mouse_genome/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf"

# Output directories
RAW_OUTPUT_DIR="./featureCounts/raw"
MATRIX_OUTPUT_DIR="./featureCounts/count_matrices"
LOG_DIR="./featureCounts/logs"

mkdir -p "$RAW_OUTPUT_DIR"
mkdir -p "$MATRIX_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Number of threads for featureCounts
THREADS=8  # Adjust this based on your system's capabilities

# Output files
COUNTS_FILE="$RAW_OUTPUT_DIR/counts_matrix.txt"
SORTED_COUNTS_FILE="$MATRIX_OUTPUT_DIR/sorted_counts_matrix.txt"
LOG_FILE="$LOG_DIR/featureCounts.log"

# Find all BAM files in the current directory (non-recursive) and sort them lexicographically
mapfile -t SORTED_BAM_FILES < <(ls "$BAM_DIR"/*.bam | sort)

echo "Processing BAM files: ${SORTED_BAM_FILES[*]}"

# Add debugging information
echo "Number of BAM files: ${#SORTED_BAM_FILES[@]}" >> "$LOG_FILE"
for file in "${SORTED_BAM_FILES[@]}"; do
    echo "File: $file" >> "$LOG_FILE"
done

# Check if any BAM files were found
if [ ${#SORTED_BAM_FILES[@]} -eq 0 ]; then
    echo "No BAM files found in $BAM_DIR" >> "$LOG_FILE"
    exit 1
fi

# Run featureCounts on all sorted BAM files in one command
featureCounts \
    -a "$ANNOTATION_FILE" \
    -o "$COUNTS_FILE" \
    -p \
    --countReadPairs \
    -B \
    -C \
    -s 2 \
    -t exon \
    -T "$THREADS" \
    "${SORTED_BAM_FILES[@]}" \
    >> "$LOG_FILE" 2>&1

if [ $? -ne 0 ]; then
    echo "featureCounts failed. Check the log file at $LOG_FILE for details."
    exit 1
fi

# Reorder columns in the output file
awk 'BEGIN {FS=OFS="\t"}
    NR==1 {next}  # Skip the first line (featureCounts summary)
    NR==2 {
        printf "Geneid";
        for (i=7; i<=NF; i++) {
            split($i, a, "/");  # Split the path
            split(a[length(a)], b, ".");  # Split the filename
            printf "\t%s", b[1];  # Print only the sample name
        }
        printf "\n";
        next;
    }
    {
        printf "%s", $1;
        for (i=7; i<=NF; i++) printf "\t%s", $i;
        printf "\n";
    }' "$COUNTS_FILE" > "$SORTED_COUNTS_FILE"

# Optional: Logging completion
echo "featureCounts analysis and count matrix creation completed." >> "$LOG_FILE"
echo "Completed processing all BAM files. Results saved to $SORTED_COUNTS_FILE"
