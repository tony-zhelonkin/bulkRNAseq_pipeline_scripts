#!/bin/bash

# Start timing the script
START_TIME=$(date +%s)

# Define paths
GENOME_DIR="/workspaces/Yasmine-retroT/0_Data/ref/ref_index100"
BATCH_1_DIR="$(pwd)"
OUTPUT_DIR="../STAR_alignment/bam"
LOG_DIR="../STAR_alignment/logs"

# Define computational resources
JOB_NUMBER="1" # number of GNU parallel jobs
GB_FREE="10G" # number of GB left free before any new process gets spawned

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Function to log and output a command
log_and_run() {
    local CMD="$1"
    local LOG_FILE="$2"
    echo "$(date '+%Y-%m-%d %H:%M:%S') Running command, logging to file" >> "$LOG_FILE"
    eval "$CMD" 2>&1 >> "$LOG_FILE"
}

# Function to process a pair of FASTQ files
process_fastq_pair() {
    local R1_FILE="$1"
    local R2_FILE="$2"
    local BASE_NAME="$3"
    local BATCH_DIR="$4"

    # Define read group (@RG) information
    local ID="${BASE_NAME}"
    local SM="Sample_${BASE_NAME}"
    local PL="ILLUMINA"
    local PU="${BATCH_DIR##*/}_${BASE_NAME}"

    # Create log file for this sample
    local LOG_FILE="${LOG_DIR}/${BASE_NAME}_log.txt"

    {
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting processing sample: $BASE_NAME"

        # STAR alignment with all original parameters
        log_and_run "STAR --runThreadN 10 \
                    --genomeDir \"$GENOME_DIR\" \
                    --readFilesIn \"$R1_FILE\" \"$R2_FILE\" \
                    --readFilesCommand zcat \
                    --sjdbOverhang 100 \
                    --limitBAMsortRAM 50000000000 \
                    --twopassMode Basic \
                    --outSAMattributes Standard \
                    --outSAMattrRGline ID:$ID LB:$ID SM:$SM PL:$PL PU:$PU \
                    --outFileNamePrefix \"$OUTPUT_DIR/${BASE_NAME}_\" \
                    --outSAMunmapped None \
                    --outSAMtype BAM Unsorted \
                    --alignEndsType EndToEnd \
                    --outMultimapperOrder Random \
                    --runRNGseed 777 \
                    --outFilterMultimapNmax 100 \
                    --winAnchorMultimapNmax 200 \
                    --outFilterMismatchNmax 10 \
                    --outSAMprimaryFlag AllBestScore \
                    --outFilterType BySJout \
                    --outFilterScoreMinOverLread 0.4 \
                    --outFilterMatchNminOverLread 0.4" "$LOG_FILE"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed processing sample: $BASE_NAME"
    } 2>&1 | tee -a "$LOG_FILE"
}

# Export function and variables for parallel execution
export -f process_fastq_pair log_and_run
export GENOME_DIR OUTPUT_DIR LOG_DIR

# Log script start
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting STAR alignment pipeline" | tee "$LOG_DIR/pipeline.log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working directory: $BATCH_1_DIR" | tee -a "$LOG_DIR/pipeline.log"

# Find and process all read pairs
(find "$BATCH_1_DIR" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | while read -r R1_FILE; do
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    if [[ -f "$R2_FILE" ]]; then
        BASE_NAME=$(basename "$R1_FILE" | sed -E 's/_R1_[^\.]*\.fastq\.gz//')
        BATCH_DIR=$(dirname "$R1_FILE")
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found pair: $BASE_NAME" >> "$LOG_DIR/pipeline.log"
        printf "%s\t%s\t%s\t%s\n" "$R1_FILE" "$R2_FILE" "$BASE_NAME" "$BATCH_DIR"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: No matching R2 file found for $R1_FILE" >> "$LOG_DIR/pipeline.log"
    fi
done) | tee "$LOG_DIR/pairs.log" | parallel --progress -j "$JOB_NUMBER" --memfree "$GB_FREE" --colsep '\t' process_fastq_pair {1} {2} {3} {4}


# End timing and calculate total duration
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

# Save total time to log file
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Total execution time: $TOTAL_TIME seconds" | tee -a "$LOG_DIR/pipeline.log" "$LOG_DIR/time.txt"
