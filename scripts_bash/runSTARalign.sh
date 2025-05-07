#!/bin/bash

# runSTARalign.sh
#
# Description:
#   Performs RNA-seq alignment using STAR aligner with support for different modes:
#   - Standard mode: Regular RNA-seq alignment with sorted BAM output
#   - TE mode: Optimized for transposable element analysis (for TEtranscripts)
#
# Usage:
#   ./runSTARalign.sh -i input_dir -o output_dir -g genome_dir [-m mode] [-j jobs] [-f mem_free]
#
# Parameters:
#   -i  Input directory containing FASTQ files (required)
#   -o  Output directory for alignment results (required)
#   -g  STAR genome index directory (required)
#   -m  Alignment mode: 'standard' or 'te' (default: standard)
#   -j  Number of parallel jobs to run (default: 1)
#   -f  Memory threshold for GNU parallel (default: 10G)
#
# Supported FASTQ file formats:
#   - File extensions: .fastq.gz, .fq.gz
#   - Paired-end naming conventions:
#     * _R1_/_R2_ (e.g., sample_R1_001.fastq.gz/sample_R2_001.fastq.gz)
#     * _R1/_R2 (e.g., sample_R1.fastq.gz/sample_R2.fastq.gz)
#     * _1/_2 (e.g., sample_1.fq.gz/sample_2.fq.gz)
#
# Examples:
#   # Standard RNA-seq alignment
#   ./runSTARalign.sh -i /path/to/fastq -o /path/to/output -g /path/to/genome
#
#   # Transposable element optimized alignment
#   ./runSTARalign.sh -i /path/to/fastq -o /path/to/output -g /path/to/genome -m te
#
# Author: Anton Zhelonkin
# Date: 2025

# Start timing the script
START_TIME=$(date +%s)

# Default parameters
MODE="standard"
JOB_NUMBER="1"
GB_FREE="10G"

# Usage message function
usage() {
    echo "Usage: $0 -i input_dir -o output_dir -g genome_dir [-m mode] [-j jobs] [-f mem_free]"
    echo "  -m  Alignment mode: 'standard' or 'te' (default: standard)"
    echo "  -f Keep free RAM, e.g.: '5G' (default: 10G)"
    exit 1
}

# Parse command-line arguments
while getopts "i:o:g:m:j:f:" opt; do
    case ${opt} in
        i )
            BATCH_DIR=$OPTARG
            ;;
        o )
            OUTPUT_BASE_DIR=$OPTARG
            ;;
        g )
            GENOME_DIR=$OPTARG
            ;;
        m )
            MODE=$OPTARG
            ;;
        j )
            JOB_NUMBER=$OPTARG
            ;;
        f )
            GB_FREE=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Ensure required arguments are provided
if [ -z "$BATCH_DIR" ] || [ -z "$OUTPUT_BASE_DIR" ] || [ -z "$GENOME_DIR" ]; then
    usage
fi

# Check if input directory exists
if [ ! -d "$BATCH_DIR" ]; then
    echo "Error: Input directory '$BATCH_DIR' does not exist."
    exit 1
fi

# Check if genome directory exists
if [ ! -d "$GENOME_DIR" ]; then
    echo "Error: Genome directory '$GENOME_DIR' does not exist."
    exit 1
fi

# Validate mode
if [[ "$MODE" != "standard" && "$MODE" != "te" ]]; then
    echo "Error: Mode must be 'standard' or 'te'"
    usage
fi

# Define output directories
OUTPUT_DIR="${OUTPUT_BASE_DIR}/bam"
LOG_DIR="${OUTPUT_BASE_DIR}/logs"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Function to log and output a command
log_and_run() {
    local CMD="$1"
    local LOG_FILE="$2"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running command and logging to $LOG_FILE"
    eval "$CMD" >> "$LOG_FILE" 2>> >(tee -a "$LOG_FILE" >&2)
}


# Function to process a pair of FASTQ files
process_fastq_pair() {
    local R1_FILE="$1"
    local R2_FILE="$2"
    local BASE_NAME="$3"
    local BATCH_DIR="$4"
    local MODE="$5"

    # Define read group (@RG) information
    local ID="${BASE_NAME}"
    local SM="Sample_${BASE_NAME}"
    local PL="ILLUMINA"
    local PU="${BATCH_DIR##*/}_${BASE_NAME}"

    # Create log file for this sample
    local LOG_FILE="${LOG_DIR}/${BASE_NAME}_log.txt"


    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: $BASE_NAME (Mode: $MODE)" >&2

    # Set STAR parameters based on mode
    if [[ "$MODE" == "standard" ]]; then
        # Standard RNA-seq alignment parameters
        log_and_run "STAR --runThreadN 10 \
                     --genomeDir \"$GENOME_DIR\" \
                     --readFilesIn \"$R1_FILE\" \"$R2_FILE\" \
                     --readFilesCommand zcat \
                     --sjdbOverhang 100 \
                     --limitBAMsortRAM 45000000000 \
                     --outFileNamePrefix \"$OUTPUT_DIR/${BASE_NAME}_\" \
                     --outSAMtype BAM SortedByCoordinate \
                     --outSAMattrRGline ID:$ID LB:$ID SM:$SM PL:$PL PU:$PU \
                     --twopassMode Basic" "$LOG_FILE"
    elif [[ "$MODE" == "te" ]]; then
        # Transposable element optimized parameters
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
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished processing sample: $BASE_NAME" >&2
}

# Export function and variables for parallel execution
export -f process_fastq_pair log_and_run
export GENOME_DIR OUTPUT_DIR LOG_DIR MODE

# Log script start
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting STAR alignment pipeline (Mode: $MODE)" | tee "$LOG_DIR/pipeline.log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working directory: $BATCH_DIR" | tee -a "$LOG_DIR/pipeline.log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUTPUT_BASE_DIR" | tee -a "$LOG_DIR/pipeline.log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Genome directory: $GENOME_DIR" | tee -a "$LOG_DIR/pipeline.log"

# Create a temporary file to store read pairs
PAIRS_FILE=$(mktemp)
PAIRS_LOG="$LOG_DIR/found_pairs.log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] List of all found read pairs:" > "$PAIRS_LOG"

# Find read pairs with various naming conventions
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Searching for paired-end FASTQ files..." | tee -a "$LOG_DIR/pipeline.log"

# ================================
# Improved file pair detection
# ================================

# Find all R1 files first
find "$BATCH_DIR" -maxdepth 1 -type f -name "*_1.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_R1.fastq.gz" -o -name "*_R1_*.fq.gz" -o -name "*_R1_*.fastq.gz" | sort > /tmp/r1_files.txt

# Process each R1 file and check for corresponding R2 file
while read -r R1_FILE; do
    FILENAME=$(basename "$R1_FILE")
    DIRNAME=$(dirname "$R1_FILE")
    
    # Construct the R2 filename based on pattern
    if [[ "$FILENAME" =~ _R1_.*\.(fastq|fq)\.gz$ ]]; then
        # _R1_ pattern
        R2_FILENAME="${FILENAME/_R1_/_R2_}"
        BASE_NAME="${FILENAME%%_R1_*}"
    elif [[ "$FILENAME" =~ _R1\.(fastq|fq)\.gz$ ]]; then
        # _R1. pattern
        R2_FILENAME="${FILENAME/_R1/_R2}"
        BASE_NAME="${FILENAME%%_R1.*}"
    elif [[ "$FILENAME" =~ _1\.(fastq|fq)\.gz$ ]]; then
        # _1. pattern - Specifically handle files ending with _1.fq.gz
        R2_FILENAME="${FILENAME%_1.fq.gz}_2.fq.gz"
        if [[ "$FILENAME" == *".fastq.gz" ]]; then
            R2_FILENAME="${FILENAME%_1.fastq.gz}_2.fastq.gz"
        fi
        BASE_NAME="${FILENAME%_1.*}"
    else
        # Skip files that don't match any expected pattern
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: Unrecognized file pattern: $R1_FILE" | tee -a "$LOG_DIR/pipeline.log"
        continue
    fi
    
    R2_FILE="$DIRNAME/$R2_FILENAME"
    
    # Check if R2 file exists
    if [[ -f "$R2_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found pair: $BASE_NAME ($FILENAME - $R2_FILENAME)" | tee -a "$LOG_DIR/pipeline.log" "$PAIRS_LOG"
        echo "$R1_FILE $R2_FILE $BASE_NAME $DIRNAME $MODE" >> "$PAIRS_FILE"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: No matching R2 file for: $R1_FILE (expected: $R2_FILE)" | tee -a "$LOG_DIR/pipeline.log"
    fi
done < /tmp/r1_files.txt

# Clean up temporary file
rm /tmp/r1_files.txt

# Check if any pairs were found
if [[ ! -s "$PAIRS_FILE" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No paired FASTQ files found in $BATCH_DIR" | tee -a "$LOG_DIR/pipeline.log"
    rm "$PAIRS_FILE"
    exit 1
fi

# Count the number of pairs
PAIR_COUNT=$(wc -l < "$PAIRS_FILE")
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found $PAIR_COUNT paired FASTQ files to process" | tee -a "$LOG_DIR/pipeline.log" "$PAIRS_LOG"

# Process all pairs using GNU parallel
cat "$PAIRS_FILE" | parallel --progress -j "$JOB_NUMBER" --memfree "$GB_FREE" --colsep ' ' process_fastq_pair {1} {2} {3} {4} {5}

# Clean up temporary file
rm "$PAIRS_FILE"

# End timing and calculate total duration
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

# Save total time to log file
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Total execution time: $TOTAL_TIME seconds" | tee -a "$LOG_DIR/pipeline.log" "$LOG_DIR/time.txt"
