#!/bin/bash

# Number of parallel FastQC jobs
JOB_NUMBER=6

# Define R1 patterns and their corresponding R2 replacements
declare -A R2_REPLACEMENTS=(
    ["_1."]="_2."
    ["_R1."]="_R2."
    ["_R1_001."]="_R2_001."
    ["_1_001."]="_2_001."
)

# Trimming parameters
LEADING=10
TRAILING=10
SLIDINGWINDOW="4:20"
MINLEN=36
ILLUMINACLIP_SETTINGS="2:30:10"
ADAPTERS="/usr/share/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# Directories setup
TRIMMED_DIR="/workspaces/Yasmine-retroT/0_Data/Preprocessing/trimmed/reads/add/trimmed"
FASTQC_OUTPUT_DIR="/workspaces/Yasmine-retroT/0_Data/Preprocessing/trimmed/reads/add/trimmed/postTrimQC/fastQC"
MULTIQC_OUTPUT_DIR="/workspaces/Yasmine-retroT/0_Data/Preprocessing/trimmed/reads/add/trimmed/postTrimQC/multiQC"
LOG_DIR="/workspaces/Yasmine-retroT/0_Data/Preprocessing/trimmed/reads/add/trimmed/logs"

mkdir -p "${TRIMMED_DIR}" "${FASTQC_OUTPUT_DIR}" "${MULTIQC_OUTPUT_DIR}" "${LOG_DIR}"

# Find all R1 files first
mapfile -t R1_FILES < <(find . -maxdepth 1 -type f \( -name "*_1.fastq.gz" -o -name "*_1.fq.gz" -o -name "*_R1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_R1_001.fastq.gz" -o -name "*_R1_001.fq.gz" \))

# Process only R1 files and their matching R2 pairs
for R1_FILE in "${R1_FILES[@]}"; do
    R1_FILE=$(basename "$R1_FILE")
    
    # Find matching R2 file
    for R1_PATTERN in "${!R2_REPLACEMENTS[@]}"; do
        if [[ $R1_FILE =~ $R1_PATTERN ]]; then
            R2_FILE="${R1_FILE/$R1_PATTERN/${R2_REPLACEMENTS[$R1_PATTERN]}}"
            
            if [[ -f "$R2_FILE" ]]; then
                # Extract sample name
                SAMPLE_NAME=$(basename "$R1_FILE" | sed -E 's/_(R)?1(_001)?\..*$//')
                
                TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_paired_trimmed.fastq.gz"
                TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_paired_trimmed.fastq.gz"
                LOG_FILE="${LOG_DIR}/${SAMPLE_NAME}_trimmomatic.log"
                
                echo "Processing ${SAMPLE_NAME} - R1: ${R1_FILE}, R2: ${R2_FILE}"
                
                trimmomatic PE -threads 8 \
                    "${R1_FILE}" "${R2_FILE}" \
                    "${TRIMMED_R1}" /dev/null \
                    "${TRIMMED_R2}" /dev/null \
                    ILLUMINACLIP:"${ADAPTERS}":"${ILLUMINACLIP_SETTINGS}" \
                    LEADING:"${LEADING}" TRAILING:"${TRAILING}" \
                    SLIDINGWINDOW:"${SLIDINGWINDOW}" MINLEN:"${MINLEN}" \
                    > "${LOG_FILE}" 2>&1
                
                echo "${TRIMMED_R1}" >> "${LOG_DIR}/fastqc_files_list.txt"
                echo "${TRIMMED_R2}" >> "${LOG_DIR}/fastqc_files_list.txt"
                break
            fi
        fi
    done
done

# Run FastQC function
run_fastqc() {
    local file=$1
    local base_name=$(basename "$file")
    local log_file="${LOG_DIR}/${base_name}_fastqc.log"
    echo "Running FastQC on $file..."
    fastqc "$file" --outdir "${FASTQC_OUTPUT_DIR}" --quiet > "$log_file" 2>&1
}

export -f run_fastqc
export FASTQC_OUTPUT_DIR LOG_DIR

# Run FastQC and MultiQC
if [[ -f "${LOG_DIR}/fastqc_files_list.txt" ]]; then
    parallel -j "$JOB_NUMBER" run_fastqc :::: "${LOG_DIR}/fastqc_files_list.txt"
    echo "Running MultiQC on FastQC results"
    multiqc "${FASTQC_OUTPUT_DIR}" -o "${MULTIQC_OUTPUT_DIR}"
else
    echo "No trimmed files found for FastQC processing"
fi
