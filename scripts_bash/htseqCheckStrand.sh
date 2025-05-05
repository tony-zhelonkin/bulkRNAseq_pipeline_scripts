#!/bin/bash

# Specify the BAM file, GTF file, and number of reads to sample
BAM_FILE=$1      # BAM file passed as the first argument
GTF_FILE=$2      # GTF annotation file passed as the second argument
N_READS_TO_SAMPLE=$3  # Number of reads to sample

# Ensure arguments are provided
if [ -z "$BAM_FILE" ] || [ -z "$GTF_FILE" ] || [ -z "$N_READS_TO_SAMPLE" ]; then
    echo "Usage: $0 <BAM file> <GTF file> <Number of reads to sample>"
    exit 1
fi

# Ensure htseq-count and samtools are installed
if ! command -v htseq-count &> /dev/null; then
    echo "htseq-count is not installed. Please install it and try again."
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "samtools is not installed. Please install it and try again."
    exit 1
fi

# Check BAM file header for SQ lines
if ! samtools view -H "$BAM_FILE" | grep -q '^@SQ'; then
    echo "Error: BAM file header is missing SQ lines. Please check the BAM file."
    exit 1
fi

# Create a temporary BAM file for sampled reads
SAMPLED_BAM=$(mktemp --suffix=.bam)

# Ensure the temporary BAM file gets deleted on exit
trap 'rm -f "$SAMPLED_BAM"' EXIT

# Sample the first N_READS_TO_SAMPLE reads from the BAM file
echo "Sampling the first $N_READS_TO_SAMPLE reads from $BAM_FILE..."
samtools view -h "$BAM_FILE" | head -n "$((N_READS_TO_SAMPLE + 1))" | samtools view -b -o "$SAMPLED_BAM"

# Check if the sampled BAM file is created
if [ ! -s "$SAMPLED_BAM" ]; then
    echo "Error: Sampled BAM file was not created or is empty."
    exit 1
fi

# Generate output file name based on the BAM file name
OUTPUT_FILE="${BAM_FILE%.bam}_strand_results.txt"

# Run htseq-count for each strandedness option and collect results
echo "Running htseq-count to infer strandedness..."
{
    echo "BAM File: $BAM_FILE"
    echo "Reference GTF File: $GTF_FILE"
    echo "Results for strandedness inference:"

    # Run for unstranded (-s no)
    echo "Testing unstranded (-s no)..."
    UNSTRANDED_COUNT=$(htseq-count -f bam -r pos -s no "$SAMPLED_BAM" "$GTF_FILE" 2>/dev/null | grep -v "__" | awk '{sum += $2} END {print sum}')
    echo "Unstranded (-s no) total reads: $UNSTRANDED_COUNT"

    # Run for stranded forward (-s yes)
    echo "Testing stranded forward (-s yes)..."
    STRANDED_YES_COUNT=$(htseq-count -f bam -r pos -s yes "$SAMPLED_BAM" "$GTF_FILE" 2>/dev/null | grep -v "__" | awk '{sum += $2} END {print sum}')
    echo "Stranded forward (-s yes) total reads: $STRANDED_YES_COUNT"

    # Run for stranded reverse (-s reverse)
    echo "Testing stranded reverse (-s reverse)..."
    STRANDED_REVERSE_COUNT=$(htseq-count -f bam -r pos -s reverse "$SAMPLED_BAM" "$GTF_FILE" 2>/dev/null | grep -v "__" | awk '{sum += $2} END {print sum}')
    echo "Stranded reverse (-s reverse) total reads: $STRANDED_REVERSE_COUNT"

    # Determine which strandedness option had the highest read count
    echo "Strandedness inference summary:"
    if [[ "$UNSTRANDED_COUNT" -ge "$STRANDED_YES_COUNT" && "$UNSTRANDED_COUNT" -ge "$STRANDED_REVERSE_COUNT" ]]; then
        echo "Inferred strandedness: Unstranded (best)"
    elif [[ "$STRANDED_YES_COUNT" -ge "$STRANDED_REVERSE_COUNT" ]]; then
        echo "Inferred strandedness: Stranded forward (best)"
    else
        echo "Inferred strandedness: Stranded reverse (best)"
    fi
} > "$OUTPUT_FILE"

echo "Strandedness inference completed. Results are stored in $OUTPUT_FILE."

