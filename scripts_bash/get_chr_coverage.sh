#! /bin/bash

# dir to indexed BAM files (bai files have to be in the same directory)
BAM_DIR="/data/star_run/sep_run/bam"

# output file to store coverage stats
OUTPUT_FILE="chr_coverage.txt"

# Initialize the output file with a header
echo -e "File\tChromosome\tCoverage"
> "$OUTPUT_FILE"

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
        #exctract filename without path for output
        FILE_NAME=$(basename "$BAM_FILE")

        # Run the samtools idxstats command
	CHRX_COVERAGE=$(samtools idxstats "$BAM_FILE" | awk '$1=="NC_000086.8" {print $3}')
	CHRY_COVERAGE=$(samtools idxstats "$BAM_FILE" | awk '$1=="NC_000087.8" {print $3}')

	echo -e "${FILE_NAME}\t${CHRX_COVERAGE}\t${CHRY_COVERAGE}" >> "$OUTPUT_FILE"

done

echo "Coverage data for specified chromosomes has been written to $OUTPUT_FILE."
