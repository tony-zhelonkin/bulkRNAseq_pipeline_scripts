#! /bin/bash

# dir with BAM files 
BAM_DIR="/data/star_run/sep_run/bam"

# Loop through BAM files in the folder 

for BAM_FILE in "$BAM_DIR"/*.bam; do
	echo "Indexing $BAM_FILE..."
	samtools index "$BAM_FILE"
	echo "$BAM_FILE indexed successfully!"
done
