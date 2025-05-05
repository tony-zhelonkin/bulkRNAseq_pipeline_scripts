#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -i input.gtf -o output.refFlat"
    exit 1
}

# Parse command-line arguments
while getopts "i:o:" opt; do
    case ${opt} in
        i )
            GTF_FILE=$OPTARG
            ;;
        o )
            REFLAT_OUTPUT=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Ensure both input and output arguments are provided
if [ -z "$GTF_FILE" ] || [ -z "$REFLAT_OUTPUT" ]; then
    usage
fi

# Extracting fields from the GTF file and formatting them into RefFlat format
awk '
BEGIN {OFS="\t"}
{
    if ($3 == "exon") {
        # Extract relevant fields from GTF
        chrom = $1;  # Chromosome
        start = $4;  # Exon start position
        end = $5;    # Exon end position
        strand = $7; # Strand (+ or -)

        # Parse gene name and transcript name (adjust based on GTF format)
        geneName = ""; transcriptName = "";
        for (i = 9; i <= NF; i++) {
            if ($i ~ /gene_id/) {
                geneName = $(i+1); gsub(/"|;/, "", geneName);
            }
            if ($i ~ /transcript_id/ && $(i+1) != "\"\"") {
                transcriptName = $(i+1); gsub(/"|;/, "", transcriptName);
            }
        }

        # Skip lines where transcript name is empty
        if (transcriptName == "") {
            next;
        }

        # Store information for each transcript
        key = geneName "\t" transcriptName "\t" chrom "\t" strand;

        # Initialize or update transcript details
        if (!(key in transcripts)) {
            transcripts[key] = geneName "\t" transcriptName "\t" chrom "\t" strand "\t" start "\t" end "\t" start "\t" end "\t1\t" start "," "\t" end ",";
        } else {
            split(transcripts[key], fields, "\t");
            txStart = (start < fields[5]) ? start : fields[5];
            txEnd = (end > fields[6]) ? end : fields[6];
            exonCount = fields[9] + 1;
            exonStarts = fields[10] "" start ",";
            exonEnds = fields[11] "" end ",";
            transcripts[key] = geneName "\t" transcriptName "\t" chrom "\t" strand "\t" txStart "\t" txEnd "\t" fields[7] "\t" fields[8] "\t" exonCount "\t" exonStarts "\t" exonEnds;
        }
    }
}

# Output in RefFlat format
END {
    for (key in transcripts) {
        print transcripts[key];
    }
}' "$GTF_FILE" > "$REFLAT_OUTPUT"

# Confirm output file creation
echo "RefFlat file created: $REFLAT_OUTPUT"

