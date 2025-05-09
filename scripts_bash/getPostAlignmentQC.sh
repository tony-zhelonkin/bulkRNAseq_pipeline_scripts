#!/usr/bin/env bash
#
#  runPostAlignmentQC.sh  – UNIFIED POST-ALIGNMENT QC PIPELINE
#  -----------------------------------------------------------
#  Performs (optionally) duplicate removal, a suite of Picard / samtools
#  metrics, and the full RSeQC battery.
#
#  WORKFLOW
#    1. (If needed) coordinate-sort BAMs
#    2. Picard MarkDuplicates   ── mark or remove dupes
#    3. Picard CollectRnaSeqMetrics  (strand-aware if requested)
#    4. samtools flagstat
#    5. RSeQC modules
#         • infer_experiment.py
#         • read_duplication.py
#         • read_distribution.py
#         • geneBody_coverage.py
#         • inner_distance.py
#         • junction_annotation.py
#         • junction_saturation.py
#    6. Aggregate everything (plus any *_STARpass1 folders) with MultiQC
#
#  REQUIRED ARGUMENTS
#    -i  DIR   directory containing BAM files
#    -o  DIR   output base directory for all QC results
#    -r  FILE  reference genome FASTA
#    -e  FILE  ribosomal.intervals file for Picard
#    -f  FILE  refFlat annotation for Picard
#    -b  FILE  BED gene model for all RSeQC tools
#
#  KEY OPTIONS
#       --remove-duplicates          physically drop duplicates (default: mark only)
#       --strand STRAND              NONE | FIRST_READ_TRANSCRIPTION_STRAND | SECOND_READ_TRANSCRIPTION_STRAND
#    -j  INT   parallel samples to process      [4]
#    -t  INT   threads per tool invocation      [8]
#    -T  DIR   temporary directory (default: /tmp or \$TMPDIR)
#
#  EXAMPLES
#
#  1) **Default** – mark duplicates, strand NONE
#     ./runPostAlignmentQC.sh \\
#          -i ./bam \\
#          -o ./QC \\
#          -r genome.fa \\
#          -e rRNA.interval_list \\
#          -f refFlat.txt \\
#          -b genes.bed
#
#  2) **Remove duplicates** and use first-strand library setting
#     ./runPostAlignmentQC.sh \\
#          -i ./bam -o ./QC -r genome.fa -e rRNA.interval_list -f refFlat.txt -b genes.bed \\
#          --remove-duplicates \\
#          --strand FIRST_READ_TRANSCRIPTION_STRAND
#
#  3) **Docker** with explicit TMPDIR (duplicates kept, strand NONE)
#     docker run --rm -it \\
#       -v \$PWD/bam:/in:ro \\
#       -v \$PWD/QC:/out \\
#       -v \$PWD/QC/tmp:/out/tmp -e TMPDIR=/out/tmp \\
#       -v \$PWD/ref:/genome:ro \\
#       myimage:latest \\
#       bash -c '/pipeline/runPostAlignmentQC.sh \\
#                   -i /in -o /out \\
#                   -r /genome/genome.fa \\
#                   -e /genome/rRNA.interval_list \\
#                   -f /genome/refFlat.txt \\
#                   -b /genome/genes.bed \\
#                   -j 4 -t 8 --strand NONE'
#
#  Author: Anton Zhelonkin
#  Last update: 2025-05-09
# ---------------------------------------------------------------------------

set -euo pipefail

# ----------------------------- defaults -------------------------------------
INPUT_DIR="."
OUTPUT_BASE_DIR="./post_align_QC"
REF_GENOME=""; RIBO_INTERVALS=""; REF_FLAT=""; BED_ANNOT=""
TMPDIR="${TMPDIR:-/tmp}"
JOBS=4
THREADS=8
REMOVE_DUPES="false"
STRAND="NONE"
# ---------------------------------------------------------------------------

usage() { grep '^#' "$0" | cut -c4-; exit 1; }

SHORT="i:o:r:e:f:b:j:t:T:h"
LONG="remove-duplicates,strand:,ribosomal-intervals:,ref-flat:,bed-annot:,help"
PARSED=$(getopt -o "$SHORT" --long "$LONG" -n "$0" -- "$@") || usage
eval set -- "$PARSED"

while true; do
  case "$1" in
    -i) INPUT_DIR=$2; shift 2 ;;
    -o) OUTPUT_BASE_DIR=$2; shift 2 ;;
    -r) REF_GENOME=$2; shift 2 ;;
    -e|--ribosomal-intervals) RIBO_INTERVALS=$2; shift 2 ;;
    -f|--ref-flat)            REF_FLAT=$2; shift 2 ;;
    -b|--bed-annot)           BED_ANNOT=$2; shift 2 ;;
    -j) JOBS=$2; shift 2 ;;
    -t) THREADS=$2; shift 2 ;;
    -T) TMPDIR=$2; shift 2 ;;
    --remove-duplicates) REMOVE_DUPES="true"; shift ;;
    --strand) STRAND=$2; shift 2 ;;
    -h|--help) usage ;;
    --) shift; break ;;
    *) usage ;;
  esac
done

# ----------------------------- validation -----------------------------------
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: INPUT_DIR not found"; exit 1; }
for f in "$REF_GENOME" "$RIBO_INTERVALS" "$REF_FLAT" "$BED_ANNOT"; do
  [[ -f "$f" ]] || { echo "ERROR: file not found → $f"; exit 1; }
done
[[ "$STRAND" =~ ^(NONE|FIRST_READ_TRANSCRIPTION_STRAND|SECOND_READ_TRANSCRIPTION_STRAND)$ ]] \
  || { echo "ERROR: invalid --strand value"; exit 1; }

mkdir -p "$OUTPUT_BASE_DIR"/{picard,rseqc,samtools,multiQC,logs,tmp}
export TMPDIR

PICARD_DIR="$OUTPUT_BASE_DIR/picard"
RSEQ_DIR="$OUTPUT_BASE_DIR/rseqc"
SAMTOOLS_DIR="$OUTPUT_BASE_DIR/samtools"
LOG_DIR="$OUTPUT_BASE_DIR/logs"

# --------------------------- helper -----------------------------------------
run_cmd() { echo "[$(date '+%F %T')] $1" >> "$2"; eval "$1" >> "$2" 2>&1; }

# --------------------------- per-sample workhorse ---------------------------
process_bam() {
  local bam="$1"
  local base=$(basename "${bam%.bam}")
  local wlog="$LOG_DIR/${base}.log"

  echo -e "\n===== ${base} =====" >> "$wlog"

  # 1) coordinate-sorted BAM
  local sorted="${bam%.bam}.sorted.bam"
  if samtools view -H "$bam" | grep -q "@HD.*SO:coordinate"; then
      cp "$bam" "$sorted"
  else
      run_cmd "samtools sort -@ ${THREADS} -o ${sorted} ${bam}" "$wlog"
  fi
  run_cmd "samtools index ${sorted}" "$wlog"

  # 2) Picard MarkDuplicates
  local mkd="${PICARD_DIR}/${base}_marked.bam"
  local dupM="${PICARD_DIR}/${base}_dup_metrics.txt"
  run_cmd "picard MarkDuplicates I=${sorted} O=${mkd} M=${dupM} \
           REMOVE_DUPLICATES=${REMOVE_DUPES} CREATE_INDEX=true TMP_DIR=${TMPDIR}" "$wlog"
  local qc_bam="$mkd"

  # 3) Picard CollectRnaSeqMetrics
  run_cmd "picard CollectRnaSeqMetrics I=${qc_bam} \
           O=${PICARD_DIR}/${base}_rna_metrics.txt \
           R=${REF_GENOME} REF_FLAT=${REF_FLAT} \
           STRAND_SPECIFICITY=${STRAND} RIBOSOMAL_INTERVALS=${RIBO_INTERVALS}" "$wlog"

  # 4) samtools flagstat
  run_cmd "samtools flagstat ${qc_bam} > ${SAMTOOLS_DIR}/${base}_flagstat.txt" "$wlog"

  # 5) RSeQC
  run_cmd "infer_experiment.py   -i ${qc_bam} -r ${BED_ANNOT} \
           > ${RSEQ_DIR}/${base}_infer_experiment.txt" "$wlog"
  run_cmd "read_duplication.py   -i ${qc_bam} -o ${RSEQ_DIR}/${base}_dup" "$wlog"
  run_cmd "read_distribution.py  -i ${qc_bam} -r ${BED_ANNOT} \
           > ${RSEQ_DIR}/${base}_read_distribution.txt" "$wlog"
  run_cmd "geneBody_coverage.py  -i ${qc_bam} -r ${BED_ANNOT} \
           -o ${RSEQ_DIR}/${base}_geneBody" "$wlog"
  run_cmd "inner_distance.py     -i ${qc_bam} -o ${RSEQ_DIR}/${base}_innerDist \
           -r ${BED_ANNOT}" "$wlog"
  run_cmd "junction_annotation.py -i ${qc_bam} -o ${RSEQ_DIR}/${base}_junctionAnnot \
           -r ${BED_ANNOT}" "$wlog"
  run_cmd "junction_saturation.py -i ${qc_bam} -o ${RSEQ_DIR}/${base}_junctionSat \
           -r ${BED_ANNOT}" "$wlog"

  echo "[$(date '+%F %T')] DONE ${base}" >> "$wlog"
}
export -f process_bam run_cmd
export PICARD_DIR RSEQ_DIR SAMTOOLS_DIR LOG_DIR REF_GENOME REF_FLAT \
       RIBO_INTERVALS BED_ANNOT TMPDIR THREADS REMOVE_DUPES STRAND

# --------------------------- run in parallel --------------------------------
find "$INPUT_DIR" -maxdepth 1 -type f -name "*.bam" | \
  parallel -j "$JOBS" process_bam {}

# --------------------------- MultiQC ----------------------------------------
STAR_DIRS=$(find "$INPUT_DIR" -maxdepth 2 -type d -name "*_STARpass1" 2>/dev/null | tr '\n' ' ')
multiqc $PICARD_DIR $RSEQ_DIR $SAMTOOLS_DIR $STAR_DIRS \
        -o "$OUTPUT_BASE_DIR/multiQC" -n "PostAlignmentQC" --force

echo -e "\n✓ All QC complete → ${OUTPUT_BASE_DIR}/multiQC/PostAlignmentQC.html"
