#!/usr/bin/env bash
#
#  runPostAlignmentQC.sh  – UNIFIED POST-ALIGNMENT QC PIPELINE
#  -----------------------------------------------------------
#  Performs (optionally) duplicate removal, a suite of Picard / samtools
#  metrics, and the full RSeQC battery
##  Description:
#    Consolidated post-alignment QC for RNA-seq BAM files.
#    ├─ Optionally coordinate-sorts BAMs that aren’t already
#    ├─ Picard MarkDuplicates    (duplicate marking **or** removal)
#    ├─ Picard CollectRnaSeqMetrics  (customisable strand)
#    ├─ Samtools flagstat
#    ├─ Core RSeQC modules
#    │     • infer_experiment.py
#    │     • read_duplication.py
#    │     • read_distribution.py
#    │     • geneBody_coverage.py
#    |     • junction_annotation.py
#    |     • inner_distance.py
#    |     • junction_saturation.py
#    └─ Aggregates everything (and STAR logs if present) with MultiQC
#
#
#  REQUIRED ARGUMENTS
#    -i  INPUT_DIR        directory containing BAM files
#    -o  OUTPUT_DIR       where QC results will be written
#    -r  REF_GENOME       reference genome FASTA
#    -ri RIBO_INTERVALS   Picard ribosomal.interval_list
#    -rf REF_FLAT         refFlat gene annotation (for Picard)
#    -b  BED_ANNOT        BED gene model (for all RSeQC modules)
#
#  KEY OPTIONS
#       --remove-duplicates   physically drop duplicates (default: mark only)
#       --strand STRAND       NONE | FIRST_READ_TRANSCRIPTION_STRAND | SECOND_READ_TRANSCRIPTION_STRAND
#    -j  JOBS                 parallel samples to process        [4]
#    -t  THREADS              threads per tool invocation        [8]
#    -T  TMPDIR               temp directory (falls back to /tmp)
#
#  EXAMPLE
#    ./runPostAlignmentQC.sh -i ./bam -o ./QC \
#         -r ref.fa -ri rRNA.interval_list -rf refFlat.txt -b genes.bed \
#         --remove-duplicates --strand FIRST_READ_TRANSCRIPTION_STRAND -j 6
#
#  Author:  Anton Zhelonkin
#  Date:    2025
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

SHORT="i:o:r:ri:rf:b:j:t:T:h"
LONG="remove-duplicates,strand:,help"
PARSED=$(getopt -o "$SHORT" --long "$LONG" -n "$0" -- "$@") || usage
eval set -- "$PARSED"

while true; do
  case "$1" in
    -i) INPUT_DIR=$2; shift 2 ;;
    -o) OUTPUT_BASE_DIR=$2; shift 2 ;;
    -r) REF_GENOME=$2; shift 2 ;;
    -ri) RIBO_INTERVALS=$2; shift 2 ;;
    -rf) REF_FLAT=$2; shift 2 ;;
    -b) BED_ANNOT=$2; shift 2 ;;
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

  # (1) coordinate-sorted BAM
  local sorted="${bam%.bam}.sorted.bam"
  if samtools view -H "$bam" | grep -q "@HD.*SO:coordinate"; then
      cp "$bam" "$sorted"
  else
      run_cmd "samtools sort -@ ${THREADS} -o ${sorted} ${bam}" "$wlog"
  fi
  run_cmd "samtools index ${sorted}" "$wlog"

  # (2) Picard MarkDuplicates
  local mkd="${PICARD_DIR}/${base}_marked.bam"
  local dupM="${PICARD_DIR}/${base}_dup_metrics.txt"
  run_cmd "picard MarkDuplicates I=${sorted} O=${mkd} M=${dupM} \
           REMOVE_DUPLICATES=${REMOVE_DUPES} CREATE_INDEX=true TMP_DIR=${TMPDIR}" "$wlog"
  local qc_bam="$mkd"

  # (3) Picard CollectRnaSeqMetrics
  run_cmd "picard CollectRnaSeqMetrics I=${qc_bam} \
           O=${PICARD_DIR}/${base}_rna_metrics.txt \
           R=${REF_GENOME} REF_FLAT=${REF_FLAT} \
           STRAND_SPECIFICITY=${STRAND} RIBOSOMAL_INTERVALS=${RIBO_INTERVALS}" "$wlog"

  # (4) samtools flagstat
  run_cmd "samtools flagstat ${qc_bam} > ${SAMTOOLS_DIR}/${base}_flagstat.txt" "$wlog"

  # (5) RSeQC – basic set
  run_cmd "infer_experiment.py   -i ${qc_bam} -r ${BED_ANNOT} \
           > ${RSEQ_DIR}/${base}_infer_experiment.txt" "$wlog"
  run_cmd "read_duplication.py   -i ${qc_bam} -o ${RSEQ_DIR}/${base}_dup" "$wlog"
  run_cmd "read_distribution.py  -i ${qc_bam} -r ${BED_ANNOT} \
           > ${RSEQ_DIR}/${base}_read_distribution.txt" "$wlog"
  run_cmd "geneBody_coverage.py  -i ${qc_bam} -r ${BED_ANNOT} \
           -o ${RSEQ_DIR}/${base}_geneBody" "$wlog"

  # (6) NEW  ► inner_distance / junction_annotation / junction_saturation
  run_cmd "inner_distance.py      -i ${qc_bam} -o ${RSEQ_DIR}/${base}_innerDist \
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
