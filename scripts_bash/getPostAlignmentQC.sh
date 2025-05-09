#!/usr/bin/env bash
#
#  runPostAlignmentQC.sh  – unified post-alignment QC pipeline
#  -----------------------------------------------------------
#  INPUT  : directory of BAMs (read-only OK)
#  OUTPUT : one MultiQC report + per-tool logs/metrics
#
#  WORKFLOW
#    1. Ensure a writable, coordinate-sorted BAM
#         ─ already sorted?  • default  : use in-place (read-only)
#                            • --copy-sorted : symlink to $TMPDIR
#         ─ unsorted?        • samtools sort  → $TMPDIR
#    2. Picard MarkDuplicates   (mark OR remove dupes)
#    3. Picard CollectRnaSeqMetrics   (strand-aware)
#    4. samtools flagstat
#    5. RSeQC modules
#         • infer_experiment.py
#         • read_duplication.py
#         • read_distribution.py
#         • geneBody_coverage.py
#         • inner_distance.py
#         • junction_annotation.py
#         • junction_saturation.py
#    6. Aggregate everything (and *_STARpass1) with MultiQC
#
#  REQUIRED
#    -i DIR   BAM directory
#    -o DIR   output directory
#    -r FILE  genome.fa
#    -e FILE  ribosomal.interval_list      (Picard)
#    -f FILE  refFlat.txt                  (Picard)
#    -b FILE  genes.bed (BED12)            (RSeQC)
#
#  COMMON OPTIONS
#       --remove-duplicates   drop duplicates (default: mark only)
#       --strand STRAND       NONE | FIRST_READ_TRANSCRIPTION_STRAND | SECOND_READ_TRANSCRIPTION_STRAND
#       --copy-sorted         copy / symlink pre-sorted BAMs to $TMPDIR instead of using in place
#    -j INT  parallel samples      [4]
#    -t INT  threads per tool      [8]
#    -T DIR  temp dir (default: /tmp or $TMPDIR)
#
#  EXAMPLE (duplicates kept, strand NONE, copy sorted BAMs)
#    runPostAlignmentQC.sh \
#         -i /in -o /out \
#         -r genome.fa -e rRNA.interval_list -f refFlat.txt -b genes.bed \
#         -j 4 -t 8 --strand NONE --copy-sorted
#
#  Author : Anton Zhelonkin · 2025-05-09
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
COPY_SORTED="false"
# ---------------------------------------------------------------------------

usage() { grep '^#' "$0" | cut -c4-; exit 1; }

SHORT="i:o:r:e:f:b:j:t:T:h"
LONG="remove-duplicates,strand:,copy-sorted,ribosomal-intervals:,ref-flat:,bed-annot:,help"
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
    --copy-sorted) COPY_SORTED="true"; shift ;;
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

  # 1) obtain writable, coordinate-sorted BAM
  local sorted
  if samtools view -H "$bam" | grep -q "@HD.*SO:coordinate"; then         # already sorted
      if [[ "$COPY_SORTED" == "true" ]]; then
          sorted="${TMPDIR}/${base}.sorted.bam"
          ln -s "$bam" "$sorted"                # lightweight copy
      else
          sorted="$bam"                         # use in place (read-only OK)
      fi
  else                                                                    # need sorting
      sorted="${TMPDIR}/${base}.sorted.bam"
      run_cmd "samtools sort -@ ${THREADS} -o ${sorted} ${bam}" "$wlog"
      run_cmd "samtools index ${sorted}" "$wlog"
  fi

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

  # 5) RSeQC suite
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
       RIBO_INTERVALS BED_ANNOT TMPDIR THREADS REMOVE_DUPES STRAND COPY_SORTED

# --------------------------- run in parallel --------------------------------
find "$INPUT_DIR" -maxdepth 1 -type f -name "*.bam" | \
  parallel -j "$JOBS" process_bam {}

# --------------------------- MultiQC ----------------------------------------
STAR_DIRS=$(find "$INPUT_DIR" -maxdepth 2 -type d -name "*_STARpass1" 2>/dev/null | tr '\n' ' ')
multiqc $PICARD_DIR $RSEQ_DIR $SAMTOOLS_DIR $STAR_DIRS \
        -o "$OUTPUT_BASE_DIR/multiQC" -n "PostAlignmentQC" --force

echo -e "\n✓ All QC complete → ${OUTPUT_BASE_DIR}/multiQC/PostAlignmentQC.html"
