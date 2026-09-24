#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_pipeline.sh reference.fasta reads.fastq [output_dir] [original_name] [design_name]
#
# Example:
#   ./run_pipeline.sh reference.fasta reads.fastq analysis

FASTA="${1:?Usage: ./run_pipeline.sh reference.fasta reads.fastq [output_dir] [original_name] [design_name]}"
FASTQ="${2:?Usage: ./run_pipeline.sh reference.fasta reads.fastq [output_dir] [original_name] [design_name]}"
OUTDIR="${3:-variant_analysis}"
ORIGINAL_NAME="${4:-}"
DESIGN_NAME="${5:-}"

NAME_ARGS=()
if [[ -n "$ORIGINAL_NAME" ]]; then
  NAME_ARGS+=(--original-name "$ORIGINAL_NAME")
fi
if [[ -n "$DESIGN_NAME" ]]; then
  NAME_ARGS+=(--degenerate-name "$DESIGN_NAME")
fi

python3 analyze_variants.py \
  --fasta "$FASTA" \
  --fastq "$FASTQ" \
  --outdir "$OUTDIR" \
  "${NAME_ARGS[@]}" \
  --min-mapq 20 \
  --min-baseq 8
