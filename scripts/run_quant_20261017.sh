#!/usr/bin/env bash
set -euo pipefail

# Activate the independently installed quantification environment before running.
PROCESSED_DIR="/mnt/sda/lxk/database/WWS/Probe_panel_test_Proj/Processed_fastq/Processed_fastq_HybCapSeq_Quant_test_20261017"
ESVIRITU_DIR="${PROCESSED_DIR}/EsViritu_analysis_HybCapSeq_Quant_test_20261017"
INPUT_DIR="${PROCESSED_DIR}/HybCapSeq_Quant_inputs"
RESULT_DIR="${PROCESSED_DIR}/HybCapSeq_Quant_results"

hybcapseq_quant \
    --samples "${INPUT_DIR}/samples_20261017.tsv" \
    --esviritu-dir "$ESVIRITU_DIR" \
    --ercc-standards "${INPUT_DIR}/ercc92_mix1.tsv" \
    --outdir "$RESULT_DIR" \
    --method both \
    --wastewater-unit copies/mL \
    "$@"
