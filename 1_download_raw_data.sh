#!/bin/bash
# =============================================================================
# STEP 1: DOWNLOAD RAW DATABASE FILES (FAST)
# =============================================================================
#
# DESCRIPTION:
# This script ONLY downloads the required raw files. It is fast and can be
# run on any machine with 'aria2c' installed. It does NOT do any slow
# processing.
#
# USAGE: ./1_download_raw_data.sh /path/to/raw_data_output
# =============================================================================
set -euo pipefail
OUTPUT_DIR="${1:-./data/raw}"
mkdir -p "$OUTPUT_DIR"
echo "Downloading raw files to: $OUTPUT_DIR"

# --- URLs ---
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
KRAKEN_VIRAL_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz"
VIRAL_GENOMES_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/viral.1.1.genomic.fna.gz"

# --- Use aria2c for fast, parallel downloads ---
aria2c -x 16 -s 16 -k 1M -d "$OUTPUT_DIR" \
    "$GRCH38_URL" \
    "$KRAKEN_VIRAL_URL" \
    "$VIRAL_GENOMES_URL"

echo "Raw file download complete." 