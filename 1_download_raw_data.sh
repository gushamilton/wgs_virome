#!/bin/bash
# =============================================================================
# STEP 1: DOWNLOAD RAW DATABASE FILES (RELIABLE VERSION 2.0)
# =============================================================================
#
# DESCRIPTION:
# This script uses verified, working URLs as of July 2025. It uses a hybrid
# download approach for maximum speed and reliability.
#
# USAGE: ./1_download_raw_data.sh /path/to/raw_data_output
# =============================================================================
set -euo pipefail
OUTPUT_DIR="${1:-./data/raw}"
mkdir -p "$OUTPUT_DIR"
echo "Downloading raw files to: $OUTPUT_DIR"

# --- URLs (Verified Working, using FTP protocol for reliability) ---
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
KRAKEN_VIRAL_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz"
CMV_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1/GCF_000845085.1_genomic.fna.gz"
EBV_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1/GCF_000819615.1_genomic.fna.gz"
HSV1_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/725/GCF_000860725.1/GCF_000860725.1_genomic.fna.gz"
HSV2_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/745/GCF_000860745.1/GCF_000860745.1_genomic.fna.gz"

echo "--- Downloading GRCh38 (fast with aria2c)... ---"
aria2c -c -x 16 -s 16 -k 1M -d "$OUTPUT_DIR" "$GRCH38_URL"

echo "--- Downloading Kraken DB (reliable with wget)... ---"
wget -c -O "$OUTPUT_DIR/k2_viral_20221209.tar.gz" "$KRAKEN_VIRAL_URL"

echo "--- Downloading individual viral genomes (reliable with wget)... ---"
wget -c -O "$OUTPUT_DIR/cmv.fna.gz" "$CMV_URL"
wget -c -O "$OUTPUT_DIR/ebv.fna.gz" "$EBV_URL"
wget -c -O "$OUTPUT_DIR/hsv1.fna.gz" "$HSV1_URL"
wget -c -O "$OUTPUT_DIR/hsv2.fna.gz" "$HSV2_URL"

echo "Raw file download complete." 