#!/bin/bash
# =============================================================================
# STEP 1: DOWNLOAD RAW DATABASE FILES (RELIABLE VERSION)
# =============================================================================
#
# DESCRIPTION:
# This script uses the best downloader for each file to ensure reliability.
# - aria2c (fast) for the large GRCh38 reference from NCBI's FTP server.
# - wget (reliable) for the Kraken DB (from S3) and the smaller viral genomes.
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
CMV_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14534/GCF_000845085.1_ViralProj14534_genomic.fna.gz"
EBV_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
HSV1_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/725/GCF_000860725.1_ViralProj15354/GCF_000860725.1_ViralProj15354_genomic.fna.gz"
HSV2_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/745/GCF_000860745.1_ViralProj15355/GCF_000860745.1_ViralProj15355_genomic.fna.gz"


echo "--- Downloading GRCh38 (fast with aria2c)... ---"
aria2c -x 16 -s 16 -k 1M -d "$OUTPUT_DIR" "$GRCH38_URL"

echo "--- Downloading Kraken DB (reliable with wget)... ---"
wget -O "$OUTPUT_DIR/k2_viral_20221209.tar.gz" "$KRAKEN_VIRAL_URL"

echo "--- Downloading individual viral genomes (reliable with wget)... ---"
wget -O "$OUTPUT_DIR/cmv.fna.gz" "$CMV_URL"
wget -O "$OUTPUT_DIR/ebv.fna.gz" "$EBV_URL"
wget -O "$OUTPUT_DIR/hsv1.fna.gz" "$HSV1_URL"
wget -O "$OUTPUT_DIR/hsv2.fna.gz" "$HSV2_URL"

echo "Raw file download complete." 