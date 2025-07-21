#!/bin/bash
# =============================================================================
# COMPLETE DATABASE DOWNLOAD SCRIPT
# =============================================================================
#
# DESCRIPTION:
# Downloads ALL required files in one script:
# 1. GRCh38 reference genome (for CRAM decompression)
# 2. Kraken viral database (for classification)
# 3. CMV genome from GenBank X17403 (for alignment)
#
# USAGE: ./download_all.sh [output_directory]
# =============================================================================

set -euo pipefail

OUTPUT_DIR="${1:-./data/raw}"
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "COMPLETE DATABASE DOWNLOAD"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"
echo "=========================================="

# --- URLs ---
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
KRAKEN_VIRAL_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz"

echo "STEP 1: Downloading GRCh38 reference (fast with aria2c)..."
aria2c -c -x 16 -s 16 -k 1M -d "$OUTPUT_DIR" "$GRCH38_URL"

echo ""
echo "STEP 2: Downloading Kraken viral database (reliable with wget)..."
wget -c -O "$OUTPUT_DIR/k2_viral_20221209.tar.gz" "$KRAKEN_VIRAL_URL"

echo ""
echo "STEP 3: Downloading CMV genome from GenBank X17403..."
wget -O "$OUTPUT_DIR/cmv.fna" \
     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=X17403&rettype=fasta&retmode=text"

# Compress CMV to match expected format
echo "Compressing CMV genome..."
gzip "$OUTPUT_DIR/cmv.fna"

echo ""
echo "=========================================="
echo "ALL DOWNLOADS COMPLETE!"
echo "=========================================="
echo "End time: $(date)"
echo ""
echo "Downloaded files:"
ls -lh "$OUTPUT_DIR"/
echo ""
echo "Next step: Run ./process_all.sh to prepare the databases." 