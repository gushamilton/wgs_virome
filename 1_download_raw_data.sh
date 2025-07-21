#!/bin/bash
# =============================================================================
# STEP 1: DOWNLOAD RAW DATABASE FILES (CMV-ONLY VERSION)
# =============================================================================
#
# DESCRIPTION:
# Downloads only the essential files: GRCh38, Kraken viral DB, and CMV genome.
# Removes the broken EBV, HSV-1, HSV-2 downloads that keep failing.
#
# USAGE: ./1_download_raw_data.sh /path/to/raw_data_output
# =============================================================================
set -euo pipefail
OUTPUT_DIR="${1:-./data/raw}"
mkdir -p "$OUTPUT_DIR"
echo "Downloading raw files to: $OUTPUT_DIR"

# --- URLs (Only the working ones) ---
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
KRAKEN_VIRAL_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz"
CMV_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14534/GCF_000845085.1_ViralProj14534_genomic.fna.gz"

echo "--- Downloading GRCh38 (fast with aria2c)... ---"
aria2c -c -x 16 -s 16 -k 1M -d "$OUTPUT_DIR" "$GRCH38_URL"

echo "--- Downloading Kraken DB (reliable with wget)... ---"
wget -c -O "$OUTPUT_DIR/k2_viral_20221209.tar.gz" "$KRAKEN_VIRAL_URL"

echo "--- Downloading CMV genome (reliable with wget)... ---"
wget -c -O "$OUTPUT_DIR/cmv.fna.gz" "$CMV_URL"

echo "Raw file download complete."
echo ""
echo "Downloaded files:"
ls -lh "$OUTPUT_DIR"/
echo ""
echo "Next step: Run ./prepare_cmv.sh to process these files." 