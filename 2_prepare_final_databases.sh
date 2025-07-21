#!/bin/bash
# =============================================================================
# STEP 2: PREPARE FINAL DATABASES (SLOW, REQUIRES DOCKER)
# =============================================================================
#
# DESCRIPTION:
# This script runs the slow, heavy processing steps (k-mer generation and
# indexing) inside the Docker container, which guarantees all tools are
# available. Run this script AFTER '1_download_raw_data.sh'.
#
# USAGE: ./2_prepare_final_databases.sh /path/to/raw_data /path/to/final_db_output
# =============================================================================
set -euo pipefail
RAW_DIR="${1:-./data/raw}"
FINAL_DIR="${2:-./data/final}"

# --- Setup Directories ---
mkdir -p "$FINAL_DIR"/{references,databases}
DOCKER_IMAGE="viral-screen:latest"
DOCKER_CMD="docker run --rm -v ${PWD}/${RAW_DIR}:/raw_data -v ${PWD}/${FINAL_DIR}:/final_data"

echo "Preparing final databases in: $FINAL_DIR"
echo "Using raw data from: $RAW_DIR"

# --- Run Preparation inside Docker ---
$DOCKER_CMD $DOCKER_IMAGE bash -c '
    set -euxo pipefail
    echo "--- STEP 2.1: Unpacking raw files ---"
    gunzip -c /raw_data/*.fna.gz > /final_data/references/GRCh38.fa
    gunzip -c /raw_data/viral.1.1.genomic.fna.gz > /final_data/references/viral_panel.fa
    
    echo "--- STEP 2.2: Creating Human K-mers (slow step)... ---"
    bbmap.sh ref=/final_data/references/GRCh38.fa out=/final_data/references/human_kmers.fa k=31

    echo "--- STEP 2.3: Preparing Kraken Database ---"
    mkdir -p /final_data/databases/viral_kraken
    tar -xzf /raw_data/k2_viral_20221209.tar.gz -C /final_data/databases/viral_kraken
    
    echo "--- STEP 2.4: Creating Minimap2 Index ---"
    minimap2 -d /final_data/references/viral_panel.mmi /final_data/references/viral_panel.fa
    
    echo "--- FINAL DATABASE PREPARATION COMPLETE ---"
'
echo "All databases are now ready in: $FINAL_DIR" 