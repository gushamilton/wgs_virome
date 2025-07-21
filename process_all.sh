#!/bin/bash
# =============================================================================
# COMPLETE DATABASE PROCESSING SCRIPT (DOCKER-BASED)
# =============================================================================
#
# DESCRIPTION:
# Processes ALL downloaded files using Docker containers:
# 1. Unpacks GRCh38, Kraken DB, and CMV
# 2. Creates human k-mer file (slow, ~30 min)
# 3. Creates minimap2 index for CMV
#
# REQUIRES:
# - Docker with viral-screen:latest image built
# - Raw files in ./data/raw/ (from download_all.sh)
#
# USAGE: ./process_all.sh
# =============================================================================

set -euo pipefail

RAW_DIR="./data/raw"
OUT_DIR="./data/final"
DOCKER_IMAGE="viral-screen:latest"

echo "=========================================="
echo "COMPLETE DATABASE PROCESSING (DOCKER)"
echo "=========================================="
echo "Raw data from: $RAW_DIR"
echo "Output to: $OUT_DIR"
echo "Start time: $(date)"
echo "=========================================="

# Check if Docker image exists
if ! docker image inspect "$DOCKER_IMAGE" >/dev/null 2>&1; then
    echo "ERROR: Docker image '$DOCKER_IMAGE' not found!"
    echo "Please build it first: docker build -t viral-screen:latest ."
    exit 1
fi

# Check if raw files exist
if [ ! -f "$RAW_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" ]; then
    echo "ERROR: GRCh38 file not found in $RAW_DIR"
    echo "Please run ./download_all.sh first"
    exit 1
fi

if [ ! -f "$RAW_DIR/cmv.fna.gz" ]; then
    echo "ERROR: CMV file not found in $RAW_DIR"
    echo "Please run ./download_all.sh first"
    exit 1
fi

if [ ! -f "$RAW_DIR/k2_viral_20221209.tar.gz" ]; then
    echo "ERROR: Kraken DB file not found in $RAW_DIR"
    echo "Please run ./download_all.sh first"
    exit 1
fi

# Create output directories
mkdir -p "$OUT_DIR"/{references,databases}

echo ""
echo "STEP 1: Unpacking all reference files..."
docker run --rm --entrypoint bash \
    -v "$PWD/$RAW_DIR":/raw \
    -v "$PWD/$OUT_DIR":/final \
    "$DOCKER_IMAGE" \
    -c "
        echo '  → Unpacking GRCh38...'
        gunzip -c /raw/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > /final/references/GRCh38.fa
        
        echo '  → Unpacking CMV genome...'
        gunzip -c /raw/cmv.fna.gz > /final/references/viral_panel.fa
        
        echo '  → Setting up Kraken DB...'
        mkdir -p /final/databases/viral_kraken
        tar -xzf /raw/k2_viral_20221209.tar.gz -C /final/databases/viral_kraken
    "

echo ""
echo "STEP 2: Creating human k-mer file (slow, ~30 min, needs 50 GB RAM)..."
docker run --rm --entrypoint bash \
    -v "$PWD/$OUT_DIR/references":/refs \
    "$DOCKER_IMAGE" \
    -c "bbmap.sh ref=/refs/GRCh38.fa out=/refs/human_kmers.fa k=31"

echo ""
echo "STEP 3: Creating minimap2 index for CMV..."
docker run --rm --entrypoint bash \
    -v "$PWD/$OUT_DIR/references":/refs \
    "$DOCKER_IMAGE" \
    -c "minimap2 -d /refs/viral_panel.mmi /refs/viral_panel.fa"

echo ""
echo "=========================================="
echo "ALL PROCESSING COMPLETE!"
echo "=========================================="
echo "End time: $(date)"
echo ""
echo "Final directory structure:"
echo "References:"
ls -la "$OUT_DIR"/references/
echo ""
echo "Databases:"
ls -la "$OUT_DIR"/databases/
echo ""
echo "Ready to run pipeline with:"
echo "docker run \\"
echo "  -v $PWD/$OUT_DIR/references:/references \\"
echo "  -v $PWD/$OUT_DIR/databases:/databases \\"
echo "  -v /path/to/crams:/crams \\"
echo "  viral-screen:latest \\"
echo "  /crams/sample.cram /crams/sample.crai /references/GRCh38.fa 16 32g false" 