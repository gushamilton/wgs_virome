#!/bin/bash
# =============================================================================
# DIRECT CMV GENOME DOWNLOAD FROM GENBANK
# =============================================================================
#
# DESCRIPTION:
# Downloads the CMV reference genome directly from GenBank record X17403
# (Human cytomegalovirus strain AD169 complete genome) using NCBI E-utilities
#
# USAGE: ./download_cmv_direct.sh [output_directory]
# =============================================================================

set -euo pipefail

OUTPUT_DIR="${1:-./data/raw}"
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "DIRECT CMV GENOME DOWNLOAD"
echo "=========================================="
echo "Downloading from: https://www.ncbi.nlm.nih.gov/nuccore/X17403"
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"
echo "=========================================="

# Download CMV genome using NCBI E-utilities (more reliable)
echo "Downloading CMV genome (X17403) from NCBI E-utilities..."
wget -O "$OUTPUT_DIR/cmv.fna" \
     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=X17403&rettype=fasta&retmode=text"

# Check if download was successful
if [ -f "$OUTPUT_DIR/cmv.fna" ]; then
    echo "✓ CMV genome downloaded successfully"
    echo "File size: $(ls -lh "$OUTPUT_DIR/cmv.fna" | awk '{print $5}')"
    
    # Compress it to match our expected format
    echo "Compressing to gzip format..."
    gzip "$OUTPUT_DIR/cmv.fna"
    
    echo "✓ CMV genome ready at: $OUTPUT_DIR/cmv.fna.gz"
else
    echo "✗ Download failed"
    exit 1
fi

echo ""
echo "CMV genome is ready at: $OUTPUT_DIR/cmv.fna.gz"
echo "You can now run: ./prepare_cmv.sh" 