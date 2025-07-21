#!/bin/bash

# =============================================================================
# MINIMAL DATABASE DOWNLOAD FOR VIRAL SCREENING PIPELINE
# =============================================================================
# 
# DESCRIPTION:
# Downloads ONLY what the pipeline actually needs - no extra processing!
#
# REQUIRED BY PIPELINE:
# 1. /references/human_kmers.fa     (for bbduk.sh host filtering)
# 2. /databases/viral_kraken/       (for krakenuniq classification)
# 3. /references/viral_panel.mmi    (for minimap2 alignment)
#
# USAGE:
# ./download_databases_minimal.sh [output_directory]
# Example: ./download_databases_minimal.sh /opt/viral_databases
#
# AUTHOR: Viral-GWAS Project
# VERSION: 1.0
# DATE: January 2025
# =============================================================================

set -euo pipefail

# Default output directory
DEFAULT_OUTPUT="/opt/viral_databases"
OUTPUT_DIR="${1:-$DEFAULT_OUTPUT}"

echo "=========================================="
echo "MINIMAL VIRAL SCREENING DATABASE DOWNLOAD"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"
echo "=========================================="

# Create directory structure
echo "Creating directories..."
mkdir -p "$OUTPUT_DIR"/{references,databases}
cd "$OUTPUT_DIR"

# =============================================================================
# DOWNLOAD 1: HUMAN K-MER DATABASE (for bbduk.sh)
# =============================================================================
echo "Downloading human k-mer database for host filtering..."
wget -O references/human_kmers.fa.gz \
     "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/human_genome.fna.gz"
gunzip references/human_kmers.fa.gz
echo "✓ Human k-mer database ready"

# =============================================================================
# DOWNLOAD 2: PRE-BUILT KRAKEN VIRAL DATABASE
# =============================================================================
echo "Downloading pre-built Kraken viral database..."
wget -O databases/kraken_viral.tar.gz \
     "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz"
echo "Extracting Kraken database..."
tar -xzf databases/kraken_viral.tar.gz -C databases/
mv databases/k2_viral_20221209 databases/viral_kraken
rm databases/kraken_viral.tar.gz
echo "✓ Kraken viral database ready"

# =============================================================================
# DOWNLOAD 3: VIRAL REFERENCE GENOMES (for minimap2)
# =============================================================================
echo "Downloading viral reference genomes..."

# Download viral genomes in parallel
(
    wget -O references/cmv.fa.gz \
         "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14534/GCF_000845085.1_ViralProj14534_genomic.fna.gz" &
    
    wget -O references/ebv.fa.gz \
         "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" &
    
    wget -O references/hsv1.fa.gz \
         "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/725/GCF_000860725.1_ViralProj15354/GCF_000860725.1_ViralProj15354_genomic.fna.gz" &
    
    wget -O references/hsv2.fa.gz \
         "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/745/GCF_000860745.1_ViralProj15355/GCF_000860745.1_ViralProj15355_genomic.fna.gz" &
    
    wait
)

# Extract and combine
echo "Creating viral panel..."
gunzip references/*.fa.gz
cat references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa > references/viral_panel.fa
rm references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa
echo "✓ Viral panel ready"

# =============================================================================
# CREATE MINIMAP2 INDEX (if minimap2 available)
# =============================================================================
if command -v minimap2 &> /dev/null; then
    echo "Creating minimap2 index..."
    minimap2 -d references/viral_panel.mmi references/viral_panel.fa
    echo "✓ Minimap2 index created"
else
    echo "⚠ minimap2 not found - index will be created when Docker runs"
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo "=========================================="
echo "MINIMAL DATABASE DOWNLOAD COMPLETE!"
echo "=========================================="
echo "End time: $(date)"
echo ""
echo "Database sizes:"
du -sh references/* databases/* 2>/dev/null || echo "Some files may not exist yet"
echo ""
echo "Ready for Docker:"
echo "docker run -v $OUTPUT_DIR/references:/references \\"
echo "           -v $OUTPUT_DIR/databases:/databases \\"
echo "           viral-screen:v1.1 \\"
echo "           input.cram input.crai 16 32g false" 