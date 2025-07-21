#!/bin/bash

# =============================================================================
# FAST DATABASE DOWNLOAD SCRIPT FOR VIRAL SCREENING PIPELINE
# =============================================================================
# 
# DESCRIPTION:
# Downloads pre-built databases for the viral screening pipeline.
# Much faster than building from scratch!
#
# DATABASES:
# 1. Pre-built human k-mer database
# 2. Pre-built Kraken viral database
# 3. Viral reference genomes
#
# USAGE:
# ./download_databases_fast.sh [output_directory]
# Example: ./download_databases_fast.sh /path/to/databases
#
# AUTHOR: Viral-GWAS Project
# VERSION: 1.0
# DATE: January 2025
# =============================================================================

set -euo pipefail

# Default output directory
DEFAULT_OUTPUT="/opt/viral_screen_databases"
OUTPUT_DIR="${1:-$DEFAULT_OUTPUT}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}==========================================${NC}"
echo -e "${GREEN}FAST VIRAL SCREENING DATABASE DOWNLOAD${NC}"
echo -e "${GREEN}==========================================${NC}"
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"
echo "=========================================="

# Create directory structure
echo -e "${YELLOW}Creating directory structure...${NC}"
mkdir -p "$OUTPUT_DIR"/{references,databases,logs}
cd "$OUTPUT_DIR"

# Function to download with progress
download_with_progress() {
    local url="$1"
    local output="$2"
    local description="$3"
    
    echo -e "${BLUE}Downloading $description...${NC}"
    wget --progress=bar:force:noscroll -O "$output" "$url" 2>&1 | \
    while read -r line; do
        if [[ $line =~ ([0-9]+%) ]]; then
            echo -ne "\r${BLUE}Progress: ${BASH_REMATCH[1]}${NC}"
        fi
    done
    echo ""
}

# =============================================================================
# PARALLEL DOWNLOADS - ALL AT ONCE!
# =============================================================================
echo -e "${YELLOW}Starting parallel downloads...${NC}"

# Download 1: Human k-mer database (pre-built)
(
    download_with_progress \
        "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/human_genome.fna.gz" \
        "references/human_kmers.fa.gz" \
        "human k-mer database"
    gunzip references/human_kmers.fa.gz
    echo -e "${GREEN}✓ Human k-mer database ready${NC}"
) &

# Download 2: Pre-built Kraken viral database
(
    download_with_progress \
        "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz" \
        "databases/kraken_viral.tar.gz" \
        "Kraken viral database"
    echo -e "${GREEN}✓ Kraken viral database downloaded${NC}"
) &

# Download 3: Viral reference genomes (parallel)
(
    # CMV
    download_with_progress \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14534/GCF_000845085.1_ViralProj14534_genomic.fna.gz" \
        "references/cmv.fa.gz" \
        "CMV genome"
    
    # EBV
    download_with_progress \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" \
        "references/ebv.fa.gz" \
        "EBV genome"
    
    # HSV-1
    download_with_progress \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/725/GCF_000860725.1_ViralProj15354/GCF_000860725.1_ViralProj15354_genomic.fna.gz" \
        "references/hsv1.fa.gz" \
        "HSV-1 genome"
    
    # HSV-2
    download_with_progress \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/745/GCF_000860745.1_ViralProj15355/GCF_000860745.1_ViralProj15355_genomic.fna.gz" \
        "references/hsv2.fa.gz" \
        "HSV-2 genome"
    
    # Extract all
    gunzip references/*.fa.gz
    
    # Combine into viral panel
    cat references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa > references/viral_panel.fa
    
    echo -e "${GREEN}✓ Viral reference panel ready${NC}"
) &

# Wait for all downloads to complete
echo -e "${YELLOW}Waiting for downloads to complete...${NC}"
wait

# =============================================================================
# POST-PROCESSING
# =============================================================================
echo -e "${YELLOW}Post-processing downloads...${NC}"

# Extract Kraken database
echo "Extracting Kraken database..."
tar -xzf databases/kraken_viral.tar.gz -C databases/
mv databases/k2_viral_20221209 databases/viral_kraken
rm databases/kraken_viral.tar.gz

# Create minimap2 index (if minimap2 is available)
if command -v minimap2 &> /dev/null; then
    echo "Creating minimap2 index..."
    minimap2 -d references/viral_panel.mmi references/viral_panel.fa
    echo -e "${GREEN}✓ Minimap2 index created${NC}"
else
    echo -e "${YELLOW}⚠ minimap2 not found - index will be created when Docker runs${NC}"
fi

# Clean up individual viral files
rm -f references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa

# =============================================================================
# SUMMARY
# =============================================================================
echo -e "${YELLOW}Creating summary...${NC}"

# Calculate sizes
echo "Database sizes:"
du -sh references/* databases/* 2>/dev/null || echo "Some files may not exist yet"

# Create summary file
cat > databases/README.md << EOF
# Viral Screening Pipeline Databases (Fast Download)

## Directory Structure
\`\`\`
$OUTPUT_DIR/
├── references/
│   ├── human_kmers.fa      # Pre-built human k-mer database
│   ├── viral_panel.fa      # Combined viral reference genomes
│   └── viral_panel.mmi     # Minimap2 index (if minimap2 available)
└── databases/
    └── viral_kraken/       # Pre-built Kraken viral database
\`\`\`

## Usage with Docker
\`\`\`bash
docker run -v $OUTPUT_DIR/references:/references \\
           -v $OUTPUT_DIR/databases:/databases \\
           viral-screen:v1.1 \\
           input.cram input.crai 16 32g false
\`\`\`

## Database Sources
- Human k-mers: Pre-built from NCBI
- Viral panel: CMV, EBV, HSV-1, HSV-2 from RefSeq
- Kraken database: Pre-built viral database from Kraken team

## Download Time
- This fast version: ~5-10 minutes
- Original build-from-scratch: ~30-60 minutes

Generated on: $(date)
EOF

echo -e "${GREEN}==========================================${NC}"
echo -e "${GREEN}FAST DATABASE DOWNLOAD COMPLETE!${NC}"
echo -e "${GREEN}==========================================${NC}"
echo "End time: $(date)"
echo ""
echo -e "${GREEN}Next steps:${NC}"
echo "1. Mount these databases when running the Docker container:"
echo "   docker run -v $OUTPUT_DIR/references:/references \\"
echo "              -v $OUTPUT_DIR/databases:/databases \\"
echo "              viral-screen:v1.1 input.cram input.crai 16 32g false"
echo ""
echo "2. For DNAnexus deployment, upload these directories to your project."
echo ""
echo -e "${GREEN}Database summary written to: databases/README.md${NC}" 