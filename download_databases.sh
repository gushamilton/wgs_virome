#!/bin/bash

# =============================================================================
# DATABASE DOWNLOAD SCRIPT FOR VIRAL SCREENING PIPELINE
# =============================================================================
# 
# DESCRIPTION:
# Downloads and sets up all required databases for the viral screening pipeline
# in the exact format expected by the Docker container.
#
# DATABASES:
# 1. Human k-mer database (/references/human_kmers.fa)
# 2. Viral Kraken database (/databases/viral_kraken)
# 3. Viral reference panel (/references/viral_panel.fa)
#
# USAGE:
# ./download_databases.sh [output_directory]
# Example: ./download_databases.sh /path/to/databases
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
NC='\033[0m' # No Color

echo -e "${GREEN}==========================================${NC}"
echo -e "${GREEN}VIRAL SCREENING DATABASE DOWNLOAD${NC}"
echo -e "${GREEN}==========================================${NC}"
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"
echo "=========================================="

# Create directory structure
echo -e "${YELLOW}Creating directory structure...${NC}"
mkdir -p "$OUTPUT_DIR"/{references,databases,logs}
cd "$OUTPUT_DIR"

# =============================================================================
# 1. HUMAN K-MER DATABASE
# =============================================================================
echo -e "${YELLOW}1. Downloading human k-mer database...${NC}"

# Download GRCh38 reference genome
echo "  Downloading GRCh38 reference genome..."
wget -O references/GRCh38.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Extract and generate k-mers
echo "  Extracting and generating k-mers..."
gunzip references/GRCh38.fa.gz

# Generate k-mers using seqtk (if available) or simple awk
if command -v seqtk &> /dev/null; then
    echo "  Using seqtk to generate k-mers..."
    seqtk seq -F 'I' references/GRCh38.fa | \
    awk 'length($0) >= 31' | \
    head -1000000 > references/human_kmers.fa
else
    echo "  Using awk to generate k-mers..."
    awk '/^>/ {next} {if (length($0) >= 31) print ">kmer_" NR "\n" substr($0, 1, 31)}' references/GRCh38.fa | \
    head -2000000 > references/human_kmers.fa
fi

echo -e "${GREEN}  ✓ Human k-mer database created: references/human_kmers.fa${NC}"

# =============================================================================
# 2. VIRAL REFERENCE PANEL
# =============================================================================
echo -e "${YELLOW}2. Downloading viral reference panel...${NC}"

# Create viral reference panel with key viruses
echo "  Downloading CMV (Human Cytomegalovirus)..."
wget -O references/cmv.fa https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14534/GCF_000845085.1_ViralProj14534_genomic.fna.gz
gunzip references/cmv.fa.gz

echo "  Downloading EBV (Epstein-Barr virus)..."
wget -O references/ebv.fa https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
gunzip references/ebv.fa.gz

echo "  Downloading HSV-1 (Herpes Simplex Virus 1)..."
wget -O references/hsv1.fa https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/725/GCF_000860725.1_ViralProj15354/GCF_000860725.1_ViralProj15354_genomic.fna.gz
gunzip references/hsv1.fa.gz

echo "  Downloading HSV-2 (Herpes Simplex Virus 2)..."
wget -O references/hsv2.fa https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/745/GCF_000860745.1_ViralProj15355/GCF_000860745.1_ViralProj15355_genomic.fna.gz
gunzip references/hsv2.fa.gz

# Combine into viral panel
echo "  Creating combined viral panel..."
cat references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa > references/viral_panel.fa

echo -e "${GREEN}  ✓ Viral reference panel created: references/viral_panel.fa${NC}"

# =============================================================================
# 3. KRAKEN VIRAL DATABASE
# =============================================================================
echo -e "${YELLOW}3. Building Kraken viral database...${NC}"

# Create Kraken database directory
mkdir -p databases/viral_kraken

# Download viral genomes from RefSeq
echo "  Downloading viral genomes from RefSeq..."
wget -O databases/viral_genomes.tar.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gunzip databases/viral_genomes.tar.gz

# Build Kraken database (this requires krakenuniq to be installed)
if command -v krakenuniq-build &> /dev/null; then
    echo "  Building KrakenUniq database..."
    krakenuniq-build --db databases/viral_kraken \
        --kmer-len 31 \
        --minimizer-spaces 6 \
        --threads $(nproc) \
        --add-to-library databases/viral_genomes.fna
    
    echo "  Finalizing database..."
    krakenuniq-build --db databases/viral_kraken --build
else
    echo -e "${RED}  ⚠ Warning: krakenuniq-build not found. Database build skipped.${NC}"
    echo "  You can build this later when the Docker container is running."
fi

echo -e "${GREEN}  ✓ Kraken database directory created: databases/viral_kraken${NC}"

# =============================================================================
# 4. CREATE MINIMAP2 INDEX
# =============================================================================
echo -e "${YELLOW}4. Creating minimap2 index...${NC}"

if command -v minimap2 &> /dev/null; then
    echo "  Building minimap2 index for viral panel..."
    minimap2 -d references/viral_panel.mmi references/viral_panel.fa
    echo -e "${GREEN}  ✓ Minimap2 index created: references/viral_panel.mmi${NC}"
else
    echo -e "${RED}  ⚠ Warning: minimap2 not found. Index build skipped.${NC}"
    echo "  You can build this later when the Docker container is running."
fi

# =============================================================================
# 5. CLEANUP AND SUMMARY
# =============================================================================
echo -e "${YELLOW}5. Cleanup and summary...${NC}"

# Remove intermediate files
rm -f references/GRCh38.fa
rm -f references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa
rm -f databases/viral_genomes.fna

# Calculate sizes
echo "Database sizes:"
du -sh references/* databases/* 2>/dev/null || echo "Some files may not exist yet"

# Create a summary file
cat > databases/README.md << EOF
# Viral Screening Pipeline Databases

## Directory Structure
\`\`\`
$OUTPUT_DIR/
├── references/
│   ├── human_kmers.fa      # Human k-mer database for host filtering
│   ├── viral_panel.fa      # Combined viral reference genomes
│   └── viral_panel.mmi     # Minimap2 index for viral panel
└── databases/
    └── viral_kraken/       # KrakenUniq viral database
\`\`\`

## Usage with Docker
\`\`\`bash
docker run -v $OUTPUT_DIR/references:/references \\
           -v $OUTPUT_DIR/databases:/databases \\
           viral-screen:v1.1 \\
           input.cram input.crai 16 32g false
\`\`\`

## Database Sources
- Human k-mers: Generated from GRCh38 reference genome
- Viral panel: CMV, EBV, HSV-1, HSV-2 from RefSeq
- Kraken database: Viral genomes from RefSeq

Generated on: $(date)
EOF

echo -e "${GREEN}==========================================${NC}"
echo -e "${GREEN}DATABASE DOWNLOAD COMPLETE${NC}"
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