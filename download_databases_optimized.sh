#!/bin/bash

# =============================================================================
# OPTIMIZED & FINAL DATABASE DOWNLOAD FOR VIRAL SCREENING PIPELINE
# =============================================================================
#
# DESCRIPTION:
# This is the fastest possible version. It uses the multi-threaded downloader
# 'aria2c' for maximum download speed and clarifies why each step is essential.
#
# REQUIRED HOST TOOLS: aria2c, wget, gzip, tar, and optionally 'bbmap.sh'.
#
# ASSETS CREATED:
# 1. /references/GRCh38.fa          (PRIMARY CRAM REFERENCE)
# 2. /references/human_kmers.fa     (for bbduk.sh host filtering, made from GRCh38.fa)
# 3. /databases/viral_kraken/       (for krakenuniq classification)
# 4. /references/viral_panel.mmi    (for minimap2 alignment)
#
# USAGE:
# ./download_databases_optimized.sh [output_directory]
# =============================================================================

set -euo pipefail

# --- Configuration ---
DEFAULT_OUTPUT="/opt/viral_screen_databases"
OUTPUT_DIR="${1:-$DEFAULT_OUTPUT}"
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
KRAKEN_VIRAL_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz"

# --- Colors for output ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}==========================================${NC}"
echo -e "${GREEN}OPTIMIZED VIRAL SCREENING DATABASE SETUP${NC}"
echo -e "${GREEN}==========================================${NC}"

if ! command -v aria2c &> /dev/null; then
    echo -e "${RED}ERROR: 'aria2c' is not installed. Please install it to use this optimized script.${NC}"
    echo "On Debian/Ubuntu: sudo apt-get update && sudo apt-get install -y aria2"
    echo "On CentOS/RHEL: sudo yum install -y aria2"
    echo "On macOS: brew install aria2"
    exit 1
fi

echo "Output directory will be: ${YELLOW}$OUTPUT_DIR${NC}"
echo "Start time: $(date)"
echo "------------------------------------------"

# --- Directory Setup ---
mkdir -p "$OUTPUT_DIR"/{references,databases}
cd "$OUTPUT_DIR"

# =============================================================================
# STEP 1: FAST DOWNLOAD OF PRIMARY GRCh38 REFERENCE
# =============================================================================
echo -e "\n${YELLOW}STEP 1: Downloading GRCh38 Reference (using multi-threaded aria2c)...${NC}"
if [ -f "references/GRCh38.fa" ]; then
    echo -e "  ${GREEN}GRCh38.fa already exists. Skipping download.${NC}"
else
    aria2c -x 16 -s 16 -k 1M -o references/GRCh38.fa.gz "$GRCH38_URL"
    gunzip references/GRCh38.fa.gz
    echo -e "  ${GREEN}✓ GRCh38.fa is ready.${NC}"
fi

# =============================================================================
# STEP 2: CREATE HUMAN K-MER FILE (ESSENTIAL FOR SPEED)
# =============================================================================
echo -e "\n${YELLOW}STEP 2: Creating Human K-mer File for Host Filtering...${NC}"
if [ -f "references/human_kmers.fa" ]; then
    echo -e "  ${GREEN}human_kmers.fa already exists. Skipping creation.${NC}"
else
    if command -v bbmap.sh &> /dev/null; then
        echo "  'bbmap.sh' found. Creating optimized k-mer file."
        bbmap.sh ref=references/GRCh38.fa out=references/human_kmers.fa k=31
        echo -e "  ${GREEN}✓ Optimized human_kmers.fa created.${NC}"
    else
        echo -e "  ${RED}WARNING: 'bbmap.sh' not found. Creating a slower, fallback k-mer file.${NC}"
        cp references/GRCh38.fa references/human_kmers.fa
        echo -e "  ${YELLOW}✓ Fallback human_kmers.fa created.${NC}"
    fi
fi

# =============================================================================
# STEP 3: FAST DOWNLOAD OF VIRAL DATABASES
# =============================================================================
echo -e "\n${YELLOW}STEP 3: Downloading Viral-Specific Databases (using multi-threaded aria2c)...${NC}"

# --- Pre-built Kraken2 Viral Database ---
if [ -d "databases/viral_kraken" ]; then
    echo "  Kraken DB exists. Skipping."
else
    echo "  Downloading pre-built Kraken viral database..."
    aria2c -x 16 -s 16 -k 1M -o databases/kraken_viral.tar.gz "$KRAKEN_VIRAL_URL"
    echo "  Extracting Kraken database..."
    tar -xzf databases/kraken_viral.tar.gz -C databases/
    mv databases/k2_viral_20221209 databases/viral_kraken
    rm databases/kraken_viral.tar.gz
    echo -e "  ${GREEN}✓ Kraken viral database is ready.${NC}"
fi

# --- Viral Reference Panel (for Minimap2) ---
if [ -f "references/viral_panel.fa" ]; then
    echo "  Viral panel exists. Skipping."
else
    echo "  Downloading viral reference genomes for alignment panel (using single-threaded wget)..."
    # Use wget for these smaller files as aria2c can be overkill
    wget -q -O references/cmv.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14534/GCF_000845085.1_ViralProj14534_genomic.fna.gz"
    wget -q -O references/ebv.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
    wget -q -O references/hsv1.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/725/GCF_000860725.1_ViralProj15354/GCF_000860725.1_ViralProj15354_genomic.fna.gz"
    wget -q -O references/hsv2.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/745/GCF_000860745.1_ViralProj15355/GCF_000860745.1_ViralProj15355_genomic.fna.gz"
    
    gunzip -f references/*.fa.gz
    cat references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa > references/viral_panel.fa
    rm references/cmv.fa references/ebv.fa references/hsv1.fa references/hsv2.fa
    echo -e "  ${GREEN}✓ Viral panel FASTA is ready.${NC}"
fi

# =============================================================================
# STEP 4: CREATE MINIMAP2 INDEX
# =============================================================================
echo -e "\n${YELLOW}STEP 4: Creating Minimap2 Index...${NC}"
if [ -f "references/viral_panel.mmi" ]; then
    echo -e "  ${GREEN}Minimap2 index already exists. Skipping.${NC}"
else
    if command -v minimap2 &> /dev/null; then
        echo "  'minimap2' found. Creating index..."
        minimap2 -d references/viral_panel.mmi references/viral_panel.fa
        echo -e "  ${GREEN}✓ Minimap2 index created.${NC}"
    else
        echo -e "  ${RED}WARNING: 'minimap2' not found. The Docker container will create it on first run.${NC}"
    fi
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo -e "\n------------------------------------------"
echo -e "${GREEN}OPTIMIZED DATABASE SETUP COMPLETE!${NC}"
echo "End time: $(date)"
echo -e "------------------------------------------\n"
echo "Final directory contents:"
ls -R 
echo -e "\nReady for Docker. Use this command, replacing paths and sample info:"
echo "docker run -v ${PWD}/references:/references \\"
echo "           -v ${PWD}/databases:/databases \\"
echo "           -v /path/to/crams:/crams \\"
echo "           viral-screen:latest \\"
echo "           /crams/sample.cram /crams/sample.crai /references/GRCh38.fa 16 32g false" 