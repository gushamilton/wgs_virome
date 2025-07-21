#!/bin/bash

# =============================================================================
# ULTRA-SCALE VIRAL SCREENING PIPELINE
# =============================================================================
# 
# DESCRIPTION:
# Stream-based pipeline for screening human WGS CRAMs for viral sequences
# with minimal storage footprint and maximum speed.
#
# WORKFLOW:
# 1. Extract candidate reads from CRAM (unmapped, low-MAPQ, soft-clipped)
# 2. Quality control with fastp
# 3. Host k-mer filtering with BBduk Bloom filter
# 4. Taxonomic classification with KrakenUniq + Bracken
# 5. Targeted virus alignment with minimap2
# 6. Summary generation and QC
#
# USAGE:
# ./viral_screen.sh <input_cram> <input_crai> <threads> <memory> <keep_viral_bam>
# Example: ./viral_screen.sh sample.cram sample.crai 16 32g false
#
# OUTPUT:
# - /results/summary.tsv: Viral detection summary
# - /results/qc.json: Quality control metrics
# - /results/viral_reads.bam: Viral reads only (optional)
# - /results/logs/: Processing logs
#
# AUTHOR: Viral-GWAS Project
# VERSION: 1.0
# DATE: January 2025
# =============================================================================

set -euo pipefail

# =============================================================================
# INPUT VALIDATION
# =============================================================================

# Check if correct number of arguments provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_cram> <input_crai> <threads> <memory> <keep_viral_bam>"
    echo "Example: $0 sample.cram sample.crai 16 32g false"
    exit 1
fi

# Parse arguments
INPUT_CRAM="$1"
INPUT_CRAI="$2"
THREADS="$3"
MEMORY="$4"
KEEP_VIRAL_BAM="$5"

# Validate inputs
if [ ! -f "$INPUT_CRAM" ]; then
    echo "Error: CRAM file not found: $INPUT_CRAM"
    exit 1
fi

if [ ! -f "$INPUT_CRAI" ]; then
    echo "Error: CRAI file not found: $INPUT_CRAI"
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [ "$THREADS" -lt 1 ]; then
    echo "Error: Threads must be a positive integer"
    exit 1
fi

if [ "$KEEP_VIRAL_BAM" != "true" ] && [ "$KEEP_VIRAL_BAM" != "false" ]; then
    echo "Error: keep_viral_bam must be 'true' or 'false'"
    exit 1
fi

# =============================================================================
# SETUP AND INITIALIZATION
# =============================================================================

# Get sample ID from CRAM filename
SAMPLE_ID=$(basename "$INPUT_CRAM" .cram)

# Create output directories
mkdir -p /results/logs /results/tmp

# Setup logging
LOG_FILE="/results/logs/viral_screen_${SAMPLE_ID}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "VIRAL SCREENING PIPELINE STARTED"
echo "=========================================="
echo "Sample: $SAMPLE_ID"
echo "CRAM: $INPUT_CRAM"
echo "CRAI: $INPUT_CRAI"
echo "Threads: $THREADS"
echo "Memory: $MEMORY"
echo "Keep viral BAM: $KEEP_VIRAL_BAM"
echo "Start time: $(date)"
echo "=========================================="

START_TIME=$(date +%s)

# =============================================================================
# STAGE A: READ EXTRACTION
# =============================================================================
echo "STAGE A: Extracting candidate reads from CRAM..."

# Extract reads that could be viral:
# - unmapped primary reads (flag 0x4)
# - low MAPQ reads (< 30)
# - supplementary/soft-clipped reads
# - reads mapping to non-autosomal regions

echo "Extracting unmapped and low-quality reads..."
samtools view -u -f 4 -F 256 "$INPUT_CRAM" | \
samtools collate -u -O - | \
samtools fastq -1 /results/tmp/unmapped_R1.fastq -2 /results/tmp/unmapped_R2.fastq -0 /dev/null -s /dev/null -n - &

echo "Extracting low MAPQ reads..."
samtools view -u -q 29 "$INPUT_CRAM" | \
samtools collate -u -O - | \
samtools fastq -1 /results/tmp/lowmapq_R1.fastq -2 /results/tmp/lowmapq_R2.fastq -0 /dev/null -s /dev/null -n - &

echo "Extracting supplementary reads..."
samtools view -u -f 2048 "$INPUT_CRAM" | \
samtools collate -u -O - | \
samtools fastq -1 /results/tmp/supplementary_R1.fastq -2 /results/tmp/supplementary_R2.fastq -0 /dev/null -s /dev/null -n - &

wait

# Combine all candidate reads
echo "Combining candidate reads..."
cat /results/tmp/unmapped_R1.fastq /results/tmp/lowmapq_R1.fastq /results/tmp/supplementary_R1.fastq > /results/tmp/candidate_R1.fastq
cat /results/tmp/unmapped_R2.fastq /results/tmp/lowmapq_R2.fastq /results/tmp/supplementary_R2.fastq > /results/tmp/candidate_R2.fastq

# Count extracted reads
CANDIDATE_READS=$(grep -c "^@" /results/tmp/candidate_R1.fastq)
echo "Extracted $CANDIDATE_READS candidate read pairs"

# =============================================================================
# STAGE B: QUALITY CONTROL
# =============================================================================
echo "STAGE B: Quality control with fastp..."

fastp \
    --in1 /results/tmp/candidate_R1.fastq \
    --in2 /results/tmp/candidate_R2.fastq \
    --out1 /results/tmp/qc_R1.fastq \
    --out2 /results/tmp/qc_R2.fastq \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --n_base_limit 5 \
    --thread $THREADS \
    --json /results/qc_fastp.json \
    --html /results/qc_fastp.html \
    --report_title "Viral Screen QC - $SAMPLE_ID"

# Count QC-passed reads
QC_READS=$(grep -c "^@" /results/tmp/qc_R1.fastq)
echo "QC passed: $QC_READS read pairs"

# =============================================================================
# STAGE C: HOST K-MER FILTERING
# =============================================================================
echo "STAGE C: Host k-mer filtering with BBduk..."

# Note: In production, the human Bloom filter would be pre-built and mounted
# For now, we'll use a simple approach with BBduk
bbduk.sh \
    in1=/results/tmp/qc_R1.fastq \
    in2=/results/tmp/qc_R2.fastq \
    out1=/results/tmp/nonhuman_R1.fastq \
    out2=/results/tmp/nonhuman_R2.fastq \
    ref=/references/human_kmers.fa \
    k=31 \
    hdist=1 \
    threads=$THREADS \
    stats=/results/tmp/bbduk_stats.txt

# Count non-human reads
NONHUMAN_READS=$(grep -c "^@" /results/tmp/nonhuman_R1.fastq)
echo "Non-human reads: $NONHUMAN_READS read pairs"

# =============================================================================
# STAGE D: TAXONOMIC CLASSIFICATION
# =============================================================================
echo "STAGE D: Taxonomic classification with KrakenUniq..."

# Run KrakenUniq on non-human reads
krakenuniq \
    --db /databases/viral_kraken \
    --paired \
    --threads $THREADS \
    --report /results/kraken_report.txt \
    --output /results/kraken_output.txt \
    /results/tmp/nonhuman_R1.fastq \
    /results/tmp/nonhuman_R2.fastq

# Run Bracken for abundance estimation
bracken \
    -d /databases/viral_kraken \
    -i /results/kraken_report.txt \
    -o /results/bracken_output.txt \
    -r 150 \
    -l S \
    -t $THREADS

# =============================================================================
# STAGE E: TARGETED VIRUS ALIGNMENT
# =============================================================================
echo "STAGE E: Targeted virus alignment with minimap2..."

# Create viral references index if not exists
if [ ! -f "/references/viral_panel.mmi" ]; then
    echo "Building viral panel index..."
    minimap2 -d /references/viral_panel.mmi /references/viral_panel.fa
fi

# Align to viral panel
minimap2 \
    -x sr \
    -t $THREADS \
    -a \
    /references/viral_panel.mmi \
    /results/tmp/nonhuman_R1.fastq \
    /results/tmp/nonhuman_R2.fastq \
    | samtools view -bS - \
    | samtools sort -@ $THREADS -o /results/tmp/viral_aligned.bam

# Index the BAM
samtools index /results/tmp/viral_aligned.bam

# Extract viral reads if requested
if [ "$KEEP_VIRAL_BAM" = "true" ]; then
    echo "Extracting viral reads to BAM..."
    samtools view -b -F 4 /results/tmp/viral_aligned.bam > /results/viral_reads.bam
    samtools index /results/viral_reads.bam
    echo "Viral BAM saved: /results/viral_reads.bam"
fi

# =============================================================================
# STAGE F: SUMMARY GENERATION
# =============================================================================
echo "STAGE F: Generating summary and QC..."

# Run summary generation script
python3 /scripts/generate_summary.py \
    --sample "$SAMPLE_ID" \
    --kraken-report /results/kraken_report.txt \
    --bracken-output /results/bracken_output.txt \
    --viral-bam /results/tmp/viral_aligned.bam \
    --qc-json /results/qc_fastp.json \
    --output /results/summary.tsv

# Run QC checks
python3 /scripts/qc_checks.py \
    --sample "$SAMPLE_ID" \
    --candidate-reads "$CANDIDATE_READS" \
    --qc-reads "$QC_READS" \
    --nonhuman-reads "$NONHUMAN_READS" \
    --qc-json /results/qc_fastp.json \
    --output /results/qc.json

# =============================================================================
# CLEANUP AND FINALIZATION
# =============================================================================
echo "Cleaning up temporary files..."
rm -rf /results/tmp/*

# Compress outputs
gzip /results/kraken_report.txt
gzip /results/bracken_output.txt
gzip /results/kraken_output.txt

# Calculate runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "=========================================="
echo "VIRAL SCREENING PIPELINE COMPLETED"
echo "=========================================="
echo "Sample: $SAMPLE_ID"
echo "Runtime: $RUNTIME seconds"
echo "End time: $(date)"
echo "Outputs:"
echo "  - /results/summary.tsv"
echo "  - /results/qc.json"
echo "  - /results/kraken_report.txt.gz"
echo "  - /results/bracken_output.txt.gz"
if [ "$KEEP_VIRAL_BAM" = "true" ]; then
    echo "  - /results/viral_reads.bam"
fi
echo "=========================================="

exit 0 