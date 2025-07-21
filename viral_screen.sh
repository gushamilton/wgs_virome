#!/bin/bash

# =============================================================================
# ULTRA-SCALE VIRAL SCREENING PIPELINE (v1.2 - Explicit GRCh38)
# =============================================================================
#
# DESCRIPTION:
# Stream-based, parallelized pipeline for screening human WGS CRAMs for
# viral sequences with minimal storage footprint and maximum speed.
#
# v1.2 Changes:
# - Added explicit GRCh38 reference FASTA as an input argument. This is
#   CRITICAL for CRAM decompression and makes the pipeline more robust
#   and portable, removing the hidden dependency on a system-wide reference.
#
# v1.1 Changes:
# - Corrected read extraction to a robust "Hybrid Strategy" (unmapped,
#   mate unmapped, or low MAPQ) using a single, efficient samtools command.
# - Implemented a true piped workflow for extraction, QC, and host filtering
#   to eliminate intermediate FASTQ files and reduce I/O.
# - Parallelized KrakenUniq and minimap2 alignment steps to run concurrently.
#
# USAGE:
# ./viral_screen.sh <input_cram> <input_crai> <grch38_fasta> <threads> <memory> <keep_viral_bam>
# Example: ./viral_screen.sh sample.cram sample.crai ref.fa 16 32g false
# =============================================================================

set -euo pipefail

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_cram> <input_crai> <grch38_fasta> <threads> <memory> <keep_viral_bam>"
    exit 1
fi

INPUT_CRAM="$1"
INPUT_CRAI="$2"
GRCH38_FASTA="$3"
THREADS="$4"
MEMORY="$5"
KEEP_VIRAL_BAM="$6"

# Check if GRCh38 reference exists
if [ ! -f "$GRCH38_FASTA" ]; then
    echo "Error: GRCh38 reference FASTA not found at '$GRCH38_FASTA'"
    exit 1
fi

# =============================================================================
# SETUP AND INITIALIZATION
# =============================================================================
SAMPLE_ID=$(basename "$INPUT_CRAM" .cram)
mkdir -p /results/logs /results/tmp

LOG_FILE="/results/logs/viral_screen_${SAMPLE_ID}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "VIRAL SCREENING PIPELINE v1.1 STARTED"
echo "Sample: $SAMPLE_ID"
echo "Threads: $THREADS"
echo "Start time: $(date)"
echo "=========================================="
START_TIME=$(date +%s)

# =============================================================================
# STAGES A, B, C: Piped Read Extraction, QC, and Host Filtering
# =============================================================================
echo "STAGES A-C: Starting piped extraction, QC, and host filtering..."

# This stream performs the first three pipeline stages without writing
# intermediate FASTQ files, significantly improving speed and reducing I/O.

# 1. samtools view: Extracts reads using the "Hybrid Strategy".
#    -T: Explicitly provides the reference FASTA for CRAM decompression.
#    -e '...' filter: gets reads that are unmapped (f4), have an unmapped
#      mate (f8), or have low mapping quality (mapq < 30).
#    -F 2304: Excludes secondary (256) and supplementary (2048) alignments.
# 2. samtools sort: Sorts reads by name (-n) so pairs are together.
# 3. samtools fastq: Converts the BAM stream to paired-end FASTQ streams.
# 4. fastp: Performs QC on the fly from piped input.
# 5. bbduk.sh: Filters out human reads from the QC'd stream.
#    Note: This uses a reference k-mer set. For maximum performance, a
#    pre-built Bloom filter can be used instead.

samtools view -u -h -T "$GRCH38_FASTA" -e 'flag.f4 || flag.f8 || mapq < 30' -F 2304 "$INPUT_CRAM" \
    | samtools sort -n -@ "$THREADS" -T "/results/tmp/collate" - \
    | samtools fastq -@ "$THREADS" -1 >(
        fastp --stdin --interleaved_in \
              --stdout \
              --qualified_quality_phred 20 \
              --length_required 50 \
              --n_base_limit 5 \
              --thread $THREADS \
              --json /results/qc_fastp.json \
              --html /results/qc_fastp.html \
              --report_title "Viral Screen QC - $SAMPLE_ID" \
        | bbduk.sh in=stdin.fastq \
                     out1=/results/tmp/nonhuman_R1.fastq \
                     out2=/results/tmp/nonhuman_R2.fastq \
                     ref=/references/human_kmers.fa \
                     k=31 hdist=1 threads=$THREADS \
                     stats=/results/tmp/bbduk_stats.txt
      ) -0 /dev/null -s /dev/null -n

NONHUMAN_READS=$(grep -c "^@" /results/tmp/nonhuman_R1.fastq)
echo "Piped workflow complete. Found $NONHUMAN_READS non-human read pairs."

# =============================================================================
# STAGE D & E: Parallel Classification and Alignment
# =============================================================================
echo "STAGE D & E: Running classification and alignment in parallel..."

# Run KrakenUniq in the background
(
    echo "Starting KrakenUniq..."
    krakenuniq \
        --db /databases/viral_kraken \
        --paired \
        --threads $THREADS \
        --report-file /results/kraken_report.txt \
        --output /results/kraken_output.txt \
        /results/tmp/nonhuman_R1.fastq \
        /results/tmp/nonhuman_R2.fastq
    echo "KrakenUniq finished."
) &
KRAKEN_PID=$!

# Run Bracken in the background, dependent on Kraken's report
(
    wait $KRAKEN_PID
    echo "Starting Bracken..."
    bracken \
        -d /databases/viral_kraken \
        -i /results/kraken_report.txt \
        -o /results/bracken_output.txt \
        -r 150 -l S -t $THREADS
    echo "Bracken finished."
) &
BRACKEN_PID=$!

# Run minimap2 alignment in the background
(
    echo "Starting minimap2 alignment..."
    minimap2 -x sr -t $THREADS -a /references/viral_panel.mmi \
        /results/tmp/nonhuman_R1.fastq \
        /results/tmp/nonhuman_R2.fastq \
        | samtools view -bS - \
        | samtools sort -@ $THREADS -o /results/tmp/viral_aligned.bam
    samtools index /results/tmp/viral_aligned.bam
    echo "minimap2 finished."
) &
MINIMAP_PID=$!

# Wait for all background jobs to complete
wait $BRACKEN_PID
wait $MINIMAP_PID
echo "All parallel analysis tasks complete."

# =============================================================================
# STAGE F: SUMMARY GENERATION
# =============================================================================
echo "STAGE F: Generating summary and QC..."

# Extract viral reads if requested
if [ "$KEEP_VIRAL_BAM" = "true" ]; then
    echo "Saving viral reads to BAM..."
    cp /results/tmp/viral_aligned.bam /results/viral_reads.bam
    cp /results/tmp/viral_aligned.bam.bai /results/viral_reads.bam.bai
fi

# Run summary generation script
python3 /scripts/generate_summary.py \
    --sample "$SAMPLE_ID" \
    --kraken-report /results/kraken_report.txt \
    --bracken-output /results/bracken_output.txt \
    --viral-bam /results/tmp/viral_aligned.bam \
    --qc-json /results/qc_fastp.json \
    --output /results/summary.tsv

# =ad hoc QC for now
CANDIDATE_READS=$(jq '.summary.before_filtering.total_reads' /results/qc_fastp.json)
QC_READS=$(jq '.summary.after_filtering.total_reads' /results/qc_fastp.json)

# Run QC checks script
python3 /scripts/qc_checks.py \
    --sample "$SAMPLE_ID" \
    --candidate-reads "$CANDIDATE_READS" \
    --qc-reads "$QC_READS" \
    --nonhuman-reads "$NONHUMAN_READS" \
    --qc-json /results/qc_fastp.json \
    --kraken-report /results/kraken_report.txt \
    --output /results/qc.json

# =============================================================================
# CLEANUP AND FINALIZATION
# =============================================================================
echo "Cleaning up temporary files..."
rm -rf /results/tmp/*

gzip /results/kraken_report.txt
gzip /results/bracken_output.txt
gzip /results/kraken_output.txt

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "=========================================="
echo "VIRAL SCREENING PIPELINE COMPLETED"
echo "Runtime: $RUNTIME seconds"
echo "=========================================="

exit 0 