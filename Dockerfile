# =============================================================================
# ULTRA-SCALE VIRAL SCREENING PIPELINE DOCKER CONTAINER
# =============================================================================
# 
# DESCRIPTION:
# Docker container for screening >500,000 human WGS CRAMs for viral sequences
# (CMV, EBV, HSV, etc.) with minimal storage footprint and maximum speed.
#
# USAGE:
# docker run -v /data:/data -v /results:/results viral-screen:latest \
#   input.cram input.crai 16 32g false
#
# AUTHOR: Viral-GWAS Project
# VERSION: 1.0
# DATE: January 2025
# =============================================================================

FROM ubuntu:22.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip \
    python3-dev \
    gfortran \
    libgsl-dev \
    libboost-all-dev \
    cmake \
    pkg-config \
    autoconf \
    automake \
    libtool \
    bzip2 \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install Mambaforge
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh" \
    && bash Mambaforge-Linux-x86_64.sh -b -p /opt/conda \
    && rm Mambaforge-Linux-x86_64.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Create conda environment and install bioinformatics tools
RUN mamba create -n viral_env python=3.10 -y \
    && mamba install -n viral_env -c bioconda -c conda-forge \
    samtools=1.19 \
    fastp=0.23.4 \
    minimap2=2.26 \
    krakenuniq=1.0.3 \
    bracken=2.8 \
    bbmap=39.06 \
    pigz=2.8 \
    parallel=20240222 \
    bwa=0.7.17 \
    blast=2.15.0 \
    bedtools=2.31.0 \
    seqtk=1.4 \
    htslib=1.19 \
    bcftools=1.19 \
    -y

# Install Python packages
RUN pip install \
    pandas>=2.0.0 \
    numpy>=1.24.0 \
    biopython>=1.81 \
    pysam>=0.21.0 \
    click>=8.1.0 \
    tqdm>=4.65.0

# Create directories
RUN mkdir -p /data /results /references /databases /scripts /tmp

# Set working directory
WORKDIR /data

# Copy pipeline scripts
COPY viral_screen.sh /scripts/viral_screen.sh
COPY generate_summary.py /scripts/generate_summary.py
COPY qc_checks.py /scripts/qc_checks.py

RUN chmod +x /scripts/viral_screen.sh \
    && chmod +x /scripts/generate_summary.py \
    && chmod +x /scripts/qc_checks.py

# Set environment variables
ENV CONDA_DEFAULT_ENV=viral_env
ENV PATH="/opt/conda/envs/viral_env/bin:$PATH"
ENV PYTHONPATH="/scripts:$PYTHONPATH"

# Set entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "viral_env", "/scripts/viral_screen.sh"] 