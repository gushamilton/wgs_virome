# Viral Screening Pipeline Setup Guide

## Overview

This guide covers setting up the ultra-scale viral screening pipeline for processing >500,000 human WGS CRAMs on DNAnexus.

## Prerequisites

### Local Development
- Docker Desktop (for local testing)
- Git
- Python 3.8+ (for local script testing)

### DNAnexus Deployment
- DNAnexus account with appropriate permissions
- dx-toolkit installed and configured

## Installation

### 1. Clone and Setup

```bash
# Clone the repository
git clone <repository-url>
cd viral-screen

# Make scripts executable
chmod +x build.sh
chmod +x viral_screen.sh
chmod +x generate_summary.py
chmod +x qc_checks.py
```

### 2. Build Docker Container

```bash
# Build the container
./build.sh viral-screen:v1.0.0

# Verify the build
docker images viral-screen:v1.0.0
```

### 3. Test Locally

```bash
# Create test directories
mkdir -p test_data test_results

# Run with test data (replace with actual CRAM/CRAI files)
docker run --rm \
  -v $(pwd)/test_data:/data \
  -v $(pwd)/test_results:/results \
  viral-screen:v1.0.0 \
  sample.cram sample.crai 4 8g false
```

## DNAnexus Deployment

### 1. Push to Container Registry

```bash
# Tag for your registry
docker tag viral-screen:v1.0.0 your-registry/viral-screen:v1.0.0

# Push to registry
docker push your-registry/viral-screen:v1.0.0
```

### 2. Create DNAnexus App

Create a `dxapp.json` file:

```json
{
  "name": "viral-screen",
  "title": "Viral Screening Pipeline",
  "summary": "Ultra-scale viral screening for human WGS data",
  "version": "1.0.0",
  "inputSpec": [
    {
      "name": "cram_file",
      "label": "Input CRAM",
      "class": "file",
      "optional": false,
      "patterns": ["*.cram"]
    },
    {
      "name": "crai_file", 
      "label": "Input CRAI",
      "class": "file",
      "optional": false,
      "patterns": ["*.crai"]
    },
    {
      "name": "threads",
      "label": "Number of threads",
      "class": "int",
      "optional": true,
      "default": 16
    },
    {
      "name": "memory",
      "label": "Memory (e.g., 32g)",
      "class": "string", 
      "optional": true,
      "default": "32g"
    },
    {
      "name": "keep_viral_bam",
      "label": "Keep viral BAM",
      "class": "boolean",
      "optional": true,
      "default": false
    }
  ],
  "outputSpec": [
    {
      "name": "summary",
      "label": "Viral detection summary",
      "class": "file",
      "patterns": ["*.tsv"]
    },
    {
      "name": "qc_report",
      "label": "QC report",
      "class": "file", 
      "patterns": ["*.json"]
    },
    {
      "name": "viral_bam",
      "label": "Viral reads BAM",
      "class": "file",
      "optional": true,
      "patterns": ["*.bam"]
    }
  ],
  "runSpec": {
    "interpreter": "bash",
    "file": "src/code.sh",
    "systemRequirements": {
      "*": {
        "instanceType": "mem3_ssd1_v2_x16"
      }
    },
    "distribution": "Ubuntu",
    "release": "20.04",
    "execDepends": [
      {"name": "docker"}
    ]
  }
}
```

### 3. Deploy to DNAnexus

```bash
# Build and deploy the app
dx build viral-screen

# Test the app
dx run viral-screen -icram_file=<file-id> -icrai_file=<file-id>
```

## Database Setup

### Required Databases

The pipeline requires the following databases to be available:

1. **Human k-mer database** (`/references/human_kmers.fa`)
   - Generated from GRCh38 reference genome
   - Used for host filtering

2. **Viral Kraken database** (`/databases/viral_kraken`)
   - KrakenUniq viral-only database
   - ~10GB compressed

3. **Viral reference panel** (`/references/viral_panel.fa`)
   - FASTA file containing target viral genomes
   - CMV, EBV, HSV, etc.

### Database Generation Scripts

```bash
# Generate human k-mers
seqtk seq -F 'I' /path/to/GRCh38.fa | \
awk 'length($0) >= 31' | \
head -1000000 > /references/human_kmers.fa

# Build Kraken viral database
krakenuniq-build --db /databases/viral_kraken --threads 16 \
  --kmer-len 31 --minimizer-spaces 6 \
  --add-to-library /path/to/viral_genomes.fa

# Create viral panel
cat /path/to/cmv.fa /path/to/ebv.fa /path/to/hsv.fa > /references/viral_panel.fa
```

## Configuration

### Environment Variables

```bash
# Set in Dockerfile or runtime
export CONDA_DEFAULT_ENV=viral_env
export PATH="/opt/conda/envs/viral_env/bin:$PATH"
export PYTHONPATH="/scripts:$PYTHONPATH"
```

### Resource Requirements

| Component | CPU | Memory | Storage |
|-----------|-----|--------|---------|
| Read extraction | 4-8 | 8GB | Minimal |
| Quality control | 8-16 | 16GB | Minimal |
| Host filtering | 8-16 | 32GB | Minimal |
| Taxonomic classification | 16-32 | 64GB | 10GB DB |
| Viral alignment | 8-16 | 16GB | Minimal |

## Testing

### Unit Tests

```bash
# Test Python scripts
python3 -m pytest tests/

# Test pipeline components
./test_pipeline.sh
```

### Integration Tests

```bash
# Test with known positive samples
./test_integration.sh --positive-samples test_data/positive/

# Test with known negative samples  
./test_integration.sh --negative-samples test_data/negative/
```

## Monitoring and Logging

### Log Files

- `/results/logs/viral_screen_<sample>.log` - Main pipeline log
- `/results/qc_fastp.json` - fastp QC metrics
- `/results/qc.json` - Overall QC summary

### Performance Metrics

- Runtime per sample: ~30 minutes (16 cores)
- Storage per sample: <25KB (summary only)
- Throughput: ~32,000 samples/day (1000 workers)

## Troubleshooting

### Common Issues

1. **Docker build fails**
   - Check Docker is running
   - Ensure sufficient disk space
   - Verify network connectivity

2. **Memory errors**
   - Increase container memory allocation
   - Reduce thread count
   - Check for memory leaks

3. **Database not found**
   - Verify database paths in Dockerfile
   - Check database files are accessible
   - Rebuild databases if corrupted

### Debug Mode

```bash
# Run with debug output
docker run --rm -it \
  -v $(pwd)/test_data:/data \
  -v $(pwd)/test_results:/results \
  viral-screen:v1.0.0 \
  bash

# Inside container
./viral_screen.sh sample.cram sample.crai 4 8g false
```

## Support

For issues and questions:
- Check the logs in `/results/logs/`
- Review QC metrics in `/results/qc.json`
- Consult the README.md for detailed documentation 