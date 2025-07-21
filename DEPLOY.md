# Viral Screening Pipeline - Cloud Deployment

## Quick Start

### 1. Build on Linux Cloud Instance

```bash
# Clone your repo on a Linux instance (AWS, GCP, HPC cluster)
git clone <your-repo>
cd viral-screen

# Build the container
./build_cloud.sh viral-screen:v1.0.0
```

### 2. Push to Registry

```bash
# Tag for your registry
docker tag viral-screen:v1.0.0 your-registry/viral-screen:v1.0.0

# Push to registry
docker push your-registry/viral-screen:v1.0.0
```

### 3. Deploy to DNAnexus

```bash
# Create dxapp.json (see SETUP.md for full config)
dx build viral-screen

# Test the app
dx run viral-screen -icram_file=<file-id> -icrai_file=<file-id>
```

## What's Included

**Core Tools:**
- samtools, fastp, minimap2, krakenuniq, bracken, bbmap
- Python: pandas, numpy, biopython, pysam

**Pipeline Scripts:**
- `viral_screen.sh` - Main pipeline
- `generate_summary.py` - Summary generation
- `qc_checks.py` - Quality control

**Build Scripts:**
- `build_cloud.sh` - Simple cloud build
- `Dockerfile` - Clean, minimal container

## File Structure

```
viral-screen/
├── Dockerfile              # Container definition
├── viral_screen.sh         # Main pipeline
├── generate_summary.py     # Summary generator
├── qc_checks.py           # QC checks
├── build_cloud.sh         # Cloud build script
├── README.md              # Overview
├── SETUP.md               # Detailed setup
└── DEPLOY.md              # This file
```

## Requirements

**Build System:**
- Linux (x86_64)
- Docker
- 8GB+ RAM
- 20GB+ disk space

**Runtime:**
- 16-32 CPU cores
- 32-64GB RAM
- 50GB+ disk space

## Next Steps

1. **Build on cloud** using `build_cloud.sh`
2. **Test locally** with sample data
3. **Deploy to DNAnexus** for production
4. **Scale up** to 500k+ samples

That's it! Clean and simple for cloud deployment. 