# Ultra-Scale Viral-Read Screening Pipeline  
### (Docker container to be executed on DNAnexus Platform)

---

## 0.  Purpose

Screen **> 500 000** human whole-genome CRAMs for **CMV, EBV and other viruses** as cheaply and quickly as possible.

* minimise **storage** (keep only KB-sized summaries ± viral-only BAM)
* run **once** per sample; if additional viruses become interesting we can re-screen the tiny residual FASTQ
* single **OCI/Docker image** → upload once to DNAnexus → run under native Docker executor  
  (DNAnexus automatically converts Docker → dx-container)

---

## 1.  High-Level Workflow

| Stage | Goal | Chosen tool (fast, parallel-ready, container-friendly) | Notes |
|-------|------|--------------------------------------------------------|-------|
| **A. Read extraction** | Pull only candidate reads from CRAM | `samtools view` (HTSlib ≥ 1.19) | Uses CRAI slicing → streams BAM blocks, keeps: <br>• unmapped primaries (flag 0x4) <br>• low-MAPQ (<30) <br>• supplementary/soft-clipped reads |
| **B. QC / trimming** | Remove adapters, Ns, very short reads | `fastp` (multithread, single binary, JSON+HTML report) | Runtime ≈ 2–3 min / sample on 16 cores |
| **C. Rapid host k-mer filter** | Drop >90 % residual human reads **before** alignment | `BBduk` Bloom filter mode **or** `Mash Screen` | 16 GB RAM Bloom of GRCh38 decoy+alt <br>Speeds up downstream by 3× |
| **D. Taxonomic screen** | Quick viral hit calling | `KrakenUniq` + `Bracken` (viral DB only, ~10 GB) | Thread-safe, scales linearly; outputs read IDs per virus |
| **E. Targeted virus alignment** | Sensitive CMV/EBV coverage stats & consensus | `minimap2 -x sr` to small panel (CMV, EBV, HSV) | < 30 s/sample; keeps BAM of viral reads only (20-100 MB **optional**) |
| **F. (Optional) Integration scan** | Host-virus chimeras for EBV | `VirusClip` (soft-clip realignment with BWA‐MEM) | Run **only** if EBV reads ≥ X to save CPU |
| **G. Summaries + QC** | TSV & json per sample | Small Python script | TSV: counts, RPM, breadth, depth, QC flags |

**Everything is a UNIX pipe**; no intermediate FASTQ written to disk.

---

## 2.  Container Layout ( `ghcr.io/yourorg/viral-screen:latest` )

* **Base**: `ubuntu:22.04`
* **Conda/Mamba** installs  
  `samtools 1.19  fastp 0.23  minimap2 2.26  krakenuniq 1.0  bracken 2.8  bbmap 39.06  mash 2.3  pigz`
* **Viral database & human Bloom** downloaded *at job start* into DNAnexus project-level cache (avoids duplicating for every job).

---

## 3.  DNAnexus Execution Sketch

```bash
dx run viral_screen \
  -iin_cram=<cram-file-id> \
  -iin_crai=<crai-id> \
  -icpu=16 -imem=32g \
  -iviral_bam=false   # skip viral-only BAMs unless needed
```

* App entrypoint is the same Docker image; no wrapper glue required.  
* Parallelisation: one **dx-job per CRAM** (embarrassingly parallel).  
  DNAnexus pools outputs into `/results/`.

---

## 4.  Storage vs Re-run Economics

| File type | Size per sample | Keep? | Rationale |
|-----------|-----------------|-------|-----------|
| Original CRAM | 90 GB | already stored by data provider |
| **Candidate FASTQ after host k-mer filter** | 50–150 MB | **No** | Can be regenerated in < 15 min from CRAM slice |
| **Kraken/Bracken reports (gz)** | 3–10 kB | **Yes** | GWAS covariates |
| Viral-only BAM | 20–100 MB | Optional (flag) | Useful for downstream variant calling; if not kept, regenerate on positives only |
| fastp JSON/HTML | 30–50 kB | Yes (optional) | QC audit trail |
| Logs | 5–10 kB | Yes | Debugging |

=> **< 25 kB mandatory storage / sample**; at 500 k samples ≈ 12 GB total.

Re-running complete screen costs ~30 CPU-minutes ≥ $0.015 (pre-emptible).  
Cheaper to **recompute** full screen than to store 50 MB FASTQ forever (>$0.01/mo).

---

## 5.  Parallelism & Throughput

* 16 CPU job → wall-time ≈ 30 min (see previous timing table).  
* 1,000 DNAnexus workers ⇒ ~32 000 samples / day.  
  → Entire 500 k cohort processed in ~16 days.

Tools chosen are multi-threaded; KrakenUniq & Mash scale almost linearly to 32 threads if nodes allow.

---

## 6.  QC Checks (minimal storage)

| Metric | Tool | Threshold |
|--------|------|-----------|
| Total reads extracted | `samtools flagstat` | > 5 k (else flag low-coverage) |
| % reads passing fastp | fastp JSON | > 80 % |
| Human contamination after Bloom + minimap2 | fraction of reads aligning to hg38 | < 5 % |
| Positive-control PhiX | KrakenUniq | < 0.1 % (lab contaminant check) |

Only the numerical values are persisted in the TSV; full JSONs optional.

---

## 7.  Alternative / Related Pipelines Reviewed

| Pipeline | Pros | Cons |
|----------|------|------|
| **GATK PathSeq v2 (Spark)** | Scales to cloud clusters; Bloom filter + BWA; well-validated | Requires Spark runtime; heavier image; DNAnexus Spark costs higher |
| **SURPI-rt** | Real-time, SNAP aligner | Not maintained; Docker image outdated |
| **Centrifuge2 + Bowtie2** | Small DB, very fast | Less sensitive than KrakenUniq for viruses |
| **WGSpipeline** (Zhang 2015) | Deep host-filter chain | Multiple slow steps; large disk footprint |

The proposed mash/bloom + KrakenUniq approach achieves ≈ PathSeq sensitivity with simpler dependencies and no Spark overhead.

---

## 8.  Next Steps

1. Build and publish `viral-screen` Docker image.  
2. Draft DNAnexus `dxapp.json` (inputs: cram, crai, cpu, mem, keep_viral_bam flag).  
3. Beta-run on 100 CRAMs; verify summary TSV versus known CMV/EBV positives.  
4. Scale-up run; monitor cost via DNAnexus job metrics.

*When satisfied we can layer in host-virus integration detection for EBV-positive subset without changing early stages.*

---
© 2025 Viral-GWAS project 