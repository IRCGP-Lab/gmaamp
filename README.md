# GenoMycAnalyzer-Ampliseq (gmaamp)

**GenoMycAnalyzer-Ampliseq (gmaamp) CLI** is a comprehensive, automated bioinformatics pipeline designed for Targeted NGS analysis of *Mycobacterium tuberculosis* (MTB).

It streamlines the entire process from raw data quality control to detailed clinical reporting, integrating multiple state-of-the-art tools and custom algorithms.

## ✨ Features

* **Quality Control**: Automated trimming (Cutadapt) and QC (FastQC, MultiQC).
* **Identification**: Species identification using **Kraken2/Bracken**.
* **Variant Calling**: Robust variant detection using BWA and Bcftools.
* **Lineage & Resistance**:
  * Accurate Lineage calling.
  * Drug Resistance prediction based on **WHO v2 Database**.
  * Integration with **TB-Profiler**.
* **Reporting**: Generates comprehensive **PDF** and **Excel** reports for clinical use.

## 🚀 Quick Start

### 1. Installation

To install gmaamp via Conda:

```bash
conda create -n gmaamp -c shyun -c conda-forge -c bioconda gmaamp

conda activate gmaamp
```

### 2. Database Preparation
Due to size constraints, the Kraken2 database is NOT included in the package. You must prepare it before running the pipeline.

```bash
# We recommend using a Kraken/Bracken database specialized for AmpliSeq.
# This specialized database is available upon request (hyun@catholic.ac.kr).
# Alternatively, you can use the Pre-built Database from the Ben Langmead Lab.

# Option A: Build Standard DB
kraken2-build --standard --threads 16 --db my_kraken_db

# Option B: Download Pre-built Index (e.g., Ben Langmead's Index)
# [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)
```

### 3. Usage Example
```bash
gmaamp \
  -i /path/to/fastq_dir \
  -o /path/to/output_dir \
  -d /path/to/kraken_db \
  -t 16 \
  --run-tbprofiler \
  --pdf
```

📚 Documentation
📖 User Manual : https://github.com/IRCGP-Lab/gmaamp/wiki : **Go here for the Full Documentation.**

📝 Citation

Kim et al. (2024). "GenoMycAnalyzer: a web-based tool for species and drug resistance prediction for Mycobacterium genomes". BMC Genomics. [Link](https://link.springer.com/article/10.1186/s12864-024-10320-3)
