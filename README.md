# superinfection

# Mosquito Viral Analysis Pipeline

A comprehensive bioinformatics pipeline for analyzing viral sequences in mosquito genomes from public databases. This tool identifies viral co-occurrence patterns, detects rare viruses, and analyzes viral exclusion relationships in mosquito species.

## Overview

This pipeline addresses the research question: "Which viruses co-occur together in mosquito genomes, which are rare, and what are the viral exclusion patterns?"

The analysis workflow includes:
1. Data Collection: Download mosquito genome assemblies from NCBI
2. Viral Detection: Use BLAST to identify viral sequences in genomes
3. Co-occurrence Analysis: Analyze which viruses appear together
4. Rare Virus Identification: Find viruses that appear in isolation
5. Statistical Analysis: Generate comprehensive reports and visualizations

## Features

### 🧬 Genome Data Collection
- Downloads high-quality mosquito genome assemblies from NCBI RefSeq
- Supports major disease vector species (Aedes, Culex, Anopheles)
- Automated handling of compressed genome files
- Rate-limited API calls respecting NCBI guidelines

### 🦠 Viral Sequence Detection
- BLAST-based viral sequence identification
- Uses NCBI RefSeq Viral database as reference
- Configurable significance thresholds
- Comprehensive hit filtering and annotation

### 📊 Co-occurrence Analysis
- Statistical analysis of viral co-occurrence patterns
- Identification of significant viral associations
- Detection of viral exclusion relationships
- Rare virus identification and isolation scoring

### 📈 Visualization & Reporting
- Interactive Streamlit dashboard
- Co-occurrence network visualizations
- Statistical summaries and reports
- Export capabilities for further analysis

## Installation

### Prerequisites

1. Python 3.8+
2. BLAST+ toolkit (required for viral detection)
   - macOS (Homebrew): brew install blast
   - Ubuntu/Debian: sudo apt-get install ncbi-blast+
   - Or download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

### Setup

1. Install Python dependencies
   - pip install -r requirements.txt

2. Verify BLAST installation
   - blastn -version

## Quick Start

### 1. Run the Complete Pipeline

```
# Basic usage (replace with your email)
python src/main_pipeline.py --email your.email@example.com

# Analyze specific species
python src/main_pipeline.py --email your.email@example.com --species aedes_aegypti culex_quinquefasciatus

# Use existing genomes (skip download)
python src/main_pipeline.py --email your.email@example.com --no-download
```

### 2. Launch Interactive Dashboard

```
# Start the Streamlit dashboard
streamlit run app/main_simple.py
```

Then open http://localhost:8501 in your browser.

## Pipeline Components

### Data Collection Module

File: src/data_collection/genome_downloader.py

- Downloads mosquito genome assemblies from NCBI
- Supports RefSeq and GenBank databases
- Handles multiple file formats (FASTA, GFF, protein sequences)

Key Species Supported:
- Aedes aegypti (Zika, dengue, yellow fever vector)
- Culex quinquefasciatus (West Nile virus vector)
- Anopheles gambiae (malaria vector)
- Aedes albopictus (Asian tiger mosquito)

### Viral Detection Pipeline

File: src/viral_detection/blast_pipeline.py

- Downloads NCBI RefSeq Viral database
- Creates BLAST databases for viral sequences
- Runs BLASTN searches against mosquito genomes
- Filters hits based on significance criteria:
  - E-value ≤ 1e-5
  - Identity ≥ 70%
  - Alignment length ≥ 50 bp
  - Query coverage ≥ 30%

### Co-occurrence Analysis

File: src/analysis/cooccurrence_analyzer.py

- Creates viral presence/absence matrices
- Calculates co-occurrence statistics:
  - Jaccard similarity coefficients
  - Phi correlation coefficients
  - Fisher's exact tests for associations
- Identifies rare viruses (prevalence < 5%)
- Detects viral exclusion patterns
- Generates co-occurrence networks

## Configuration

### Command Line Options

```
python src/main_pipeline.py --help
```

Key Parameters:
- --email: NCBI email (required)
- --api-key: NCBI API key (optional, recommended)
- --species: Target species list
- --max-assemblies: Max assemblies per species (default: 3)
- --threads: BLAST threads (default: 4)
- --work-dir: Output directory

### Configuration File

Create a config.json file:

```
{
  "work_dir": "mosquito_viral_analysis",
  "ncbi_email": "your.email@example.com",
  "ncbi_api_key": "your_api_key_here",
  "download_genomes": true,
  "target_species": ["aedes_aegypti", "culex_quinquefasciatus"],
  "max_assemblies_per_species": 5,
  "blast_threads": 8
}
```

Then run:
```
python src/main_pipeline.py --config config.json
```

## Output Files

The pipeline generates comprehensive output in the working directory:

### Main Results
- final_analysis_report.json: Complete analysis results
- analysis_summary.txt: Human-readable summary
- pipeline.log: Detailed execution log

### Genome Data
- genomes/: Downloaded genome assemblies organized by species
- genome_download_summary.json: Download statistics

### Viral Detection
- viral_detection/results/: BLAST results for each genome
- viral_detection/databases/: Viral reference databases
- viral_detection_summary.csv: Summary of viral hits

### Co-occurrence Analysis
- cooccurrence_analysis/viral_presence_matrix.csv: Binary presence/absence matrix
- cooccurrence_analysis/viral_cooccurrence_matrix.csv: Co-occurrence counts
- cooccurrence_analysis/significant_viral_associations.csv: Statistical associations
- cooccurrence_analysis/rare_virus_analysis.csv: Rare virus details
- cooccurrence_analysis/viral_cooccurrence_network.json: Network data

## Research Applications

### Vector Biology Research
- Understand viral diversity in disease vectors
- Identify potential viral interference patterns
- Study host-virus evolutionary relationships

### Public Health
- Monitor viral co-circulation patterns
- Identify emerging viral threats
- Assess vector competence factors

### Comparative Genomics
- Cross-species viral comparison
- Phylogenetic analysis of viral sequences
- Host range determination
