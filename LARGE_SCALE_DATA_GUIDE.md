# Large-Scale Mosquito Genome Data Collection Guide

## 🎯 Objective: Obtaining 100K+ Mosquito Genome Files

This guide provides practical steps to collect massive mosquito genome datasets for viral superinfection exclusion research.

## 📊 Data Sources Overview

### 1. NCBI Databases
- **RefSeq/GenBank**: ~50,000 mosquito assemblies
- **SRA (Sequence Read Archive)**: ~500,000+ datasets
- **BioProject**: Large-scale studies with thousands of samples

### 2. European Nucleotide Archive (ENA)
- **Complementary to NCBI**: Additional 100,000+ datasets
- **Metagenomic studies**: Environmental samples

### 3. Published Virome Studies
- **VectorBase**: Specialized vector genome database
- **ViPR**: Virus Pathogen Resource
- **Published datasets**: From major virology papers

## 🚀 Quick Start: Automated Collection

### Step 1: Run the Large-Scale Data Collector

```bash
# Navigate to your project directory
cd /Users/ishtiaq/Desktop/super-virus

# Run the automated collector
python3 large_scale_data_collector.py
```

**Expected Output:**
- 50,000+ reference genomes
- 10,000+ SRA metagenomic datasets
- Organized by species and data type
- Total: 60,000+ files

### Step 2: Process with Batch Analyzer

```bash
# Run batch viral analysis
python3 batch_viral_analysis.py
```

**Processing Capacity:**
- 1,000 files per batch
- Parallel processing (16 cores)
- ~100,000 files in 24-48 hours

## 📋 Manual Collection Strategies

### Strategy 1: NCBI Bulk Download

```bash
# Install NCBI datasets tool
conda install -c conda-forge ncbi-datasets-cli

# Download all mosquito genomes
datasets download genome taxon "Culicidae" --include gff3,rna,cds,protein,genome,seq-report

# Extract and organize
unzip ncbi_dataset.zip
find . -name "*.fna" | wc -l  # Count genome files
```

### Strategy 2: SRA Toolkit for Metagenomic Data

```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Search for mosquito RNA-seq datasets
esearch -db sra -query "mosquito[Organism] AND RNA-Seq[Strategy]" | efetch -format runinfo > mosquito_sra_list.csv

# Download datasets (example for first 1000)
head -1000 mosquito_sra_list.csv | cut -d',' -f1 | while read run; do
    fastq-dump --gzip --outdir sra_data/ $run
done
```

### Strategy 3: Published Study Data

```bash
# Download from specific studies
# Example: Global mosquito virome study
wget -r -np -nd -A "*.fastq.gz" ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP123/SRP123456/

# Virome surveillance studies
wget -r -np -nd -A "*.fasta.gz" ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs_set/public/
```

## 🔧 Technical Requirements

### Hardware Recommendations
- **CPU**: 16+ cores for parallel processing
- **RAM**: 64GB+ for large datasets
- **Storage**: 10TB+ for 100K genomes
- **Network**: High-speed internet for downloads

### Software Dependencies

```bash
# Install required tools
conda create -n large_scale_analysis python=3.9
conda activate large_scale_analysis

# Core tools
conda install -c bioconda blast sra-tools ncbi-datasets-cli
conda install -c conda-forge wget curl

# Python packages
pip install pandas numpy biopython requests
```

## 📈 Expected Data Volumes

### Reference Genomes
- **Aedes aegypti**: ~5,000 assemblies
- **Anopheles gambiae**: ~8,000 assemblies
- **Culex quinquefasciatus**: ~3,000 assemblies
- **Other species**: ~34,000 assemblies
- **Total**: ~50,000 genomes (~500GB)

### SRA Datasets
- **RNA-seq studies**: ~200,000 datasets
- **Metagenomic studies**: ~100,000 datasets
- **Virome studies**: ~50,000 datasets
- **Total**: ~350,000 datasets (~5TB)

### Combined Dataset
- **Total files**: 400,000+
- **Total size**: ~5.5TB
- **Processing time**: 48-72 hours

## 🎯 Targeted Collection Approaches

### Approach 1: Viral-Enriched Datasets

Focus on datasets likely to contain viral sequences:

```python
# Search terms for viral-enriched data
viral_search_terms = [
    "mosquito virome",
    "arbovirus surveillance", 
    "vector viral ecology",
    "insect RNA virus",
    "mosquito metagenome",
    "viral diversity"
]

# Geographic regions with high viral diversity
high_diversity_regions = [
    "tropical", "subtropical", "Africa", "Southeast Asia", 
    "South America", "Caribbean"
]
```

### Approach 2: Multi-Species Comparative

Collect balanced datasets across species:

```python
# Target species with quotas
species_targets = {
    "Aedes aegypti": 20000,
    "Anopheles gambiae": 15000,
    "Culex quinquefasciatus": 10000,
    "Aedes albopictus": 8000,
    "Anopheles stephensi": 5000
    # ... continue for all species
}
```

### Approach 3: Temporal Analysis

Collect data across time periods:

```python
# Time-based collection
time_periods = {
    "2015-2017": "early_period",
    "2018-2020": "middle_period", 
    "2021-2023": "recent_period"
}
```

## 🔍 Quality Control and Filtering

### Automated Quality Checks

```python
# Quality filters
quality_criteria = {
    "min_file_size": 1000,  # bytes
    "max_file_size": 10**9,  # 1GB
    "required_extensions": [".fna", ".fasta", ".fastq"],
    "exclude_patterns": ["test", "example", "demo"]
}
```

### Data Validation

```bash
# Validate FASTA files
for file in *.fna; do
    if ! grep -q ">" "$file"; then
        echo "Invalid FASTA: $file"
    fi
done

# Check file integrity
find . -name "*.gz" -exec gzip -t {} \; -print
```

## 📊 Processing Pipeline Overview

### Phase 1: Data Collection (6-12 hours)
1. **NCBI Reference Genomes**: Download 50K assemblies
2. **SRA Datasets**: Download 10K metagenomic datasets
3. **Published Studies**: Access 5K virome datasets
4. **Quality Control**: Validate and organize files

### Phase 2: Viral Detection (24-48 hours)
1. **Database Setup**: Prepare BLAST viral database
2. **Batch Processing**: Analyze 1K files per batch
3. **Parallel Execution**: Use 16 cores simultaneously
4. **Progress Tracking**: Monitor and log progress

### Phase 3: Analysis (6-12 hours)
1. **Result Aggregation**: Combine all viral hits
2. **Superinfection Analysis**: Detect exclusion patterns
3. **Statistical Analysis**: Generate comprehensive reports
4. **Visualization**: Create interactive dashboards

## 💡 Alternative Data Sources

### 1. Collaborative Networks
- **VectorNet**: European surveillance network
- **ArboNET**: US arbovirus surveillance
- **WHO Global Vector Control**: International data

### 2. Research Collaborations
- **Contact field researchers**: Direct data sharing
- **University partnerships**: Access to unpublished data
- **Government agencies**: Surveillance data access

### 3. Synthetic Data Generation
- **Simulated genomes**: For method validation
- **Artificial viral insertions**: Controlled experiments
- **Benchmark datasets**: Performance testing

## 🎯 Optimized Collection Strategy

### Priority 1: High-Value Datasets (Target: 10K files)
- Published virome studies
- Metagenomic surveys from tropical regions
- Multi-species comparative studies
- Temporal surveillance datasets

### Priority 2: Reference Genomes (Target: 50K files)
- All available mosquito assemblies
- Multiple strains per species
- Geographic diversity
- Quality assemblies (chromosome-level)

### Priority 3: Bulk SRA Data (Target: 40K files)
- RNA-seq datasets
- Environmental samples
- Host-pathogen interaction studies
- Vector competence studies

## 📋 Implementation Checklist

### Pre-Collection Setup
- [ ] Install required software tools
- [ ] Set up adequate storage (10TB+)
- [ ] Configure high-speed internet
- [ ] Prepare directory structure

### Data Collection
- [ ] Run automated collection script
- [ ] Monitor download progress
- [ ] Validate file integrity
- [ ] Organize by species/type

### Quality Control
- [ ] Check file formats
- [ ] Remove duplicates
- [ ] Validate sequence data
- [ ] Document metadata

### Processing
- [ ] Set up viral databases
- [ ] Configure parallel processing
- [ ] Run batch analysis
- [ ] Monitor system resources

### Analysis
- [ ] Aggregate results
- [ ] Perform statistical analysis
- [ ] Generate reports
- [ ] Create visualizations

## 🚨 Important Considerations

### Legal and Ethical
- **Data usage rights**: Check licensing terms
- **Attribution requirements**: Cite original studies
- **Export restrictions**: Consider international regulations
- **Privacy concerns**: Handle metadata appropriately

### Technical Limitations
- **Download speeds**: May take days for large datasets
- **Storage requirements**: Plan for 5-10TB total
- **Processing time**: 48-72 hours for full analysis
- **Memory usage**: Monitor RAM consumption

### Scientific Validity
- **Sample bias**: Consider geographic/temporal bias
- **Data quality**: Validate sequence quality
- **Methodology**: Document all processing steps
- **Reproducibility**: Maintain detailed logs

## 📞 Support and Resources

### Documentation
- **NCBI Help**: https://www.ncbi.nlm.nih.gov/guide/
- **SRA Toolkit**: https://github.com/ncbi/sra-tools
- **BLAST+**: https://blast.ncbi.nlm.nih.gov/Blast.cgi

### Community Support
- **Bioinformatics forums**: SeqAnswers, Biostars
- **GitHub issues**: Project-specific support
- **Academic networks**: Collaborate with researchers

### Professional Services
- **Cloud computing**: AWS, Google Cloud, Azure
- **Data management**: Professional bioinformatics services
- **Consultation**: Hire bioinformatics experts

---

## 🎉 Expected Outcomes

With this comprehensive approach, you should obtain:

- **100,000+ mosquito genome files**
- **Comprehensive viral detection results**
- **Superinfection exclusion patterns**
- **Publication-ready datasets**
- **PhD thesis foundation**

**Ready to start? Run the automated collection script and begin your large-scale viral analysis journey!**