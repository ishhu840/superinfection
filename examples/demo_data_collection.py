#!/usr/bin/env python3
"""
Demo Data Collection Script

Demonstrates how to collect mosquito genome data starting small and scaling up.
This script shows the process with a manageable dataset before attempting 100K+ files.
"""

import os
import sys
import json
import time
import logging
import requests
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class DemoDataCollector:
    """Demonstrates data collection process with a small, manageable dataset."""
    
    def __init__(self, output_dir: str = "demo_mosquito_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        self.genomes_dir = self.output_dir / "genomes"
        self.sra_dir = self.output_dir / "sra_samples"
        self.metadata_dir = self.output_dir / "metadata"
        
        for dir_path in [self.genomes_dir, self.sra_dir, self.metadata_dir]:
            dir_path.mkdir(exist_ok=True)
        
        # Demo targets (small numbers for demonstration)
        self.demo_species = [
            "Aedes aegypti",
            "Anopheles gambiae", 
            "Culex quinquefasciatus"
        ]
        
        self.collection_stats = {
            "start_time": None,
            "genomes_collected": 0,
            "sra_samples_collected": 0,
            "species_covered": set(),
            "total_size_mb": 0
        }
    
    def get_sample_ncbi_genomes(self, max_per_species: int = 5) -> List[Dict]:
        """Get a small sample of NCBI genomes for demonstration."""
        logger.info(f"Collecting sample genomes ({max_per_species} per species)...")
        
        sample_genomes = []
        
        for species in self.demo_species:
            try:
                # Search for assemblies
                search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                search_params = {
                    "db": "assembly",
                    "term": f'"{species}"[Organism] AND "latest refseq"[filter]',
                    "retmax": max_per_species,
                    "retmode": "json"
                }
                
                response = requests.get(search_url, params=search_params)
                if response.status_code == 200:
                    data = response.json()
                    assembly_ids = data.get("esearchresult", {}).get("idlist", [])
                    
                    logger.info(f"Found {len(assembly_ids)} assemblies for {species}")
                    
                    # Get assembly details
                    if assembly_ids:
                        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                        summary_params = {
                            "db": "assembly",
                            "id": ",".join(assembly_ids),
                            "retmode": "json"
                        }
                        
                        summary_response = requests.get(summary_url, params=summary_params)
                        if summary_response.status_code == 200:
                            summary_data = summary_response.json()
                            
                            for assembly_id in assembly_ids:
                                if assembly_id in summary_data.get("result", {}):
                                    assembly_info = summary_data["result"][assembly_id]
                                    sample_genomes.append({
                                        "species": species,
                                        "assembly_id": assembly_id,
                                        "accession": assembly_info.get("assemblyaccession", ""),
                                        "name": assembly_info.get("assemblyname", ""),
                                        "ftp_path": assembly_info.get("ftppath_genbank", ""),
                                        "size_mb": assembly_info.get("total_length", 0) / 1000000
                                    })
                
                time.sleep(1)  # Rate limiting
                
            except Exception as e:
                logger.error(f"Error collecting genomes for {species}: {e}")
        
        logger.info(f"Collected {len(sample_genomes)} sample genomes")
        return sample_genomes
    
    def get_sample_sra_datasets(self, max_per_species: int = 3) -> List[Dict]:
        """Get a small sample of SRA datasets for demonstration."""
        logger.info(f"Collecting sample SRA datasets ({max_per_species} per species)...")
        
        sample_sra = []
        
        for species in self.demo_species:
            try:
                # Search for RNA-seq datasets
                search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                search_params = {
                    "db": "sra",
                    "term": f'"{species}"[Organism] AND "RNA-Seq"[Strategy]',
                    "retmax": max_per_species,
                    "retmode": "json"
                }
                
                response = requests.get(search_url, params=search_params)
                if response.status_code == 200:
                    data = response.json()
                    sra_ids = data.get("esearchresult", {}).get("idlist", [])
                    
                    logger.info(f"Found {len(sra_ids)} SRA datasets for {species}")
                    
                    for sra_id in sra_ids:
                        sample_sra.append({
                            "species": species,
                            "sra_id": sra_id,
                            "type": "RNA-Seq"
                        })
                
                time.sleep(1)  # Rate limiting
                
            except Exception as e:
                logger.error(f"Error collecting SRA for {species}: {e}")
        
        logger.info(f"Collected {len(sample_sra)} sample SRA datasets")
        return sample_sra
    
    def download_sample_genomes(self, genomes: List[Dict]) -> None:
        """Download sample genome files."""
        logger.info("Downloading sample genome files...")
        
        for genome in genomes[:10]:  # Limit to 10 for demo
            try:
                species_dir = self.genomes_dir / genome["species"].replace(" ", "_")
                species_dir.mkdir(exist_ok=True)
                
                accession = genome["accession"]
                ftp_path = genome["ftp_path"]
                
                if ftp_path:
                    # Download genomic FASTA
                    fasta_url = f"{ftp_path}/{accession}_genomic.fna.gz"
                    fasta_file = species_dir / f"{accession}_genomic.fna.gz"
                    
                    if not fasta_file.exists():
                        logger.info(f"Downloading {accession}...")
                        subprocess.run(["wget", "-q", fasta_url, "-O", str(fasta_file)], 
                                     check=True, timeout=120)
                        
                        # Update stats
                        if fasta_file.exists():
                            size_mb = fasta_file.stat().st_size / (1024 * 1024)
                            self.collection_stats["genomes_collected"] += 1
                            self.collection_stats["total_size_mb"] += size_mb
                            self.collection_stats["species_covered"].add(genome["species"])
                            
                            logger.info(f"Downloaded {accession} ({size_mb:.1f} MB)")
                    else:
                        logger.info(f"Already exists: {accession}")
                        
            except Exception as e:
                logger.error(f"Failed to download {genome.get('accession', 'unknown')}: {e}")
    
    def create_demo_sra_files(self, sra_datasets: List[Dict]) -> None:
        """Create demo SRA files (placeholder for actual download)."""
        logger.info("Creating demo SRA files...")
        
        # For demo purposes, create small placeholder files
        # In real implementation, you would use fastq-dump
        
        for sra in sra_datasets[:5]:  # Limit to 5 for demo
            try:
                species_dir = self.sra_dir / sra["species"].replace(" ", "_")
                species_dir.mkdir(exist_ok=True)
                
                sra_id = sra["sra_id"]
                demo_file = species_dir / f"{sra_id}_demo.fastq"
                
                if not demo_file.exists():
                    # Create a small demo FASTQ file
                    demo_content = f"""@{sra_id}_read_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@{sra_id}_read_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
                    
                    with open(demo_file, 'w') as f:
                        f.write(demo_content)
                    
                    self.collection_stats["sra_samples_collected"] += 1
                    logger.info(f"Created demo file: {sra_id}")
                    
            except Exception as e:
                logger.error(f"Failed to create demo file for {sra.get('sra_id', 'unknown')}: {e}")
    
    def generate_demo_summary(self) -> None:
        """Generate summary of demo collection."""
        logger.info("Generating demo collection summary...")
        
        # Update final stats
        self.collection_stats["species_covered"] = list(self.collection_stats["species_covered"])
        self.collection_stats["end_time"] = datetime.now().isoformat()
        
        # Save summary
        summary_file = self.metadata_dir / "demo_collection_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(self.collection_stats, f, indent=2)
        
        # Create README
        readme_content = f"""
# Demo Mosquito Genome Data Collection

## Collection Summary
- **Collection Date**: {self.collection_stats['start_time']}
- **Genomes Collected**: {self.collection_stats['genomes_collected']}
- **SRA Samples**: {self.collection_stats['sra_samples_collected']}
- **Species Covered**: {len(self.collection_stats['species_covered'])}
- **Total Size**: {self.collection_stats['total_size_mb']:.1f} MB

## Species Analyzed
{chr(10).join(f'- {species}' for species in self.collection_stats['species_covered'])}

## Directory Structure
```
{self.output_dir.name}/
├── genomes/           # Reference genome assemblies
├── sra_samples/       # SRA dataset samples
└── metadata/          # Collection metadata
```

## Next Steps

### 1. Run Viral Analysis
```bash
# Analyze the demo dataset
python3 analyze_downloaded_genomes.py
```

### 2. Scale Up Collection
```bash
# Collect larger dataset (1000+ files)
python3 large_scale_data_collector.py
```

### 3. Batch Processing
```bash
# Process large datasets efficiently
python3 batch_viral_analysis.py
```

## Demo vs. Production

| Aspect | Demo | Production (100K+) |
|--------|------|--------------------|
| Files | ~15 | 100,000+ |
| Size | ~100 MB | ~5 TB |
| Time | 5 minutes | 24-48 hours |
| Species | 3 | 15+ |
| Purpose | Learning | Research |

## Scaling Strategy

1. **Start Small**: Use this demo (15 files)
2. **Medium Scale**: Collect 1,000 files
3. **Large Scale**: Collect 10,000 files  
4. **Full Scale**: Collect 100,000+ files

Each step validates your pipeline before scaling up!
"""
        
        readme_file = self.output_dir / "README.md"
        with open(readme_file, 'w') as f:
            f.write(readme_content)
        
        logger.info(f"Demo summary saved to: {summary_file}")
        logger.info(f"README saved to: {readme_file}")
    
    def run_demo_collection(self) -> None:
        """Run the complete demo collection process."""
        logger.info("Starting demo mosquito genome data collection...")
        self.collection_stats["start_time"] = datetime.now().isoformat()
        
        try:
            # 1. Get sample genomes
            logger.info("=== PHASE 1: Sample Reference Genomes ===")
            sample_genomes = self.get_sample_ncbi_genomes(max_per_species=5)
            
            # 2. Download sample genomes
            if sample_genomes:
                self.download_sample_genomes(sample_genomes)
            
            # 3. Get sample SRA datasets
            logger.info("=== PHASE 2: Sample SRA Datasets ===")
            sample_sra = self.get_sample_sra_datasets(max_per_species=3)
            
            # 4. Create demo SRA files
            if sample_sra:
                self.create_demo_sra_files(sample_sra)
            
            # 5. Generate summary
            logger.info("=== PHASE 3: Summary Generation ===")
            self.generate_demo_summary()
            
            # 6. Show results
            logger.info("=" * 50)
            logger.info("DEMO COLLECTION COMPLETED")
            logger.info(f"Genomes: {self.collection_stats['genomes_collected']}")
            logger.info(f"SRA Samples: {self.collection_stats['sra_samples_collected']}")
            logger.info(f"Species: {len(self.collection_stats['species_covered'])}")
            logger.info(f"Total Size: {self.collection_stats['total_size_mb']:.1f} MB")
            logger.info(f"Output: {self.output_dir}")
            logger.info("=" * 50)
            
        except Exception as e:
            logger.error(f"Error during demo collection: {e}")
            raise

def show_scaling_roadmap():
    """Show the roadmap for scaling from demo to 100K+ files."""
    print("""
🗺️  SCALING ROADMAP: Demo → 100K+ Files

📊 Phase 1: Demo Collection (Current)
   ├── Files: ~15
   ├── Time: 5 minutes
   ├── Purpose: Learn the process
   └── Next: Run viral analysis on demo data

📈 Phase 2: Small Scale (1K files)
   ├── Files: ~1,000
   ├── Time: 1-2 hours
   ├── Purpose: Validate pipeline
   └── Command: python3 large_scale_data_collector.py (set max to 1000)

📊 Phase 3: Medium Scale (10K files)
   ├── Files: ~10,000
   ├── Time: 6-12 hours
   ├── Purpose: Test batch processing
   └── Command: python3 batch_viral_analysis.py

🚀 Phase 4: Large Scale (100K+ files)
   ├── Files: 100,000+
   ├── Time: 24-48 hours
   ├── Purpose: Full research dataset
   └── Command: Full automated pipeline

💡 Key Success Factors:
   ✅ Start with demo to understand the process
   ✅ Validate each phase before scaling up
   ✅ Monitor storage and processing capacity
   ✅ Document everything for reproducibility

🎯 Your PhD Research Goals:
   • Viral diversity analysis across species
   • Superinfection exclusion pattern detection
   • Comparative virome studies
   • Publication-ready datasets
""")

def main():
    """Main function for demo data collection."""
    print("🦟 Demo Mosquito Genome Data Collection")
    print("=" * 50)
    print("This demo shows how to collect mosquito genome data.")
    print("We'll start with a small, manageable dataset (~15 files).")
    print("")
    
    # Show scaling roadmap
    show_scaling_roadmap()
    
    # Get user confirmation
    proceed = input("\nProceed with demo collection? (y/n): ").strip().lower()
    if proceed != 'y':
        print("Demo cancelled. Run again when ready!")
        return
    
    # Run demo collection
    print("\nStarting demo collection...")
    collector = DemoDataCollector()
    collector.run_demo_collection()
    
    print("\n✅ Demo collection completed!")
    print(f"Check {collector.output_dir}/ for your demo dataset.")
    print("\n🔬 Next Steps:")
    print("1. Run viral analysis: python3 analyze_downloaded_genomes.py")
    print("2. View results in dashboard: Check http://localhost:8503")
    print("3. Scale up: python3 large_scale_data_collector.py")
    print("\n📚 Read LARGE_SCALE_DATA_GUIDE.md for complete instructions!")

if __name__ == "__main__":
    main()