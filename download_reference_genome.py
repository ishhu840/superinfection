#!/usr/bin/env python3
"""
High-Quality Reference Genome Downloader

Downloads the latest high-quality Anopheles gambiae reference genome
(GCA_943734735.2) from NCBI for use in viral detection pipeline.
"""

import os
import sys
import requests
import gzip
import shutil
from pathlib import Path
import logging
from typing import Optional, Dict

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ReferenceGenomeDownloader:
    """
    Downloads high-quality reference genomes for mosquito viral analysis.
    """
    
    def __init__(self, output_dir: str = "reference_genomes"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # High-quality reference genomes
        self.reference_genomes = {
            'anopheles_gambiae_ifakara': {
                'accession': 'GCA_943734735.2',
                'name': 'Anopheles_gambiae_Ifakara_strain',
                'description': 'High-quality chromosomal reference genome (264 Mb, 33 gaps)',
                'url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/943/734/735/GCA_943734735.2_idAnoGambNW_F1_1/GCA_943734735.2_idAnoGambNW_F1_1_genomic.fna.gz',
                'species': 'Anopheles gambiae',
                'strain': 'Ifakara',
                'quality': 'chromosomal',
                'size_mb': 264,
                'gaps': 33
            },
            'anopheles_gambiae_pest': {
                'accession': 'GCA_000005575.2',
                'name': 'Anopheles_gambiae_PEST_strain',
                'description': 'Classic PEST reference genome (AgamP4)',
                'url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/575/GCA_000005575.2_AgamP4/GCA_000005575.2_AgamP4_genomic.fna.gz',
                'species': 'Anopheles gambiae',
                'strain': 'PEST',
                'quality': 'standard',
                'size_mb': 225,
                'gaps': 6000
            }
        }
    
    def download_file(self, url: str, output_path: Path, chunk_size: int = 8192) -> bool:
        """
        Download a file from URL with progress tracking.
        """
        try:
            logger.info(f"Downloading {url}")
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        
                        if total_size > 0:
                            progress = (downloaded / total_size) * 100
                            print(f"\rProgress: {progress:.1f}%", end='', flush=True)
            
            print()  # New line after progress
            logger.info(f"Downloaded: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to download {url}: {e}")
            return False
    
    def extract_genome(self, gz_file: Path, output_file: Path) -> bool:
        """
        Extract gzipped genome file.
        """
        try:
            logger.info(f"Extracting {gz_file}")
            with gzip.open(gz_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            logger.info(f"Extracted: {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to extract {gz_file}: {e}")
            return False
    
    def download_reference_genome(self, genome_key: str) -> Optional[Path]:
        """
        Download a specific reference genome.
        """
        if genome_key not in self.reference_genomes:
            logger.error(f"Unknown genome: {genome_key}")
            return None
        
        genome_info = self.reference_genomes[genome_key]
        
        # Create species directory
        species_dir = self.output_dir / genome_info['species'].replace(' ', '_')
        species_dir.mkdir(exist_ok=True)
        
        # File paths
        gz_file = species_dir / f"{genome_info['name']}.fna.gz"
        fasta_file = species_dir / f"{genome_info['name']}.fna"
        
        # Skip if already exists
        if fasta_file.exists():
            logger.info(f"Reference genome already exists: {fasta_file}")
            return fasta_file
        
        # Download compressed file
        if not self.download_file(genome_info['url'], gz_file):
            return None
        
        # Extract genome
        if not self.extract_genome(gz_file, fasta_file):
            return None
        
        # Clean up compressed file
        gz_file.unlink()
        
        # Log genome info
        logger.info(f"Successfully downloaded {genome_info['description']}")
        logger.info(f"Size: {genome_info['size_mb']} Mb, Gaps: {genome_info['gaps']}")
        logger.info(f"Location: {fasta_file}")
        
        return fasta_file
    
    def download_all_references(self) -> Dict[str, Path]:
        """
        Download all available reference genomes.
        """
        results = {}
        
        for genome_key in self.reference_genomes:
            logger.info(f"\n=== Downloading {genome_key} ===")
            result = self.download_reference_genome(genome_key)
            if result:
                results[genome_key] = result
        
        return results
    
    def get_best_reference(self) -> Optional[Path]:
        """
        Download the highest quality reference genome (Ifakara strain).
        """
        logger.info("Downloading highest quality Anopheles gambiae reference genome...")
        return self.download_reference_genome('anopheles_gambiae_ifakara')

def main():
    """
    Main function to download reference genomes.
    """
    print("🧬 High-Quality Reference Genome Downloader")
    print("=" * 50)
    
    downloader = ReferenceGenomeDownloader()
    
    # Download the best reference genome
    best_ref = downloader.get_best_reference()
    
    if best_ref:
        print(f"\n✅ Successfully downloaded reference genome: {best_ref}")
        print("\nThis high-quality reference can be used for:")
        print("- Improved viral detection sensitivity")
        print("- Better sequence alignment quality")
        print("- Species-specific viral analysis")
        print("- Contamination filtering")
        
        # Integration suggestions
        print("\n🔧 Integration with existing pipeline:")
        print("1. Update blast_pipeline.py to use this reference")
        print("2. Add host genome filtering step")
        print("3. Improve viral hit validation")
        
    else:
        print("❌ Failed to download reference genome")
        sys.exit(1)

if __name__ == "__main__":
    main()