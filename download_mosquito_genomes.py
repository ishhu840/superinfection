#!/usr/bin/env python3
"""
Mosquito Genome Download Script for PhD Research

This script downloads real mosquito genomes from NCBI for superinfection exclusion research.
Focuses on key vector species: Aedes aegypti, Anopheles gambiae, Culex quinquefasciatus.
"""

import os
import sys
from pathlib import Path

# Add src to path
src_dir = Path(__file__).parent / "src"
sys.path.insert(0, str(src_dir))

from data_collection.genome_downloader import MosquitoGenomeDownloader

def main():
    """
    Download mosquito genomes for superinfection exclusion research.
    """
    print("=" * 60)
    print("MOSQUITO GENOME DOWNLOAD FOR PhD RESEARCH")
    print("Superinfection Exclusion Analysis")
    print("=" * 60)
    
    # Get email for NCBI access
    email = input("Enter your email for NCBI access: ").strip()
    if not email or '@' not in email:
        print("Error: Valid email required for NCBI access")
        return
    
    # Initialize downloader
    output_dir = "research_genomes"
    downloader = MosquitoGenomeDownloader(
        email=email,
        output_dir=output_dir
    )
    
    print(f"\nDownloading genomes to: {output_dir}/")
    print("This may take several minutes depending on your connection...\n")
    
    # Key species for superinfection exclusion research
    target_species = [
        'aedes_aegypti',        # Primary DENV/ZIKV vector
        'anopheles_gambiae',    # Malaria vector
        'culex_quinquefasciatus' # WNV vector
    ]
    
    results = {}
    
    for species in target_species:
        print(f"\n{'='*50}")
        print(f"Downloading {species.replace('_', ' ').title()} genomes...")
        print(f"{'='*50}")
        
        try:
            species_results = downloader.download_species_genomes(
                species, 
                max_assemblies=2  # Download 2 assemblies per species
            )
            results[species] = species_results
            print(f"✓ Successfully downloaded {len(species_results)} assemblies for {species}")
            
        except Exception as e:
            print(f"✗ Error downloading {species}: {str(e)}")
            results[species] = []
    
    # Summary
    print("\n" + "="*60)
    print("DOWNLOAD SUMMARY")
    print("="*60)
    
    total_files = 0
    summary = downloader.get_downloaded_genomes_summary()
    
    for species, files in summary.items():
        print(f"\n{species.replace('_', ' ').title()}:")
        print(f"  - {len(files)} files downloaded")
        total_files += len(files)
        
        # Show first few files
        for i, file in enumerate(files[:3]):
            file_size = ""
            file_path = Path(output_dir) / species / file
            if file_path.exists():
                size_mb = file_path.stat().st_size / (1024 * 1024)
                file_size = f" ({size_mb:.1f} MB)"
            print(f"    {i+1}. {file}{file_size}")
        
        if len(files) > 3:
            print(f"    ... and {len(files) - 3} more files")
    
    print(f"\nTotal files downloaded: {total_files}")
    print(f"Output directory: {Path(output_dir).absolute()}")
    
    if total_files > 0:
        print("\n✓ Genome download completed successfully!")
        print("\nNext steps for your PhD research:")
        print("1. Run viral detection pipeline on downloaded genomes")
        print("2. Analyze superinfection exclusion patterns")
        print("3. Use the comprehensive dashboard for visualization")
        
        # Suggest running the main pipeline
        print(f"\nTo analyze these genomes, run:")
        print(f"python3 src/main_pipeline.py --genome_dir {output_dir}")
    else:
        print("\n✗ No genomes were downloaded. Please check your internet connection and try again.")

if __name__ == "__main__":
    main()