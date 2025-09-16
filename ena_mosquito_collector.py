#!/usr/bin/env python3
"""
European Nucleotide Archive (ENA) Mosquito Genome Collector
Collects additional mosquito genome data from ENA to maximize dataset size.
"""

import requests
import json
import time
from pathlib import Path
from typing import Dict, List, Optional
import logging
from datetime import datetime

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('ena_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ENAMosquitoCollector:
    """Collector for mosquito genome data from European Nucleotide Archive."""
    
    def __init__(self, output_dir: str = "massive_mosquito_genomes/ena_datasets"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # ENA API endpoints
        self.ena_search_url = "https://www.ebi.ac.uk/ena/portal/api/search"
        self.ena_download_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
        
        # Mosquito species to search for
        self.mosquito_species = [
            "Aedes aegypti",
            "Aedes albopictus", 
            "Aedes japonicus",
            "Anopheles gambiae",
            "Anopheles stephensi",
            "Anopheles funestus",
            "Anopheles arabiensis",
            "Anopheles coluzzii",
            "Culex pipiens",
            "Culex quinquefasciatus",
            "Culex tarsalis",
            "Ochlerotatus triseriatus"
        ]
        
        # Search terms for comprehensive collection
        self.search_terms = [
            "mosquito genome",
            "mosquito transcriptome", 
            "mosquito metagenome",
            "mosquito virome",
            "arbovirus surveillance",
            "vector microbiome",
            "insect virome",
            "mosquito RNA-seq",
            "vector pathogen",
            "mosquito metagenomic"
        ]
        
    def search_ena_datasets(self, species: str, search_term: str, max_results: int = 5000) -> List[Dict]:
        """Search ENA for datasets matching species and search term."""
        
        # Simplified query format that works with ENA API
        species_clean = species.replace(" ", "%20")
        search_clean = search_term.replace(" ", "%20")
        query = f"scientific_name={species_clean}"
        
        params = {
            'result': 'read_run',
            'query': query,
            'format': 'json',
            'limit': max_results,
            'fields': 'run_accession,sample_accession,experiment_accession,study_accession,scientific_name,library_strategy,library_source,instrument_platform,base_count,read_count'
        }
        
        try:
            logger.info(f"Searching ENA for '{search_term}' + '{species}'...")
            response = requests.get(self.ena_search_url, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                logger.info(f"Found {len(data)} ENA datasets for '{search_term}' + '{species}'")
                return data
            else:
                logger.warning(f"ENA search failed for {species} + {search_term}: {response.status_code}")
                return []
                
        except Exception as e:
            logger.error(f"Error searching ENA for {species} + {search_term}: {e}")
            return []
    
    def download_ena_metadata(self, datasets: List[Dict], species: str, search_term: str) -> None:
        """Download and save ENA dataset metadata."""
        
        if not datasets:
            return
            
        # Create species directory
        species_dir = self.output_dir / species.replace(" ", "_")
        species_dir.mkdir(exist_ok=True)
        
        # Save metadata
        metadata_file = species_dir / f"{search_term.replace(' ', '_')}_metadata.json"
        
        with open(metadata_file, 'w') as f:
            json.dump(datasets, f, indent=2)
            
        logger.info(f"Saved {len(datasets)} ENA metadata records to {metadata_file}")
        
        # Create download URLs file for large datasets
        if len(datasets) > 100:
            urls_file = species_dir / f"{search_term.replace(' ', '_')}_download_urls.txt"
            
            with open(urls_file, 'w') as f:
                for dataset in datasets:
                    if 'fastq_ftp' in dataset and dataset['fastq_ftp']:
                        urls = dataset['fastq_ftp'].split(';')
                        for url in urls:
                            if url.strip():
                                f.write(f"ftp://{url.strip()}\n")
                                
            logger.info(f"Created download URLs file: {urls_file}")
    
    def collect_ena_data(self, max_datasets_per_search: int = 5000) -> Dict:
        """Collect comprehensive ENA mosquito datasets."""
        
        logger.info("Starting ENA mosquito genome collection...")
        logger.info(f"Target: {len(self.mosquito_species)} species x {len(self.search_terms)} search terms")
        logger.info(f"Maximum datasets per search: {max_datasets_per_search:,}")
        
        collection_summary = {
            'collection_start': datetime.now().isoformat(),
            'total_datasets_found': 0,
            'species_results': {},
            'search_term_results': {}
        }
        
        total_found = 0
        
        for species in self.mosquito_species:
            logger.info(f"\n=== Processing {species} ===")
            species_total = 0
            
            for search_term in self.search_terms:
                # Search ENA
                datasets = self.search_ena_datasets(species, search_term, max_datasets_per_search)
                
                if datasets:
                    # Download metadata
                    self.download_ena_metadata(datasets, species, search_term)
                    
                    species_total += len(datasets)
                    total_found += len(datasets)
                    
                    # Update summary
                    if search_term not in collection_summary['search_term_results']:
                        collection_summary['search_term_results'][search_term] = 0
                    collection_summary['search_term_results'][search_term] += len(datasets)
                
                # Rate limiting
                time.sleep(2)
            
            collection_summary['species_results'][species] = species_total
            logger.info(f"Total datasets found for {species}: {species_total:,}")
        
        # Final summary
        collection_summary['total_datasets_found'] = total_found
        collection_summary['collection_end'] = datetime.now().isoformat()
        
        # Save summary
        summary_file = self.output_dir / "ena_collection_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(collection_summary, f, indent=2)
            
        logger.info(f"\n=== ENA COLLECTION COMPLETE ===")
        logger.info(f"Total ENA datasets found: {total_found:,}")
        logger.info(f"Summary saved to: {summary_file}")
        
        return collection_summary

def main():
    """Main execution function."""
    
    logger.info("=== ENA Mosquito Genome Collector ===")
    logger.info("Collecting maximum mosquito datasets from European Nucleotide Archive")
    
    # Initialize collector
    collector = ENAMosquitoCollector()
    
    # Start collection
    summary = collector.collect_ena_data(max_datasets_per_search=10000)
    
    # Print final results
    print("\n" + "="*50)
    print("ENA COLLECTION SUMMARY")
    print("="*50)
    print(f"Total datasets found: {summary['total_datasets_found']:,}")
    print(f"Species processed: {len(summary['species_results'])}")
    print(f"Search terms used: {len(summary['search_term_results'])}")
    
    print("\nTop species by dataset count:")
    sorted_species = sorted(summary['species_results'].items(), key=lambda x: x[1], reverse=True)
    for species, count in sorted_species[:5]:
        print(f"  {species}: {count:,} datasets")
    
    print("\nTop search terms by dataset count:")
    sorted_terms = sorted(summary['search_term_results'].items(), key=lambda x: x[1], reverse=True)
    for term, count in sorted_terms[:5]:
        print(f"  {term}: {count:,} datasets")
    
    print(f"\nCollection completed! Check 'massive_mosquito_genomes/ena_datasets/' for results.")
    print(f"Summary file: massive_mosquito_genomes/ena_datasets/ena_collection_summary.json")

if __name__ == "__main__":
    main()