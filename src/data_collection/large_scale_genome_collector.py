#!/usr/bin/env python3
"""
Large-Scale Mosquito Genome Data Collector

This module implements comprehensive data collection strategies for obtaining
extensive, high-quality mosquito genome datasets from NCBI and other databases
suitable for robust scientific analysis and publication.

Authors: Mosquito Viral Analysis Pipeline Team
Date: 2025
"""

import os
import sys
import time
import requests
import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from pathlib import Path
import json
from datetime import datetime, timedelta
import hashlib
from typing import List, Dict, Optional, Tuple
import xml.etree.ElementTree as ET
from urllib.parse import urlencode
import gzip
import shutil
from dataclasses import dataclass

@dataclass
class GenomeMetadata:
    """Structured metadata for genome assemblies."""
    accession: str
    species: str
    strain: str
    assembly_level: str
    genome_size: int
    contig_count: int
    n50: int
    collection_date: str
    geographic_location: str
    bioproject: str
    biosample: str
    submitter: str
    release_date: str
    last_update: str
    assembly_stats: Dict
    quality_score: float

class LargeScaleGenomeCollector:
    """
    Comprehensive genome data collector for large-scale mosquito viral analysis.
    
    This class implements strategies to collect thousands of high-quality
    mosquito genome assemblies from multiple databases, with proper metadata
    tracking and quality control suitable for scientific publication.
    """
    
    # Major mosquito species for comprehensive analysis
    TARGET_SPECIES = {
        'Aedes aegypti': 'aedes_aegypti',
        'Aedes albopictus': 'aedes_albopictus', 
        'Anopheles gambiae': 'anopheles_gambiae',
        'Anopheles stephensi': 'anopheles_stephensi',
        'Anopheles funestus': 'anopheles_funestus',
        'Culex quinquefasciatus': 'culex_quinquefasciatus',
        'Culex pipiens': 'culex_pipiens',
        'Anopheles coluzzii': 'anopheles_coluzzii',
        'Anopheles arabiensis': 'anopheles_arabiensis',
        'Aedes japonicus': 'aedes_japonicus'
    }
    
    # Quality thresholds for publication-grade data
    QUALITY_THRESHOLDS = {
        'min_genome_size': 100000,  # 100kb minimum
        'max_genome_size': 2000000000,  # 2Gb maximum
        'max_n_content': 0.05,  # 5% N content maximum
        'min_n50': 1000,  # 1kb minimum N50
        'max_contigs': 50000,  # Maximum number of contigs
        'min_coverage': 10,  # Minimum sequencing coverage
        'assembly_levels': ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
    }
    
    def __init__(self, email: str, api_key: Optional[str] = None, 
                 output_dir: str = "large_scale_genomes", max_workers: int = 8):
        """
        Initialize the large-scale genome collector.
        
        Args:
            email: Email for NCBI API access
            api_key: NCBI API key for increased rate limits
            output_dir: Directory for storing downloaded genomes
            max_workers: Number of parallel download workers
        """
        self.email = email
        self.api_key = api_key
        self.output_dir = Path(output_dir)
        self.max_workers = max_workers
        
        # Set up NCBI Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            
        # Create output directories
        self.output_dir.mkdir(exist_ok=True)
        self.metadata_dir = self.output_dir / "metadata"
        self.metadata_dir.mkdir(exist_ok=True)
        
        # Set up logging
        self.setup_logging()
        
        # Initialize tracking
        self.download_stats = {
            'total_attempted': 0,
            'successful_downloads': 0,
            'failed_downloads': 0,
            'quality_filtered': 0,
            'duplicate_filtered': 0
        }
        
        self.collected_genomes = []
        self.failed_downloads = []
        
    def setup_logging(self):
        """Set up comprehensive logging for data collection."""
        log_file = self.output_dir / f"collection_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        
        self.logger = logging.getLogger(__name__)
        
    def search_comprehensive_assemblies(self, species: str, max_results: int = 10000) -> List[str]:
        """
        Search for all available genome assemblies for a species.
        
        Args:
            species: Species name (e.g., 'Aedes aegypti')
            max_results: Maximum number of results to retrieve
            
        Returns:
            List of assembly accession numbers
        """
        self.logger.info(f"Searching for {species} genome assemblies...")
        
        # Comprehensive search terms
        search_terms = [
            f'"{species}"[Organism] AND ("latest refseq"[filter] OR "latest genbank"[filter])',
            f'"{species}"[Organism] AND "complete genome"[All Fields]',
            f'"{species}"[Organism] AND "chromosome"[All Fields]',
            f'"{species}"[Organism] AND "scaffold"[All Fields]',
            f'"{species}"[Organism] AND "contig"[All Fields]'
        ]
        
        all_assemblies = set()
        
        for search_term in search_terms:
            try:
                # Search assembly database
                handle = Entrez.esearch(
                    db="assembly",
                    term=search_term,
                    retmax=max_results,
                    sort="relevance"
                )
                search_results = Entrez.read(handle)
                handle.close()
                
                assembly_ids = search_results['IdList']
                self.logger.info(f"Found {len(assembly_ids)} assemblies with search term: {search_term[:50]}...")
                
                # Get assembly accessions
                if assembly_ids:
                    handle = Entrez.esummary(db="assembly", id=",".join(assembly_ids))
                    summaries = Entrez.read(handle)
                    handle.close()
                    
                    for summary in summaries['DocumentSummarySet']['DocumentSummary']:
                        accession = summary.get('AssemblyAccession', '')
                        if accession:
                            all_assemblies.add(accession)
                            
                time.sleep(0.1)  # Rate limiting
                
            except Exception as e:
                self.logger.error(f"Error searching with term '{search_term}': {e}")
                continue
        
        assembly_list = list(all_assemblies)
        self.logger.info(f"Total unique assemblies found for {species}: {len(assembly_list)}")
        
        return assembly_list
    
    def get_detailed_assembly_metadata(self, accession: str) -> Optional[GenomeMetadata]:
        """
        Retrieve comprehensive metadata for an assembly.
        
        Args:
            accession: Assembly accession number
            
        Returns:
            GenomeMetadata object or None if failed
        """
        try:
            # Get assembly summary
            handle = Entrez.esummary(db="assembly", id=accession)
            summary = Entrez.read(handle)
            handle.close()
            
            if not summary['DocumentSummarySet']['DocumentSummary']:
                return None
                
            assembly_data = summary['DocumentSummarySet']['DocumentSummary'][0]
            
            # Extract metadata
            metadata = GenomeMetadata(
                accession=accession,
                species=assembly_data.get('SpeciesName', ''),
                strain=assembly_data.get('Biosource', {}).get('InfraspecificName', ''),
                assembly_level=assembly_data.get('AssemblyStatus', ''),
                genome_size=int(assembly_data.get('TotalLength', 0)),
                contig_count=int(assembly_data.get('ContigN50', 0)),
                n50=int(assembly_data.get('ScaffoldN50', 0)),
                collection_date=assembly_data.get('AsmReleaseDate_GenBank', ''),
                geographic_location=assembly_data.get('Biosource', {}).get('GeographicLocation', ''),
                bioproject=assembly_data.get('GB_BioProjects', [{}])[0].get('BioprojectAccn', '') if assembly_data.get('GB_BioProjects') else '',
                biosample=assembly_data.get('BiosampleAccn', ''),
                submitter=assembly_data.get('SubmitterOrganization', ''),
                release_date=assembly_data.get('AsmReleaseDate_GenBank', ''),
                last_update=assembly_data.get('LastUpdateDate', ''),
                assembly_stats=assembly_data.get('Meta', {}),
                quality_score=0.0  # Will be calculated
            )
            
            # Calculate quality score
            metadata.quality_score = self.calculate_quality_score(metadata)
            
            return metadata
            
        except Exception as e:
            self.logger.error(f"Error getting metadata for {accession}: {e}")
            return None
    
    def calculate_quality_score(self, metadata: GenomeMetadata) -> float:
        """
        Calculate a quality score for genome assembly.
        
        Args:
            metadata: GenomeMetadata object
            
        Returns:
            Quality score (0-100)
        """
        score = 0.0
        
        # Assembly level scoring
        level_scores = {
            'Complete Genome': 40,
            'Chromosome': 30,
            'Scaffold': 20,
            'Contig': 10
        }
        score += level_scores.get(metadata.assembly_level, 0)
        
        # Genome size scoring (optimal range)
        if 200000000 <= metadata.genome_size <= 800000000:  # 200Mb - 800Mb
            score += 20
        elif 100000000 <= metadata.genome_size <= 1000000000:  # 100Mb - 1Gb
            score += 15
        elif metadata.genome_size > 50000000:  # > 50Mb
            score += 10
        
        # N50 scoring
        if metadata.n50 > 1000000:  # > 1Mb
            score += 20
        elif metadata.n50 > 100000:  # > 100kb
            score += 15
        elif metadata.n50 > 10000:  # > 10kb
            score += 10
        
        # Contig count scoring (fewer is better)
        if metadata.contig_count < 100:
            score += 10
        elif metadata.contig_count < 1000:
            score += 5
        
        # Recent assembly bonus
        try:
            release_year = int(metadata.release_date[:4]) if metadata.release_date else 2000
            current_year = datetime.now().year
            if current_year - release_year <= 2:
                score += 10
            elif current_year - release_year <= 5:
                score += 5
        except:
            pass
        
        return min(score, 100.0)
    
    def filter_high_quality_assemblies(self, assemblies_metadata: List[GenomeMetadata], 
                                     min_quality_score: float = 50.0) -> List[GenomeMetadata]:
        """
        Filter assemblies based on quality criteria.
        
        Args:
            assemblies_metadata: List of GenomeMetadata objects
            min_quality_score: Minimum quality score threshold
            
        Returns:
            Filtered list of high-quality assemblies
        """
        self.logger.info(f"Filtering {len(assemblies_metadata)} assemblies for quality...")
        
        high_quality = []
        
        for metadata in assemblies_metadata:
            # Quality score filter
            if metadata.quality_score < min_quality_score:
                continue
                
            # Size filters
            if not (self.QUALITY_THRESHOLDS['min_genome_size'] <= 
                   metadata.genome_size <= self.QUALITY_THRESHOLDS['max_genome_size']):
                continue
                
            # Assembly level filter
            if metadata.assembly_level not in self.QUALITY_THRESHOLDS['assembly_levels']:
                continue
                
            # N50 filter
            if metadata.n50 < self.QUALITY_THRESHOLDS['min_n50']:
                continue
                
            # Contig count filter
            if metadata.contig_count > self.QUALITY_THRESHOLDS['max_contigs']:
                continue
                
            high_quality.append(metadata)
        
        self.logger.info(f"Selected {len(high_quality)} high-quality assemblies")
        return high_quality
    
    def download_genome_assembly(self, metadata: GenomeMetadata) -> Tuple[bool, str]:
        """
        Download a genome assembly with comprehensive error handling.
        
        Args:
            metadata: GenomeMetadata object
            
        Returns:
            Tuple of (success, file_path or error_message)
        """
        species_dir = self.output_dir / "genomes" / self.TARGET_SPECIES.get(metadata.species, metadata.species.replace(' ', '_').lower())
        species_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = species_dir / f"{metadata.accession}.fna.gz"
        
        # Skip if already downloaded
        if output_file.exists():
            self.logger.info(f"Assembly {metadata.accession} already exists, skipping")
            return True, str(output_file)
        
        try:
            # Get FTP URLs for assembly
            handle = Entrez.esummary(db="assembly", id=metadata.accession)
            summary = Entrez.read(handle)
            handle.close()
            
            if not summary['DocumentSummarySet']['DocumentSummary']:
                return False, "No assembly summary found"
                
            assembly_data = summary['DocumentSummarySet']['DocumentSummary'][0]
            ftp_path = assembly_data.get('FtpPath_GenBank', '') or assembly_data.get('FtpPath_RefSeq', '')
            
            if not ftp_path:
                return False, "No FTP path available"
            
            # Construct download URL
            assembly_name = ftp_path.split('/')[-1]
            download_url = f"{ftp_path}/{assembly_name}_genomic.fna.gz"
            
            # Download with retries
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    self.logger.info(f"Downloading {metadata.accession} (attempt {attempt + 1}/{max_retries})")
                    
                    response = requests.get(download_url, stream=True, timeout=300)
                    response.raise_for_status()
                    
                    with open(output_file, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)
                    
                    # Verify download
                    if output_file.stat().st_size > 1000:  # At least 1KB
                        self.logger.info(f"Successfully downloaded {metadata.accession}")
                        return True, str(output_file)
                    else:
                        output_file.unlink(missing_ok=True)
                        return False, "Downloaded file too small"
                        
                except Exception as e:
                    self.logger.warning(f"Download attempt {attempt + 1} failed for {metadata.accession}: {e}")
                    if attempt < max_retries - 1:
                        time.sleep(2 ** attempt)  # Exponential backoff
                    continue
            
            return False, f"All {max_retries} download attempts failed"
            
        except Exception as e:
            return False, f"Download error: {e}"
    
    def collect_comprehensive_dataset(self, target_samples_per_species: int = 500,
                                    min_quality_score: float = 60.0) -> Dict:
        """
        Collect comprehensive dataset for all target mosquito species.
        
        Args:
            target_samples_per_species: Target number of samples per species
            min_quality_score: Minimum quality score for inclusion
            
        Returns:
            Collection summary dictionary
        """
        self.logger.info(f"Starting comprehensive data collection for {len(self.TARGET_SPECIES)} species")
        self.logger.info(f"Target: {target_samples_per_species} high-quality genomes per species")
        
        collection_summary = {
            'species_results': {},
            'total_collected': 0,
            'total_attempted': 0,
            'collection_start': datetime.now().isoformat(),
            'parameters': {
                'target_samples_per_species': target_samples_per_species,
                'min_quality_score': min_quality_score
            }
        }
        
        for species_name, species_code in self.TARGET_SPECIES.items():
            self.logger.info(f"\n=== Processing {species_name} ===")
            
            # Search for assemblies
            assemblies = self.search_comprehensive_assemblies(species_name, max_results=5000)
            
            if not assemblies:
                self.logger.warning(f"No assemblies found for {species_name}")
                collection_summary['species_results'][species_name] = {
                    'assemblies_found': 0,
                    'high_quality': 0,
                    'downloaded': 0,
                    'failed': 0
                }
                continue
            
            # Get metadata for all assemblies
            self.logger.info(f"Retrieving metadata for {len(assemblies)} assemblies...")
            assemblies_metadata = []
            
            for i, accession in enumerate(assemblies):
                if i % 100 == 0:
                    self.logger.info(f"Processed {i}/{len(assemblies)} metadata records")
                    
                metadata = self.get_detailed_assembly_metadata(accession)
                if metadata:
                    assemblies_metadata.append(metadata)
                    
                time.sleep(0.1)  # Rate limiting
            
            # Filter for high quality
            high_quality = self.filter_high_quality_assemblies(
                assemblies_metadata, min_quality_score
            )
            
            # Sort by quality score and select top samples
            high_quality.sort(key=lambda x: x.quality_score, reverse=True)
            selected_assemblies = high_quality[:target_samples_per_species]
            
            self.logger.info(f"Selected {len(selected_assemblies)} assemblies for download")
            
            # Download assemblies in parallel
            downloaded = 0
            failed = 0
            
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_metadata = {
                    executor.submit(self.download_genome_assembly, metadata): metadata 
                    for metadata in selected_assemblies
                }
                
                for future in as_completed(future_to_metadata):
                    metadata = future_to_metadata[future]
                    try:
                        success, result = future.result()
                        if success:
                            downloaded += 1
                            self.collected_genomes.append(metadata)
                        else:
                            failed += 1
                            self.failed_downloads.append((metadata.accession, result))
                            
                    except Exception as e:
                        failed += 1
                        self.logger.error(f"Download failed for {metadata.accession}: {e}")
                        self.failed_downloads.append((metadata.accession, str(e)))
            
            # Save species metadata
            species_metadata_file = self.metadata_dir / f"{species_code}_metadata.json"
            with open(species_metadata_file, 'w') as f:
                metadata_dict = [
                    {
                        'accession': m.accession,
                        'species': m.species,
                        'strain': m.strain,
                        'assembly_level': m.assembly_level,
                        'genome_size': m.genome_size,
                        'quality_score': m.quality_score,
                        'collection_date': m.collection_date,
                        'geographic_location': m.geographic_location
                    } for m in selected_assemblies
                ]
                json.dump(metadata_dict, f, indent=2)
            
            # Update summary
            collection_summary['species_results'][species_name] = {
                'assemblies_found': len(assemblies),
                'high_quality': len(high_quality),
                'selected': len(selected_assemblies),
                'downloaded': downloaded,
                'failed': failed,
                'success_rate': downloaded / len(selected_assemblies) if selected_assemblies else 0
            }
            
            collection_summary['total_collected'] += downloaded
            collection_summary['total_attempted'] += len(selected_assemblies)
            
            self.logger.info(f"Species {species_name} complete: {downloaded}/{len(selected_assemblies)} downloaded")
        
        # Save comprehensive summary
        collection_summary['collection_end'] = datetime.now().isoformat()
        collection_summary['total_species'] = len(self.TARGET_SPECIES)
        collection_summary['overall_success_rate'] = (
            collection_summary['total_collected'] / collection_summary['total_attempted']
            if collection_summary['total_attempted'] > 0 else 0
        )
        
        summary_file = self.metadata_dir / "collection_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(collection_summary, f, indent=2)
        
        self.logger.info(f"\n=== COLLECTION COMPLETE ===")
        self.logger.info(f"Total genomes collected: {collection_summary['total_collected']}")
        self.logger.info(f"Overall success rate: {collection_summary['overall_success_rate']:.2%}")
        
        return collection_summary
    
    def generate_data_quality_report(self) -> Dict:
        """
        Generate comprehensive data quality report for publication.
        
        Returns:
            Quality report dictionary
        """
        self.logger.info("Generating data quality report...")
        
        if not self.collected_genomes:
            return {'error': 'No genomes collected yet'}
        
        # Calculate quality metrics
        quality_scores = [g.quality_score for g in self.collected_genomes]
        genome_sizes = [g.genome_size for g in self.collected_genomes]
        n50_values = [g.n50 for g in self.collected_genomes]
        
        # Species distribution
        species_counts = {}
        for genome in self.collected_genomes:
            species_counts[genome.species] = species_counts.get(genome.species, 0) + 1
        
        # Assembly level distribution
        assembly_levels = {}
        for genome in self.collected_genomes:
            assembly_levels[genome.assembly_level] = assembly_levels.get(genome.assembly_level, 0) + 1
        
        # Geographic distribution
        geographic_locations = {}
        for genome in self.collected_genomes:
            if genome.geographic_location:
                geographic_locations[genome.geographic_location] = geographic_locations.get(genome.geographic_location, 0) + 1
        
        quality_report = {
            'dataset_overview': {
                'total_genomes': len(self.collected_genomes),
                'species_count': len(species_counts),
                'geographic_locations': len(geographic_locations)
            },
            'quality_metrics': {
                'mean_quality_score': np.mean(quality_scores),
                'median_quality_score': np.median(quality_scores),
                'min_quality_score': np.min(quality_scores),
                'max_quality_score': np.max(quality_scores)
            },
            'genome_size_stats': {
                'mean_size': np.mean(genome_sizes),
                'median_size': np.median(genome_sizes),
                'min_size': np.min(genome_sizes),
                'max_size': np.max(genome_sizes)
            },
            'n50_stats': {
                'mean_n50': np.mean(n50_values),
                'median_n50': np.median(n50_values),
                'min_n50': np.min(n50_values),
                'max_n50': np.max(n50_values)
            },
            'species_distribution': species_counts,
            'assembly_level_distribution': assembly_levels,
            'geographic_distribution': dict(list(geographic_locations.items())[:20])  # Top 20 locations
        }
        
        # Save quality report
        report_file = self.metadata_dir / "data_quality_report.json"
        with open(report_file, 'w') as f:
            json.dump(quality_report, f, indent=2)
        
        return quality_report

def main():
    """
    Example usage for large-scale data collection.
    """
    print("Large-Scale Mosquito Genome Data Collector")
    print("===========================================")
    
    # Initialize collector
    collector = LargeScaleGenomeCollector(
        email="your.email@institution.edu",
        api_key="your_ncbi_api_key",  # Optional but recommended
        output_dir="publication_dataset",
        max_workers=6
    )
    
    print("\nRecommended data collection strategy for publication:")
    print("- Target: 500+ high-quality genomes per major species")
    print("- Quality threshold: 60+ (out of 100)")
    print("- Total expected: 5,000+ genomes across 10 species")
    print("- Estimated time: 6-12 hours depending on network")
    print("- Storage required: ~500GB - 2TB")
    
    # Start comprehensive collection
    print("\nStarting comprehensive data collection...")
    summary = collector.collect_comprehensive_dataset(
        target_samples_per_species=500,
        min_quality_score=60.0
    )
    
    # Generate quality report
    quality_report = collector.generate_data_quality_report()
    
    print("\n=== FINAL SUMMARY ===")
    print(f"Total genomes collected: {summary['total_collected']}")
    print(f"Success rate: {summary['overall_success_rate']:.2%}")
    print(f"Mean quality score: {quality_report['quality_metrics']['mean_quality_score']:.1f}")
    print(f"Species represented: {quality_report['dataset_overview']['species_count']}")
    
    print("\nDataset is ready for publication-quality viral co-occurrence analysis!")

if __name__ == "__main__":
    main()