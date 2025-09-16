#!/usr/bin/env python3
"""
Mosquito Genome Data Collection Module

This module provides functionality to download mosquito genome assemblies
from NCBI databases using the Entrez API and BioPython.

Based on research findings:
- NCBI has chromosome-scale assemblies for major mosquito species
- Aedes aegypti (AaegL5), Culex quinquefasciatus, Anopheles gambiae available
- RefSeq provides high-quality reference genomes
"""

import os
import urllib.request
import gzip
import logging
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import time
import json
import xml.etree.ElementTree as ET
from urllib.parse import urljoin

try:
    from Bio import SeqIO, Entrez
except ImportError:
    SeqIO = None
    Entrez = None
    logging.warning("BioPython not available. Some features may be limited.")

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Comprehensive mosquito genome database with real NCBI accessions
MOSQUITO_GENOMES = {
    "aedes_aegypti": {
        "species_name": "Aedes aegypti",
        "common_name": "Yellow fever mosquito",
        "ncbi_assembly_accession": "GCF_002204515.2",
        "refseq_accession": "GCF_002204515.2",
        "genbank_accession": "GCA_002204515.2",
        "assembly_name": "AaegL5.0",
        "chromosomes": ["NC_035107.1", "NC_035108.1", "NC_035109.1"],
        "genome_size_mb": 1282,
        "vector_competence": ["DENV", "ZIKV", "CHIKV", "YFV"]
    },
    "anopheles_gambiae": {
        "species_name": "Anopheles gambiae",
        "common_name": "African malaria mosquito",
        "ncbi_assembly_accession": "GCF_000005575.2",
        "refseq_accession": "GCF_000005575.2",
        "genbank_accession": "GCA_000005575.2",
        "assembly_name": "AgamP3",
        "chromosomes": ["NC_004818.2", "NC_004819.2", "NC_004820.2"],
        "genome_size_mb": 265,
        "vector_competence": ["Plasmodium", "O'nyong-nyong virus"]
    },
    "culex_quinquefasciatus": {
        "species_name": "Culex quinquefasciatus",
        "common_name": "Southern house mosquito",
        "ncbi_assembly_accession": "GCF_015732765.1",
        "refseq_accession": "GCF_015732765.1",
        "genbank_accession": "GCA_015732765.1",
        "assembly_name": "VPISU_Cqui_1.0_pri",
        "chromosomes": ["NC_052092.1", "NC_052093.1", "NC_052094.1"],
        "genome_size_mb": 579,
        "vector_competence": ["WNV", "SLEV", "EEEV", "Rift Valley fever virus"]
    },
    "aedes_albopictus": {
        "species_name": "Aedes albopictus",
        "common_name": "Asian tiger mosquito",
        "ncbi_assembly_accession": "GCF_006496715.1",
        "refseq_accession": "GCF_006496715.1",
        "genbank_accession": "GCA_006496715.1",
        "assembly_name": "Aalbo_primary.1",
        "chromosomes": ["NC_044393.1", "NC_044394.1", "NC_044395.1"],
        "genome_size_mb": 1967,
        "vector_competence": ["DENV", "CHIKV", "ZIKV"]
    },
    "anopheles_stephensi": {
        "species_name": "Anopheles stephensi",
        "common_name": "Asian malaria mosquito",
        "ncbi_assembly_accession": "GCF_013141755.1",
        "refseq_accession": "GCF_013141755.1",
        "genbank_accession": "GCA_013141755.1",
        "assembly_name": "UCI_ANSTEP_V1.0",
        "chromosomes": ["NC_048211.1", "NC_048212.1", "NC_048213.1"],
        "genome_size_mb": 221,
        "vector_competence": ["Plasmodium", "Dengue virus"]
    }
}

class MosquitoGenomeDownloader:
    """
    Downloads mosquito genome assemblies from NCBI and other public databases.
    
    Key mosquito species targeted:
    - Aedes aegypti (yellow fever, dengue, Zika vector)
    - Culex quinquefasciatus (West Nile virus vector)
    - Anopheles gambiae (malaria vector)
    - Aedes albopictus (Asian tiger mosquito)
    """
    
    # Major mosquito species of interest
    MOSQUITO_SPECIES = {
        'aedes_aegypti': 'Aedes aegypti',
        'culex_quinquefasciatus': 'Culex quinquefasciatus', 
        'anopheles_gambiae': 'Anopheles gambiae',
        'aedes_albopictus': 'Aedes albopictus',
        'anopheles_coluzzii': 'Anopheles coluzzii',
        'culex_pipiens': 'Culex pipiens'
    }
    
    def __init__(self, email: str, api_key: Optional[str] = None, output_dir: str = "genomes"):
        """
        Initialize the genome downloader.
        
        Args:
            email: Required by NCBI for API access
            api_key: Optional NCBI API key for higher rate limits
            output_dir: Directory to save downloaded genomes
        """
        if not email:
            raise ValueError("Email is required for NCBI API access")
            
        if Entrez:
            Entrez.email = email
            if api_key:
                Entrez.api_key = api_key
            
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.mosquito_db = MOSQUITO_GENOMES
        
        # NCBI FTP base URLs
        self.ncbi_ftp_base = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        self.ncbi_api_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        logger.info(f"Initialized MosquitoGenomeDownloader with output directory: {self.output_dir}")
    
    def get_assembly_summary(self, assembly_id: str) -> Dict:
        """
        Get assembly summary information from NCBI.
        
        Args:
            assembly_id: NCBI assembly ID
            
        Returns:
            Dictionary containing assembly summary information
        """
        try:
            handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
            summary = Entrez.read(handle)
            handle.close()
            return summary
        except Exception as e:
            logger.error(f"Error getting assembly summary for {assembly_id}: {e}")
            return {}
    
    def search_mosquito_assemblies(self, species_name: str, max_results: int = 50) -> List[str]:
        """
        Search for mosquito genome assemblies in NCBI.
        
        Args:
            species_name: Scientific name of mosquito species
            max_results: Maximum number of results to return
            
        Returns:
            List of assembly IDs
        """
        try:
            # Search for assemblies with preference for RefSeq and complete genomes
            search_term = f'{species_name}[Organism] AND ("latest refseq"[filter] OR "complete genome"[filter])'
            
            handle = Entrez.esearch(
                db="assembly", 
                term=search_term, 
                retmax=max_results,
                sort="relevance"
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            assembly_ids = search_results['IdList']
            logger.info(f"Found {len(assembly_ids)} assemblies for {species_name}")
            
            return assembly_ids
            
        except Exception as e:
            logger.error(f"Error searching assemblies for {species_name}: {e}")
            return []
    
    def get_assembly_download_info(self, assembly_id: str) -> Optional[Dict]:
        """
        Get download information for a specific assembly.
        
        Args:
            assembly_id: NCBI assembly ID
            
        Returns:
            Dictionary with download information or None if error
        """
        summary = self.get_assembly_summary(assembly_id)
        
        if not summary or 'DocumentSummarySet' not in summary:
            return None
            
        doc_summary = summary['DocumentSummarySet']['DocumentSummary'][0]
        
        # Prefer RefSeq over GenBank
        ftp_path = doc_summary.get('FtpPath_RefSeq') or doc_summary.get('FtpPath_GenBank')
        
        if not ftp_path:
            logger.warning(f"No FTP path found for assembly {assembly_id}")
            return None
            
        assembly_name = os.path.basename(ftp_path)
        organism = doc_summary.get('Organism', 'Unknown')
        assembly_level = doc_summary.get('AssemblyStatus', 'Unknown')
        
        return {
            'assembly_id': assembly_id,
            'assembly_name': assembly_name,
            'organism': organism,
            'assembly_level': assembly_level,
            'ftp_path': ftp_path,
            'genomic_fna_url': f"{ftp_path}/{assembly_name}_genomic.fna.gz",
            'gff_url': f"{ftp_path}/{assembly_name}_genomic.gff.gz",
            'protein_faa_url': f"{ftp_path}/{assembly_name}_protein.faa.gz"
        }
    
    def download_file(self, url: str, output_path: Path, max_retries: int = 3) -> bool:
        """
        Download a file from URL with retry logic.
        
        Args:
            url: URL to download
            output_path: Local path to save file
            max_retries: Maximum number of retry attempts
            
        Returns:
            True if successful, False otherwise
        """
        for attempt in range(max_retries):
            try:
                logger.info(f"Downloading {url} to {output_path}")
                urllib.request.urlretrieve(url, output_path)
                
                # Verify file was downloaded and has content
                if output_path.exists() and output_path.stat().st_size > 0:
                    logger.info(f"Successfully downloaded {output_path.name}")
                    return True
                else:
                    logger.warning(f"Downloaded file {output_path.name} is empty or missing")
                    
            except Exception as e:
                logger.error(f"Attempt {attempt + 1} failed for {url}: {e}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff
                    
        return False
    
    def download_assembly(self, assembly_info: Dict, download_types: List[str] = None) -> Dict[str, bool]:
        """
        Download files for a specific assembly.
        
        Args:
            assembly_info: Assembly information dictionary
            download_types: List of file types to download ['genomic', 'gff', 'protein']
            
        Returns:
            Dictionary indicating success/failure for each file type
        """
        if download_types is None:
            download_types = ['genomic', 'gff']  # Default to genome and annotation
            
        assembly_name = assembly_info['assembly_name']
        organism_dir = self.output_dir / assembly_info['organism'].replace(' ', '_')
        organism_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        # Download genomic FASTA
        if 'genomic' in download_types:
            genomic_path = organism_dir / f"{assembly_name}_genomic.fna.gz"
            results['genomic'] = self.download_file(assembly_info['genomic_fna_url'], genomic_path)
        
        # Download GFF annotation
        if 'gff' in download_types:
            gff_path = organism_dir / f"{assembly_name}_genomic.gff.gz"
            results['gff'] = self.download_file(assembly_info['gff_url'], gff_path)
        
        # Download protein sequences
        if 'protein' in download_types:
            protein_path = organism_dir / f"{assembly_name}_protein.faa.gz"
            results['protein'] = self.download_file(assembly_info['protein_faa_url'], protein_path)
        
        return results
    
    def download_species_genomes(self, species_key: str, max_assemblies: int = 5) -> List[Dict]:
        """
        Download genome assemblies for a specific mosquito species.
        
        Args:
            species_key: Key from MOSQUITO_SPECIES dictionary
            max_assemblies: Maximum number of assemblies to download
            
        Returns:
            List of download results
        """
        if species_key not in self.MOSQUITO_SPECIES:
            raise ValueError(f"Unknown species key: {species_key}. Available: {list(self.MOSQUITO_SPECIES.keys())}")
        
        species_name = self.MOSQUITO_SPECIES[species_key]
        logger.info(f"Starting download for {species_name}")
        
        # Search for assemblies
        assembly_ids = self.search_mosquito_assemblies(species_name, max_assemblies * 2)
        
        if not assembly_ids:
            logger.warning(f"No assemblies found for {species_name}")
            return []
        
        download_results = []
        downloaded_count = 0
        
        for assembly_id in assembly_ids:
            if downloaded_count >= max_assemblies:
                break
                
            # Get assembly information
            assembly_info = self.get_assembly_download_info(assembly_id)
            if not assembly_info:
                continue
            
            logger.info(f"Processing {assembly_info['assembly_name']} ({assembly_info['assembly_level']})")
            
            # Download assembly files
            download_result = self.download_assembly(assembly_info)
            download_result['assembly_info'] = assembly_info
            download_results.append(download_result)
            
            if any(download_result.values()):
                downloaded_count += 1
            
            # Rate limiting - respect NCBI guidelines
            time.sleep(0.5)
        
        logger.info(f"Completed download for {species_name}: {downloaded_count} assemblies")
        return download_results
    
    def download_all_mosquito_species(self, max_assemblies_per_species: int = 3) -> Dict[str, List[Dict]]:
        """
        Download genome assemblies for all major mosquito species.
        
        Args:
            max_assemblies_per_species: Maximum assemblies per species
            
        Returns:
            Dictionary mapping species to download results
        """
        all_results = {}
        
        for species_key in self.MOSQUITO_SPECIES:
            try:
                results = self.download_species_genomes(species_key, max_assemblies_per_species)
                all_results[species_key] = results
                
                # Longer pause between species to be respectful to NCBI
                time.sleep(2)
                
            except Exception as e:
                logger.error(f"Error downloading {species_key}: {e}")
                all_results[species_key] = []
        
        return all_results
    
    def download_real_mosquito_genomes(self, species_list: List[str] = None, max_assemblies: int = 3) -> Dict[str, List[str]]:
        """
        Download real mosquito genome assemblies from NCBI using verified accessions.
        
        Args:
            species_list: List of species keys to download (default: all in database)
            max_assemblies: Maximum number of assemblies per species
            
        Returns:
            Dictionary mapping species to list of downloaded file paths
        """
        if species_list is None:
            species_list = list(self.mosquito_db.keys())
            
        downloaded_files = {}
        
        for species_key in species_list:
            if species_key not in self.mosquito_db:
                logger.warning(f"Species {species_key} not found in database")
                continue
                
            species_info = self.mosquito_db[species_key]
            species_name = species_info["species_name"]
            
            logger.info(f"Downloading real genome for {species_name} ({species_key})")
            
            try:
                species_files = []
                species_dir = self.output_dir / species_key
                species_dir.mkdir(parents=True, exist_ok=True)
                
                # Download main assembly
                assembly_acc = species_info["ncbi_assembly_accession"]
                assembly_name = species_info["assembly_name"]
                
                # Construct NCBI FTP URL
                acc_parts = assembly_acc.split("_")
                if len(acc_parts) >= 2:
                    # Format: GCF_000005575.2 -> GCF/000/005/575/GCF_000005575.2_AssemblyName
                    prefix = acc_parts[0]
                    number = acc_parts[1].split(".")[0]
                    
                    # Create directory structure
                    dir_path = f"{prefix}/{number[:3]}/{number[3:6]}/{number[6:9]}"
                    filename_base = f"{assembly_acc}_{assembly_name}"
                    
                    # Try different file extensions
                    file_extensions = [
                        "_genomic.fna.gz",
                        "_genomic.gbff.gz",
                        "_cds_from_genomic.fna.gz"
                    ]
                    
                    for ext in file_extensions:
                        download_url = f"{self.ncbi_ftp_base}{dir_path}/{filename_base}/{filename_base}{ext}"
                        output_file = species_dir / f"{species_key}_genome{ext}"
                        
                        if self._download_file_with_retry(download_url, output_file):
                            # Decompress if needed
                            if output_file.suffix == ".gz":
                                decompressed_file = output_file.with_suffix("")
                                try:
                                    with gzip.open(output_file, 'rb') as f_in:
                                        with open(decompressed_file, 'wb') as f_out:
                                            import shutil
                                            shutil.copyfileobj(f_in, f_out)
                                    output_file.unlink()  # Remove compressed file
                                    species_files.append(str(decompressed_file))
                                    logger.info(f"Downloaded and decompressed {ext} for {species_name}")
                                except Exception as e:
                                    logger.error(f"Error decompressing {output_file}: {e}")
                            else:
                                species_files.append(str(output_file))
                                logger.info(f"Downloaded {ext} for {species_name}")
                            
                            break  # Successfully downloaded one format
                    
                    # Also try to download individual chromosomes
                    for i, chr_acc in enumerate(species_info.get("chromosomes", [])[:max_assemblies]):
                        chr_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chr_acc}&rettype=fasta&retmode=text"
                        chr_file = species_dir / f"{species_key}_chr{i+1}.fna"
                        
                        if self._download_file_with_retry(chr_url, chr_file):
                            species_files.append(str(chr_file))
                            logger.info(f"Downloaded chromosome {chr_acc} for {species_name}")
                
                downloaded_files[species_key] = species_files
                
                # Save metadata
                metadata_file = species_dir / "genome_metadata.json"
                with open(metadata_file, 'w') as f:
                    json.dump(species_info, f, indent=2)
                
            except Exception as e:
                logger.error(f"Error downloading genomes for {species_name}: {e}")
                downloaded_files[species_key] = []
                
        return downloaded_files
    
    def _download_file_with_retry(self, url: str, output_path: Path, max_retries: int = 3) -> bool:
        """
        Download a file with retry logic and better error handling.
        
        Args:
            url: URL to download from
            output_path: Local path to save file
            max_retries: Maximum number of retry attempts
            
        Returns:
            True if download successful, False otherwise
        """
        for attempt in range(max_retries):
            try:
                logger.info(f"Downloading {url} (attempt {attempt + 1}/{max_retries})")
                
                headers = {
                    'User-Agent': 'Mozilla/5.0 (compatible; GenomeDownloader/1.0)'
                }
                
                import requests
                response = requests.get(url, stream=True, timeout=60, headers=headers)
                response.raise_for_status()
                
                # Check if we got actual content
                content_length = response.headers.get('content-length')
                if content_length and int(content_length) < 100:
                    logger.warning(f"File seems too small: {content_length} bytes")
                    continue
                
                with open(output_path, 'wb') as f:
                    downloaded = 0
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            
                # Verify file was downloaded
                if output_path.exists() and output_path.stat().st_size > 0:
                    logger.info(f"Successfully downloaded {output_path} ({output_path.stat().st_size} bytes)")
                    return True
                else:
                    logger.warning(f"Downloaded file is empty or missing")
                    
            except Exception as e:
                if hasattr(e, 'response') and hasattr(e.response, 'status_code') and e.response.status_code == 404:
                    logger.warning(f"File not found (404): {url}")
                    return False  # Don't retry 404s
                else:
                    logger.warning(f"Download attempt {attempt + 1} failed: {e}")
                
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
                    
        logger.error(f"Failed to download {url} after {max_retries} attempts")
        return False
    
    def get_downloaded_genomes_summary(self) -> Dict[str, List[str]]:
        """
        Get summary of downloaded genome files.
        
        Returns:
            Dictionary mapping species to list of downloaded files
        """
        summary = {}
        
        for species_dir in self.output_dir.iterdir():
            if species_dir.is_dir():
                files = [f.name for f in species_dir.iterdir() if f.is_file()]
                summary[species_dir.name] = files
        
        return summary


def main():
    """
    Example usage of the MosquitoGenomeDownloader.
    """
    # Initialize downloader (replace with your email)
    downloader = MosquitoGenomeDownloader(
        email="your.email@example.com",  # Replace with actual email
        output_dir="mosquito_genomes"
    )
    
    # Download genomes for Aedes aegypti (primary Zika/dengue vector)
    print("Downloading Aedes aegypti genomes...")
    aedes_results = downloader.download_species_genomes('aedes_aegypti', max_assemblies=2)
    
    # Download genomes for Culex quinquefasciatus (West Nile vector)
    print("Downloading Culex quinquefasciatus genomes...")
    culex_results = downloader.download_species_genomes('culex_quinquefasciatus', max_assemblies=2)
    
    # Get summary of downloaded files
    summary = downloader.get_downloaded_genomes_summary()
    print("\nDownloaded genomes summary:")
    for species, files in summary.items():
        print(f"{species}: {len(files)} files")
        for file in files[:3]:  # Show first 3 files
            print(f"  - {file}")
        if len(files) > 3:
            print(f"  ... and {len(files) - 3} more")


if __name__ == "__main__":
    main()