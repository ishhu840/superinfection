#!/usr/bin/env python3
"""
Viral Detection Pipeline for Mosquito Genomes

This module implements viral sequence detection in mosquito genomes using BLAST
against viral reference databases. Based on research findings:
- NCBI RefSeq Viral database provides comprehensive viral sequences
- BLAST/Diamond can efficiently search for viral homologs
- Metagenomic analysis approaches for virus identification
"""

import os
import subprocess
import gzip
import logging
from typing import List, Dict, Optional, Tuple, Set
from pathlib import Path
import pandas as pd
import tempfile
from dataclasses import dataclass
from collections import defaultdict

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class BlastHit:
    """Represents a BLAST hit result."""
    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    query_coverage: float
    subject_coverage: float
    
    @property
    def is_significant(self) -> bool:
        """Check if hit meets significance criteria for viral detection."""
        return (
            self.evalue <= 1e-5 and
            self.identity >= 70.0 and
            self.alignment_length >= 50 and
            self.query_coverage >= 30.0
        )

@dataclass
class ViralDetectionResult:
    """Results from viral detection analysis."""
    genome_file: str
    species: str
    total_contigs: int
    viral_hits: List[BlastHit]
    unique_viruses: Set[str]
    virus_families: Dict[str, int]
    
    @property
    def viral_prevalence(self) -> float:
        """Calculate percentage of contigs with viral hits."""
        if self.total_contigs == 0:
            return 0.0
        viral_contigs = len(set(hit.query_id for hit in self.viral_hits))
        return (viral_contigs / self.total_contigs) * 100

class ViralDetectionPipeline:
    """
    Pipeline for detecting viral sequences in mosquito genomes using BLAST.
    
    This pipeline:
    1. Downloads and prepares viral reference databases
    2. Processes mosquito genome assemblies
    3. Runs BLAST searches against viral databases
    4. Filters and annotates viral hits
    5. Generates detection reports
    """
    
    # Viral database sources
    VIRAL_DATABASES = {
        'refseq_viral': {
            'url': 'https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz',
            'description': 'NCBI RefSeq Viral genomes'
        },
        'refseq_viral_protein': {
            'url': 'https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz',
            'description': 'NCBI RefSeq Viral proteins'
        }
    }
    
    # Common mosquito-associated virus families
    MOSQUITO_VIRUS_FAMILIES = {
        'Flaviviridae': ['Zika', 'Dengue', 'Yellow fever', 'West Nile'],
        'Togaviridae': ['Chikungunya', 'Eastern equine encephalitis'],
        'Bunyaviridae': ['La Crosse', 'Rift Valley fever'],
        'Reoviridae': ['Bluetongue'],
        'Rhabdoviridae': ['Vesicular stomatitis'],
        'Picornaviridae': ['Hepatitis A'],
        'Phenuiviridae': ['Rift Valley fever'],
        'Peribunyaviridae': ['La Crosse encephalitis']
    }
    
    def __init__(self, work_dir: str = "viral_detection", blast_threads: int = 4):
        """
        Initialize the viral detection pipeline.
        
        Args:
            work_dir: Working directory for databases and results
            blast_threads: Number of threads for BLAST searches
        """
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        self.db_dir = self.work_dir / "databases"
        self.db_dir.mkdir(parents=True, exist_ok=True)
        
        self.results_dir = self.work_dir / "results"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        self.blast_threads = blast_threads
        
        logger.info(f"Initialized ViralDetectionPipeline in {self.work_dir}")
    
    def check_blast_installation(self) -> bool:
        """
        Check if BLAST+ is installed and available.
        
        Returns:
            True if BLAST+ is available, False otherwise
        """
        try:
            result = subprocess.run(['blastn', '-version'], 
                                  capture_output=True, text=True, check=True)
            logger.info(f"BLAST+ found: {result.stdout.split()[1]}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error("BLAST+ not found. Please install BLAST+ toolkit.")
            logger.error("Install instructions: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download")
            return False
    
    def download_viral_database(self, db_name: str) -> Optional[Path]:
        """
        Download viral reference database.
        
        Args:
            db_name: Name of database to download
            
        Returns:
            Path to downloaded database file or None if failed
        """
        if db_name not in self.VIRAL_DATABASES:
            logger.error(f"Unknown database: {db_name}")
            return None
        
        db_info = self.VIRAL_DATABASES[db_name]
        db_file = self.db_dir / f"{db_name}.fna.gz"
        
        if db_file.exists():
            logger.info(f"Database {db_name} already exists")
            return db_file
        
        try:
            logger.info(f"Downloading {db_info['description']}...")
            import urllib.request
            urllib.request.urlretrieve(db_info['url'], db_file)
            
            if db_file.exists() and db_file.stat().st_size > 0:
                logger.info(f"Successfully downloaded {db_name}")
                return db_file
            else:
                logger.error(f"Failed to download {db_name}")
                return None
                
        except Exception as e:
            logger.error(f"Error downloading {db_name}: {e}")
            return None
    
    def prepare_blast_database(self, fasta_file: Path, db_type: str = "nucl") -> Optional[Path]:
        """
        Create BLAST database from FASTA file.
        
        Args:
            fasta_file: Path to FASTA file
            db_type: Database type ('nucl' or 'prot')
            
        Returns:
            Path to BLAST database or None if failed
        """
        if not fasta_file.exists():
            logger.error(f"FASTA file not found: {fasta_file}")
            return None
        
        # Decompress if needed
        if fasta_file.suffix == '.gz':
            decompressed_file = fasta_file.with_suffix('')
            if not decompressed_file.exists():
                logger.info(f"Decompressing {fasta_file}...")
                with gzip.open(fasta_file, 'rb') as f_in:
                    with open(decompressed_file, 'wb') as f_out:
                        f_out.write(f_in.read())
            fasta_file = decompressed_file
        
        db_path = fasta_file.with_suffix('')
        
        # Check if database already exists
        db_files = list(db_path.parent.glob(f"{db_path.name}.*"))
        if any(f.suffix in ['.nhr', '.nin', '.nsq', '.phr', '.pin', '.psq'] for f in db_files):
            logger.info(f"BLAST database already exists for {fasta_file.name}")
            return db_path
        
        try:
            cmd = [
                'makeblastdb',
                '-in', str(fasta_file),
                '-dbtype', db_type,
                '-out', str(db_path),
                '-title', f"Viral database from {fasta_file.name}"
            ]
            
            logger.info(f"Creating BLAST database: {db_path.name}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            if result.returncode == 0:
                logger.info(f"Successfully created BLAST database: {db_path.name}")
                return db_path
            else:
                logger.error(f"makeblastdb failed: {result.stderr}")
                return None
                
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating BLAST database: {e.stderr}")
            return None
    
    def run_blast_search(self, query_file: Path, database: Path, 
                        output_file: Path, blast_type: str = "blastn") -> bool:
        """
        Run BLAST search against viral database.
        
        Args:
            query_file: Path to query sequences (mosquito genome)
            database: Path to BLAST database
            output_file: Path for BLAST output
            blast_type: Type of BLAST search ('blastn', 'blastx', 'tblastn')
            
        Returns:
            True if successful, False otherwise
        """
        try:
            cmd = [
                blast_type,
                '-query', str(query_file),
                '-db', str(database),
                '-out', str(output_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs scovs',
                '-evalue', '1e-5',
                '-num_threads', str(self.blast_threads),
                '-max_target_seqs', '10'
            ]
            
            logger.info(f"Running {blast_type} search: {query_file.name} vs {database.name}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            if result.returncode == 0 and output_file.exists():
                logger.info(f"BLAST search completed: {output_file.name}")
                return True
            else:
                logger.error(f"BLAST search failed: {result.stderr}")
                return False
                
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running BLAST: {e.stderr}")
            return False
    
    def parse_blast_results(self, blast_output: Path) -> List[BlastHit]:
        """
        Parse BLAST output file into BlastHit objects.
        
        Args:
            blast_output: Path to BLAST output file
            
        Returns:
            List of BlastHit objects
        """
        hits = []
        
        if not blast_output.exists() or blast_output.stat().st_size == 0:
            logger.warning(f"BLAST output file is empty or missing: {blast_output}")
            return hits
        
        try:
            with open(blast_output, 'r') as f:
                for line in f:
                    if line.strip():
                        fields = line.strip().split('\t')
                        if len(fields) >= 12:
                            # Handle both 13 and 14 field formats
                            query_coverage = float(fields[12]) if len(fields) > 12 else 0.0
                            subject_coverage = float(fields[13]) if len(fields) > 13 else 0.0
                            
                            hit = BlastHit(
                                query_id=fields[0],
                                subject_id=fields[1],
                                identity=float(fields[2]),
                                alignment_length=int(fields[3]),
                                mismatches=int(fields[4]),
                                gap_opens=int(fields[5]),
                                query_start=int(fields[6]),
                                query_end=int(fields[7]),
                                subject_start=int(fields[8]),
                                subject_end=int(fields[9]),
                                evalue=float(fields[10]),
                                bit_score=float(fields[11]),
                                query_coverage=query_coverage,
                                subject_coverage=subject_coverage
                            )
                            hits.append(hit)
            
            logger.info(f"Parsed {len(hits)} BLAST hits from {blast_output.name}")
            
        except Exception as e:
            logger.error(f"Error parsing BLAST results: {e}")
        
        return hits
    
    def filter_significant_hits(self, hits: List[BlastHit]) -> List[BlastHit]:
        """
        Filter BLAST hits for significant viral matches.
        
        Args:
            hits: List of all BLAST hits
            
        Returns:
            List of significant viral hits
        """
        significant_hits = [hit for hit in hits if hit.is_significant]
        
        logger.info(f"Filtered {len(significant_hits)} significant hits from {len(hits)} total hits")
        
        return significant_hits
    
    def annotate_viral_hits(self, hits: List[BlastHit]) -> Dict[str, List[BlastHit]]:
        """
        Annotate viral hits with virus family information.
        
        Args:
            hits: List of significant viral hits
            
        Returns:
            Dictionary mapping virus families to hits
        """
        family_hits = defaultdict(list)
        
        for hit in hits:
            # Extract virus information from subject ID
            subject_lower = hit.subject_id.lower()
            
            # Try to match known virus families
            matched_family = None
            for family, viruses in self.MOSQUITO_VIRUS_FAMILIES.items():
                for virus in viruses:
                    if virus.lower() in subject_lower:
                        matched_family = family
                        break
                if matched_family:
                    break
            
            if not matched_family:
                # Generic classification based on common viral terms
                if any(term in subject_lower for term in ['virus', 'viral', 'phage']):
                    matched_family = 'Other_viruses'
                else:
                    matched_family = 'Unknown'
            
            family_hits[matched_family].append(hit)
        
        return dict(family_hits)
    
    def count_contigs_in_genome(self, genome_file: Path) -> int:
        """
        Count the number of contigs/sequences in a genome file.
        
        Args:
            genome_file: Path to genome FASTA file
            
        Returns:
            Number of sequences/contigs
        """
        count = 0
        
        try:
            if genome_file.suffix == '.gz':
                opener = gzip.open
                mode = 'rt'
            else:
                opener = open
                mode = 'r'
            
            with opener(genome_file, mode) as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
        
        except Exception as e:
            logger.error(f"Error counting contigs in {genome_file}: {e}")
        
        return count
    
    def analyze_genome(self, genome_file: Path, species: str, 
                      viral_db: Path) -> ViralDetectionResult:
        """
        Analyze a single mosquito genome for viral sequences.
        
        Args:
            genome_file: Path to mosquito genome file
            species: Species name
            viral_db: Path to viral BLAST database
            
        Returns:
            ViralDetectionResult object
        """
        logger.info(f"Analyzing {genome_file.name} for viral sequences")
        
        # Count total contigs
        total_contigs = self.count_contigs_in_genome(genome_file)
        
        # Prepare output file
        blast_output = self.results_dir / f"{genome_file.stem}_viral_hits.txt"
        
        # Run BLAST search
        blast_success = self.run_blast_search(genome_file, viral_db, blast_output)
        
        if not blast_success:
            return ViralDetectionResult(
                genome_file=str(genome_file),
                species=species,
                total_contigs=total_contigs,
                viral_hits=[],
                unique_viruses=set(),
                virus_families={}
            )
        
        # Parse and filter results
        all_hits = self.parse_blast_results(blast_output)
        significant_hits = self.filter_significant_hits(all_hits)
        
        # Annotate hits
        family_hits = self.annotate_viral_hits(significant_hits)
        
        # Extract unique viruses
        unique_viruses = set(hit.subject_id for hit in significant_hits)
        
        # Count hits per family
        virus_families = {family: len(hits) for family, hits in family_hits.items()}
        
        result = ViralDetectionResult(
            genome_file=str(genome_file),
            species=species,
            total_contigs=total_contigs,
            viral_hits=significant_hits,
            unique_viruses=unique_viruses,
            virus_families=virus_families
        )
        
        logger.info(f"Found {len(significant_hits)} viral hits in {genome_file.name} ({result.viral_prevalence:.1f}% prevalence)")
        
        return result
    
    def setup_viral_databases(self) -> Optional[Path]:
        """
        Download and prepare viral databases for BLAST searches.
        
        Returns:
            Path to prepared BLAST database or None if failed
        """
        if not self.check_blast_installation():
            return None
        
        # Download viral database
        viral_fasta = self.download_viral_database('refseq_viral')
        if not viral_fasta:
            return None
        
        # Create BLAST database
        viral_db = self.prepare_blast_database(viral_fasta, 'nucl')
        
        return viral_db
    
    def analyze_mosquito_genomes(self, genome_dir: Path) -> List[ViralDetectionResult]:
        """
        Analyze all mosquito genomes in a directory for viral sequences.
        
        Args:
            genome_dir: Directory containing mosquito genome files
            
        Returns:
            List of ViralDetectionResult objects
        """
        # Setup viral databases
        viral_db = self.setup_viral_databases()
        if not viral_db:
            logger.error("Failed to setup viral databases")
            return []
        
        results = []
        
        # Find all genome files
        genome_files = []
        for pattern in ['*.fna', '*.fna.gz', '*.fa', '*.fa.gz', '*.fasta', '*.fasta.gz']:
            genome_files.extend(genome_dir.rglob(pattern))
        
        logger.info(f"Found {len(genome_files)} genome files to analyze")
        
        for genome_file in genome_files:
            try:
                # Extract species from path
                species = genome_file.parent.name.replace('_', ' ')
                
                # Analyze genome
                result = self.analyze_genome(genome_file, species, viral_db)
                results.append(result)
                
            except Exception as e:
                logger.error(f"Error analyzing {genome_file}: {e}")
        
        return results
    
    def generate_summary_report(self, results: List[ViralDetectionResult]) -> pd.DataFrame:
        """
        Generate summary report of viral detection results.
        
        Args:
            results: List of ViralDetectionResult objects
            
        Returns:
            DataFrame with summary statistics
        """
        summary_data = []
        
        for result in results:
            summary_data.append({
                'Species': result.species,
                'Genome_File': Path(result.genome_file).name,
                'Total_Contigs': result.total_contigs,
                'Viral_Hits': len(result.viral_hits),
                'Viral_Prevalence_%': result.viral_prevalence,
                'Unique_Viruses': len(result.unique_viruses),
                'Virus_Families': len(result.virus_families),
                'Top_Family': max(result.virus_families.items(), key=lambda x: x[1])[0] if result.virus_families else 'None'
            })
        
        summary_df = pd.DataFrame(summary_data)
        
        # Save summary report
        summary_file = self.results_dir / "viral_detection_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Summary report saved to {summary_file}")
        
        return summary_df


def main():
    """
    Example usage of the ViralDetectionPipeline.
    """
    # Initialize pipeline
    pipeline = ViralDetectionPipeline(work_dir="viral_analysis")
    
    # Analyze genomes (assuming they were downloaded by genome_downloader)
    genome_dir = Path("mosquito_genomes")
    
    if genome_dir.exists():
        print("Starting viral detection analysis...")
        results = pipeline.analyze_mosquito_genomes(genome_dir)
        
        if results:
            print(f"\nAnalyzed {len(results)} genomes")
            
            # Generate summary report
            summary = pipeline.generate_summary_report(results)
            print("\nSummary Report:")
            print(summary.to_string(index=False))
            
            # Show some interesting findings
            total_viral_hits = sum(len(r.viral_hits) for r in results)
            species_with_viruses = len([r for r in results if r.viral_hits])
            
            print(f"\nKey Findings:")
            print(f"- Total viral hits found: {total_viral_hits}")
            print(f"- Species with viral sequences: {species_with_viruses}/{len(results)}")
            
            if results:
                avg_prevalence = sum(r.viral_prevalence for r in results) / len(results)
                print(f"- Average viral prevalence: {avg_prevalence:.1f}%")
        
    else:
        print(f"Genome directory not found: {genome_dir}")
        print("Please run genome_downloader.py first to download mosquito genomes")


if __name__ == "__main__":
    main()