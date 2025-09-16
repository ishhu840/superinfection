"""Main data processing pipeline for mosquito sequencing data."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..utils import get_logger, log_processing_step, get_config_value


class ViralAnalysisPipeline:
    """Main pipeline for processing mosquito sequencing data to detect viral interference."""
    
    def __init__(self, config_path: str = None):
        """Initialize the viral analysis pipeline.
        
        Args:
            config_path: Path to configuration file
        """
        self.logger = get_logger('data_processing')
        self.config_path = config_path
        
        # Load configuration parameters
        self.min_read_length = get_config_value('data_processing.min_read_length', 50)
        self.min_quality_score = get_config_value('data_processing.min_quality_score', 20)
        self.min_contig_length = get_config_value('data_processing.min_contig_length', 200)
        self.min_abundance_threshold = get_config_value('data_processing.min_abundance_threshold', 1.0)
        self.normalization_method = get_config_value('data_processing.normalization_method', 'rpm')
        
        self.logger.info("Viral analysis pipeline initialized")
    
    def process_sample(self, 
                      sample_id: str,
                      fastq_files: List[str],
                      metadata: Dict,
                      output_dir: str) -> Dict:
        """Process a single mosquito sample through the complete pipeline.
        
        Args:
            sample_id: Unique identifier for the sample
            fastq_files: List of FASTQ file paths (single-end or paired-end)
            metadata: Sample metadata dictionary
            output_dir: Output directory for results
            
        Returns:
            Dictionary containing processing results and viral abundance data
        """
        self.logger.info(f"Processing sample {sample_id}")
        
        results = {
            'sample_id': sample_id,
            'metadata': metadata,
            'processing_stats': {},
            'viral_abundance': {},
            'viral_contigs': [],
            'quality_metrics': {}
        }
        
        output_path = Path(output_dir) / sample_id
        output_path.mkdir(parents=True, exist_ok=True)
        
        try:
            # Step 1: Quality Control
            with log_processing_step("Quality Control", self.logger):
                qc_files, qc_stats = self._quality_control(fastq_files, output_path)
                results['processing_stats']['quality_control'] = qc_stats
            
            # Step 2: Host Removal
            with log_processing_step("Host Removal", self.logger):
                non_host_files, host_stats = self._remove_host_reads(qc_files, output_path)
                results['processing_stats']['host_removal'] = host_stats
            
            # Step 3: Viral Classification
            with log_processing_step("Viral Classification", self.logger):
                viral_reads, classification_stats = self._classify_viral_reads(non_host_files, output_path)
                results['processing_stats']['viral_classification'] = classification_stats
            
            # Step 4: Assembly
            with log_processing_step("Viral Assembly", self.logger):
                contigs, assembly_stats = self._assemble_viral_contigs(viral_reads, output_path)
                results['viral_contigs'] = contigs
                results['processing_stats']['assembly'] = assembly_stats
            
            # Step 5: Quantification
            with log_processing_step("Abundance Quantification", self.logger):
                abundance_data = self._quantify_viral_abundance(viral_reads, contigs, output_path)
                results['viral_abundance'] = abundance_data
            
            # Step 6: Quality Metrics
            with log_processing_step("Quality Assessment", self.logger):
                quality_metrics = self._assess_quality(results)
                results['quality_metrics'] = quality_metrics
            
            self.logger.info(f"Successfully processed sample {sample_id}")
            
        except Exception as e:
            self.logger.error(f"Error processing sample {sample_id}: {str(e)}")
            results['error'] = str(e)
            raise
        
        return results
    
    def _quality_control(self, fastq_files: List[str], output_path: Path) -> Tuple[List[str], Dict]:
        """Perform quality control on raw sequencing reads.
        
        Args:
            fastq_files: List of input FASTQ files
            output_path: Output directory
            
        Returns:
            Tuple of (processed_files, statistics)
        """
        stats = {
            'total_reads_input': 0,
            'total_reads_output': 0,
            'reads_filtered': 0,
            'mean_quality_before': 0,
            'mean_quality_after': 0
        }
        
        processed_files = []
        
        for i, fastq_file in enumerate(fastq_files):
            input_path = Path(fastq_file)
            output_file = output_path / f"qc_{input_path.name}"
            
            # Simulate quality control (in real implementation, would use fastp or similar)
            reads_in, reads_out = self._simulate_quality_control(input_path, output_file)
            
            stats['total_reads_input'] += reads_in
            stats['total_reads_output'] += reads_out
            stats['reads_filtered'] += (reads_in - reads_out)
            
            processed_files.append(str(output_file))
        
        stats['retention_rate'] = stats['total_reads_output'] / max(stats['total_reads_input'], 1)
        
        return processed_files, stats
    
    def _simulate_quality_control(self, input_file: Path, output_file: Path) -> Tuple[int, int]:
        """Simulate quality control step (placeholder for real implementation).
        
        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            
        Returns:
            Tuple of (reads_input, reads_output)
        """
        # In a real implementation, this would use tools like fastp
        # For now, we'll simulate the process
        
        reads_input = np.random.randint(10000, 100000)  # Simulate read count
        retention_rate = np.random.uniform(0.8, 0.95)  # Simulate QC retention
        reads_output = int(reads_input * retention_rate)
        
        # Create a dummy output file
        output_file.touch()
        
        return reads_input, reads_output
    
    def _remove_host_reads(self, fastq_files: List[str], output_path: Path) -> Tuple[List[str], Dict]:
        """Remove mosquito host reads from the data.
        
        Args:
            fastq_files: Quality-controlled FASTQ files
            output_path: Output directory
            
        Returns:
            Tuple of (non_host_files, statistics)
        """
        stats = {
            'total_reads_input': 0,
            'host_reads_removed': 0,
            'non_host_reads': 0
        }
        
        non_host_files = []
        
        for fastq_file in fastq_files:
            input_path = Path(fastq_file)
            output_file = output_path / f"non_host_{input_path.name}"
            
            # Simulate host removal
            reads_in, non_host_reads = self._simulate_host_removal(input_path, output_file)
            
            stats['total_reads_input'] += reads_in
            stats['non_host_reads'] += non_host_reads
            stats['host_reads_removed'] += (reads_in - non_host_reads)
            
            non_host_files.append(str(output_file))
        
        stats['host_removal_rate'] = stats['host_reads_removed'] / max(stats['total_reads_input'], 1)
        
        return non_host_files, stats
    
    def _simulate_host_removal(self, input_file: Path, output_file: Path) -> Tuple[int, int]:
        """Simulate host read removal (placeholder for real implementation).
        
        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            
        Returns:
            Tuple of (reads_input, non_host_reads)
        """
        # In real implementation, would use bowtie2 or minimap2 against mosquito genome
        reads_input = np.random.randint(8000, 80000)
        host_removal_rate = np.random.uniform(0.85, 0.95)  # Most reads are host
        non_host_reads = int(reads_input * (1 - host_removal_rate))
        
        output_file.touch()
        
        return reads_input, non_host_reads
    
    def _classify_viral_reads(self, fastq_files: List[str], output_path: Path) -> Tuple[List[str], Dict]:
        """Classify reads as viral using taxonomic classification.
        
        Args:
            fastq_files: Non-host FASTQ files
            output_path: Output directory
            
        Returns:
            Tuple of (viral_read_files, statistics)
        """
        stats = {
            'total_reads_input': 0,
            'viral_reads_identified': 0,
            'classification_rate': 0
        }
        
        viral_files = []
        
        for fastq_file in fastq_files:
            input_path = Path(fastq_file)
            output_file = output_path / f"viral_{input_path.name}"
            
            # Simulate viral classification
            reads_in, viral_reads = self._simulate_viral_classification(input_path, output_file)
            
            stats['total_reads_input'] += reads_in
            stats['viral_reads_identified'] += viral_reads
            
            viral_files.append(str(output_file))
        
        stats['classification_rate'] = stats['viral_reads_identified'] / max(stats['total_reads_input'], 1)
        
        return viral_files, stats
    
    def _simulate_viral_classification(self, input_file: Path, output_file: Path) -> Tuple[int, int]:
        """Simulate viral read classification (placeholder for real implementation).
        
        Args:
            input_file: Input FASTQ file
            output_file: Output FASTQ file
            
        Returns:
            Tuple of (reads_input, viral_reads)
        """
        # In real implementation, would use Kraken2, Diamond, or BLAST
        reads_input = np.random.randint(1000, 10000)
        viral_rate = np.random.uniform(0.01, 0.1)  # Small fraction are viral
        viral_reads = int(reads_input * viral_rate)
        
        output_file.touch()
        
        return reads_input, viral_reads
    
    def _assemble_viral_contigs(self, viral_files: List[str], output_path: Path) -> Tuple[List[Dict], Dict]:
        """Assemble viral reads into contigs.
        
        Args:
            viral_files: Viral read files
            output_path: Output directory
            
        Returns:
            Tuple of (contig_list, assembly_statistics)
        """
        stats = {
            'total_contigs': 0,
            'viral_contigs': 0,
            'mean_contig_length': 0,
            'n50': 0
        }
        
        # Simulate viral contig assembly
        contigs = self._simulate_viral_assembly(output_path)
        
        stats['total_contigs'] = len(contigs)
        stats['viral_contigs'] = len([c for c in contigs if c['classification'] == 'viral'])
        
        if contigs:
            lengths = [c['length'] for c in contigs]
            stats['mean_contig_length'] = np.mean(lengths)
            stats['n50'] = self._calculate_n50(lengths)
        
        return contigs, stats
    
    def _simulate_viral_assembly(self, output_path: Path) -> List[Dict]:
        """Simulate viral contig assembly (placeholder for real implementation).
        
        Args:
            output_path: Output directory
            
        Returns:
            List of contig dictionaries
        """
        # In real implementation, would use MEGAHIT, metaSPAdes, etc.
        num_contigs = np.random.randint(5, 50)
        contigs = []
        
        # Simulate some known viruses that might be found in mosquitoes
        virus_names = [
            'Aedes aegypti densovirus',
            'Culex flavivirus',
            'Aedes albopictus C6/36 virus',
            'Phasi Charoen-like virus',
            'Hubei mosquito virus 2',
            'Menghai flavivirus',
            'Aedes anphevirus'
        ]
        
        for i in range(num_contigs):
            length = np.random.randint(200, 5000)
            is_viral = np.random.random() < 0.3  # 30% chance of being viral
            
            contig = {
                'contig_id': f'contig_{i+1}',
                'length': length,
                'coverage': np.random.uniform(1.0, 100.0),
                'classification': 'viral' if is_viral else 'unknown',
                'virus_name': np.random.choice(virus_names) if is_viral else 'Unknown',
                'taxonomy': 'Viruses' if is_viral else 'Unclassified',
                'sequence': 'N' * length  # Placeholder sequence
            }
            
            contigs.append(contig)
        
        return contigs
    
    def _calculate_n50(self, lengths: List[int]) -> int:
        """Calculate N50 statistic for contig lengths.
        
        Args:
            lengths: List of contig lengths
            
        Returns:
            N50 value
        """
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        target_length = total_length / 2
        
        cumulative_length = 0
        for length in sorted_lengths:
            cumulative_length += length
            if cumulative_length >= target_length:
                return length
        
        return 0
    
    def _quantify_viral_abundance(self, viral_files: List[str], contigs: List[Dict], output_path: Path) -> Dict:
        """Quantify viral abundance in the sample.
        
        Args:
            viral_files: Viral read files
            contigs: Assembled viral contigs
            output_path: Output directory
            
        Returns:
            Dictionary of viral abundance data
        """
        abundance_data = {}
        
        # Simulate abundance quantification
        viral_contigs = [c for c in contigs if c['classification'] == 'viral']
        
        total_reads = sum([np.random.randint(100, 10000) for _ in viral_files])
        
        for contig in viral_contigs:
            # Simulate read mapping to contigs
            mapped_reads = np.random.randint(1, 1000)
            
            # Calculate abundance metrics
            rpm = (mapped_reads / total_reads) * 1e6 if total_reads > 0 else 0
            coverage = mapped_reads * 150 / contig['length']  # Assuming 150bp reads
            
            abundance_data[contig['virus_name']] = {
                'contig_id': contig['contig_id'],
                'mapped_reads': mapped_reads,
                'rpm': rpm,
                'coverage': coverage,
                'length': contig['length'],
                'present': rpm >= self.min_abundance_threshold
            }
        
        return abundance_data
    
    def _assess_quality(self, results: Dict) -> Dict:
        """Assess overall quality of the sample processing.
        
        Args:
            results: Processing results dictionary
            
        Returns:
            Quality metrics dictionary
        """
        metrics = {
            'overall_quality': 'good',
            'sequencing_depth': 'adequate',
            'viral_diversity': 'normal',
            'contamination_risk': 'low',
            'reliability_score': 0.8
        }
        
        # Calculate quality metrics based on processing stats
        stats = results['processing_stats']
        
        # Check sequencing depth
        total_input_reads = stats.get('quality_control', {}).get('total_reads_input', 0)
        if total_input_reads < 10000:
            metrics['sequencing_depth'] = 'low'
            metrics['reliability_score'] -= 0.2
        elif total_input_reads > 100000:
            metrics['sequencing_depth'] = 'high'
        
        # Check viral detection
        viral_abundance = results['viral_abundance']
        num_viruses = len([v for v in viral_abundance.values() if v['present']])
        
        if num_viruses == 0:
            metrics['viral_diversity'] = 'none_detected'
            metrics['reliability_score'] -= 0.3
        elif num_viruses > 10:
            metrics['viral_diversity'] = 'high'
            metrics['contamination_risk'] = 'medium'
            metrics['reliability_score'] -= 0.1
        
        # Overall quality assessment
        if metrics['reliability_score'] >= 0.8:
            metrics['overall_quality'] = 'excellent'
        elif metrics['reliability_score'] >= 0.6:
            metrics['overall_quality'] = 'good'
        elif metrics['reliability_score'] >= 0.4:
            metrics['overall_quality'] = 'fair'
        else:
            metrics['overall_quality'] = 'poor'
        
        return metrics
    
    def process_batch(self, 
                     sample_manifest: str,
                     output_dir: str,
                     max_parallel: int = 4) -> pd.DataFrame:
        """Process multiple samples in batch.
        
        Args:
            sample_manifest: Path to CSV file with sample information
            output_dir: Output directory for all results
            max_parallel: Maximum number of parallel processes
            
        Returns:
            DataFrame with processing results for all samples
        """
        self.logger.info(f"Starting batch processing with manifest: {sample_manifest}")
        
        # Load sample manifest
        manifest_df = pd.read_csv(sample_manifest)
        
        results_list = []
        
        for _, row in manifest_df.iterrows():
            try:
                sample_id = row['sample_id']
                fastq_files = [row['fastq_1']]
                if 'fastq_2' in row and pd.notna(row['fastq_2']):
                    fastq_files.append(row['fastq_2'])
                
                metadata = row.to_dict()
                
                # Process single sample
                result = self.process_sample(sample_id, fastq_files, metadata, output_dir)
                results_list.append(result)
                
            except Exception as e:
                self.logger.error(f"Failed to process sample {row.get('sample_id', 'unknown')}: {e}")
                continue
        
        # Convert results to DataFrame
        results_df = self._results_to_dataframe(results_list)
        
        # Save batch results
        output_path = Path(output_dir)
        results_df.to_csv(output_path / 'batch_processing_results.csv', index=False)
        
        self.logger.info(f"Batch processing completed. Processed {len(results_list)} samples.")
        
        return results_df
    
    def _results_to_dataframe(self, results_list: List[Dict]) -> pd.DataFrame:
        """Convert processing results to a structured DataFrame.
        
        Args:
            results_list: List of processing result dictionaries
            
        Returns:
            DataFrame with structured results
        """
        rows = []
        
        for result in results_list:
            # Extract basic information
            row = {
                'sample_id': result['sample_id'],
                'overall_quality': result['quality_metrics'].get('overall_quality', 'unknown'),
                'reliability_score': result['quality_metrics'].get('reliability_score', 0),
                'total_input_reads': result['processing_stats'].get('quality_control', {}).get('total_reads_input', 0),
                'viral_reads': result['processing_stats'].get('viral_classification', {}).get('viral_reads_identified', 0),
                'num_viral_contigs': result['processing_stats'].get('assembly', {}).get('viral_contigs', 0),
                'num_viruses_detected': len([v for v in result['viral_abundance'].values() if v['present']])
            }
            
            # Add metadata
            if 'metadata' in result:
                for key, value in result['metadata'].items():
                    if key not in row:  # Avoid overwriting existing columns
                        row[f'metadata_{key}'] = value
            
            rows.append(row)
        
        return pd.DataFrame(rows)