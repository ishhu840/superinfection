#!/usr/bin/env python3
"""
Batch Viral Analysis for Large-Scale Mosquito Genome Data

Processes 100K+ mosquito genome files for viral detection and superinfection analysis.
Optimized for high-throughput processing with parallel execution and progress tracking.
"""

import os
import sys
import json
import time
import logging
import multiprocessing as mp
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Optional, Tuple
import pandas as pd
import numpy as np
from dataclasses import asdict

# Add project root to path
sys.path.append(str(Path(__file__).parent))

from src.viral_detection.blast_pipeline import ViralDetectionPipeline
from src.analysis.superinfection_analysis import SuperinfectionAnalyzer
from src.utils.file_utils import get_genome_files

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('batch_viral_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class BatchViralAnalyzer:
    """High-throughput viral analysis for massive mosquito genome datasets."""
    
    def __init__(self, 
                 input_dir: str = "massive_mosquito_data",
                 output_dir: str = "massive_viral_results",
                 max_workers: int = None):
        
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Create output subdirectories
        self.results_dir = self.output_dir / "results"
        self.summaries_dir = self.output_dir / "summaries"
        self.databases_dir = self.output_dir / "databases"
        self.logs_dir = self.output_dir / "logs"
        
        for dir_path in [self.results_dir, self.summaries_dir, self.databases_dir, self.logs_dir]:
            dir_path.mkdir(exist_ok=True)
        
        # Set up parallel processing
        self.max_workers = max_workers or min(mp.cpu_count(), 16)  # Limit to 16 cores max
        
        # Initialize viral detection pipeline
        self.pipeline = ViralDetectionPipeline()
        
        # Initialize superinfection analyzer
        self.superinfection_analyzer = SuperinfectionAnalyzer()
        
        # Batch processing stats
        self.processing_stats = {
            "start_time": None,
            "end_time": None,
            "total_files": 0,
            "processed_files": 0,
            "failed_files": 0,
            "total_viral_hits": 0,
            "species_analyzed": set(),
            "viral_families_detected": set(),
            "superinfection_pairs": 0
        }
    
    def discover_genome_files(self) -> Dict[str, List[Path]]:
        """Discover all genome files organized by species."""
        logger.info(f"Discovering genome files in {self.input_dir}...")
        
        genome_files = {
            "reference_genomes": [],
            "sra_datasets": [],
            "total_count": 0
        }
        
        # Find reference genomes
        ref_dir = self.input_dir / "reference_genomes"
        if ref_dir.exists():
            ref_files = list(ref_dir.rglob("*.fna.gz")) + list(ref_dir.rglob("*.fna"))
            genome_files["reference_genomes"] = ref_files
            logger.info(f"Found {len(ref_files):,} reference genome files")
        
        # Find SRA datasets
        sra_dir = self.input_dir / "sra_datasets"
        if sra_dir.exists():
            sra_files = list(sra_dir.rglob("*.fastq.gz")) + list(sra_dir.rglob("*.fastq"))
            genome_files["sra_datasets"] = sra_files
            logger.info(f"Found {len(sra_files):,} SRA dataset files")
        
        total_files = len(genome_files["reference_genomes"]) + len(genome_files["sra_datasets"])
        genome_files["total_count"] = total_files
        
        logger.info(f"Total genome files discovered: {total_files:,}")
        return genome_files
    
    def process_single_genome(self, file_info: Tuple[Path, str]) -> Optional[Dict]:
        """Process a single genome file for viral detection."""
        genome_file, file_type = file_info
        
        try:
            # Extract species from file path
            species = self._extract_species_from_path(genome_file)
            
            # Run viral detection
            if file_type == "reference":
                result = self.pipeline.analyze_genome(str(genome_file))
            else:  # SRA dataset
                # For SRA datasets, we might need different processing
                result = self._process_sra_dataset(genome_file)
            
            if result:
                # Convert result to dictionary for JSON serialization
                result_dict = asdict(result) if hasattr(result, '__dict__') else result
                result_dict["species"] = species
                result_dict["file_type"] = file_type
                result_dict["file_path"] = str(genome_file)
                result_dict["processing_time"] = datetime.now().isoformat()
                
                # Save individual result
                result_file = self.results_dir / f"{genome_file.stem}_viral_result.json"
                with open(result_file, 'w') as f:
                    json.dump(result_dict, f, indent=2)
                
                return result_dict
            
        except Exception as e:
            logger.error(f"Error processing {genome_file}: {e}")
            return None
    
    def _extract_species_from_path(self, file_path: Path) -> str:
        """Extract species name from file path."""
        # Try to extract from parent directory name
        parent_name = file_path.parent.name
        if "_" in parent_name:
            return parent_name.replace("_", " ")
        
        # Try to extract from filename
        filename = file_path.stem
        if "_" in filename:
            parts = filename.split("_")
            if len(parts) >= 2:
                return f"{parts[0]} {parts[1]}"
        
        return "Unknown species"
    
    def _process_sra_dataset(self, sra_file: Path) -> Optional[Dict]:
        """Process SRA dataset (FASTQ) for viral sequences."""
        # For SRA datasets, we might need to:
        # 1. Convert FASTQ to FASTA
        # 2. Run viral detection on reads
        # 3. Assemble contigs if needed
        
        # Simplified approach: convert FASTQ to FASTA and analyze
        try:
            fasta_file = sra_file.with_suffix('.fasta')
            
            # Convert FASTQ to FASTA (simplified)
            if not fasta_file.exists():
                self._fastq_to_fasta(sra_file, fasta_file)
            
            # Run viral detection on converted file
            if fasta_file.exists():
                result = self.pipeline.analyze_genome(str(fasta_file))
                return result
            
        except Exception as e:
            logger.error(f"Error processing SRA dataset {sra_file}: {e}")
            return None
    
    def _fastq_to_fasta(self, fastq_file: Path, fasta_file: Path) -> None:
        """Convert FASTQ to FASTA format."""
        import gzip
        
        opener = gzip.open if fastq_file.suffix == '.gz' else open
        
        with opener(fastq_file, 'rt') as infile, open(fasta_file, 'w') as outfile:
            line_count = 0
            for line in infile:
                line_count += 1
                if line_count % 4 == 1:  # Header line
                    outfile.write(f">{line[1:]}")
                elif line_count % 4 == 2:  # Sequence line
                    outfile.write(line)
                # Skip quality lines (3 and 4)
    
    def process_batch(self, genome_files: Dict[str, List[Path]], 
                     batch_size: int = 1000) -> List[Dict]:
        """Process genome files in batches with parallel execution."""
        logger.info(f"Starting batch processing with {self.max_workers} workers...")
        
        all_results = []
        
        # Prepare file list with types
        file_list = []
        for ref_file in genome_files["reference_genomes"]:
            file_list.append((ref_file, "reference"))
        for sra_file in genome_files["sra_datasets"]:
            file_list.append((sra_file, "sra"))
        
        total_files = len(file_list)
        self.processing_stats["total_files"] = total_files
        self.processing_stats["start_time"] = datetime.now().isoformat()
        
        logger.info(f"Processing {total_files:,} files in batches of {batch_size}...")
        
        # Process in batches
        for i in range(0, total_files, batch_size):
            batch = file_list[i:i + batch_size]
            batch_num = i // batch_size + 1
            total_batches = (total_files + batch_size - 1) // batch_size
            
            logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} files)...")
            
            # Process batch in parallel
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                futures = [executor.submit(self.process_single_genome, file_info) 
                          for file_info in batch]
                
                batch_results = []
                for future in as_completed(futures):
                    try:
                        result = future.result(timeout=300)  # 5 minute timeout per file
                        if result:
                            batch_results.append(result)
                            self.processing_stats["processed_files"] += 1
                            
                            # Update stats
                            if "viral_hits" in result:
                                self.processing_stats["total_viral_hits"] += len(result["viral_hits"])
                            if "species" in result:
                                self.processing_stats["species_analyzed"].add(result["species"])
                            if "viral_families" in result:
                                self.processing_stats["viral_families_detected"].update(result["viral_families"])
                        else:
                            self.processing_stats["failed_files"] += 1
                            
                    except Exception as e:
                        logger.error(f"Batch processing error: {e}")
                        self.processing_stats["failed_files"] += 1
                
                all_results.extend(batch_results)
                
                # Log progress
                processed = self.processing_stats["processed_files"]
                failed = self.processing_stats["failed_files"]
                logger.info(f"Batch {batch_num} complete. Total processed: {processed:,}, Failed: {failed:,}")
                
                # Save intermediate results
                if batch_num % 10 == 0:  # Save every 10 batches
                    self._save_intermediate_results(all_results, batch_num)
        
        self.processing_stats["end_time"] = datetime.now().isoformat()
        logger.info(f"Batch processing completed. Total results: {len(all_results):,}")
        
        return all_results
    
    def _save_intermediate_results(self, results: List[Dict], batch_num: int) -> None:
        """Save intermediate results to prevent data loss."""
        intermediate_file = self.summaries_dir / f"intermediate_results_batch_{batch_num}.json"
        with open(intermediate_file, 'w') as f:
            json.dump({
                "batch_number": batch_num,
                "results_count": len(results),
                "processing_stats": self._serialize_stats(),
                "results": results
            }, f, indent=2)
        
        logger.info(f"Intermediate results saved: {intermediate_file}")
    
    def _serialize_stats(self) -> Dict:
        """Serialize processing stats for JSON output."""
        stats = self.processing_stats.copy()
        stats["species_analyzed"] = list(stats["species_analyzed"])
        stats["viral_families_detected"] = list(stats["viral_families_detected"])
        return stats
    
    def analyze_superinfection_patterns(self, results: List[Dict]) -> Dict:
        """Analyze superinfection exclusion patterns in the results."""
        logger.info("Analyzing superinfection exclusion patterns...")
        
        # Group results by species
        species_results = {}
        for result in results:
            species = result.get("species", "Unknown")
            if species not in species_results:
                species_results[species] = []
            species_results[species].append(result)
        
        superinfection_analysis = {
            "total_species": len(species_results),
            "species_analysis": {},
            "global_patterns": {},
            "exclusion_pairs": []
        }
        
        # Analyze each species
        for species, species_data in species_results.items():
            viral_cooccurrences = []
            
            for result in species_data:
                if "viral_hits" in result and len(result["viral_hits"]) > 1:
                    # Multiple viruses found - potential superinfection
                    viral_families = [hit.get("viral_family", "Unknown") 
                                    for hit in result["viral_hits"]]
                    viral_cooccurrences.append(viral_families)
            
            # Analyze cooccurrence patterns
            if viral_cooccurrences:
                exclusion_pairs = self.superinfection_analyzer.find_exclusion_pairs(viral_cooccurrences)
                superinfection_analysis["species_analysis"][species] = {
                    "samples_with_multiple_viruses": len(viral_cooccurrences),
                    "exclusion_pairs": exclusion_pairs,
                    "cooccurrence_patterns": self._analyze_cooccurrence_patterns(viral_cooccurrences)
                }
                superinfection_analysis["exclusion_pairs"].extend(exclusion_pairs)
        
        # Global analysis
        all_exclusion_pairs = superinfection_analysis["exclusion_pairs"]
        superinfection_analysis["global_patterns"] = {
            "total_exclusion_pairs": len(all_exclusion_pairs),
            "most_common_exclusions": self._get_most_common_exclusions(all_exclusion_pairs),
            "cross_species_patterns": self._analyze_cross_species_patterns(species_results)
        }
        
        self.processing_stats["superinfection_pairs"] = len(all_exclusion_pairs)
        
        return superinfection_analysis
    
    def _analyze_cooccurrence_patterns(self, cooccurrences: List[List[str]]) -> Dict:
        """Analyze viral cooccurrence patterns."""
        from collections import Counter
        
        # Count pairwise cooccurrences
        pairs = []
        for viruses in cooccurrences:
            for i in range(len(viruses)):
                for j in range(i + 1, len(viruses)):
                    pair = tuple(sorted([viruses[i], viruses[j]]))
                    pairs.append(pair)
        
        pair_counts = Counter(pairs)
        
        return {
            "total_cooccurrences": len(cooccurrences),
            "unique_pairs": len(pair_counts),
            "most_common_pairs": dict(pair_counts.most_common(10))
        }
    
    def _get_most_common_exclusions(self, exclusion_pairs: List[Tuple]) -> List[Dict]:
        """Get most common exclusion pairs."""
        from collections import Counter
        
        pair_counts = Counter(exclusion_pairs)
        return [{
            "pair": list(pair),
            "count": count
        } for pair, count in pair_counts.most_common(10)]
    
    def _analyze_cross_species_patterns(self, species_results: Dict) -> Dict:
        """Analyze patterns across species."""
        # This would involve comparing exclusion patterns between species
        # Simplified implementation
        return {
            "species_count": len(species_results),
            "analysis": "Cross-species analysis would be implemented here"
        }
    
    def generate_comprehensive_report(self, results: List[Dict], 
                                    superinfection_analysis: Dict) -> None:
        """Generate comprehensive analysis report."""
        logger.info("Generating comprehensive analysis report...")
        
        # Create summary statistics
        summary = {
            "analysis_metadata": {
                "analysis_date": datetime.now().isoformat(),
                "total_files_processed": len(results),
                "processing_stats": self._serialize_stats()
            },
            "viral_detection_summary": self._create_viral_summary(results),
            "superinfection_analysis": superinfection_analysis,
            "species_breakdown": self._create_species_breakdown(results),
            "recommendations": self._generate_recommendations(results, superinfection_analysis)
        }
        
        # Save comprehensive summary
        summary_file = self.summaries_dir / "comprehensive_analysis_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Save all results
        all_results_file = self.summaries_dir / "all_viral_detection_results.json"
        with open(all_results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Generate CSV for easy analysis
        self._generate_csv_reports(results)
        
        # Generate markdown report
        self._generate_markdown_report(summary)
        
        logger.info(f"Comprehensive report saved to: {summary_file}")
    
    def _create_viral_summary(self, results: List[Dict]) -> Dict:
        """Create viral detection summary."""
        total_viral_hits = sum(len(r.get("viral_hits", [])) for r in results)
        files_with_viruses = sum(1 for r in results if len(r.get("viral_hits", [])) > 0)
        
        return {
            "total_viral_hits": total_viral_hits,
            "files_with_viruses": files_with_viruses,
            "viral_prevalence": files_with_viruses / len(results) if results else 0,
            "average_viruses_per_positive_sample": total_viral_hits / files_with_viruses if files_with_viruses > 0 else 0
        }
    
    def _create_species_breakdown(self, results: List[Dict]) -> Dict:
        """Create species-wise breakdown."""
        species_data = {}
        
        for result in results:
            species = result.get("species", "Unknown")
            if species not in species_data:
                species_data[species] = {
                    "total_samples": 0,
                    "samples_with_viruses": 0,
                    "total_viral_hits": 0,
                    "viral_families": set()
                }
            
            species_data[species]["total_samples"] += 1
            viral_hits = result.get("viral_hits", [])
            
            if viral_hits:
                species_data[species]["samples_with_viruses"] += 1
                species_data[species]["total_viral_hits"] += len(viral_hits)
                
                for hit in viral_hits:
                    family = hit.get("viral_family", "Unknown")
                    species_data[species]["viral_families"].add(family)
        
        # Convert sets to lists for JSON serialization
        for species in species_data:
            species_data[species]["viral_families"] = list(species_data[species]["viral_families"])
            species_data[species]["viral_prevalence"] = (
                species_data[species]["samples_with_viruses"] / 
                species_data[species]["total_samples"]
            )
        
        return species_data
    
    def _generate_recommendations(self, results: List[Dict], 
                                superinfection_analysis: Dict) -> List[str]:
        """Generate research recommendations based on results."""
        recommendations = []
        
        total_hits = sum(len(r.get("viral_hits", [])) for r in results)
        
        if total_hits > 1000:
            recommendations.append("Rich viral diversity detected - proceed with detailed phylogenetic analysis")
        
        if superinfection_analysis["global_patterns"]["total_exclusion_pairs"] > 10:
            recommendations.append("Significant superinfection exclusion patterns detected - investigate mechanisms")
        
        if len(self.processing_stats["species_analyzed"]) > 5:
            recommendations.append("Multi-species dataset suitable for comparative virome analysis")
        
        recommendations.extend([
            "Consider temporal analysis if sample collection dates are available",
            "Investigate geographic patterns if location data is available",
            "Validate findings with experimental superinfection studies",
            "Prepare results for publication in virology or vector biology journals"
        ])
        
        return recommendations
    
    def _generate_csv_reports(self, results: List[Dict]) -> None:
        """Generate CSV reports for easy analysis."""
        # Flatten results for CSV
        flattened_data = []
        
        for result in results:
            base_data = {
                "species": result.get("species", "Unknown"),
                "file_type": result.get("file_type", "Unknown"),
                "file_path": result.get("file_path", ""),
                "total_contigs": result.get("total_contigs", 0),
                "viral_hits_count": len(result.get("viral_hits", []))
            }
            
            viral_hits = result.get("viral_hits", [])
            if viral_hits:
                for i, hit in enumerate(viral_hits):
                    row_data = base_data.copy()
                    row_data.update({
                        "viral_hit_index": i,
                        "viral_family": hit.get("viral_family", "Unknown"),
                        "virus_name": hit.get("virus_name", "Unknown"),
                        "e_value": hit.get("e_value", ""),
                        "bit_score": hit.get("bit_score", ""),
                        "identity": hit.get("identity", "")
                    })
                    flattened_data.append(row_data)
            else:
                flattened_data.append(base_data)
        
        # Save to CSV
        df = pd.DataFrame(flattened_data)
        csv_file = self.summaries_dir / "viral_detection_results.csv"
        df.to_csv(csv_file, index=False)
        
        logger.info(f"CSV report saved to: {csv_file}")
    
    def _generate_markdown_report(self, summary: Dict) -> None:
        """Generate markdown report."""
        stats = summary["processing_stats"]
        viral_summary = summary["viral_detection_summary"]
        
        report_content = f"""
# Large-Scale Mosquito Viral Analysis Report

## Analysis Overview
- **Analysis Date**: {summary['analysis_metadata']['analysis_date']}
- **Total Files Processed**: {stats['processed_files']:,}
- **Failed Files**: {stats['failed_files']:,}
- **Processing Success Rate**: {(stats['processed_files'] / stats['total_files'] * 100):.1f}%

## Viral Detection Results
- **Total Viral Hits**: {viral_summary['total_viral_hits']:,}
- **Files with Viruses**: {viral_summary['files_with_viruses']:,}
- **Viral Prevalence**: {viral_summary['viral_prevalence']:.2%}
- **Average Viruses per Positive Sample**: {viral_summary['average_viruses_per_positive_sample']:.1f}

## Species Analysis
- **Species Analyzed**: {len(stats['species_analyzed'])}
- **Viral Families Detected**: {len(stats['viral_families_detected'])}

## Superinfection Exclusion Analysis
- **Exclusion Pairs Detected**: {summary['superinfection_analysis']['global_patterns']['total_exclusion_pairs']}
- **Species with Exclusion Patterns**: {len([s for s in summary['superinfection_analysis']['species_analysis'] if summary['superinfection_analysis']['species_analysis'][s]['exclusion_pairs']])}

## Key Findings
{chr(10).join(f'- {rec}' for rec in summary['recommendations'])}

## Next Steps
1. Validate exclusion patterns with experimental studies
2. Investigate molecular mechanisms of exclusion
3. Prepare findings for publication
4. Apply insights to vector control strategies
"""
        
        report_file = self.summaries_dir / "analysis_report.md"
        with open(report_file, 'w') as f:
            f.write(report_content)
        
        logger.info(f"Markdown report saved to: {report_file}")
    
    def run_full_analysis(self, batch_size: int = 1000) -> None:
        """Run complete large-scale viral analysis pipeline."""
        logger.info("Starting large-scale viral analysis pipeline...")
        
        try:
            # 1. Discover genome files
            genome_files = self.discover_genome_files()
            
            if genome_files["total_count"] == 0:
                logger.error("No genome files found. Please run large_scale_data_collector.py first.")
                return
            
            # 2. Setup viral databases
            logger.info("Setting up viral databases...")
            db_path = self.pipeline.setup_viral_databases()
            logger.info(f"Viral database ready at: {db_path}")
            
            # 3. Process all files
            results = self.process_batch(genome_files, batch_size)
            
            # 4. Analyze superinfection patterns
            superinfection_analysis = self.analyze_superinfection_patterns(results)
            
            # 5. Generate comprehensive report
            self.generate_comprehensive_report(results, superinfection_analysis)
            
            # 6. Final summary
            logger.info("=" * 60)
            logger.info("LARGE-SCALE ANALYSIS COMPLETED")
            logger.info(f"Files Processed: {self.processing_stats['processed_files']:,}")
            logger.info(f"Viral Hits: {self.processing_stats['total_viral_hits']:,}")
            logger.info(f"Species: {len(self.processing_stats['species_analyzed'])}")
            logger.info(f"Exclusion Pairs: {self.processing_stats['superinfection_pairs']}")
            logger.info(f"Results saved to: {self.output_dir}")
            logger.info("=" * 60)
            
        except Exception as e:
            logger.error(f"Error in full analysis pipeline: {e}")
            raise

def main():
    """Main function for batch viral analysis."""
    print("🦟 Batch Viral Analysis for Large-Scale Data")
    print("=" * 50)
    
    # Get user input
    input_dir = input("Input directory with genome data (default 'massive_mosquito_data'): ").strip()
    input_dir = input_dir if input_dir else "massive_mosquito_data"
    
    output_dir = input("Output directory for results (default 'massive_viral_results'): ").strip()
    output_dir = output_dir if output_dir else "massive_viral_results"
    
    batch_size = input("Batch size for processing (default 1000): ").strip()
    batch_size = int(batch_size) if batch_size.isdigit() else 1000
    
    max_workers = input(f"Max parallel workers (default {min(mp.cpu_count(), 16)}): ").strip()
    max_workers = int(max_workers) if max_workers.isdigit() else None
    
    print(f"\nStarting analysis:")
    print(f"- Input: {input_dir}")
    print(f"- Output: {output_dir}")
    print(f"- Batch size: {batch_size:,}")
    print(f"- Workers: {max_workers or min(mp.cpu_count(), 16)}")
    print("\nThis may take several hours for large datasets...")
    
    # Initialize and run analyzer
    analyzer = BatchViralAnalyzer(input_dir, output_dir, max_workers)
    analyzer.run_full_analysis(batch_size)
    
    print("\n✅ Large-scale viral analysis completed!")
    print(f"Check {output_dir}/ for comprehensive results and reports.")

if __name__ == "__main__":
    main()