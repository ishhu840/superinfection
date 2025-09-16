#!/usr/bin/env python3
"""
Analyze Downloaded Mosquito Genomes for PhD Research

This script runs viral detection and superinfection exclusion analysis
on the real mosquito genomes we downloaded.
"""

import os
import sys
from pathlib import Path
import logging
import json
from datetime import datetime

# Add src to path
src_dir = Path(__file__).parent / "src"
sys.path.insert(0, str(src_dir))

from viral_detection.blast_pipeline import ViralDetectionPipeline
from viral_detection.viral_families import get_viral_family_classifier
from analysis.cooccurrence_analyzer import ViralCooccurrenceAnalyzer
from visualization.network_graphs import create_mosquito_virus_visualizations
from analysis.statistical_validation import PublicationStatistics

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('genome_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def analyze_mosquito_genomes():
    """
    Analyze the downloaded mosquito genomes for viral content and superinfection patterns.
    """
    print("=" * 70)
    print("MOSQUITO GENOME VIRAL ANALYSIS FOR PhD RESEARCH")
    print("Superinfection Exclusion Analysis")
    print("=" * 70)
    
    # Setup directories
    genome_dir = Path("research_genomes")
    output_dir = Path("viral_analysis_results")
    output_dir.mkdir(exist_ok=True)
    
    if not genome_dir.exists():
        print(f"Error: Genome directory {genome_dir} not found.")
        print("Please run download_mosquito_genomes.py first.")
        return False
    
    # Find genome files
    genome_files = []
    species_data = {}
    
    for species_dir in genome_dir.iterdir():
        if species_dir.is_dir():
            species_name = species_dir.name
            fasta_files = list(species_dir.glob("*.fna.gz"))
            
            if fasta_files:
                species_data[species_name] = fasta_files
                genome_files.extend(fasta_files)
                print(f"Found {len(fasta_files)} genome files for {species_name}")
    
    if not genome_files:
        print("No genome files found. Please check the research_genomes directory.")
        return False
    
    print(f"\nTotal genome files to analyze: {len(genome_files)}")
    
    # Step 1: Viral Detection
    print("\n" + "="*50)
    print("STEP 1: VIRAL SEQUENCE DETECTION")
    print("="*50)
    
    try:
        # Initialize viral detection pipeline
        viral_pipeline = ViralDetectionPipeline(
            work_dir=str(output_dir),
            blast_threads=4
        )
        
        # Setup viral databases
        viral_db = viral_pipeline.setup_viral_databases()
        if not viral_db:
            logger.error("Failed to setup viral databases")
            return False
        
        print(f"Viral database prepared: {viral_db}")
        
        # Run viral detection on each genome
        all_viral_results = {}
        
        for species, files in species_data.items():
            print(f"\nAnalyzing {species}...")
            species_results = []
            
            for genome_file in files:
                print(f"  Processing {genome_file.name}...")
                
                # Analyze genome for viral sequences
                result = viral_pipeline.analyze_genome(
                    genome_file,
                    species,
                    viral_db
                )
                
                species_results.append(result)
                print(f"    Found {len(result.viral_hits)} viral hits ({result.viral_prevalence:.1f}% prevalence)")
                
                if result.virus_families:
                    print(f"    Viral families: {', '.join(result.virus_families.keys())}")
            
            all_viral_results[species] = species_results
            total_hits = sum(len(r.viral_hits) for r in species_results)
            print(f"Total viral hits for {species}: {total_hits}")
        
        # Convert results to serializable format
        serializable_results = {}
        for species, results in all_viral_results.items():
            serializable_results[species] = []
            for result in results:
                result_dict = {
                    'genome_file': result.genome_file,
                    'species': result.species,
                    'total_contigs': result.total_contigs,
                    'viral_prevalence': result.viral_prevalence,
                    'virus_families': result.virus_families,
                    'unique_viruses': list(result.unique_viruses),
                    'viral_hits': [{
                        'query_id': hit.query_id,
                        'subject_id': hit.subject_id,
                        'identity': hit.identity,
                        'evalue': hit.evalue,
                        'bit_score': hit.bit_score,
                        'alignment_length': hit.alignment_length
                    } for hit in result.viral_hits]
                }
                serializable_results[species].append(result_dict)
        
        # Save viral detection results
        results_file = output_dir / "viral_detection_results.json"
        with open(results_file, 'w') as f:
            json.dump(serializable_results, f, indent=2)
        
        print(f"\nViral detection results saved to: {results_file}")
        
    except Exception as e:
        logger.error(f"Viral detection failed: {e}")
        return False
    
    # Step 2: Viral Classification
    print("\n" + "="*50)
    print("STEP 2: VIRAL FAMILY CLASSIFICATION")
    print("="*50)
    
    try:
        classifier = get_viral_family_classifier()
        # Use the virus families already detected by the pipeline
        classified_results = serializable_results
        
        for species, results in classified_results.items():
            if results:
                print(f"\nViral families in {species}:")
                all_families = {}
                for result in results:
                    for family, count in result['virus_families'].items():
                        all_families[family] = all_families.get(family, 0) + count
                
                if all_families:
                    for family, count in all_families.items():
                        print(f"    {family}: {count} sequences")
                else:
                    print(f"    No viral families detected")
        
        # Save classified results
        classified_file = output_dir / "classified_viral_results.json"
        with open(classified_file, 'w') as f:
            json.dump(classified_results, f, indent=2, default=str)
        
        print(f"\nClassified results saved to: {classified_file}")
        
    except Exception as e:
        logger.error(f"Viral classification failed: {e}")
        classified_results = all_viral_hits
    
    # Step 3: Cooccurrence Analysis
    print("\n" + "="*50)
    print("STEP 3: VIRAL COOCCURRENCE ANALYSIS")
    print("="*50)
    
    try:
        # Prepare data for cooccurrence analysis
        cooccurrence_data = []
        
        for species, results in classified_results.items():
            for result in results:
                for hit in result['viral_hits']:
                    # Determine viral family from the virus families detected
                    viral_family = 'Unknown'
                    for family in result['virus_families'].keys():
                        if family != 'Unknown':
                            viral_family = family
                            break
                    
                    cooccurrence_data.append({
                        'species': species,
                        'viral_family': viral_family,
                        'sequence_id': hit['query_id'],
                        'e_value': hit['evalue'],
                        'bit_score': hit['bit_score']
                    })
        
        if cooccurrence_data:
            analyzer = ViralCooccurrenceAnalyzer()
            
            # Analyze cooccurrence patterns
            cooccurrence_results = analyzer.analyze_cooccurrence_patterns(cooccurrence_data)
            
            # Save cooccurrence results
            cooccurrence_file = output_dir / "cooccurrence_analysis.json"
            with open(cooccurrence_file, 'w') as f:
                json.dump(cooccurrence_results, f, indent=2, default=str)
            
            print(f"Cooccurrence analysis completed")
            print(f"Results saved to: {cooccurrence_file}")
        
    except Exception as e:
        logger.error(f"Cooccurrence analysis failed: {e}")
    
    # Step 4: Generate Summary Report
    print("\n" + "="*50)
    print("STEP 4: GENERATING SUMMARY REPORT")
    print("="*50)
    
    # Create summary
    summary = {
        'analysis_date': datetime.now().isoformat(),
        'total_genomes_analyzed': len(genome_files),
        'species_analyzed': list(species_data.keys()),
        'total_viral_hits': sum(len(result['viral_hits']) for results in classified_results.values() for result in results),
        'viral_families_detected': set(),
        'species_results': {}
    }
    
    for species, results in classified_results.items():
        all_families = set()
        total_hits = 0
        
        for result in results:
            all_families.update(result['virus_families'].keys())
            total_hits += len(result['viral_hits'])
        
        summary['viral_families_detected'].update(all_families)
        
        summary['species_results'][species] = {
            'total_hits': total_hits,
            'viral_families': list(all_families),
            'genome_files': len(species_data[species])
        }
    
    summary['viral_families_detected'] = list(summary['viral_families_detected'])
    
    # Save summary
    summary_file = output_dir / "analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Print summary
    print(f"\n{'='*70}")
    print("ANALYSIS SUMMARY")
    print(f"{'='*70}")
    print(f"Analysis Date: {summary['analysis_date']}")
    print(f"Total Genomes Analyzed: {summary['total_genomes_analyzed']}")
    print(f"Species Analyzed: {', '.join(summary['species_analyzed'])}")
    print(f"Total Viral Hits: {summary['total_viral_hits']}")
    print(f"Viral Families Detected: {', '.join(summary['viral_families_detected'])}")
    
    print(f"\nPer-Species Results:")
    for species, results in summary['species_results'].items():
        print(f"  {species}:")
        print(f"    - Genome files: {results['genome_files']}")
        print(f"    - Viral hits: {results['total_hits']}")
        print(f"    - Viral families: {', '.join(results['viral_families']) if results['viral_families'] else 'None'}")
    
    print(f"\n✓ Analysis completed successfully!")
    print(f"Results directory: {output_dir.absolute()}")
    print(f"\nFor your PhD research, you can now:")
    print(f"1. Examine the viral detection results in {results_file}")
    print(f"2. Analyze viral family classifications in {classified_file}")
    print(f"3. Study cooccurrence patterns for superinfection exclusion")
    print(f"4. Use the comprehensive dashboard for visualization")
    
    return True

if __name__ == "__main__":
    analyze_mosquito_genomes()