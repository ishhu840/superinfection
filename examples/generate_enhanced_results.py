#!/usr/bin/env python3

from enhanced_viral_pipeline import EnhancedViralDetectionPipeline
from pathlib import Path
import json
import logging
from collections import defaultdict, Counter
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO)

def generate_enhanced_dashboard_data():
    """Generate enhanced results data for the dashboard"""
    
    # Initialize enhanced pipeline
    pipeline = EnhancedViralDetectionPipeline()
    
    # Analyze the mosquito collection
    print("Running enhanced viral detection on mosquito collection...")
    results = pipeline.analyze_mosquito_collection(Path('real_mosquito_data'))
    
    # Generate enhanced report
    report = pipeline.generate_enhanced_report(results)
    
    # Calculate additional metrics for dashboard
    all_unique_viruses = set()
    species_distribution = defaultdict(int)
    viral_detection_by_species = defaultdict(lambda: {'samples': 0, 'with_viruses': 0})
    
    viral_results = []
    samples_with_viruses = []
    
    for result in results:
        all_unique_viruses.update(result.unique_viruses)
        species_distribution[result.species] += 1
        
        species_stats = viral_detection_by_species[result.species]
        species_stats['samples'] += 1
        
        if result.filtered_viral_hits:
            species_stats['with_viruses'] += 1
            samples_with_viruses.append(result.genome_file)
            
            for hit in result.filtered_viral_hits:
                viral_results.append({
                    'sample_id': Path(result.genome_file).stem,
                    'species': result.species,
                    'virus_id': hit.subject_id,
                    'virus_family': 'Unknown',  # Would need virus family mapping
                    'identity': hit.identity,
                    'evalue': hit.evalue,
                    'alignment_length': hit.alignment_length
                })
    
    # Calculate cooccurrence patterns (simplified for now)
    cooccurrence_pairs = {}
    cooccurrence_probabilities = {}
    exclusion_pairs = {}
    
    # For samples with multiple viruses, find cooccurrence
    for result in results:
        if len(result.unique_viruses) > 1:
            viruses = sorted(list(result.unique_viruses))
            for i in range(len(viruses)):
                for j in range(i+1, len(viruses)):
                    pair_key = f"{viruses[i]} + {viruses[j]}"
                    cooccurrence_pairs[pair_key] = cooccurrence_pairs.get(pair_key, 0) + 1
    
    # Calculate probabilities
    total_samples_with_viruses = len(samples_with_viruses)
    for pair, count in cooccurrence_pairs.items():
        cooccurrence_probabilities[pair] = count / total_samples_with_viruses if total_samples_with_viruses > 0 else 0
    
    # Calculate superinfection rate
    samples_with_multiple_viruses = len([r for r in results if len(r.unique_viruses) > 1])
    superinfection_rate = samples_with_multiple_viruses / total_samples_with_viruses if total_samples_with_viruses > 0 else 0
    
    # Create dashboard data structure
    dashboard_data = {
        'analysis_metadata': {
            'analysis_date': datetime.now().isoformat(),
            'pipeline_version': 'Enhanced v2.0 with Reference Genome Integration',
            'total_samples': len(results),
            'total_viral_hits': sum(len(r.filtered_viral_hits) for r in results),
            'samples_with_viruses': len(samples_with_viruses),
            'unique_viruses_detected': len(all_unique_viruses),
            'viral_detection_rate': len(samples_with_viruses) / len(results) * 100
        },
        'sample_summary': {
            'species_distribution': dict(species_distribution),
            'viral_detection_by_species': {
                species: {
                    'total_samples': stats['samples'],
                    'samples_with_viruses': stats['with_viruses'],
                    'detection_rate': stats['with_viruses'] / stats['samples'] * 100 if stats['samples'] > 0 else 0
                }
                for species, stats in viral_detection_by_species.items()
            }
        },
        'viral_detection_results': viral_results,
        'cooccurrence_analysis': {
            'superinfection_rate': superinfection_rate,
            'samples_with_multiple_viruses': samples_with_multiple_viruses,
            'cooccurrence_pairs': cooccurrence_pairs,
            'cooccurrence_probabilities': cooccurrence_probabilities,
            'exclusion_pairs': exclusion_pairs
        },
        'enhanced_metrics': {
            'average_viral_diversity': sum(r.viral_diversity_index for r in results) / len(results),
            'average_host_contamination_rate': sum(r.host_contamination_rate for r in results) / len(results),
            'pipeline_improvement_factor': report['summary']['pipeline_improvement_factor'],
            'unique_viruses_list': sorted(list(all_unique_viruses))
        },
        'detailed_results': report['detailed_results']
    }
    
    return dashboard_data

def main():
    # Generate enhanced results
    print("Generating enhanced dashboard data...")
    dashboard_data = generate_enhanced_dashboard_data()
    
    # Save results
    output_dir = Path("real_viral_results")
    output_dir.mkdir(exist_ok=True)
    
    output_file = output_dir / "enhanced_viral_analysis_results.json"
    with open(output_file, 'w') as f:
        json.dump(dashboard_data, f, indent=2, default=str)
    
    print(f"\nEnhanced results saved to: {output_file}")
    
    # Print summary
    print("\nEnhanced Analysis Summary:")
    print("-" * 40)
    print(f"Total samples: {dashboard_data['analysis_metadata']['total_samples']}")
    print(f"Samples with viral hits: {dashboard_data['analysis_metadata']['samples_with_viruses']}")
    print(f"Viral detection rate: {dashboard_data['analysis_metadata']['viral_detection_rate']:.1f}%")
    print(f"Total viral hits: {dashboard_data['analysis_metadata']['total_viral_hits']}")
    print(f"Unique viruses detected: {dashboard_data['analysis_metadata']['unique_viruses_detected']}")
    print(f"Superinfection rate: {dashboard_data['cooccurrence_analysis']['superinfection_rate']:.1%}")
    print(f"Average viral diversity: {dashboard_data['enhanced_metrics']['average_viral_diversity']:.3f}")
    
if __name__ == "__main__":
    main()