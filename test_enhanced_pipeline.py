#!/usr/bin/env python3

from enhanced_viral_pipeline import EnhancedViralDetectionPipeline
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

def main():
    # Initialize enhanced pipeline
    pipeline = EnhancedViralDetectionPipeline()
    
    # Analyze the mosquito collection
    print("Running enhanced viral detection on mosquito collection...")
    results = pipeline.analyze_mosquito_collection(Path('real_mosquito_data'))
    
    print(f"\nAnalyzed {len(results)} samples")
    print("\nSample Results:")
    print("-" * 80)
    
    for result in results:
        print(f"Sample: {Path(result.genome_file).name}")
        print(f"  Species: {result.species}")
        print(f"  Viral hits: {len(result.filtered_viral_hits)}")
        print(f"  Unique viruses: {len(result.unique_viruses)}")
        print(f"  Viral diversity: {result.viral_diversity_index:.3f}")
        print(f"  Host contamination rate: {result.host_contamination_rate:.1f}%")
        if result.unique_viruses:
            print(f"  Detected viruses: {', '.join(result.unique_viruses)}")
        print()
    
    # Generate enhanced report
    report = pipeline.generate_enhanced_report(results)
    
    # Calculate unique viruses across all samples
    all_unique_viruses = set()
    total_diversity = 0
    samples_with_diversity = 0
    
    for result in results:
        all_unique_viruses.update(result.unique_viruses)
        if result.viral_diversity_index > 0:
            total_diversity += result.viral_diversity_index
            samples_with_diversity += 1
    
    avg_diversity = total_diversity / samples_with_diversity if samples_with_diversity > 0 else 0
    
    print("Enhanced Analysis Summary:")
    print("-" * 40)
    print(f"Total samples: {report['summary']['total_samples']}")
    print(f"Samples with viral hits: {report['summary']['samples_with_viruses']}")
    print(f"Viral detection rate: {report['summary']['viral_detection_rate']:.1f}%")
    print(f"Total viral hits: {report['summary']['total_viral_hits']}")
    print(f"Total raw hits: {report['summary']['total_raw_hits']}")
    print(f"Unique viruses detected: {len(all_unique_viruses)}")
    print(f"Average viral diversity: {avg_diversity:.3f}")
    print(f"Pipeline improvement factor: {report['summary']['pipeline_improvement_factor']:.2f}")
    print(f"Average host contamination rate: {report['summary']['average_host_contamination_rate']:.1f}%")
    
    if all_unique_viruses:
        print(f"\nDetected viruses: {', '.join(sorted(all_unique_viruses))}")
    
if __name__ == "__main__":
    main()