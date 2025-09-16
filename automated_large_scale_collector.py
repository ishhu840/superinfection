#!/usr/bin/env python3
"""
Automated Large-Scale Mosquito Genome Data Collector

This script automatically collects the maximum number of mosquito genome files
from multiple public sources without requiring user interaction.
"""

import os
import sys
import json
import time
import logging
import requests
import subprocess
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Optional

# Import the main collector class
from large_scale_data_collector import LargeScaleDataCollector

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('automated_large_scale_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def main():
    """
    Automated main function to run large-scale data collection.
    Downloads maximum available mosquito genomes without user interaction.
    """
    print("🦟 Automated Large-Scale Mosquito Genome Data Collector")
    print("=" * 60)
    print("Downloading MAXIMUM available mosquito genome files for viral analysis.")
    print("")
    
    # Set maximum parameters for comprehensive collection
    max_assemblies = 50000  # Maximum NCBI assemblies
    max_sra = 20000        # Increased SRA datasets
    output_dir = "massive_mosquito_genomes"
    
    print(f"Collection Parameters:")
    print(f"- Reference genomes: {max_assemblies:,}")
    print(f"- SRA datasets: {max_sra:,}")
    print(f"- Output directory: {output_dir}")
    print(f"- Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\nThis will take several hours to complete...")
    print("Collection will run automatically without further input.")
    print("=" * 60)
    
    try:
        # Initialize collector
        collector = LargeScaleDataCollector(output_dir)
        
        # Start automated collection
        logger.info("Starting automated large-scale mosquito genome collection")
        collector.collect_all_data(max_assemblies, max_sra)
        
        print("\n✅ COLLECTION COMPLETED SUCCESSFULLY!")
        print(f"Check {output_dir}/ for your massive mosquito genome dataset.")
        print(f"Total collection time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("You can now run viral detection analysis on this large-scale dataset.")
        
        # Generate final summary
        summary_file = Path(output_dir) / "COLLECTION_SUMMARY.md"
        with open(summary_file, 'w') as f:
            f.write(f"# Automated Large-Scale Mosquito Genome Collection Summary\n\n")
            f.write(f"**Collection Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Target Assemblies:** {max_assemblies:,}\n")
            f.write(f"**Target SRA Datasets:** {max_sra:,}\n")
            f.write(f"**Output Directory:** {output_dir}\n\n")
            f.write(f"## Collection Status\n")
            f.write(f"- ✅ NCBI RefSeq/GenBank assemblies\n")
            f.write(f"- ✅ SRA metagenomic datasets\n")
            f.write(f"- ✅ Published virome studies\n")
            f.write(f"- ✅ Quality control and organization\n\n")
            f.write(f"## Next Steps\n")
            f.write(f"1. Run viral detection analysis: `python3 batch_viral_analysis.py`\n")
            f.write(f"2. Update dashboard with new statistics\n")
            f.write(f"3. Begin large-scale viral superinfection analysis\n")
        
        logger.info(f"Collection summary saved to: {summary_file}")
        
    except Exception as e:
        logger.error(f"Error during automated collection: {e}")
        print(f"\n❌ Collection failed: {e}")
        print("Check the log file for detailed error information.")
        raise

if __name__ == "__main__":
    main()