#!/usr/bin/env python3
"""Test script for the viral interference research pipeline."""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import our modules
from src.data_processing import ViralAnalysisPipeline
from src.statistical_analysis import ViralInteractionAnalyzer
from src.machine_learning import ViralInteractionPredictor
from src.visualization import ViralVisualization
from src.utils import setup_logging, get_logger

def test_complete_pipeline():
    """Test the complete viral interference analysis pipeline."""
    
    # Setup logging
    setup_logging()
    logger = get_logger('test_pipeline')
    
    logger.info("Starting viral interference pipeline test")
    
    # Load sample data
    data_dir = Path('data/sample')
    
    print("\n" + "="*60)
    print("VIRAL INTERFERENCE RESEARCH PIPELINE TEST")
    print("="*60)
    
    print("\n1. Loading sample data...")
    try:
        metadata = pd.read_csv(data_dir / 'sample_metadata.csv')
        abundance_matrix = pd.read_csv(data_dir / 'viral_abundance_matrix.csv', index_col=0)
        
        print(f"   ✓ Loaded {len(metadata)} samples")
        print(f"   ✓ Loaded abundance data for {abundance_matrix.shape[1]} viruses")
        print(f"   ✓ Date range: {metadata['collection_date'].min()} to {metadata['collection_date'].max()}")
        print(f"   ✓ Countries: {', '.join(metadata['country'].unique())}")
        
    except Exception as e:
        print(f"   ✗ Error loading data: {e}")
        return False
    
    # Test data processing pipeline
    print("\n2. Testing data processing pipeline...")
    try:
        pipeline = ViralAnalysisPipeline()
        
        # Test quality control
        qc_results = pipeline.quality_control(
            abundance_matrix, 
            min_prevalence=0.01,
            min_abundance=10
        )
        
        print(f"   ✓ Quality control completed")
        print(f"   ✓ Filtered to {qc_results['filtered_matrix'].shape[1]} viruses")
        print(f"   ✓ Removed {qc_results['n_low_prevalence']} low-prevalence viruses")
        print(f"   ✓ Removed {qc_results['n_low_abundance']} low-abundance viruses")
        
        filtered_abundance = qc_results['filtered_matrix']
        
    except Exception as e:
        print(f"   ✗ Error in data processing: {e}")
        return False
    
    # Test statistical analysis
    print("\n3. Testing statistical analysis...")
    try:
        analyzer = ViralInteractionAnalyzer()
        
        # Convert to presence/absence
        presence_matrix = analyzer.convert_to_presence_absence(filtered_abundance)
        print(f"   ✓ Converted to presence/absence matrix")
        
        # Pairwise analysis
        pairwise_results = analyzer.pairwise_fisher_test(presence_matrix)
        print(f"   ✓ Performed {len(pairwise_results)} pairwise tests")
        
        # Multiple testing correction
        corrected_results = analyzer.multiple_testing_correction(
            pairwise_results, method='fdr_bh'
        )
        print(f"   ✓ Applied multiple testing correction")
        
        # Find exclusion candidates
        exclusion_candidates = analyzer.identify_exclusion_candidates(
            corrected_results, significance_threshold=0.05
        )
        print(f"   ✓ Identified {len(exclusion_candidates)} exclusion candidates")
        
        # Build interaction network
        network = analyzer.build_interaction_network(
            corrected_results, significance_threshold=0.05
        )
        print(f"   ✓ Built interaction network with {network.number_of_nodes()} nodes and {network.number_of_edges()} edges")
        
        if len(exclusion_candidates) > 0:
            print(f"   ✓ Top exclusion pair: {exclusion_candidates[0]['virus_pair']} (strength: {exclusion_candidates[0]['exclusion_strength']:.3f})")
        
    except Exception as e:
        print(f"   ✗ Error in statistical analysis: {e}")
        return False
    
    # Test machine learning
    print("\n4. Testing machine learning models...")
    try:
        predictor = ViralInteractionPredictor()
        
        # Prepare features
        features = predictor.prepare_interaction_features(
            filtered_abundance, metadata
        )
        print(f"   ✓ Prepared {features.shape[1]} features for {features.shape[0]} samples")
        
        # Clustering analysis
        clustering_results = predictor.cluster_viral_profiles(
            filtered_abundance, n_clusters=3
        )
        print(f"   ✓ Performed clustering analysis")
        print(f"   ✓ Identified {len(set(clustering_results['cluster_labels']))} clusters")
        
        # Train exclusion classifier (if we have enough exclusion pairs)
        if len(exclusion_candidates) >= 5:
            classifier_results = predictor.train_exclusion_classifier(
                filtered_abundance, exclusion_candidates[:5]
            )
            print(f"   ✓ Trained exclusion classifier")
            print(f"   ✓ Cross-validation accuracy: {classifier_results['cv_scores'].mean():.3f} ± {classifier_results['cv_scores'].std():.3f}")
        else:
            print(f"   ⚠ Skipped classifier training (insufficient exclusion pairs)")
        
    except Exception as e:
        print(f"   ✗ Error in machine learning: {e}")
        return False
    
    # Test visualization
    print("\n5. Testing visualization...")
    try:
        viz = ViralVisualization()
        
        # Create abundance heatmap
        heatmap_fig = viz.plot_abundance_heatmap(
            filtered_abundance.iloc[:50, :10],  # Subset for speed
            interactive=False
        )
        print(f"   ✓ Created abundance heatmap")
        
        # Create co-occurrence matrix
        if len(corrected_results) > 0:
            cooccurrence_fig = viz.plot_cooccurrence_matrix(
                corrected_results, interactive=False
            )
            print(f"   ✓ Created co-occurrence matrix")
        
        # Create interaction network
        if network.number_of_nodes() > 0:
            network_fig = viz.plot_interaction_network(
                network, interactive=False
            )
            print(f"   ✓ Created interaction network plot")
        
        # Create prevalence comparison
        prevalence_fig = viz.plot_prevalence_comparison(
            filtered_abundance, metadata, group_by='mosquito_species', interactive=False
        )
        print(f"   ✓ Created prevalence comparison")
        
        # Create clustering visualization
        if 'dimensionality_reduction' in clustering_results:
            clustering_fig = viz.plot_clustering_results(
                clustering_results, filtered_abundance, interactive=False
            )
            print(f"   ✓ Created clustering visualization")
        
    except Exception as e:
        print(f"   ✗ Error in visualization: {e}")
        return False
    
    # Summary statistics
    print("\n6. Analysis Summary")
    print("   " + "-"*40)
    
    # Viral prevalence
    prevalence = (filtered_abundance > 0).mean().sort_values(ascending=False)
    print(f"   Most prevalent viruses:")
    for virus, prev in prevalence.head(5).items():
        print(f"     • {virus}: {prev:.1%}")
    
    # Significant interactions
    significant_interactions = corrected_results[corrected_results['p_adj'] < 0.05]
    if len(significant_interactions) > 0:
        print(f"   \n   Significant viral interactions: {len(significant_interactions)}")
        
        # Exclusions
        exclusions = significant_interactions[significant_interactions['odds_ratio'] < 1]
        if len(exclusions) > 0:
            print(f"     • Exclusions: {len(exclusions)}")
            strongest_exclusion = exclusions.loc[exclusions['odds_ratio'].idxmin()]
            print(f"     • Strongest: {strongest_exclusion['virus_a']} ⟷ {strongest_exclusion['virus_b']} (OR: {strongest_exclusion['odds_ratio']:.3f})")
        
        # Co-occurrences
        cooccurrences = significant_interactions[significant_interactions['odds_ratio'] > 1]
        if len(cooccurrences) > 0:
            print(f"     • Co-occurrences: {len(cooccurrences)}")
            strongest_cooccurrence = cooccurrences.loc[cooccurrences['odds_ratio'].idxmax()]
            print(f"     • Strongest: {strongest_cooccurrence['virus_a']} ⟷ {strongest_cooccurrence['virus_b']} (OR: {strongest_cooccurrence['odds_ratio']:.3f})")
    else:
        print(f"   \n   No significant viral interactions found")
    
    # Sample diversity
    viral_richness = (filtered_abundance > 0).sum(axis=1)
    print(f"   \n   Sample viral diversity:")
    print(f"     • Mean richness: {viral_richness.mean():.1f} viruses/sample")
    print(f"     • Range: {viral_richness.min()}-{viral_richness.max()} viruses/sample")
    
    # Geographic patterns
    country_prevalence = {}
    for country in metadata['country'].unique():
        country_samples = metadata[metadata['country'] == country]['sample_id']
        country_abundance = filtered_abundance.loc[filtered_abundance.index.isin(country_samples)]
        country_richness = (country_abundance > 0).sum(axis=1).mean()
        country_prevalence[country] = country_richness
    
    print(f"   \n   Geographic patterns (mean viral richness):")
    for country, richness in sorted(country_prevalence.items(), key=lambda x: x[1], reverse=True):
        print(f"     • {country}: {richness:.1f} viruses/sample")
    
    print("\n" + "="*60)
    print("PIPELINE TEST COMPLETED SUCCESSFULLY! ✓")
    print("="*60)
    
    print("\nNext steps:")
    print("1. Run the Streamlit app: streamlit run app/main.py")
    print("2. Upload your own mosquito virome data")
    print("3. Explore interactive visualizations")
    print("4. Export results for publication")
    
    logger.info("Pipeline test completed successfully")
    return True

def test_streamlit_compatibility():
    """Test that all components work with Streamlit."""
    print("\n" + "="*60)
    print("STREAMLIT COMPATIBILITY TEST")
    print("="*60)
    
    try:
        import streamlit as st
        print("   ✓ Streamlit import successful")
    except ImportError:
        print("   ⚠ Streamlit not installed (install with: pip install streamlit)")
        return False
    
    try:
        # Test plotly imports
        import plotly.express as px
        import plotly.graph_objects as go
        print("   ✓ Plotly imports successful")
        
        # Test other required packages
        import networkx as nx
        import seaborn as sns
        import matplotlib.pyplot as plt
        print("   ✓ All visualization dependencies available")
        
        return True
        
    except ImportError as e:
        print(f"   ✗ Missing dependency: {e}")
        return False

if __name__ == '__main__':
    print("Viral Interference Research System - Pipeline Test")
    print("This script tests all major components of the analysis pipeline.\n")
    
    # Run main pipeline test
    success = test_complete_pipeline()
    
    if success:
        # Test Streamlit compatibility
        streamlit_ok = test_streamlit_compatibility()
        
        if streamlit_ok:
            print("\n🎉 All tests passed! The system is ready for use.")
        else:
            print("\n⚠ Pipeline works but some Streamlit dependencies may be missing.")
            print("Install missing packages with: pip install -r requirements.txt")
    else:
        print("\n❌ Pipeline test failed. Please check the error messages above.")
        sys.exit(1)