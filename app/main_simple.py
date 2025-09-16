"""Simplified Streamlit app for viral interference research - works with basic packages only."""

import streamlit as st
import pandas as pd
import numpy as np
import os
import sys
import json
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent.parent / 'src'))

def load_sample_data():
    """Load sample data if available."""
    try:
        data_dir = Path(__file__).parent.parent / 'data' / 'sample'
        
        # Load metadata
        metadata_path = data_dir / 'sample_metadata.csv'
        abundance_path = data_dir / 'viral_abundance_matrix.csv'
        summary_path = data_dir / 'dataset_summary.json'
        
        if metadata_path.exists() and abundance_path.exists():
            metadata = pd.read_csv(metadata_path)
            abundance = pd.read_csv(abundance_path)
            
            # Set sample_id as index for both dataframes
            metadata = metadata.set_index('sample_id')
            abundance = abundance.set_index('sample_id')
            
            # Find common sample IDs
            common_samples = metadata.index.intersection(abundance.index)
            metadata = metadata.loc[common_samples]
            abundance = abundance.loc[common_samples]
            
            summary = {}
            if summary_path.exists():
                with open(summary_path, 'r') as f:
                    summary = json.load(f)
            
            return metadata, abundance, summary
    except Exception as e:
        st.error(f"Error loading sample data: {e}")
    
    return None, None, {}

def show_home_page():
    """Display home page."""
    st.title("🦟 Viral Interference Research Platform")
    
    st.markdown("""
    ## Welcome to the Mosquito Virome Analysis Tool
    
    This platform helps researchers analyze viral interactions and interference patterns 
    in mosquito populations using statistical and machine learning approaches.
    
    ### Key Features:
    - **Data Upload & Processing**: Import mosquito virome datasets
    - **Statistical Analysis**: Analyze viral co-occurrence and exclusion patterns
    - **Visualization**: Interactive plots and network visualizations
    - **Machine Learning**: Predict viral interactions and clustering
    
    ### Getting Started:
    1. Navigate to **Data Upload** to load your dataset
    2. Use **Statistical Analysis** to explore viral patterns
    3. View **Visualization** for interactive plots
    
    ### Sample Data Available:
    We've generated sample mosquito virome data for demonstration purposes.
    """)
    
    # Show sample data info if available
    metadata, abundance, summary = load_sample_data()
    if metadata is not None:
        st.success("✅ Sample data loaded successfully!")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Samples", len(metadata))
        with col2:
            st.metric("Viral Species", abundance.shape[1] if abundance is not None else 0)
        with col3:
            st.metric("Countries", len(metadata['country'].unique()) if 'country' in metadata.columns else 0)
        
        if summary:
            st.json(summary)
    else:
        st.warning("⚠️ No sample data found. Generate sample data first.")

def show_data_page():
    """Display data upload and exploration page."""
    st.title("📊 Data Upload & Exploration")
    
    # Load sample data
    metadata, abundance, summary = load_sample_data()
    
    if metadata is not None and abundance is not None:
        st.success("Sample data loaded successfully!")
        
        # Show data overview
        st.subheader("Dataset Overview")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**Metadata Sample:**")
            st.dataframe(metadata.head())
            
            st.write("**Metadata Info:**")
            st.write(f"- Samples: {len(metadata)}")
            st.write(f"- Countries: {', '.join(metadata['country'].unique())}")
            st.write(f"- Species: {', '.join(metadata['mosquito_species'].unique())}")
            
        with col2:
            st.write("**Viral Abundance Matrix:**")
            st.dataframe(abundance.head())
            
            st.write("**Abundance Info:**")
            st.write(f"- Viral species: {abundance.shape[1]}")
            st.write(f"- Non-zero entries: {(abundance > 0).sum().sum()}")
            st.write(f"- Sparsity: {((abundance == 0).sum().sum() / abundance.size * 100):.1f}%")
        
        # Basic statistics
        st.subheader("Basic Statistics")
        
        # Viral prevalence
        prevalence = (abundance > 0).mean().sort_values(ascending=False)
        st.write("**Top 10 Most Prevalent Viruses:**")
        st.bar_chart(prevalence.head(10))
        
        # Sample viral richness
        richness = (abundance > 0).sum(axis=1)
        st.write("**Viral Richness Distribution:**")
        richness_counts = richness.value_counts().sort_index()
        st.bar_chart(richness_counts)
        
    else:
        st.warning("No data available. Please generate sample data first.")
        
        if st.button("Generate Sample Data"):
            try:
                # Run sample data generation
                import subprocess
                result = subprocess.run(
                    ['python3', 'data/sample/generate_sample_data.py'],
                    cwd=Path(__file__).parent.parent,
                    capture_output=True,
                    text=True
                )
                
                if result.returncode == 0:
                    st.success("Sample data generated successfully!")
                    st.rerun()
                else:
                    st.error(f"Error generating sample data: {result.stderr}")
            except Exception as e:
                st.error(f"Error: {e}")

def show_analysis_page():
    """Display statistical analysis page."""
    st.title("📈 Statistical Analysis")
    
    metadata, abundance, summary = load_sample_data()
    
    if metadata is not None and abundance is not None:
        st.subheader("Viral Co-occurrence Analysis")
        
        # Convert to presence/absence
        presence = (abundance > 0).astype(int)
        
        # Calculate co-occurrence matrix
        cooccurrence = presence.T.dot(presence)
        np.fill_diagonal(cooccurrence.values, 0)  # Remove self-interactions
        
        st.write("**Co-occurrence Matrix (Top 10x10):**")
        top_viruses = presence.sum().nlargest(10).index
        cooc_subset = cooccurrence.loc[top_viruses, top_viruses]
        st.dataframe(cooc_subset)
        
        # Simple exclusion analysis
        st.subheader("Potential Viral Exclusions")
        
        exclusions = []
        for i, virus1 in enumerate(top_viruses):
            for virus2 in top_viruses[i+1:]:
                # Count co-occurrences
                both = ((presence[virus1] == 1) & (presence[virus2] == 1)).sum()
                virus1_only = ((presence[virus1] == 1) & (presence[virus2] == 0)).sum()
                virus2_only = ((presence[virus1] == 0) & (presence[virus2] == 1)).sum()
                neither = ((presence[virus1] == 0) & (presence[virus2] == 0)).sum()
                
                # Simple exclusion score (low co-occurrence relative to individual prevalence)
                prev1 = presence[virus1].mean()
                prev2 = presence[virus2].mean()
                expected_cooc = prev1 * prev2 * len(presence)
                observed_cooc = both
                
                if expected_cooc > 0:
                    exclusion_score = 1 - (observed_cooc / expected_cooc)
                    if exclusion_score > 0.5:  # Threshold for potential exclusion
                        exclusions.append({
                            'Virus 1': virus1,
                            'Virus 2': virus2,
                            'Expected Co-occurrence': f"{expected_cooc:.1f}",
                            'Observed Co-occurrence': observed_cooc,
                            'Exclusion Score': f"{exclusion_score:.3f}"
                        })
        
        if exclusions:
            exclusions_df = pd.DataFrame(exclusions)
            st.dataframe(exclusions_df)
        else:
            st.info("No strong exclusion patterns detected in top viruses.")
        
        # Geographic patterns
        st.subheader("Geographic Patterns")
        
        if 'country' in metadata.columns:
            # Calculate richness for each sample
            sample_richness = presence.sum(axis=1)
            # Merge with metadata using proper index alignment
            metadata_with_richness = metadata.copy()
            metadata_with_richness = metadata_with_richness.join(sample_richness.rename('richness'))
            country_richness = metadata_with_richness.groupby('country')['richness'].mean()
            st.write("**Average Viral Richness by Country:**")
            st.bar_chart(country_richness)
        
        # Species patterns
        if 'mosquito_species' in metadata.columns:
            # Calculate richness for each sample
            sample_richness = presence.sum(axis=1)
            # Merge with metadata using proper index alignment
            metadata_with_richness = metadata.copy()
            metadata_with_richness = metadata_with_richness.join(sample_richness.rename('richness'))
            species_richness = metadata_with_richness.groupby('mosquito_species')['richness'].mean()
            st.write("**Average Viral Richness by Mosquito Species:**")
            st.bar_chart(species_richness)
    
    else:
        st.warning("No data available for analysis. Please load data first.")

def show_visualization_page():
    """Display visualization page."""
    st.title("📊 Data Visualization")
    
    metadata, abundance, summary = load_sample_data()
    
    if metadata is not None and abundance is not None:
        # Viral prevalence chart
        st.subheader("Viral Prevalence")
        prevalence = (abundance > 0).mean().sort_values(ascending=False)
        st.bar_chart(prevalence.head(15))
        
        # Sample richness distribution
        st.subheader("Sample Viral Richness Distribution")
        richness = (abundance > 0).sum(axis=1)
        richness_counts = richness.value_counts().sort_index()
        st.bar_chart(richness_counts)
        
        # Correlation heatmap (simplified)
        st.subheader("Viral Co-occurrence Heatmap")
        presence = (abundance > 0).astype(int)
        top_viruses = presence.sum().nlargest(10).index
        corr_matrix = presence[top_viruses].corr()
        
        # Display as a simple table since we don't have seaborn
        st.write("**Correlation Matrix (Top 10 Viruses):**")
        st.dataframe(corr_matrix.round(3))
        
        # Geographic distribution
        if 'country' in metadata.columns:
            st.subheader("Geographic Distribution")
            country_counts = metadata['country'].value_counts()
            st.bar_chart(country_counts)
        
        # Species distribution
        if 'mosquito_species' in metadata.columns:
            st.subheader("Mosquito Species Distribution")
            species_counts = metadata['mosquito_species'].value_counts()
            st.bar_chart(species_counts)
    
    else:
        st.warning("No data available for visualization. Please load data first.")

def main():
    """Main Streamlit application."""
    st.set_page_config(
        page_title="Viral Interference Research",
        page_icon="🦟",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Sidebar navigation
    st.sidebar.title("Navigation")
    page = st.sidebar.selectbox(
        "Choose a page:",
        ["Home", "Data Upload", "Statistical Analysis", "Visualization"]
    )
    
    # Display selected page
    if page == "Home":
        show_home_page()
    elif page == "Data Upload":
        show_data_page()
    elif page == "Statistical Analysis":
        show_analysis_page()
    elif page == "Visualization":
        show_visualization_page()
    
    # Footer
    st.sidebar.markdown("---")
    st.sidebar.markdown(
        "**Viral Interference Research Platform**\n\n"
        "A tool for analyzing mosquito virome data and viral interactions."
    )

if __name__ == "__main__":
    main()