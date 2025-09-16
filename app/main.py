"""Main Streamlit application for viral interference research."""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

# Add src to path for imports
sys.path.append(str(Path(__file__).parent.parent / 'src'))

from data_processing import ViralAnalysisPipeline
from statistical_analysis import ViralInteractionAnalyzer
from utils import ConfigLoader, get_logger

# Page configuration
st.set_page_config(
    page_title="Viral Interference Research Platform",
    page_icon="🦟",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
.main-header {
    font-size: 2.5rem;
    color: #1f77b4;
    text-align: center;
    margin-bottom: 2rem;
}
.sub-header {
    font-size: 1.5rem;
    color: #ff7f0e;
    margin-top: 1rem;
    margin-bottom: 1rem;
}
.metric-card {
    background-color: #f0f2f6;
    padding: 1rem;
    border-radius: 0.5rem;
    border-left: 4px solid #1f77b4;
}
.warning-box {
    background-color: #fff3cd;
    border: 1px solid #ffeaa7;
    border-radius: 0.5rem;
    padding: 1rem;
    margin: 1rem 0;
}
</style>
""", unsafe_allow_html=True)

def initialize_session_state():
    """Initialize session state variables."""
    if 'data_loaded' not in st.session_state:
        st.session_state.data_loaded = False
    if 'abundance_matrix' not in st.session_state:
        st.session_state.abundance_matrix = None
    if 'metadata' not in st.session_state:
        st.session_state.metadata = None
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = None
    if 'pipeline' not in st.session_state:
        st.session_state.pipeline = None
    if 'analyzer' not in st.session_state:
        st.session_state.analyzer = None

def load_sample_data():
    """Load sample mosquito virome data."""
    np.random.seed(42)
    
    # Sample virus names
    viruses = [
        'Aedes_aegypti_densovirus', 'Culex_flavivirus', 'Anopheles_gambiae_densovirus',
        'Mosquito_nodavirus', 'Aedes_anphevirus', 'Culex_thogotovirus',
        'Anopheles_cypovirus', 'Mosquito_iflavirus', 'Aedes_totivirus',
        'Culex_rhabdovirus', 'Anopheles_bunyavirus', 'Mosquito_picornavirus'
    ]
    
    # Sample metadata
    n_samples = 150
    species = np.random.choice(['Aedes_aegypti', 'Culex_pipiens', 'Anopheles_gambiae'], n_samples)
    locations = np.random.choice(['Site_A', 'Site_B', 'Site_C', 'Site_D'], n_samples)
    dates = pd.date_range('2023-01-01', periods=n_samples, freq='D')
    
    metadata = pd.DataFrame({
        'sample_id': [f'Sample_{i:03d}' for i in range(n_samples)],
        'mosquito_species': species,
        'location': locations,
        'collection_date': dates,
        'latitude': np.random.uniform(-10, 10, n_samples),
        'longitude': np.random.uniform(-10, 10, n_samples)
    })
    
    # Generate abundance matrix with some realistic patterns
    abundance_matrix = pd.DataFrame(
        index=metadata['sample_id'],
        columns=viruses
    )
    
    for virus in viruses:
        # Base prevalence varies by virus
        base_prevalence = np.random.uniform(0.1, 0.4)
        
        # Species-specific effects
        species_effect = {
            'Aedes_aegypti': 1.0 if 'Aedes' in virus else 0.3,
            'Culex_pipiens': 1.0 if 'Culex' in virus else 0.3,
            'Anopheles_gambiae': 1.0 if 'Anopheles' in virus else 0.3
        }
        
        # Generate abundances
        abundances = []
        for i, row in metadata.iterrows():
            prob = base_prevalence * species_effect[row['mosquito_species']]
            if np.random.random() < prob:
                # Log-normal distribution for positive detections
                abundance = np.random.lognormal(mean=2, sigma=1)
            else:
                abundance = 0
            abundances.append(abundance)
        
        abundance_matrix[virus] = abundances
    
    # Introduce some exclusion patterns
    # Aedes densovirus and Aedes anphevirus rarely co-occur
    cooccur_mask = (abundance_matrix['Aedes_aegypti_densovirus'] > 0) & (abundance_matrix['Aedes_anphevirus'] > 0)
    exclude_indices = np.random.choice(abundance_matrix[cooccur_mask].index, 
                                     size=int(0.8 * cooccur_mask.sum()), 
                                     replace=False)
    abundance_matrix.loc[exclude_indices, 'Aedes_anphevirus'] = 0
    
    # Culex viruses tend to co-occur
    culex_virus_mask = abundance_matrix['Culex_flavivirus'] > 0
    cooccur_indices = abundance_matrix[culex_virus_mask].index
    boost_indices = np.random.choice(cooccur_indices, 
                                   size=int(0.6 * len(cooccur_indices)), 
                                   replace=False)
    abundance_matrix.loc[boost_indices, 'Culex_thogotovirus'] = np.random.lognormal(2, 1, len(boost_indices))
    
    return abundance_matrix, metadata

def main_page():
    """Main dashboard page."""
    st.markdown('<h1 class="main-header">🦟 Viral Interference Research Platform</h1>', unsafe_allow_html=True)
    
    # Ethical disclaimer
    st.markdown("""
    <div class="warning-box">
        <h4>⚠️ Ethical Research Notice</h4>
        <p>This platform is designed for <strong>research and educational purposes only</strong>. 
        All analyses should be conducted in accordance with institutional ethics guidelines and 
        biosafety protocols. Results should be interpreted by qualified researchers and validated 
        through appropriate experimental methods.</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    ## About This Platform
    
    This application provides computational tools for studying **viral interference** and 
    **superinfection exclusion** in mosquito populations. The platform integrates:
    
    - 📊 **Data Processing**: Quality control and viral classification pipelines
    - 📈 **Statistical Analysis**: Co-occurrence patterns and exclusion detection
    - 🤖 **Machine Learning**: Predictive models for viral interactions
    - 🌐 **Interactive Visualization**: Network analysis and data exploration
    
    ### Key Features
    
    1. **Viral Co-occurrence Analysis**: Statistical tests for virus pair associations
    2. **Exclusion Pattern Detection**: Identification of potential superinfection exclusion
    3. **Confounding Factor Control**: Analysis accounting for mosquito species, location, etc.
    4. **Interactive Networks**: Visualization of viral interaction networks
    5. **Machine Learning Predictions**: Models for viral interaction prediction
    
    ### Getting Started
    
    1. **Load Data**: Upload your viral abundance data or use sample data
    2. **Explore**: Use the sidebar to navigate between analysis modules
    3. **Analyze**: Run statistical tests and machine learning models
    4. **Visualize**: Create interactive plots and network diagrams
    5. **Export**: Download results for further analysis
    
    ### Scientific Background
    
    **Superinfection exclusion** is a phenomenon where infection with one virus prevents 
    or reduces infection with a related virus. This mechanism has important implications for:
    
    - Vector control strategies
    - Viral evolution and diversity
    - Disease transmission dynamics
    - Biocontrol applications
    
    Use the navigation menu to explore different analysis modules.
    """, unsafe_allow_html=True)
    
    # Quick stats if data is loaded
    if st.session_state.data_loaded:
        st.markdown('<h3 class="sub-header">📊 Data Overview</h3>', unsafe_allow_html=True)
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric(
                "Samples", 
                len(st.session_state.abundance_matrix),
                help="Total number of mosquito samples"
            )
        
        with col2:
            st.metric(
                "Viruses", 
                len(st.session_state.abundance_matrix.columns),
                help="Number of viral taxa detected"
            )
        
        with col3:
            total_detections = (st.session_state.abundance_matrix > 0).sum().sum()
            st.metric(
                "Detections", 
                total_detections,
                help="Total viral detections across all samples"
            )
        
        with col4:
            if st.session_state.metadata is not None:
                n_species = st.session_state.metadata['mosquito_species'].nunique()
                st.metric(
                    "Species", 
                    n_species,
                    help="Number of mosquito species"
                )

def data_upload_page():
    """Data upload and management page."""
    st.markdown('<h2 class="section-header">📁 Data Upload & Management</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["Upload Data", "Sample Data", "Data Preview"])
    
    with tab1:
        st.markdown("""
        ### Upload Your Data
        
        Upload viral abundance data and sample metadata for analysis.
        
        **Required Format:**
        - **Abundance Matrix**: CSV with samples as rows, viruses as columns
        - **Metadata**: CSV with sample information (species, location, date, etc.)
        """, unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Viral Abundance Data")
            abundance_file = st.file_uploader(
                "Choose abundance matrix CSV",
                type=['csv'],
                help="CSV file with samples as rows and viruses as columns"
            )
            
            if abundance_file is not None:
                try:
                    abundance_df = pd.read_csv(abundance_file, index_col=0)
                    st.success(f"Loaded {abundance_df.shape[0]} samples × {abundance_df.shape[1]} viruses")
                    st.session_state.abundance_matrix = abundance_df
                except Exception as e:
                    st.error(f"Error loading abundance data: {e}")
        
        with col2:
            st.subheader("Sample Metadata")
            metadata_file = st.file_uploader(
                "Choose metadata CSV",
                type=['csv'],
                help="CSV file with sample information"
            )
            
            if metadata_file is not None:
                try:
                    metadata_df = pd.read_csv(metadata_file)
                    st.success(f"Loaded metadata for {len(metadata_df)} samples")
                    st.session_state.metadata = metadata_df
                except Exception as e:
                    st.error(f"Error loading metadata: {e}")
        
        if st.session_state.abundance_matrix is not None and st.session_state.metadata is not None:
            if st.button("Confirm Data Loading"):
                st.session_state.data_loaded = True
                st.success("✅ Data successfully loaded!")
                st.rerun()
    
    with tab2:
        st.markdown("""
        ### Use Sample Data
        
        Load synthetic mosquito virome data for testing and demonstration.
        
        **Sample Dataset Includes:**
        - 150 mosquito samples
        - 12 viral taxa
        - 3 mosquito species
        - Realistic abundance patterns
        - Some exclusion and co-occurrence patterns
        """, unsafe_allow_html=True)
        
        if st.button("Load Sample Data", type="primary"):
            with st.spinner("Generating sample data..."):
                abundance_matrix, metadata = load_sample_data()
                st.session_state.abundance_matrix = abundance_matrix
                st.session_state.metadata = metadata
                st.session_state.data_loaded = True
                st.success("✅ Sample data loaded successfully!")
                st.rerun()
    
    with tab3:
        if st.session_state.data_loaded:
            st.subheader("Abundance Matrix Preview")
            st.dataframe(
                st.session_state.abundance_matrix.head(10),
                use_container_width=True
            )
            
            st.subheader("Metadata Preview")
            st.dataframe(
                st.session_state.metadata.head(10),
                use_container_width=True
            )
            
            # Basic statistics
            st.subheader("Data Summary")
            col1, col2 = st.columns(2)
            
            with col1:
                st.write("**Abundance Matrix Statistics**")
                stats_df = st.session_state.abundance_matrix.describe()
                st.dataframe(stats_df)
            
            with col2:
                st.write("**Viral Prevalence**")
                prevalence = (st.session_state.abundance_matrix > 0).mean().sort_values(ascending=False)
                st.bar_chart(prevalence)
        else:
            st.info("Please load data first to see preview.")

def statistical_analysis_page():
    """Statistical analysis page."""
    st.markdown('<h2 class="section-header">📊 Statistical Analysis</h2>', unsafe_allow_html=True)
    
    if not st.session_state.data_loaded:
        st.warning("Please load data first from the Data Upload page.")
        return
    
    # Initialize analyzer
    if st.session_state.analyzer is None:
        st.session_state.analyzer = ViralInteractionAnalyzer()
    
    tab1, tab2, tab3 = st.tabs(["Co-occurrence Analysis", "Exclusion Detection", "Confounder Analysis"])
    
    with tab1:
        st.subheader("Viral Co-occurrence Analysis")
        
        st.markdown("""
        Analyze patterns of viral co-occurrence using Fisher's exact tests.
        This analysis identifies virus pairs that occur together more or less 
        frequently than expected by chance.
        """)
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.write("**Analysis Parameters**")
            
            presence_threshold = st.number_input(
                "Presence Threshold (RPM)",
                min_value=0.0,
                max_value=10.0,
                value=1.0,
                step=0.1,
                help="Minimum abundance to consider virus present"
            )
            
            alpha_level = st.selectbox(
                "Significance Level",
                [0.05, 0.01, 0.001],
                index=0,
                help="Alpha level for statistical tests"
            )
            
            correction_method = st.selectbox(
                "Multiple Testing Correction",
                ['fdr_bh', 'bonferroni', 'holm'],
                index=0,
                help="Method for correcting multiple comparisons"
            )
            
            if st.button("Run Co-occurrence Analysis", type="primary"):
                with st.spinner("Running statistical analysis..."):
                    # Update analyzer parameters
                    st.session_state.analyzer.presence_threshold = presence_threshold
                    st.session_state.analyzer.fisher_alpha = alpha_level
                    st.session_state.analyzer.correction_method = correction_method
                    
                    # Run analysis
                    results = st.session_state.analyzer.analyze_cooccurrence(
                        st.session_state.abundance_matrix,
                        st.session_state.metadata
                    )
                    st.session_state.analysis_results = results
                    st.success("Analysis completed!")
        
        with col2:
            if st.session_state.analysis_results is not None:
                results = st.session_state.analysis_results
                
                # Summary metrics
                st.write("**Analysis Summary**")
                summary = results['summary_stats']
                
                metric_col1, metric_col2, metric_col3 = st.columns(3)
                with metric_col1:
                    st.metric("Virus Pairs Tested", summary['num_virus_pairs_tested'])
                with metric_col2:
                    st.metric("Significant Associations", summary['num_significant_associations'])
                with metric_col3:
                    st.metric("Exclusion Candidates", summary['num_exclusion_candidates'])
                
                # Results table
                if results['pairwise_tests'] is not None:
                    st.write("**Pairwise Test Results**")
                    
                    # Filter options
                    show_significant_only = st.checkbox("Show significant results only", value=True)
                    
                    display_df = results['pairwise_tests'].copy()
                    if show_significant_only:
                        display_df = display_df[display_df['significant']]
                    
                    # Sort by adjusted p-value
                    display_df = display_df.sort_values('p_adj')
                    
                    st.dataframe(
                        display_df[['virus_a', 'virus_b', 'odds_ratio', 'p_value', 'p_adj', 
                                  'association_type', 'both_present']],
                        use_container_width=True
                    )
    
    with tab2:
        st.subheader("Exclusion Pattern Detection")
        
        if st.session_state.analysis_results is not None:
            results = st.session_state.analysis_results
            
            if results['exclusion_candidates']:
                st.write("**Potential Exclusion Pairs**")
                
                exclusion_df = pd.DataFrame(results['exclusion_candidates'])
                
                # Display exclusion candidates
                for i, candidate in enumerate(results['exclusion_candidates'][:10]):
                    with st.expander(f"Exclusion Pair {i+1}: {candidate['virus_pair'][0]} ↔ {candidate['virus_pair'][1]}"):
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.metric("Odds Ratio", f"{candidate['odds_ratio']:.3f}")
                        with col2:
                            st.metric("Adjusted P-value", f"{candidate['p_adj']:.2e}")
                        with col3:
                            st.metric("Evidence Level", candidate['evidence_level'])
                        
                        st.write(f"**Exclusion Strength:** {candidate['exclusion_strength']:.3f}")
                        st.write(f"**Observed Co-occurrence:** {candidate['both_present']} samples")
                        st.write(f"**Expected Co-occurrence:** {candidate['expected_cooccurrence']:.1f} samples")
            else:
                st.info("No significant exclusion patterns detected.")
        else:
            st.info("Run co-occurrence analysis first to detect exclusion patterns.")
    
    with tab3:
        st.subheader("Confounder Analysis")
        
        st.markdown("""
        Control for confounding factors like mosquito species, location, and collection date
        when analyzing viral associations.
        """)
        
        if st.session_state.metadata is not None:
            # Select target virus
            viruses = st.session_state.abundance_matrix.columns.tolist()
            target_virus = st.selectbox(
                "Select Target Virus",
                viruses,
                help="Virus to analyze for associations"
            )
            
            # Select confounders
            potential_confounders = [col for col in st.session_state.metadata.columns 
                                   if col not in ['sample_id']]
            
            confounders = st.multiselect(
                "Select Confounding Variables",
                potential_confounders,
                default=['mosquito_species', 'location'] if 'mosquito_species' in potential_confounders else [],
                help="Variables to control for in the analysis"
            )
            
            if st.button("Run Confounder Analysis"):
                with st.spinner("Running confounder analysis..."):
                    confounder_results = st.session_state.analyzer.analyze_confounders(
                        st.session_state.abundance_matrix,
                        st.session_state.metadata,
                        target_virus,
                        confounders
                    )
                    
                    # Display results
                    st.write(f"**Analysis Results for {target_virus}**")
                    
                    if 'error' not in confounder_results:
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("**Viral Associations (Controlled)**")
                            if confounder_results['viral_associations']:
                                viral_assoc_df = pd.DataFrame(confounder_results['viral_associations']).T
                                st.dataframe(viral_assoc_df)
                            else:
                                st.info("No significant viral associations found.")
                        
                        with col2:
                            st.write("**Confounder Effects**")
                            if confounder_results['confounder_effects']:
                                confounder_df = pd.DataFrame(confounder_results['confounder_effects']).T
                                st.dataframe(confounder_df)
                            else:
                                st.info("No significant confounder effects found.")
                        
                        # Model performance
                        st.write("**Model Performance**")
                        perf = confounder_results['model_performance']
                        st.write(f"Accuracy: {perf['accuracy']:.3f}")
                        st.write(f"Number of features: {perf['n_features']}")
                    else:
                        st.error(f"Analysis failed: {confounder_results['error']}")
        else:
            st.warning("Metadata required for confounder analysis.")

def visualization_page():
    """Visualization page."""
    st.markdown('<h2 class="section-header">📈 Data Visualization</h2>', unsafe_allow_html=True)
    
    if not st.session_state.data_loaded:
        st.warning("Please load data first from the Data Upload page.")
        return
    
    tab1, tab2, tab3, tab4 = st.tabs(["Abundance Heatmap", "Interaction Network", "Species Comparison", "Temporal Patterns"])
    
    with tab1:
        st.subheader("Viral Abundance Heatmap")
        
        # Heatmap options
        col1, col2 = st.columns([1, 3])
        
        with col1:
            log_transform = st.checkbox("Log Transform", value=True)
            cluster_samples = st.checkbox("Cluster Samples", value=False)
            cluster_viruses = st.checkbox("Cluster Viruses", value=True)
            
            # Color scale
            color_scale = st.selectbox(
                "Color Scale",
                ['Viridis', 'Plasma', 'Blues', 'Reds', 'RdBu'],
                index=0
            )
        
        with col2:
            # Prepare data for heatmap
            heatmap_data = st.session_state.abundance_matrix.copy()
            
            if log_transform:
                heatmap_data = np.log10(heatmap_data + 1)
            
            # Create heatmap
            fig = px.imshow(
                heatmap_data.T,  # Transpose to have viruses on y-axis
                aspect='auto',
                color_continuous_scale=color_scale.lower(),
                title="Viral Abundance Across Samples",
                labels={'x': 'Samples', 'y': 'Viruses', 'color': 'Log10(Abundance + 1)' if log_transform else 'Abundance'}
            )
            
            fig.update_layout(
                height=600,
                xaxis_title="Samples",
                yaxis_title="Viruses"
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        st.subheader("Viral Interaction Network")
        
        if st.session_state.analysis_results is not None and st.session_state.analysis_results['network'] is not None:
            network = st.session_state.analysis_results['network']
            
            # Network visualization options
            col1, col2 = st.columns([1, 3])
            
            with col1:
                layout_algorithm = st.selectbox(
                    "Layout Algorithm",
                    ['spring', 'circular', 'kamada_kawai'],
                    index=0
                )
                
                show_edge_labels = st.checkbox("Show Edge Labels", value=False)
                node_size_factor = st.slider("Node Size", 10, 50, 20)
                
                # Filter by association type
                association_filter = st.multiselect(
                    "Show Associations",
                    ['positive', 'negative'],
                    default=['positive', 'negative']
                )
            
            with col2:
                # Filter network
                filtered_edges = [(u, v, d) for u, v, d in network.edges(data=True) 
                                if d['association_type'] in association_filter]
                
                filtered_network = nx.Graph()
                filtered_network.add_nodes_from(network.nodes())
                filtered_network.add_edges_from([(u, v, d) for u, v, d in filtered_edges])
                
                if len(filtered_network.edges()) > 0:
                    # Get layout positions
                    if layout_algorithm == 'spring':
                        pos = nx.spring_layout(filtered_network, k=1, iterations=50)
                    elif layout_algorithm == 'circular':
                        pos = nx.circular_layout(filtered_network)
                    else:
                        pos = nx.kamada_kawai_layout(filtered_network)
                    
                    # Create plotly figure
                    fig = go.Figure()
                    
                    # Add edges
                    for edge in filtered_network.edges(data=True):
                        x0, y0 = pos[edge[0]]
                        x1, y1 = pos[edge[1]]
                        
                        color = 'red' if edge[2]['association_type'] == 'negative' else 'blue'
                        width = min(abs(np.log(edge[2]['odds_ratio'])) * 2, 10)
                        
                        fig.add_trace(go.Scatter(
                            x=[x0, x1, None],
                            y=[y0, y1, None],
                            mode='lines',
                            line=dict(color=color, width=width),
                            showlegend=False,
                            hoverinfo='none'
                        ))
                    
                    # Add nodes
                    node_x = [pos[node][0] for node in filtered_network.nodes()]
                    node_y = [pos[node][1] for node in filtered_network.nodes()]
                    node_text = list(filtered_network.nodes())
                    
                    fig.add_trace(go.Scatter(
                        x=node_x,
                        y=node_y,
                        mode='markers+text',
                        marker=dict(
                            size=node_size_factor,
                            color='lightblue',
                            line=dict(width=2, color='black')
                        ),
                        text=node_text,
                        textposition='middle center',
                        textfont=dict(size=10),
                        showlegend=False
                    ))
                    
                    fig.update_layout(
                        title="Viral Interaction Network",
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        annotations=[
                            dict(
                                text="Blue edges: positive associations (co-occurrence)<br>Red edges: negative associations (exclusion)",
                                showarrow=False,
                                xref="paper", yref="paper",
                                x=0.005, y=-0.002,
                                xanchor='left', yanchor='bottom',
                                font=dict(size=12)
                            )
                        ],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        height=600
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("No interactions to display with current filters.")
        else:
            st.info("Run statistical analysis first to generate interaction network.")
    
    with tab3:
        st.subheader("Species Comparison")
        
        if 'mosquito_species' in st.session_state.metadata.columns:
            # Viral prevalence by species
            species_prevalence = []
            
            for species in st.session_state.metadata['mosquito_species'].unique():
                species_samples = st.session_state.metadata[st.session_state.metadata['mosquito_species'] == species]['sample_id']
                species_abundance = st.session_state.abundance_matrix.loc[species_samples]
                prevalence = (species_abundance > 0).mean()
                
                for virus in prevalence.index:
                    species_prevalence.append({
                        'Species': species,
                        'Virus': virus,
                        'Prevalence': prevalence[virus]
                    })
            
            prevalence_df = pd.DataFrame(species_prevalence)
            
            # Create grouped bar chart
            fig = px.bar(
                prevalence_df,
                x='Virus',
                y='Prevalence',
                color='Species',
                title='Viral Prevalence by Mosquito Species',
                barmode='group'
            )
            
            fig.update_layout(
                xaxis_tickangle=-45,
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Species information not available in metadata.")
    
    with tab4:
        st.subheader("Temporal Patterns")
        
        if 'collection_date' in st.session_state.metadata.columns:
            # Convert to datetime if needed
            metadata_temp = st.session_state.metadata.copy()
            metadata_temp['collection_date'] = pd.to_datetime(metadata_temp['collection_date'])
            
            # Select virus for temporal analysis
            selected_virus = st.selectbox(
                "Select Virus for Temporal Analysis",
                st.session_state.abundance_matrix.columns.tolist()
            )
            
            # Aggregate by date
            temporal_data = []
            
            for date in metadata_temp['collection_date'].dt.date.unique():
                date_samples = metadata_temp[metadata_temp['collection_date'].dt.date == date]['sample_id']
                date_abundance = st.session_state.abundance_matrix.loc[date_samples, selected_virus]
                
                temporal_data.append({
                    'Date': date,
                    'Prevalence': (date_abundance > 0).mean(),
                    'Mean_Abundance': date_abundance.mean(),
                    'Sample_Count': len(date_samples)
                })
            
            temporal_df = pd.DataFrame(temporal_data).sort_values('Date')
            
            # Create temporal plot
            fig = make_subplots(
                rows=2, cols=1,
                subplot_titles=['Prevalence Over Time', 'Mean Abundance Over Time'],
                vertical_spacing=0.1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=temporal_df['Date'],
                    y=temporal_df['Prevalence'],
                    mode='lines+markers',
                    name='Prevalence',
                    line=dict(color='blue')
                ),
                row=1, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=temporal_df['Date'],
                    y=temporal_df['Mean_Abundance'],
                    mode='lines+markers',
                    name='Mean Abundance',
                    line=dict(color='red')
                ),
                row=2, col=1
            )
            
            fig.update_layout(
                title=f'Temporal Patterns for {selected_virus}',
                height=600,
                showlegend=False
            )
            
            fig.update_xaxes(title_text="Date", row=2, col=1)
            fig.update_yaxes(title_text="Prevalence", row=1, col=1)
            fig.update_yaxes(title_text="Mean Abundance", row=2, col=1)
            
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Temporal information not available in metadata.")

def main():
    """Main application function."""
    initialize_session_state()
    
    # Sidebar navigation
    st.sidebar.title("🦟 Navigation")
    
    pages = {
        "🏠 Home": main_page,
        "📁 Data Upload": data_upload_page,
        "📊 Statistical Analysis": statistical_analysis_page,
        "📈 Visualization": visualization_page
    }
    
    selected_page = st.sidebar.selectbox(
        "Choose a page",
        list(pages.keys())
    )
    
    # Data status indicator
    if st.session_state.data_loaded:
        st.sidebar.success("✅ Data Loaded")
        
        # Quick data info
        st.sidebar.markdown("**Data Summary:**")
        st.sidebar.write(f"Samples: {len(st.session_state.abundance_matrix)}")
        st.sidebar.write(f"Viruses: {len(st.session_state.abundance_matrix.columns)}")
    else:
        st.sidebar.warning("⚠️ No Data Loaded")
    
    # Analysis status
    if st.session_state.analysis_results is not None:
        st.sidebar.success("✅ Analysis Complete")
    
    # Run selected page
    pages[selected_page]()
    
    # Footer
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    **About:** This platform is for research and educational purposes only.
    Always follow institutional ethics guidelines.
    """)

if __name__ == "__main__":
    main()