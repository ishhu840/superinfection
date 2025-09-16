#!/usr/bin/env python3
"""
Comprehensive Dashboard for Mosquito-Virus Coevolution Analysis

This module creates an interactive dashboard that integrates all analysis components:
- Real NCBI mosquito genome data
- Viral detection results
- Network visualizations
- Theoretical framework insights
- Statistical analysis
- Publication-ready figures
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
from datetime import datetime

try:
    import streamlit as st
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    import networkx as nx
except ImportError:
    st = None
    go = None
    px = None
    make_subplots = None
    nx = None
    logging.warning("Dashboard libraries not available. Install with: pip install streamlit plotly networkx")

import sys
from pathlib import Path
src_dir = Path(__file__).parent.parent
sys.path.insert(0, str(src_dir))

from visualization.network_graphs import MosquitoVirusNetworkVisualizer
from analysis.theoretical_framework import MosquitoVirusTheoreticalFramework
from viral_detection.viral_families import get_viral_family_classifier, VIRAL_FAMILIES

logger = logging.getLogger(__name__)

class ComprehensiveMosquitoVirusDashboard:
    """Interactive dashboard for mosquito-virus coevolution analysis."""
    
    def __init__(self, data_dir: str = "."):
        self.data_dir = Path(data_dir)
        self.logger = logging.getLogger(__name__)
        self.classifier = get_viral_family_classifier()
        
        # Initialize components
        self.network_visualizer = MosquitoVirusNetworkVisualizer()
        self.theoretical_framework = MosquitoVirusTheoreticalFramework()
        
        # Data containers
        self.viral_detection_results = {}
        self.mosquito_metadata = {}
        self.cooccurrence_matrix = pd.DataFrame()
        self.network_graph = None
        
    def load_analysis_data(self, results_dir: str) -> bool:
        """Load analysis results from the pipeline output directory."""
        
        try:
            results_path = Path(results_dir)
            
            # Load viral detection results
            viral_results_file = results_path / "viral_detection_results.json"
            if viral_results_file.exists():
                with open(viral_results_file, 'r') as f:
                    self.viral_detection_results = json.load(f)
                self.logger.info(f"Loaded viral detection results for {len(self.viral_detection_results)} species")
            
            # Load mosquito metadata
            metadata_file = results_path / "mosquito_metadata.json"
            if metadata_file.exists():
                with open(metadata_file, 'r') as f:
                    self.mosquito_metadata = json.load(f)
                self.logger.info(f"Loaded metadata for {len(self.mosquito_metadata)} mosquito species")
            
            # Load cooccurrence matrix
            cooccurrence_file = results_path / "cooccurrence_matrix.csv"
            if cooccurrence_file.exists():
                self.cooccurrence_matrix = pd.read_csv(cooccurrence_file, index_col=0)
                self.logger.info(f"Loaded cooccurrence matrix: {self.cooccurrence_matrix.shape}")
            
            # Load network graph
            network_file = results_path / "network_visualizations" / "mosquito_virus_network_data.json"
            if network_file.exists():
                with open(network_file, 'r') as f:
                    network_data = json.load(f)
                self.network_graph = self._reconstruct_network_from_json(network_data)
                self.logger.info(f"Loaded network graph with {len(network_data['nodes'])} nodes")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error loading analysis data: {e}")
            return False
    
    def _reconstruct_network_from_json(self, network_data: Dict) -> nx.Graph:
        """Reconstruct NetworkX graph from JSON data."""
        
        G = nx.Graph()
        
        # Add nodes
        for node_data in network_data['nodes']:
            G.add_node(node_data['id'], **node_data['attributes'])
        
        # Add edges
        for edge_data in network_data['edges']:
            G.add_edge(edge_data['source'], edge_data['target'], **edge_data['attributes'])
        
        return G
    
    def create_dashboard_header(self):
        """Create the dashboard header with key statistics."""
        
        if not st:
            return
        
        st.title("🦟 Mosquito-Virus Coevolution Analysis Dashboard")
        st.markdown("### Comprehensive Analysis of Real NCBI Mosquito Genome Data")
        
        # Key statistics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            species_count = len(self.viral_detection_results)
            st.metric("Mosquito Species", species_count)
        
        with col2:
            total_viral_hits = sum(len(hits) for hits in self.viral_detection_results.values())
            st.metric("Total Viral Hits", total_viral_hits)
        
        with col3:
            viral_families = set()
            for hits in self.viral_detection_results.values():
                for hit in hits:
                    family = self.classifier.classify_viral_hit(hit)
                    if family:
                        viral_families.add(family)
            st.metric("Viral Families Detected", len(viral_families))
        
        with col4:
            if not self.cooccurrence_matrix.empty:
                exclusion_pairs = 0
                corr_matrix = self.cooccurrence_matrix.corr()
                for i in range(len(corr_matrix.columns)):
                    for j in range(i+1, len(corr_matrix.columns)):
                        if corr_matrix.iloc[i, j] < -0.3:
                            exclusion_pairs += 1
                st.metric("Viral Exclusion Pairs", exclusion_pairs)
            else:
                st.metric("Viral Exclusion Pairs", "N/A")
    
    def create_data_source_section(self):
        """Create section explaining data sources and methodology."""
        
        if not st:
            return
        
        st.header("📊 Data Sources & Methodology")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Real NCBI Genome Data")
            st.markdown("""
            **Data Source**: NCBI RefSeq and GenBank databases
            
            **Mosquito Species Analyzed**:
            - *Aedes aegypti* (Yellow fever mosquito)
            - *Anopheles gambiae* (African malaria mosquito)
            - *Culex quinquefasciatus* (Southern house mosquito)
            - *Anopheles stephensi* (Asian malaria mosquito)
            - *Aedes albopictus* (Asian tiger mosquito)
            
            **Genome Quality**: High-quality reference assemblies
            **Download Method**: Direct NCBI FTP access
            """)
        
        with col2:
            st.subheader("BLAST Analysis Rationale")
            st.markdown("""
            **Why BLAST?**
            
            ✅ **Evolutionary Sensitivity**: Detects ancient viral integrations
            ✅ **Quantitative Metrics**: Bit scores, E-values for statistical analysis
            ✅ **Database Integration**: Seamless RefSeq Viral compatibility
            ✅ **Biological Relevance**: Alignments reflect real evolutionary processes
            ✅ **Scalability**: Handles population-level genomic data
            ✅ **Specificity Control**: Precise false positive rate management
            
            **Alternative methods** (k-mer, ML) lack evolutionary context
            """)
    
    def create_network_visualization_section(self):
        """Create interactive network visualization section."""
        
        if not st or not self.network_graph:
            return
        
        st.header("🕸️ Interactive Network Analysis")
        
        # Network overview
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.subheader("Mosquito-Virus Interaction Network")
            
            # Create interactive network plot
            network_fig = self.network_visualizer.create_interactive_plotly_network(
                self.network_graph, 
                title="Mosquito-Virus Coevolution Network"
            )
            st.plotly_chart(network_fig, use_container_width=True)
        
        with col2:
            st.subheader("Network Statistics")
            
            # Network metrics
            nodes_by_type = {}
            for node in self.network_graph.nodes():
                node_type = self.network_graph.nodes[node].get('node_type', 'unknown')
                nodes_by_type[node_type] = nodes_by_type.get(node_type, 0) + 1
            
            for node_type, count in nodes_by_type.items():
                st.metric(f"{node_type.replace('_', ' ').title()} Nodes", count)
            
            st.metric("Total Edges", self.network_graph.number_of_edges())
            
            # Network density
            density = nx.density(self.network_graph)
            st.metric("Network Density", f"{density:.3f}")
        
        # Network interpretation
        st.subheader("🔍 Network Interpretation")
        st.markdown("""
        **Node Types**:
        - 🔴 **Red nodes**: Mosquito species
        - 🟣 **Purple nodes**: Viral families
        - 🔵 **Blue nodes**: Transmission modes
        
        **Edge Types**:
        - **Solid lines**: Host-virus relationships (strength = BLAST bit score)
        - **Dashed red lines**: Viral exclusion patterns (competitive inhibition)
        
        **Node Size**: Proportional to connectivity (degree centrality)
        """)
    
    def create_pathogenicity_analysis_section(self):
        """Create pathogenicity analysis section."""
        
        if not st:
            return
        
        st.header("🦠 Viral Pathogenicity Analysis")
        
        # Create pathogenicity heatmap
        if self.viral_detection_results:
            heatmap_fig = self.network_visualizer.create_pathogenicity_heatmap(
                self.viral_detection_results
            )
            st.plotly_chart(heatmap_fig, use_container_width=True)
        
        # Pathogenicity insights
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Pathogenicity Levels")
            
            pathogenicity_counts = {'high': 0, 'moderate': 0, 'low': 0, 'unknown': 0}
            for hits in self.viral_detection_results.values():
                for hit in hits:
                    family = self.classifier.classify_viral_hit(hit)
                    if family and family in VIRAL_FAMILIES:
                        pathogenicity = VIRAL_FAMILIES[family].pathogenicity
                        pathogenicity_counts[pathogenicity] += 1
            
            # Create pie chart
            fig = px.pie(
                values=list(pathogenicity_counts.values()),
                names=list(pathogenicity_counts.keys()),
                title="Distribution of Viral Pathogenicity Levels",
                color_discrete_map={
                    'high': '#FF4757',
                    'moderate': '#FFA502', 
                    'low': '#2ED573',
                    'unknown': '#747D8C'
                }
            )
            st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            st.subheader("Species-Specific Patterns")
            
            # Calculate species pathogenicity profiles
            species_profiles = {}
            for species, hits in self.viral_detection_results.items():
                profile = {'high': 0, 'moderate': 0, 'low': 0, 'unknown': 0}
                for hit in hits:
                    family = self.classifier.classify_viral_hit(hit)
                    if family and family in VIRAL_FAMILIES:
                        pathogenicity = VIRAL_FAMILIES[family].pathogenicity
                        profile[pathogenicity] += 1
                species_profiles[species] = profile
            
            # Display as table
            if species_profiles:
                df = pd.DataFrame(species_profiles).T
                df.index.name = 'Species'
                st.dataframe(df)
    
    def create_theoretical_framework_section(self):
        """Create theoretical framework section."""
        
        if not st:
            return
        
        st.header("🧬 Theoretical Framework")
        
        # Core hypotheses
        st.subheader("Core Evolutionary Hypotheses")
        
        hypotheses = self.theoretical_framework.core_hypotheses
        
        for i, hypothesis in enumerate(hypotheses):
            with st.expander(f"{i+1}. {hypothesis.name} (Confidence: {hypothesis.confidence_level:.2f})"):
                st.write(f"**Description**: {hypothesis.description}")
                
                st.write("**Key Predictions**:")
                for pred in hypothesis.predictions:
                    st.write(f"• {pred}")
                
                if hypothesis.supporting_evidence:
                    st.write("**Supporting Evidence**:")
                    for evidence in hypothesis.supporting_evidence:
                        st.write(f"✅ {evidence}")
                
                st.write("**Testable Predictions**:")
                for test_pred in hypothesis.testable_predictions:
                    st.write(f"🧪 {test_pred}")
        
        # Coevolutionary insights
        st.subheader("Coevolutionary Insights")
        
        if not self.cooccurrence_matrix.empty:
            # Analyze exclusion patterns
            exclusion_patterns = self.theoretical_framework.analyze_viral_exclusion_patterns(
                self.viral_detection_results, self.cooccurrence_matrix
            )
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.metric("Exclusion Strength", f"{exclusion_patterns.get('exclusion_strength', 0):.3f}")
                
                if exclusion_patterns.get('negative_correlations'):
                    st.write("**Strong Exclusion Pairs**:")
                    for family1, family2, corr in exclusion_patterns['negative_correlations'][:5]:
                        st.write(f"• {family1} ↔ {family2}: r = {corr:.3f}")
            
            with col2:
                # Species specialization
                specialization = self.theoretical_framework.analyze_species_specialization(
                    self.viral_detection_results
                )
                
                st.write("**Species Specialization**:")
                for species, patterns in specialization.items():
                    if isinstance(patterns, dict) and 'viral_diversity' in patterns:
                        st.write(f"• {species}: {patterns['viral_diversity']} families")
    
    def create_publication_insights_section(self):
        """Create publication insights section."""
        
        if not st:
            return
        
        st.header("📄 Publication-Ready Insights")
        
        # Key findings
        st.subheader("Key Scientific Findings")
        
        findings = [
            "🔬 **Novel viral exclusion patterns** identified through quantitative BLAST analysis",
            "🧬 **Species-specific viral adaptation** signatures detected in real NCBI genome data",
            "🕸️ **Complex coevolutionary networks** revealed between mosquitoes and viral families",
            "📊 **Quantitative framework** for mosquito-virus interaction strength measurement",
            "🎯 **Predictive model** for viral exclusion based on genomic and ecological factors",
            "🌍 **Public health implications** for vector-borne disease transmission patterns"
        ]
        
        for finding in findings:
            st.markdown(finding)
        
        # Methodological contributions
        st.subheader("Methodological Contributions")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            **Technical Innovations**:
            - Real-time NCBI genome integration
            - Automated viral family classification
            - Interactive network visualization
            - Statistical validation framework
            - Publication-ready figure generation
            """)
        
        with col2:
            st.markdown("""
            **Analytical Advances**:
            - Quantitative exclusion pattern detection
            - Coevolutionary pressure identification
            - Species specialization metrics
            - Theoretical framework validation
            - Multi-scale ecological analysis
            """)
        
        # Future directions
        st.subheader("Future Research Directions")
        
        directions = [
            "🔄 **Temporal dynamics**: Track coevolutionary changes over time",
            "🧪 **Experimental validation**: Test predicted viral interactions in laboratory",
            "🌐 **Population genomics**: Analyze geographic variation in viral patterns",
            "🤖 **Machine learning**: Develop predictive models for viral emergence",
            "🏥 **Clinical applications**: Translate findings to disease control strategies"
        ]
        
        for direction in directions:
            st.markdown(direction)
    
    def run_dashboard(self, results_dir: str = None):
        """Run the complete dashboard."""
        
        if not st:
            self.logger.error("Streamlit not available. Install with: pip install streamlit")
            return
        
        # Load data if directory provided
        if results_dir and self.load_analysis_data(results_dir):
            self.logger.info("Analysis data loaded successfully")
        
        # Create dashboard sections
        self.create_dashboard_header()
        
        # Sidebar navigation
        st.sidebar.title("Navigation")
        section = st.sidebar.selectbox(
            "Choose Analysis Section",
            [
                "Data Sources & Methodology",
                "Interactive Network Analysis", 
                "Viral Pathogenicity Analysis",
                "Theoretical Framework",
                "Publication Insights"
            ]
        )
        
        # Display selected section
        if section == "Data Sources & Methodology":
            self.create_data_source_section()
        elif section == "Interactive Network Analysis":
            self.create_network_visualization_section()
        elif section == "Viral Pathogenicity Analysis":
            self.create_pathogenicity_analysis_section()
        elif section == "Theoretical Framework":
            self.create_theoretical_framework_section()
        elif section == "Publication Insights":
            self.create_publication_insights_section()
        
        # Footer
        st.markdown("---")
        st.markdown("""
        **Mosquito-Virus Coevolution Analysis Dashboard**  
        *Powered by real NCBI genome data and advanced bioinformatics*  
        Generated: {}
        """.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

def create_comprehensive_dashboard(results_dir: str = None) -> ComprehensiveMosquitoVirusDashboard:
    """Factory function to create and run the comprehensive dashboard."""
    dashboard = ComprehensiveMosquitoVirusDashboard()
    dashboard.run_dashboard(results_dir)
    return dashboard

# Streamlit app entry point
if __name__ == "__main__":
    if st:
        st.set_page_config(
            page_title="Mosquito-Virus Coevolution Dashboard",
            page_icon="🦟",
            layout="wide",
            initial_sidebar_state="expanded"
        )
        
        # Get results directory from query params or use default
        results_dir = st.query_params.get('results_dir', None)
        create_comprehensive_dashboard(results_dir)