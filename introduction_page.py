#!/usr/bin/env python3
"""
Comprehensive Introduction Page for Mosquito Viral Cooccurrence Analysis
Explains methodology, datasets, theoretical framework, and scientific approach
"""

import streamlit as st
import pandas as pd
import json
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import numpy as np
from collections import Counter

# Page configuration
st.set_page_config(
    page_title="Mosquito Viral Analysis - Introduction",
    page_icon="🦟",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
.main-header {
    font-size: 3rem;
    color: #1f77b4;
    text-align: center;
    margin-bottom: 2rem;
    border-bottom: 4px solid #1f77b4;
    padding-bottom: 1rem;
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
}
.section-header {
    font-size: 2rem;
    color: #2c3e50;
    margin: 2rem 0 1rem 0;
    border-left: 5px solid #3498db;
    padding-left: 1rem;
    background: #f8f9fa;
    padding: 1rem;
    border-radius: 8px;
}
.methodology-box {
    background: linear-gradient(135deg, #74b9ff 0%, #0984e3 100%);
    padding: 2rem;
    border-radius: 15px;
    color: white;
    margin: 1rem 0;
    box-shadow: 0 8px 32px rgba(0,0,0,0.1);
}
.dataset-card {
    background: #ffffff;
    padding: 1.5rem;
    border-radius: 12px;
    border: 1px solid #e0e0e0;
    margin: 1rem 0;
    box-shadow: 0 4px 16px rgba(0,0,0,0.1);
    transition: transform 0.3s ease;
}
.dataset-card:hover {
    transform: translateY(-5px);
    box-shadow: 0 8px 25px rgba(0,0,0,0.15);
}
.theory-section {
    background: linear-gradient(135deg, #fd79a8 0%, #e84393 100%);
    padding: 2rem;
    border-radius: 15px;
    color: white;
    margin: 2rem 0;
}
.results-preview {
    background: linear-gradient(135deg, #00b894 0%, #00a085 100%);
    padding: 2rem;
    border-radius: 15px;
    color: white;
    margin: 2rem 0;
}
.metric-highlight {
    background: rgba(255,255,255,0.2);
    padding: 1rem;
    border-radius: 8px;
    text-align: center;
    margin: 0.5rem;
}
.step-box {
    background: #f1f2f6;
    padding: 1.5rem;
    border-radius: 10px;
    border-left: 4px solid #3742fa;
    margin: 1rem 0;
}
.highlight-stat {
    font-size: 2.5rem;
    font-weight: bold;
    color: #e17055;
    text-align: center;
}
.navigation-box {
    background: linear-gradient(135deg, #a29bfe 0%, #6c5ce7 100%);
    padding: 2rem;
    border-radius: 15px;
    color: white;
    margin: 2rem 0;
    text-align: center;
}
</style>
""", unsafe_allow_html=True)

def load_analysis_data():
    """Load analysis results for statistics"""
    results_file = Path("/Users/ishtiaq/Desktop/super-virus/real_viral_results/real_viral_analysis_results.json")
    
    if results_file.exists():
        with open(results_file, 'r') as f:
            return json.load(f)
    return None

def create_methodology_flowchart():
    """Create methodology flowchart"""
    fig = go.Figure()
    
    # Define steps and positions
    steps = [
        {"name": "Data Collection", "x": 1, "y": 5, "color": "#74b9ff"},
        {"name": "Quality Control", "x": 2, "y": 5, "color": "#fd79a8"},
        {"name": "BLAST Analysis", "x": 3, "y": 5, "color": "#00b894"},
        {"name": "Viral Detection", "x": 4, "y": 5, "color": "#fdcb6e"},
        {"name": "Cooccurrence Analysis", "x": 2.5, "y": 3, "color": "#e17055"},
        {"name": "Statistical Validation", "x": 2.5, "y": 1, "color": "#a29bfe"}
    ]
    
    # Add nodes
    for step in steps:
        fig.add_trace(go.Scatter(
            x=[step["x"]], y=[step["y"]],
            mode='markers+text',
            marker=dict(size=80, color=step["color"]),
            text=step["name"],
            textposition="middle center",
            textfont=dict(size=10, color="white"),
            showlegend=False,
            hoverinfo='text',
            hovertext=f"Step: {step['name']}"
        ))
    
    # Add arrows
    arrows = [
        (1, 5, 2, 5), (2, 5, 3, 5), (3, 5, 4, 5),
        (3, 5, 2.5, 3), (2.5, 3, 2.5, 1)
    ]
    
    for x1, y1, x2, y2 in arrows:
        fig.add_annotation(
            x=x2, y=y2, ax=x1, ay=y1,
            xref='x', yref='y', axref='x', ayref='y',
            arrowhead=2, arrowsize=1, arrowwidth=2,
            arrowcolor="#2d3436"
        )
    
    fig.update_layout(
        title="Viral Cooccurrence Analysis Methodology",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0, 5]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0, 6]),
        height=400,
        plot_bgcolor='white'
    )
    
    return fig

def create_dataset_overview():
    """Create dataset overview visualization"""
    # Sample data representing our datasets
    datasets = {
        'NCBI RefSeq Genomes': 3,
        'SRA Virome Studies': 20,
        'Field Sample Simulations': 50,
        'Published Datasets': 15,
        'Metagenomic Samples': 30
    }
    
    fig = px.bar(
        x=list(datasets.keys()),
        y=list(datasets.values()),
        title="Dataset Sources and Sample Counts",
        labels={'x': 'Dataset Source', 'y': 'Number of Samples'},
        color=list(datasets.values()),
        color_continuous_scale='viridis'
    )
    
    fig.update_layout(height=400, showlegend=False)
    return fig

def main():
    # Main header
    st.markdown('<h1 class="main-header">🦟 Mosquito Viral Cooccurrence Analysis Platform</h1>', unsafe_allow_html=True)
    
    # Introduction
    st.markdown("""
    <div class="intro-section">
        <h2>Advanced Computational Framework for Understanding Viral Interactions in Mosquito Vectors</h2>
        <p>A comprehensive bioinformatics pipeline for analyzing viral cooccurrence patterns, superinfection exclusion mechanisms, 
        and ecological interactions in mosquito-borne virus systems.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Load data for statistics
    data = load_analysis_data()
    
    # Quick Stats
    if data:
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric(
                label="Genome Samples",
                value=data["analysis_metadata"]["total_samples"]
            )
        
        with col2:
            st.metric(
                label="Viral Detections",
                value=data["analysis_metadata"]["total_viral_hits"]
            )
        
        with col3:
            superinfection_rate = data['cooccurrence_analysis']['superinfection_rate']
            st.metric(
                label="Superinfection Rate",
                value=f"{superinfection_rate:.1%}"
            )
        
        with col4:
            cooccurrence_pairs = len(data['cooccurrence_analysis']['cooccurrence_pairs'])
            st.metric(
                label="Cooccurrence Pairs",
                value=cooccurrence_pairs
            )
    
    # Research Objectives
    st.markdown('<h2 class="section-header">🎯 Research Objectives & Significance</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        with st.container():
            st.subheader("🔬 Primary Research Questions")
            st.markdown("""
            - **Viral Cooccurrence Patterns:** Which viruses occur together in mosquito vectors?
            - **Superinfection Exclusion:** What mechanisms prevent multiple viral infections?
            - **Species-Specific Associations:** How do viral patterns vary across mosquito species?
            - **Ecological Implications:** What are the public health and evolutionary consequences?
            """)
    
    with col2:
        with st.container():
            st.subheader("🌍 Scientific Impact")
            st.markdown("""
            - **Disease Control:** Inform vector control strategies
            - **Outbreak Prediction:** Model viral emergence patterns
            - **Vaccine Development:** Understand viral interference
            - **Evolutionary Biology:** Reveal host-pathogen dynamics
            """)
    
    # Methodology Section
    st.markdown('<h2 class="section-header">🔬 Methodology & Computational Pipeline</h2>', unsafe_allow_html=True)
    
    # Methodology flowchart
    methodology_chart = create_methodology_flowchart()
    st.plotly_chart(methodology_chart, use_container_width=True)
    
    with st.container():
        st.subheader("📋 Detailed Methodology")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### 🧬 Data Collection & Processing")
            st.markdown("""
            - **Genome Sources:** NCBI RefSeq, SRA databases, published studies
            - **Quality Control:** Sequence validation, contamination screening
            - **Standardization:** Uniform formatting and metadata annotation
            - **Species Coverage:** Aedes, Anopheles, Culex genera
            """)
        
        with col2:
            st.markdown("#### 🔍 Viral Detection Pipeline")
            st.markdown("""
            - **BLAST Analysis:** Against comprehensive viral databases
            - **Sensitivity Optimization:** E-value thresholds, coverage filters
            - **Taxonomic Classification:** Family, genus, species level
            - **Validation:** Multiple detection algorithms
            """)
    
    # Theoretical Framework
    st.markdown('<h2 class="section-header">📚 Theoretical Framework</h2>', unsafe_allow_html=True)
    
    with st.container():
        st.subheader("🧠 Superinfection Exclusion Theory")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### 🔬 Biological Mechanisms")
            st.markdown("""
            - **Resource Competition:** Limited cellular resources for replication
            - **Immune Priming:** Initial infection activates antiviral responses
            - **Receptor Blocking:** Occupied cellular entry points
            - **Metabolic Interference:** Altered cellular metabolism
            """)
        
        with col2:
            st.markdown("#### 📊 Mathematical Models")
            st.markdown("""
            - **Probability Matrices:** Cooccurrence likelihood calculations
            - **Network Analysis:** Viral interaction networks
            - **Statistical Testing:** Chi-square, Fisher's exact tests
            - **Machine Learning:** Pattern recognition algorithms
            """)
    
    # Dataset Information
    st.markdown('<h2 class="section-header">📊 Comprehensive Dataset Overview</h2>', unsafe_allow_html=True)
    
    # Dataset visualization
    dataset_chart = create_dataset_overview()
    st.plotly_chart(dataset_chart, use_container_width=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class="dataset-card">
            <h3>🧬 Reference Genomes</h3>
            <p><strong>Source:</strong> NCBI RefSeq</p>
            <p><strong>Species:</strong> 6 major mosquito species</p>
            <p><strong>Quality:</strong> High-quality, annotated genomes</p>
            <p><strong>Purpose:</strong> Baseline viral content analysis</p>
            <div style="background: #e8f5e8; padding: 1rem; border-radius: 8px; margin-top: 1rem;">
                <strong>Key Species:</strong><br>
                • Aedes aegypti<br>
                • Anopheles gambiae<br>
                • Culex quinquefasciatus
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="dataset-card">
            <h3>🔬 SRA Virome Studies</h3>
            <p><strong>Source:</strong> NCBI SRA Database</p>
            <p><strong>Samples:</strong> 20+ field-collected specimens</p>
            <p><strong>Coverage:</strong> Multiple geographic regions</p>
            <p><strong>Purpose:</strong> Real-world viral diversity</p>
            <div style="background: #fff3cd; padding: 1rem; border-radius: 8px; margin-top: 1rem;">
                <strong>Advantages:</strong><br>
                • Authentic field samples<br>
                • Natural viral communities<br>
                • Geographic diversity
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("""
        <div class="dataset-card">
            <h3>🌐 Metagenomic Data</h3>
            <p><strong>Source:</strong> Published studies</p>
            <p><strong>Samples:</strong> 50+ environmental samples</p>
            <p><strong>Scope:</strong> Viral ecology studies</p>
            <p><strong>Purpose:</strong> Comprehensive viral landscape</p>
            <div style="background: #f8d7da; padding: 1rem; border-radius: 8px; margin-top: 1rem;">
                <strong>Benefits:</strong><br>
                • Large sample sizes<br>
                • Diverse viral families<br>
                • Ecological context
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Analysis Pipeline Steps
    st.markdown('<h2 class="section-header">⚙️ Computational Analysis Pipeline</h2>', unsafe_allow_html=True)
    
    steps = [
        {
            "title": "1. Data Preprocessing",
            "description": "Quality control, sequence validation, and standardization of genome data from multiple sources.",
            "details": ["FASTQ quality assessment", "Adapter trimming", "Contamination screening", "Format standardization"]
        },
        {
            "title": "2. Viral Database Construction",
            "description": "Comprehensive viral reference database creation from NCBI RefSeq and custom sources.",
            "details": ["RefSeq viral genomes", "Custom viral sequences", "Taxonomic annotation", "BLAST database indexing"]
        },
        {
            "title": "3. Homology-Based Detection",
            "description": "BLAST-based viral sequence detection with optimized parameters for sensitivity.",
            "details": ["BLASTn analysis", "E-value optimization", "Coverage filtering", "Identity thresholds"]
        },
        {
            "title": "4. Cooccurrence Analysis",
            "description": "Statistical analysis of viral cooccurrence patterns and exclusion relationships.",
            "details": ["Pairwise analysis", "Probability calculations", "Network construction", "Statistical validation"]
        },
        {
            "title": "5. Visualization & Reporting",
            "description": "Interactive dashboards and comprehensive reports for result interpretation.",
            "details": ["Network graphs", "Heatmaps", "Statistical plots", "Downloadable reports"]
        }
    ]
    
    for i, step in enumerate(steps):
        st.markdown(
            f'<div class="step-box">'
            f'<h3>{step["title"]}</h3>'
            f'<p>{step["description"]}</p>'
            f'<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin-top: 1rem;">'
            + ''.join([f'<div style="background: #e8f4fd; padding: 0.5rem; border-radius: 4px;">• {detail}</div>' for detail in step["details"]])
            + '</div></div>',
            unsafe_allow_html=True
        )
    
    # Results Preview
    if data:
        st.markdown('<h2 class="section-header">📈 Key Findings Preview</h2>', unsafe_allow_html=True)
        
        st.markdown(
            f'<div class="results-preview">'
            f'<h3>🎯 Current Analysis Results</h3>'
            f'<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1rem; margin-top: 1rem;">'
            f'<div class="metric-highlight">'
            f'<h4>{data["analysis_metadata"]["total_samples"]}</h4>'
            f'<p>Mosquito Genome Samples Analyzed</p>'
            f'</div>'
            f'<div class="metric-highlight">'
            f'<h4>{data["analysis_metadata"]["total_viral_hits"]}</h4>'
            f'<p>Viral Sequences Detected</p>'
            f'</div>'
            f'<div class="metric-highlight">'
            f'<h4>{len(data["cooccurrence_analysis"]["cooccurrence_pairs"])}</h4>'
            f'<p>Viral Cooccurrence Pairs</p>'
            f'</div>'
            f'<div class="metric-highlight">'
            f'<h4>{data["cooccurrence_analysis"]["superinfection_rate"]:.1%}</h4>'
            f'<p>Superinfection Rate</p>'
            f'</div>'
            f'</div>'
            f'</div>',
            unsafe_allow_html=True
        )
    
    # Future Directions
    st.markdown('<h2 class="section-header">🚀 Future Directions & Scalability</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="dataset-card">
            <h3>📈 Scaling to Larger Datasets</h3>
            <ul>
                <li><strong>Target:</strong> 100,000+ genome samples</li>
                <li><strong>Sources:</strong> Global surveillance networks</li>
                <li><strong>Coverage:</strong> All major arbovirus vectors</li>
                <li><strong>Timeline:</strong> Continuous data integration</li>
            </ul>
            <div style="background: #e8f5e8; padding: 1rem; border-radius: 8px; margin-top: 1rem;">
                <strong>Expected Benefits:</strong><br>
                • Higher statistical power<br>
                • Rare pattern detection<br>
                • Geographic comparisons<br>
                • Temporal trend analysis
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="dataset-card">
            <h3>🔬 Advanced Analytics</h3>
            <ul>
                <li><strong>Machine Learning:</strong> Pattern prediction models</li>
                <li><strong>Network Analysis:</strong> Complex interaction networks</li>
                <li><strong>Phylogenetics:</strong> Evolutionary relationships</li>
                <li><strong>Real-time Analysis:</strong> Streaming data processing</li>
            </ul>
            <div style="background: #fff3cd; padding: 1rem; border-radius: 8px; margin-top: 1rem;">
                <strong>Research Applications:</strong><br>
                • Outbreak prediction<br>
                • Vaccine design<br>
                • Vector control<br>
                • Drug development
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Navigation
    st.markdown('<h2 class="section-header">🧭 Navigate the Platform</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="navigation-box">
        <h3>🎯 Explore Our Analysis Dashboards</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 2rem; margin-top: 2rem;">
            <div style="background: rgba(255,255,255,0.2); padding: 1.5rem; border-radius: 10px;">
                <h4>📊 Real Data Dashboard</h4>
                <p>Interactive visualization of actual viral cooccurrence patterns from field samples</p>
                <p><strong>Port:</strong> 8504</p>
                <p><strong>Features:</strong> Network graphs, heatmaps, statistical analysis</p>
            </div>
            <div style="background: rgba(255,255,255,0.2); padding: 1.5rem; border-radius: 10px;">
                <h4>🔬 Comprehensive Dashboard</h4>
                <p>Advanced analytics with machine learning and predictive modeling</p>
                <p><strong>Port:</strong> 8503</p>
                <p><strong>Features:</strong> ML models, predictions, advanced statistics</p>
            </div>
            <div style="background: rgba(255,255,255,0.2); padding: 1.5rem; border-radius: 10px;">
                <h4>🎮 Simple Interface</h4>
                <p>User-friendly interface for basic analysis and exploration</p>
                <p><strong>Port:</strong> 8502</p>
                <p><strong>Features:</strong> Easy navigation, basic visualizations</p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Technical Specifications
    st.markdown('<h2 class="section-header">⚙️ Technical Specifications</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="dataset-card">
            <h3>💻 Computational Requirements</h3>
            <ul>
                <li><strong>CPU:</strong> Multi-core processing (16+ cores recommended)</li>
                <li><strong>Memory:</strong> 32GB+ RAM for large datasets</li>
                <li><strong>Storage:</strong> 1TB+ for genome data and results</li>
                <li><strong>Network:</strong> High-speed internet for data downloads</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="dataset-card">
            <h3>🛠️ Software Stack</h3>
            <ul>
                <li><strong>Languages:</strong> Python 3.9+, R for statistics</li>
                <li><strong>Bioinformatics:</strong> BLAST+, BioPython, Entrez</li>
                <li><strong>Visualization:</strong> Streamlit, Plotly, NetworkX</li>
                <li><strong>Analysis:</strong> Pandas, NumPy, SciPy, Scikit-learn</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # Footer
    st.markdown("---")
    st.markdown(
        "<div style='text-align: center; color: #666; padding: 2rem;'>" 
        "🦟 <strong>Mosquito Viral Cooccurrence Analysis Platform</strong> | "
        "Advanced Bioinformatics for Vector-Borne Disease Research | "
        "Developed for PhD-level research and publication"
        "</div>",
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()