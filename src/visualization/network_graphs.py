#!/usr/bin/env python3
"""
Interactive network visualization for mosquito-virus relationships.
Creates comprehensive network graphs showing viral families, mosquito species,
and their complex ecological relationships.
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
import pandas as pd
import numpy as np

try:
    import networkx as nx
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.offline as pyo
except ImportError:
    nx = None
    go = None
    px = None
    make_subplots = None
    pyo = None
    logging.warning("Network visualization libraries not available. Install with: pip install networkx plotly")

import sys
from pathlib import Path
src_dir = Path(__file__).parent.parent
sys.path.insert(0, str(src_dir))

from viral_detection.viral_families import get_viral_family_classifier, VIRAL_FAMILIES

logger = logging.getLogger(__name__)

class MosquitoVirusNetworkVisualizer:
    """Creates interactive network visualizations of mosquito-virus relationships."""
    
    def __init__(self, output_dir: str = "."):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.classifier = get_viral_family_classifier()
        self.logger = logging.getLogger(__name__)
        
        # Color schemes for different node types
        self.colors = {
            "mosquito": "#FF6B6B",  # Red for mosquitoes
            "virus_high": "#FF4757",  # Dark red for high pathogenicity
            "virus_moderate": "#FFA502",  # Orange for moderate pathogenicity
            "virus_low": "#2ED573",  # Green for low pathogenicity
            "virus_unknown": "#747D8C",  # Gray for unknown
            "family": "#5352ED",  # Purple for viral families
            "transmission": "#70A1FF"  # Light blue for transmission modes
        }
        
    def create_comprehensive_network(self, viral_detection_results: Dict[str, List[Dict]], 
                                   mosquito_metadata: Dict[str, Dict]) -> nx.Graph:
        """Create a comprehensive network graph of mosquito-virus relationships."""
        
        if not nx:
            raise ImportError("NetworkX is required for network visualization")
            
        G = nx.Graph()
        
        # Add mosquito species nodes
        for species_key, metadata in mosquito_metadata.items():
            G.add_node(
                species_key,
                node_type="mosquito",
                species_name=metadata.get("species_name", species_key),
                common_name=metadata.get("common_name", ""),
                genome_size=metadata.get("genome_size_mb", 0),
                vector_competence=metadata.get("vector_competence", []),
                color=self.colors["mosquito"]
            )
        
        # Add viral family nodes and connections
        detected_families = set()
        for species_key, viral_hits in viral_detection_results.items():
            for hit in viral_hits:
                family = self.classifier.classify_viral_hit(hit)
                if family and family in VIRAL_FAMILIES:
                    detected_families.add(family)
                    
                    # Add viral family node if not exists
                    if not G.has_node(family):
                        family_info = VIRAL_FAMILIES[family]
                        color = self.colors.get(f"virus_{family_info.pathogenicity}", 
                                              self.colors["virus_unknown"])
                        
                        G.add_node(
                            family,
                            node_type="viral_family",
                            genome_type=family_info.genome_type.value,
                            pathogenicity=family_info.pathogenicity,
                            genome_size_range=family_info.genome_size_range,
                            representative_species=family_info.representative_species,
                            color=color
                        )
                    
                    # Add edge between mosquito and viral family
                    if G.has_node(species_key):
                        edge_weight = hit.get('bit_score', 1) / 100  # Normalize bit score
                        G.add_edge(
                            species_key, 
                            family,
                            relationship="hosts",
                            weight=edge_weight,
                            evalue=hit.get('evalue', 1.0),
                            identity=hit.get('identity', 0)
                        )
        
        # Add transmission mode nodes and connections
        transmission_modes = set()
        for family_name in detected_families:
            family_info = VIRAL_FAMILIES[family_name]
            for mode in family_info.transmission_modes:
                mode_str = mode.value
                transmission_modes.add(mode_str)
                
                # Add transmission mode node
                if not G.has_node(mode_str):
                    G.add_node(
                        mode_str,
                        node_type="transmission",
                        color=self.colors["transmission"]
                    )
                
                # Connect viral family to transmission mode
                G.add_edge(family_name, mode_str, relationship="transmitted_by")
        
        # Add exclusion relationships as negative edges
        exclusions = self.classifier.get_exclusion_patterns()
        for family1, excluded_families in exclusions.items():
            if family1 in detected_families:
                for family2 in excluded_families:
                    if family2 in detected_families:
                        G.add_edge(
                            family1, 
                            family2, 
                            relationship="excludes",
                            weight=-1.0,
                            edge_type="exclusion"
                        )
        
        self.logger.info(f"Created network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        return G
    
    def create_interactive_plotly_network(self, G: nx.Graph, title: str = "Mosquito-Virus Interaction Network") -> go.Figure:
        """Create an interactive Plotly network visualization."""
        
        if not go:
            raise ImportError("Plotly is required for interactive visualization")
        
        # Calculate layout
        pos = nx.spring_layout(G, k=3, iterations=50, seed=42)
        
        # Prepare node traces by type
        node_traces = {}
        node_types = set(G.nodes[node].get('node_type', 'unknown') for node in G.nodes())
        
        for node_type in node_types:
            nodes_of_type = [node for node in G.nodes() 
                           if G.nodes[node].get('node_type') == node_type]
            
            if not nodes_of_type:
                continue
                
            x_coords = [pos[node][0] for node in nodes_of_type]
            y_coords = [pos[node][1] for node in nodes_of_type]
            
            # Create hover text
            hover_text = []
            for node in nodes_of_type:
                node_data = G.nodes[node]
                if node_type == "mosquito":
                    text = f"<b>{node_data.get('species_name', node)}</b><br>"
                    text += f"Common name: {node_data.get('common_name', 'N/A')}<br>"
                    text += f"Genome size: {node_data.get('genome_size', 'N/A')} Mb<br>"
                    text += f"Vector competence: {', '.join(node_data.get('vector_competence', []))}"
                elif node_type == "viral_family":
                    text = f"<b>{node}</b><br>"
                    text += f"Genome type: {node_data.get('genome_type', 'N/A')}<br>"
                    text += f"Pathogenicity: {node_data.get('pathogenicity', 'N/A')}<br>"
                    text += f"Size range: {node_data.get('genome_size_range', 'N/A')} nt<br>"
                    species_list = node_data.get('representative_species', [])
                    text += f"Representative species: {', '.join(species_list[:3])}"
                    if len(species_list) > 3:
                        text += f" (+{len(species_list)-3} more)"
                else:
                    text = f"<b>{node}</b><br>Type: {node_type}"
                
                hover_text.append(text)
            
            # Get color for this node type
            if nodes_of_type:
                sample_node = nodes_of_type[0]
                color = G.nodes[sample_node].get('color', '#747D8C')
            else:
                color = '#747D8C'
            
            # Create trace
            trace = go.Scatter(
                x=x_coords,
                y=y_coords,
                mode='markers+text',
                marker=dict(
                    size=[15 + G.degree(node) * 2 for node in nodes_of_type],
                    color=color,
                    line=dict(width=2, color='white')
                ),
                text=[node.replace('_', ' ').title() for node in nodes_of_type],
                textposition="middle center",
                textfont=dict(size=10, color='white'),
                hovertext=hover_text,
                hoverinfo='text',
                name=node_type.replace('_', ' ').title()
            )
            
            node_traces[node_type] = trace
        
        # Create edge traces
        edge_traces = []
        
        # Regular edges
        regular_edges = [(u, v) for u, v, d in G.edges(data=True) 
                        if d.get('relationship') != 'excludes']
        
        if regular_edges:
            edge_x, edge_y = [], []
            for edge in regular_edges:
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
            
            edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                line=dict(width=2, color='#888'),
                hoverinfo='none',
                mode='lines',
                name='Relationships'
            )
            edge_traces.append(edge_trace)
        
        # Exclusion edges (dashed red lines)
        exclusion_edges = [(u, v) for u, v, d in G.edges(data=True) 
                          if d.get('relationship') == 'excludes']
        
        if exclusion_edges:
            excl_x, excl_y = [], []
            for edge in exclusion_edges:
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                excl_x.extend([x0, x1, None])
                excl_y.extend([y0, y1, None])
            
            exclusion_trace = go.Scatter(
                x=excl_x, y=excl_y,
                line=dict(width=3, color='red', dash='dash'),
                hoverinfo='none',
                mode='lines',
                name='Exclusion Patterns'
            )
            edge_traces.append(exclusion_trace)
        
        # Create figure
        fig = go.Figure(data=edge_traces + list(node_traces.values()))
        
        fig.update_layout(
            title=dict(
                text=title,
                x=0.5,
                font=dict(size=20)
            ),
            showlegend=True,
            hovermode='closest',
            margin=dict(b=20,l=5,r=5,t=40),
            annotations=[
                dict(
                    text="Interactive network showing mosquito species, viral families, and transmission modes.<br>" +
                         "Node size indicates connectivity. Red dashed lines show viral exclusion patterns.",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002,
                    xanchor='left', yanchor='bottom',
                    font=dict(size=12)
                )
            ],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white'
        )
        
        return fig
    
    def create_pathogenicity_heatmap(self, viral_detection_results: Dict[str, List[Dict]]) -> go.Figure:
        """Create a heatmap showing pathogenicity levels across mosquito species."""
        
        if not go:
            raise ImportError("Plotly is required for heatmap visualization")
        
        # Prepare data matrix
        species_list = list(viral_detection_results.keys())
        family_list = list(VIRAL_FAMILIES.keys())
        
        # Create pathogenicity matrix
        pathogenicity_scores = {'high': 3, 'moderate': 2, 'low': 1, 'unknown': 0}
        matrix = np.zeros((len(species_list), len(family_list)))
        
        for i, species in enumerate(species_list):
            for hit in viral_detection_results[species]:
                family = self.classifier.classify_viral_hit(hit)
                if family and family in VIRAL_FAMILIES:
                    j = family_list.index(family)
                    family_info = VIRAL_FAMILIES[family]
                    score = pathogenicity_scores.get(family_info.pathogenicity, 0)
                    matrix[i, j] = max(matrix[i, j], score)  # Take highest score if multiple hits
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=family_list,
            y=[s.replace('_', ' ').title() for s in species_list],
            colorscale=[
                [0, '#F8F9FA'],     # White for no detection
                [0.33, '#2ED573'],  # Green for low pathogenicity
                [0.66, '#FFA502'],  # Orange for moderate pathogenicity
                [1, '#FF4757']      # Red for high pathogenicity
            ],
            colorbar=dict(
                title="Pathogenicity Level",
                tickvals=[0, 1, 2, 3],
                ticktext=["None", "Low", "Moderate", "High"]
            ),
            hovertemplate='<b>%{y}</b><br>Family: %{x}<br>Pathogenicity: %{z}<extra></extra>'
        ))
        
        fig.update_layout(
            title="Viral Pathogenicity Heatmap Across Mosquito Species",
            xaxis_title="Viral Family",
            yaxis_title="Mosquito Species",
            font=dict(size=12)
        )
        
        return fig
    
    def save_visualizations(self, viral_detection_results: Dict[str, List[Dict]], 
                          mosquito_metadata: Dict[str, Dict], 
                          output_prefix: str = "mosquito_virus_network"):
        """Save all visualizations to HTML files."""
        
        try:
            # Create network graph
            G = self.create_comprehensive_network(viral_detection_results, mosquito_metadata)
            
            # Create interactive network visualization
            network_fig = self.create_interactive_plotly_network(G)
            network_file = self.output_dir / f"{output_prefix}_network.html"
            network_fig.write_html(str(network_file))
            self.logger.info(f"Saved network visualization to {network_file}")
            
            # Create pathogenicity heatmap
            heatmap_fig = self.create_pathogenicity_heatmap(viral_detection_results)
            heatmap_file = self.output_dir / f"{output_prefix}_pathogenicity.html"
            heatmap_fig.write_html(str(heatmap_file))
            self.logger.info(f"Saved pathogenicity heatmap to {heatmap_file}")
            
            # Save network data as JSON for further analysis
            network_data = {
                "nodes": [{
                    "id": node,
                    "type": G.nodes[node].get('node_type', 'unknown'),
                    "attributes": dict(G.nodes[node])
                } for node in G.nodes()],
                "edges": [{
                    "source": edge[0],
                    "target": edge[1],
                    "attributes": dict(G.edges[edge])
                } for edge in G.edges()]
            }
            
            json_file = self.output_dir / f"{output_prefix}_data.json"
            with open(json_file, 'w') as f:
                json.dump(network_data, f, indent=2, default=str)
            self.logger.info(f"Saved network data to {json_file}")
            
            return {
                "network_html": str(network_file),
                "heatmap_html": str(heatmap_file),
                "network_json": str(json_file),
                "graph": G
            }
            
        except Exception as e:
            self.logger.error(f"Error creating visualizations: {e}")
            return None

def create_mosquito_virus_visualizations(viral_detection_results: Dict[str, List[Dict]], 
                                       mosquito_metadata: Dict[str, Dict],
                                       output_dir: str = ".") -> Optional[Dict]:
    """Factory function to create mosquito-virus visualizations."""
    visualizer = MosquitoVirusNetworkVisualizer(output_dir)
    return visualizer.save_visualizations(viral_detection_results, mosquito_metadata)