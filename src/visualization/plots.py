"""Visualization functions for viral interference research."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import networkx as nx
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')
from typing import Dict, List, Tuple, Optional, Union

import sys
from pathlib import Path
src_dir = Path(__file__).parent.parent
sys.path.insert(0, str(src_dir))

from utils import get_logger, get_config_value


class ViralVisualization:
    """Visualization tools for viral interference analysis."""
    
    def __init__(self):
        """Initialize the visualization module."""
        self.logger = get_logger('visualization')
        
        # Load configuration
        self.figure_width = get_config_value('visualization.figure_width', 12)
        self.figure_height = get_config_value('visualization.figure_height', 8)
        self.color_palette = get_config_value('visualization.color_palette', 'viridis')
        self.dpi = get_config_value('visualization.dpi', 300)
        
        # Set style
        plt.style.use('seaborn-v0_8')
        sns.set_palette(self.color_palette)
        
        self.logger.info("Viral visualization module initialized")
    
    def plot_abundance_heatmap(self, 
                             abundance_matrix: pd.DataFrame,
                             metadata: pd.DataFrame = None,
                             log_transform: bool = True,
                             cluster_samples: bool = False,
                             cluster_viruses: bool = True,
                             interactive: bool = True) -> Union[plt.Figure, go.Figure]:
        """Create abundance heatmap.
        
        Args:
            abundance_matrix: Viral abundance matrix
            metadata: Sample metadata
            log_transform: Whether to log-transform abundances
            cluster_samples: Whether to cluster samples
            cluster_viruses: Whether to cluster viruses
            interactive: Whether to create interactive plot
            
        Returns:
            Matplotlib or Plotly figure
        """
        self.logger.info("Creating abundance heatmap")
        
        # Prepare data
        plot_data = abundance_matrix.copy()
        
        if log_transform:
            plot_data = np.log10(plot_data + 1)
            title_suffix = " (Log10 Transformed)"
            colorbar_title = "Log10(Abundance + 1)"
        else:
            title_suffix = ""
            colorbar_title = "Abundance"
        
        # Clustering
        if cluster_viruses:
            virus_linkage = linkage(plot_data.T, method='ward')
            virus_order = dendrogram(virus_linkage, no_plot=True)['leaves']
            plot_data = plot_data.iloc[:, virus_order]
        
        if cluster_samples:
            sample_linkage = linkage(plot_data, method='ward')
            sample_order = dendrogram(sample_linkage, no_plot=True)['leaves']
            plot_data = plot_data.iloc[sample_order, :]
        
        if interactive:
            # Plotly heatmap
            fig = px.imshow(
                plot_data.T,  # Transpose for viruses on y-axis
                aspect='auto',
                color_continuous_scale=self.color_palette,
                title=f"Viral Abundance Heatmap{title_suffix}",
                labels={'x': 'Samples', 'y': 'Viruses', 'color': colorbar_title}
            )
            
            fig.update_layout(
                width=self.figure_width * 100,
                height=self.figure_height * 100,
                xaxis_title="Samples",
                yaxis_title="Viruses"
            )
            
            return fig
        
        else:
            # Matplotlib heatmap
            fig, ax = plt.subplots(figsize=(self.figure_width, self.figure_height))
            
            sns.heatmap(
                plot_data.T,
                cmap=self.color_palette,
                cbar_kws={'label': colorbar_title},
                ax=ax
            )
            
            ax.set_title(f"Viral Abundance Heatmap{title_suffix}")
            ax.set_xlabel("Samples")
            ax.set_ylabel("Viruses")
            
            plt.tight_layout()
            return fig
    
    def plot_cooccurrence_matrix(self, 
                               pairwise_results: pd.DataFrame,
                               significance_threshold: float = 0.05,
                               interactive: bool = True) -> Union[plt.Figure, go.Figure]:
        """Create co-occurrence matrix visualization.
        
        Args:
            pairwise_results: Results from pairwise statistical tests
            significance_threshold: P-value threshold for significance
            interactive: Whether to create interactive plot
            
        Returns:
            Matplotlib or Plotly figure
        """
        self.logger.info("Creating co-occurrence matrix")
        
        # Get unique viruses
        viruses = sorted(set(pairwise_results['virus_a'].tolist() + pairwise_results['virus_b'].tolist()))
        
        # Create matrices
        odds_ratio_matrix = pd.DataFrame(1.0, index=viruses, columns=viruses)
        pvalue_matrix = pd.DataFrame(1.0, index=viruses, columns=viruses)
        significance_matrix = pd.DataFrame(False, index=viruses, columns=viruses)
        
        # Fill matrices
        for _, row in pairwise_results.iterrows():
            virus_a, virus_b = row['virus_a'], row['virus_b']
            odds_ratio = row['odds_ratio']
            p_adj = row['p_adj']
            significant = p_adj < significance_threshold
            
            # Symmetric matrix
            odds_ratio_matrix.loc[virus_a, virus_b] = odds_ratio
            odds_ratio_matrix.loc[virus_b, virus_a] = odds_ratio
            pvalue_matrix.loc[virus_a, virus_b] = p_adj
            pvalue_matrix.loc[virus_b, virus_a] = p_adj
            significance_matrix.loc[virus_a, virus_b] = significant
            significance_matrix.loc[virus_b, virus_a] = significant
        
        if interactive:
            # Log-transform odds ratios for better visualization
            log_odds_matrix = np.log2(odds_ratio_matrix)
            
            # Create custom colorscale (blue for exclusion, red for co-occurrence)
            colorscale = [
                [0, 'blue'],
                [0.5, 'white'],
                [1, 'red']
            ]
            
            # Create heatmap
            fig = go.Figure(data=go.Heatmap(
                z=log_odds_matrix.values,
                x=log_odds_matrix.columns,
                y=log_odds_matrix.index,
                colorscale=colorscale,
                zmid=0,
                colorbar=dict(title="Log2(Odds Ratio)"),
                hovertemplate='Virus A: %{y}<br>Virus B: %{x}<br>Log2(OR): %{z:.2f}<extra></extra>'
            ))
            
            # Add significance markers
            sig_x, sig_y = np.where(significance_matrix.values)
            if len(sig_x) > 0:
                fig.add_trace(go.Scatter(
                    x=[viruses[j] for j in sig_y],
                    y=[viruses[i] for i in sig_x],
                    mode='markers',
                    marker=dict(symbol='star', size=8, color='black'),
                    name='Significant',
                    showlegend=True
                ))
            
            fig.update_layout(
                title="Viral Co-occurrence Matrix (Log2 Odds Ratios)",
                xaxis_title="Virus B",
                yaxis_title="Virus A",
                width=self.figure_width * 80,
                height=self.figure_height * 80
            )
            
            return fig
        
        else:
            # Matplotlib version
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(self.figure_width * 2, self.figure_height))
            
            # Odds ratio heatmap
            log_odds_matrix = np.log2(odds_ratio_matrix)
            sns.heatmap(
                log_odds_matrix,
                cmap='RdBu_r',
                center=0,
                annot=significance_matrix,
                fmt='',
                cbar_kws={'label': 'Log2(Odds Ratio)'},
                ax=ax1
            )
            ax1.set_title("Log2(Odds Ratios)")
            
            # P-value heatmap
            sns.heatmap(
                -np.log10(pvalue_matrix),
                cmap='viridis',
                cbar_kws={'label': '-Log10(P-value)'},
                ax=ax2
            )
            ax2.set_title("-Log10(Adjusted P-values)")
            
            plt.tight_layout()
            return fig
    
    def plot_interaction_network(self, 
                               network: nx.Graph,
                               layout: str = 'spring',
                               node_size_attr: str = 'prevalence',
                               edge_width_attr: str = 'weight',
                               color_by_type: bool = True,
                               interactive: bool = True) -> Union[plt.Figure, go.Figure]:
        """Visualize viral interaction network.
        
        Args:
            network: NetworkX graph
            layout: Layout algorithm ('spring', 'circular', 'kamada_kawai')
            node_size_attr: Node attribute for sizing
            edge_width_attr: Edge attribute for width
            color_by_type: Whether to color edges by interaction type
            interactive: Whether to create interactive plot
            
        Returns:
            Matplotlib or Plotly figure
        """
        self.logger.info("Creating interaction network visualization")
        
        if len(network.nodes()) == 0:
            self.logger.warning("Empty network provided")
            return None
        
        # Get layout positions
        if layout == 'spring':
            pos = nx.spring_layout(network, k=1, iterations=50)
        elif layout == 'circular':
            pos = nx.circular_layout(network)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(network)
        else:
            pos = nx.spring_layout(network)
        
        if interactive:
            # Plotly network
            fig = go.Figure()
            
            # Add edges
            edge_colors = []
            edge_widths = []
            
            for edge in network.edges(data=True):
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                
                # Edge color based on interaction type
                if color_by_type and 'association_type' in edge[2]:
                    if edge[2]['association_type'] == 'negative':
                        color = 'red'
                    elif edge[2]['association_type'] == 'positive':
                        color = 'blue'
                    else:
                        color = 'gray'
                else:
                    color = 'gray'
                
                # Edge width
                if edge_width_attr in edge[2]:
                    width = min(edge[2][edge_width_attr] * 5, 10)
                else:
                    width = 2
                
                fig.add_trace(go.Scatter(
                    x=[x0, x1, None],
                    y=[y0, y1, None],
                    mode='lines',
                    line=dict(color=color, width=width),
                    showlegend=False,
                    hoverinfo='none'
                ))
            
            # Add nodes
            node_x = [pos[node][0] for node in network.nodes()]
            node_y = [pos[node][1] for node in network.nodes()]
            node_text = list(network.nodes())
            
            # Node sizes
            if node_size_attr and any(node_size_attr in network.nodes[node] for node in network.nodes()):
                node_sizes = [network.nodes[node].get(node_size_attr, 0.1) * 100 for node in network.nodes()]
            else:
                node_sizes = [20] * len(network.nodes())
            
            # Node hover text
            hover_text = []
            for node in network.nodes():
                attrs = network.nodes[node]
                hover_info = f"Virus: {node}<br>"
                for attr, value in attrs.items():
                    if isinstance(value, float):
                        hover_info += f"{attr}: {value:.3f}<br>"
                    else:
                        hover_info += f"{attr}: {value}<br>"
                hover_text.append(hover_info)
            
            fig.add_trace(go.Scatter(
                x=node_x,
                y=node_y,
                mode='markers+text',
                marker=dict(
                    size=node_sizes,
                    color='lightblue',
                    line=dict(width=2, color='black')
                ),
                text=node_text,
                textposition='middle center',
                textfont=dict(size=10),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False
            ))
            
            fig.update_layout(
                title="Viral Interaction Network",
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                annotations=[
                    dict(
                        text="Node size: prevalence | Blue edges: co-occurrence | Red edges: exclusion",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002,
                        xanchor='left', yanchor='bottom',
                        font=dict(size=12)
                    )
                ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                width=self.figure_width * 80,
                height=self.figure_height * 80
            )
            
            return fig
        
        else:
            # Matplotlib network
            fig, ax = plt.subplots(figsize=(self.figure_width, self.figure_height))
            
            # Node sizes
            if node_size_attr and any(node_size_attr in network.nodes[node] for node in network.nodes()):
                node_sizes = [network.nodes[node].get(node_size_attr, 0.1) * 1000 for node in network.nodes()]
            else:
                node_sizes = 300
            
            # Edge colors and widths
            edge_colors = []
            edge_widths = []
            
            for edge in network.edges(data=True):
                if color_by_type and 'association_type' in edge[2]:
                    if edge[2]['association_type'] == 'negative':
                        edge_colors.append('red')
                    elif edge[2]['association_type'] == 'positive':
                        edge_colors.append('blue')
                    else:
                        edge_colors.append('gray')
                else:
                    edge_colors.append('gray')
                
                if edge_width_attr in edge[2]:
                    edge_widths.append(min(edge[2][edge_width_attr] * 3, 5))
                else:
                    edge_widths.append(1)
            
            # Draw network
            nx.draw_networkx_nodes(
                network, pos, 
                node_size=node_sizes,
                node_color='lightblue',
                alpha=0.7,
                ax=ax
            )
            
            nx.draw_networkx_edges(
                network, pos,
                edge_color=edge_colors,
                width=edge_widths,
                alpha=0.6,
                ax=ax
            )
            
            nx.draw_networkx_labels(
                network, pos,
                font_size=8,
                ax=ax
            )
            
            ax.set_title("Viral Interaction Network")
            ax.axis('off')
            
            # Legend
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color='blue', lw=2, label='Co-occurrence'),
                Line2D([0], [0], color='red', lw=2, label='Exclusion'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='lightblue', 
                      markersize=10, label='Virus')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            
            plt.tight_layout()
            return fig
    
    def plot_prevalence_comparison(self, 
                                 abundance_matrix: pd.DataFrame,
                                 metadata: pd.DataFrame,
                                 group_by: str,
                                 interactive: bool = True) -> Union[plt.Figure, go.Figure]:
        """Compare viral prevalence across groups.
        
        Args:
            abundance_matrix: Viral abundance matrix
            metadata: Sample metadata
            group_by: Column name to group by
            interactive: Whether to create interactive plot
            
        Returns:
            Matplotlib or Plotly figure
        """
        self.logger.info(f"Creating prevalence comparison by {group_by}")
        
        if group_by not in metadata.columns:
            raise ValueError(f"Column '{group_by}' not found in metadata")
        
        # Calculate prevalence by group
        prevalence_data = []
        
        for group in metadata[group_by].unique():
            group_samples = metadata[metadata[group_by] == group]['sample_id']
            if hasattr(group_samples, 'values'):
                group_samples = group_samples.values
            
            # Filter abundance matrix for group samples
            group_abundance = abundance_matrix.loc[abundance_matrix.index.isin(group_samples)]
            prevalence = (group_abundance > 0).mean()
            
            for virus in prevalence.index:
                prevalence_data.append({
                    'Group': group,
                    'Virus': virus,
                    'Prevalence': prevalence[virus]
                })
        
        prevalence_df = pd.DataFrame(prevalence_data)
        
        if interactive:
            # Plotly grouped bar chart
            fig = px.bar(
                prevalence_df,
                x='Virus',
                y='Prevalence',
                color='Group',
                title=f'Viral Prevalence by {group_by}',
                barmode='group'
            )
            
            fig.update_layout(
                xaxis_tickangle=-45,
                width=self.figure_width * 80,
                height=self.figure_height * 80,
                yaxis_title="Prevalence",
                xaxis_title="Virus"
            )
            
            return fig
        
        else:
            # Matplotlib grouped bar chart
            fig, ax = plt.subplots(figsize=(self.figure_width, self.figure_height))
            
            # Pivot for grouped bar chart
            pivot_df = prevalence_df.pivot(index='Virus', columns='Group', values='Prevalence')
            
            pivot_df.plot(kind='bar', ax=ax)
            ax.set_title(f'Viral Prevalence by {group_by}')
            ax.set_ylabel('Prevalence')
            ax.set_xlabel('Virus')
            ax.legend(title=group_by)
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            return fig
    
    def plot_temporal_trends(self, 
                           abundance_matrix: pd.DataFrame,
                           metadata: pd.DataFrame,
                           date_column: str = 'collection_date',
                           viruses: List[str] = None,
                           interactive: bool = True) -> Union[plt.Figure, go.Figure]:
        """Plot temporal trends in viral prevalence.
        
        Args:
            abundance_matrix: Viral abundance matrix
            metadata: Sample metadata
            date_column: Column name with dates
            viruses: List of viruses to plot (default: all)
            interactive: Whether to create interactive plot
            
        Returns:
            Matplotlib or Plotly figure
        """
        self.logger.info("Creating temporal trends plot")
        
        if date_column not in metadata.columns:
            raise ValueError(f"Column '{date_column}' not found in metadata")
        
        # Convert to datetime
        metadata_temp = metadata.copy()
        metadata_temp[date_column] = pd.to_datetime(metadata_temp[date_column])
        
        if viruses is None:
            viruses = abundance_matrix.columns.tolist()[:5]  # Limit to first 5 for readability
        
        # Aggregate by date
        temporal_data = []
        
        for date in metadata_temp[date_column].dt.date.unique():
            date_samples = metadata_temp[metadata_temp[date_column].dt.date == date]['sample_id']
            if hasattr(date_samples, 'values'):
                date_samples = date_samples.values
            
            date_abundance = abundance_matrix.loc[abundance_matrix.index.isin(date_samples)]
            
            for virus in viruses:
                if virus in date_abundance.columns:
                    prevalence = (date_abundance[virus] > 0).mean()
                    mean_abundance = date_abundance[virus].mean()
                    
                    temporal_data.append({
                        'Date': date,
                        'Virus': virus,
                        'Prevalence': prevalence,
                        'Mean_Abundance': mean_abundance,
                        'Sample_Count': len(date_abundance)
                    })
        
        temporal_df = pd.DataFrame(temporal_data).sort_values('Date')
        
        if interactive:
            # Plotly line plot
            fig = px.line(
                temporal_df,
                x='Date',
                y='Prevalence',
                color='Virus',
                title='Temporal Trends in Viral Prevalence',
                hover_data=['Mean_Abundance', 'Sample_Count']
            )
            
            fig.update_layout(
                width=self.figure_width * 80,
                height=self.figure_height * 80,
                yaxis_title="Prevalence",
                xaxis_title="Date"
            )
            
            return fig
        
        else:
            # Matplotlib line plot
            fig, ax = plt.subplots(figsize=(self.figure_width, self.figure_height))
            
            for virus in viruses:
                virus_data = temporal_df[temporal_df['Virus'] == virus]
                ax.plot(virus_data['Date'], virus_data['Prevalence'], 
                       marker='o', label=virus, linewidth=2)
            
            ax.set_title('Temporal Trends in Viral Prevalence')
            ax.set_ylabel('Prevalence')
            ax.set_xlabel('Date')
            ax.legend()
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            return fig
    
    def plot_clustering_results(self, 
                              clustering_results: Dict,
                              abundance_matrix: pd.DataFrame,
                              interactive: bool = True) -> Union[plt.Figure, go.Figure]:
        """Visualize clustering results.
        
        Args:
            clustering_results: Results from clustering analysis
            abundance_matrix: Original abundance matrix
            interactive: Whether to create interactive plot
            
        Returns:
            Matplotlib or Plotly figure
        """
        self.logger.info("Creating clustering visualization")
        
        if 'dimensionality_reduction' not in clustering_results:
            raise ValueError("Dimensionality reduction results not found")
        
        cluster_labels = clustering_results['cluster_labels']
        
        if interactive:
            # Use PCA coordinates if available
            if 'pca' in clustering_results['dimensionality_reduction']:
                coords = clustering_results['dimensionality_reduction']['pca']['coordinates']
                x_label = f"PC1 ({clustering_results['dimensionality_reduction']['pca']['explained_variance_ratio'][0]:.1%})"
                y_label = f"PC2 ({clustering_results['dimensionality_reduction']['pca']['explained_variance_ratio'][1]:.1%})"
                title = "Sample Clustering (PCA)"
            elif 'tsne' in clustering_results['dimensionality_reduction']:
                coords = clustering_results['dimensionality_reduction']['tsne']['coordinates']
                x_label = "t-SNE 1"
                y_label = "t-SNE 2"
                title = "Sample Clustering (t-SNE)"
            else:
                raise ValueError("No suitable coordinates found for visualization")
            
            # Create scatter plot
            plot_df = pd.DataFrame({
                'x': coords[:, 0],
                'y': coords[:, 1],
                'cluster': [f'Cluster {c}' if c != -1 else 'Noise' for c in cluster_labels],
                'sample_id': abundance_matrix.index
            })
            
            fig = px.scatter(
                plot_df,
                x='x',
                y='y',
                color='cluster',
                hover_data=['sample_id'],
                title=title
            )
            
            fig.update_layout(
                xaxis_title=x_label,
                yaxis_title=y_label,
                width=self.figure_width * 80,
                height=self.figure_height * 80
            )
            
            return fig
        
        else:
            # Matplotlib version
            fig, ax = plt.subplots(figsize=(self.figure_width, self.figure_height))
            
            # Use PCA coordinates
            if 'pca' in clustering_results['dimensionality_reduction']:
                coords = clustering_results['dimensionality_reduction']['pca']['coordinates']
                x_label = f"PC1 ({clustering_results['dimensionality_reduction']['pca']['explained_variance_ratio'][0]:.1%})"
                y_label = f"PC2 ({clustering_results['dimensionality_reduction']['pca']['explained_variance_ratio'][1]:.1%})"
                title = "Sample Clustering (PCA)"
            else:
                raise ValueError("PCA coordinates not found")
            
            # Create scatter plot
            unique_clusters = set(cluster_labels)
            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_clusters)))
            
            for i, cluster in enumerate(unique_clusters):
                mask = cluster_labels == cluster
                label = f'Cluster {cluster}' if cluster != -1 else 'Noise'
                ax.scatter(coords[mask, 0], coords[mask, 1], 
                          c=[colors[i]], label=label, alpha=0.7)
            
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(title)
            ax.legend()
            plt.tight_layout()
            
            return fig
    
    def create_summary_dashboard(self, 
                               abundance_matrix: pd.DataFrame,
                               analysis_results: Dict,
                               metadata: pd.DataFrame = None) -> go.Figure:
        """Create a comprehensive summary dashboard.
        
        Args:
            abundance_matrix: Viral abundance matrix
            analysis_results: Results from statistical analysis
            metadata: Sample metadata
            
        Returns:
            Plotly dashboard figure
        """
        self.logger.info("Creating summary dashboard")
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                'Viral Prevalence',
                'Co-occurrence Network',
                'Exclusion Candidates',
                'Sample Diversity'
            ],
            specs=[
                [{'type': 'bar'}, {'type': 'scatter'}],
                [{'type': 'bar'}, {'type': 'histogram'}]
            ]
        )
        
        # 1. Viral prevalence
        prevalence = (abundance_matrix > 0).mean().sort_values(ascending=True)
        fig.add_trace(
            go.Bar(
                x=prevalence.values,
                y=prevalence.index,
                orientation='h',
                name='Prevalence'
            ),
            row=1, col=1
        )
        
        # 2. Network (simplified)
        if 'network' in analysis_results and analysis_results['network'] is not None:
            network = analysis_results['network']
            if len(network.edges()) > 0:
                # Simple network representation
                pos = nx.spring_layout(network)
                
                # Add edges
                for edge in network.edges(data=True):
                    x0, y0 = pos[edge[0]]
                    x1, y1 = pos[edge[1]]
                    
                    color = 'red' if edge[2].get('association_type') == 'negative' else 'blue'
                    
                    fig.add_trace(
                        go.Scatter(
                            x=[x0, x1, None],
                            y=[y0, y1, None],
                            mode='lines',
                            line=dict(color=color, width=2),
                            showlegend=False
                        ),
                        row=1, col=2
                    )
                
                # Add nodes
                node_x = [pos[node][0] for node in network.nodes()]
                node_y = [pos[node][1] for node in network.nodes()]
                
                fig.add_trace(
                    go.Scatter(
                        x=node_x,
                        y=node_y,
                        mode='markers+text',
                        text=list(network.nodes()),
                        textposition='middle center',
                        marker=dict(size=10, color='lightblue'),
                        showlegend=False
                    ),
                    row=1, col=2
                )
        
        # 3. Exclusion candidates
        if 'exclusion_candidates' in analysis_results and analysis_results['exclusion_candidates']:
            exclusion_data = analysis_results['exclusion_candidates'][:10]  # Top 10
            
            pair_names = [f"{pair['virus_pair'][0]} vs {pair['virus_pair'][1]}" 
                         for pair in exclusion_data]
            exclusion_strengths = [pair['exclusion_strength'] for pair in exclusion_data]
            
            fig.add_trace(
                go.Bar(
                    x=exclusion_strengths,
                    y=pair_names,
                    orientation='h',
                    name='Exclusion Strength'
                ),
                row=2, col=1
            )
        
        # 4. Sample diversity
        viral_richness = (abundance_matrix > 0).sum(axis=1)
        fig.add_trace(
            go.Histogram(
                x=viral_richness,
                name='Viral Richness',
                nbinsx=20
            ),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            height=800,
            title_text="Viral Interference Analysis Dashboard",
            showlegend=False
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="Prevalence", row=1, col=1)
        fig.update_yaxes(title_text="Virus", row=1, col=1)
        fig.update_xaxes(title_text="Exclusion Strength", row=2, col=1)
        fig.update_yaxes(title_text="Virus Pair", row=2, col=1)
        fig.update_xaxes(title_text="Viral Richness", row=2, col=2)
        fig.update_yaxes(title_text="Count", row=2, col=2)
        
        return fig
    
    def save_figure(self, fig: Union[plt.Figure, go.Figure], 
                   filename: str, format: str = 'png') -> None:
        """Save figure to file.
        
        Args:
            fig: Figure to save
            filename: Output filename
            format: Output format ('png', 'pdf', 'svg', 'html')
        """
        if isinstance(fig, plt.Figure):
            fig.savefig(filename, format=format, dpi=self.dpi, bbox_inches='tight')
        elif isinstance(fig, go.Figure):
            if format == 'html':
                fig.write_html(filename)
            else:
                fig.write_image(filename, format=format)
        
        self.logger.info(f"Figure saved to {filename}")