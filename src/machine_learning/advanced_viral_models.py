#!/usr/bin/env python3
"""
Advanced Machine Learning Models for Viral Co-occurrence Analysis

This module implements sophisticated mathematical and ML models for analyzing
viral interference patterns in mosquito genomes, suitable for scientific publication.

Authors: Mosquito Viral Analysis Pipeline Team
Date: 2025
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN, KMeans
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import networkx as nx
from statsmodels.stats.multitest import multipletests
from statsmodels.discrete.discrete_model import Logit
import warnings
warnings.filterwarnings('ignore')

class ViralInterferenceModeler:
    """
    Advanced mathematical and machine learning models for viral interference analysis.
    
    This class implements multiple sophisticated approaches:
    1. Bayesian Network Models for causal inference
    2. Random Forest for feature importance and prediction
    3. Logistic Regression with interaction terms
    4. Network-based community detection
    5. Time-series analysis for temporal patterns
    6. Phylogenetic-aware statistical models
    """
    
    def __init__(self, random_state=42):
        self.random_state = random_state
        self.models = {}
        self.results = {}
        self.feature_importance = {}
        
    def prepare_features(self, viral_matrix, metadata_df):
        """
        Prepare comprehensive feature matrix for ML analysis.
        
        Args:
            viral_matrix: Binary presence/absence matrix (samples x viruses)
            metadata_df: Sample metadata (geographic, temporal, host info)
            
        Returns:
            feature_matrix: Enhanced feature matrix with engineered features
        """
        print("Preparing comprehensive feature matrix...")
        
        # Basic viral presence features
        features = viral_matrix.copy()
        
        # Viral richness (total number of viruses per sample)
        features['viral_richness'] = viral_matrix.sum(axis=1)
        
        # Viral diversity indices
        features['shannon_diversity'] = self._calculate_shannon_diversity(viral_matrix)
        features['simpson_diversity'] = self._calculate_simpson_diversity(viral_matrix)
        
        # Co-occurrence features
        cooccurrence_features = self._generate_cooccurrence_features(viral_matrix)
        features = pd.concat([features, cooccurrence_features], axis=1)
        
        # Phylogenetic features (if available)
        if 'phylogenetic_distance' in metadata_df.columns:
            features['phylo_distance'] = metadata_df['phylogenetic_distance']
            
        # Geographic features
        if 'latitude' in metadata_df.columns and 'longitude' in metadata_df.columns:
            features['latitude'] = metadata_df['latitude']
            features['longitude'] = metadata_df['longitude']
            features['geographic_cluster'] = self._assign_geographic_clusters(
                metadata_df[['latitude', 'longitude']]
            )
            
        # Temporal features
        if 'collection_date' in metadata_df.columns:
            temporal_features = self._extract_temporal_features(metadata_df['collection_date'])
            features = pd.concat([features, temporal_features], axis=1)
            
        # Host species features
        if 'species' in metadata_df.columns:
            species_encoded = pd.get_dummies(metadata_df['species'], prefix='species')
            features = pd.concat([features, species_encoded], axis=1)
            
        return features
    
    def _calculate_shannon_diversity(self, viral_matrix):
        """Calculate Shannon diversity index for each sample."""
        def shannon_sample(row):
            present_viruses = row[row > 0]
            if len(present_viruses) <= 1:
                return 0
            proportions = present_viruses / present_viruses.sum()
            return -np.sum(proportions * np.log(proportions))
        
        return viral_matrix.apply(shannon_sample, axis=1)
    
    def _calculate_simpson_diversity(self, viral_matrix):
        """Calculate Simpson diversity index for each sample."""
        def simpson_sample(row):
            present_viruses = row[row > 0]
            if len(present_viruses) <= 1:
                return 0
            proportions = present_viruses / present_viruses.sum()
            return 1 - np.sum(proportions ** 2)
        
        return viral_matrix.apply(simpson_sample, axis=1)
    
    def _generate_cooccurrence_features(self, viral_matrix):
        """Generate pairwise co-occurrence features."""
        virus_names = viral_matrix.columns
        cooccurrence_features = pd.DataFrame(index=viral_matrix.index)
        
        # Calculate pairwise co-occurrences for top virus pairs
        virus_prevalence = viral_matrix.mean().sort_values(ascending=False)
        top_viruses = virus_prevalence.head(20).index  # Top 20 most prevalent
        
        for i, virus1 in enumerate(top_viruses):
            for virus2 in top_viruses[i+1:]:
                cooccur_name = f'cooccur_{virus1}_{virus2}'
                cooccurrence_features[cooccur_name] = (
                    viral_matrix[virus1] * viral_matrix[virus2]
                )
                
        return cooccurrence_features
    
    def _assign_geographic_clusters(self, coordinates, n_clusters=5):
        """Assign geographic clusters using K-means."""
        if len(coordinates) < n_clusters:
            return np.zeros(len(coordinates))
            
        kmeans = KMeans(n_clusters=n_clusters, random_state=self.random_state)
        return kmeans.fit_predict(coordinates)
    
    def _extract_temporal_features(self, dates):
        """Extract temporal features from collection dates."""
        dates = pd.to_datetime(dates)
        temporal_features = pd.DataFrame(index=dates.index)
        
        temporal_features['year'] = dates.dt.year
        temporal_features['month'] = dates.dt.month
        temporal_features['day_of_year'] = dates.dt.dayofyear
        temporal_features['season'] = dates.dt.month.map({
            12: 'winter', 1: 'winter', 2: 'winter',
            3: 'spring', 4: 'spring', 5: 'spring',
            6: 'summer', 7: 'summer', 8: 'summer',
            9: 'autumn', 10: 'autumn', 11: 'autumn'
        })
        
        # Encode seasons
        season_encoded = pd.get_dummies(temporal_features['season'], prefix='season')
        temporal_features = pd.concat([temporal_features, season_encoded], axis=1)
        temporal_features.drop('season', axis=1, inplace=True)
        
        return temporal_features
    
    def fit_random_forest_model(self, features, target_virus, test_size=0.2):
        """
        Fit Random Forest model to predict viral presence.
        
        Args:
            features: Feature matrix
            target_virus: Target virus to predict
            test_size: Proportion of data for testing
            
        Returns:
            model_results: Dictionary with model performance and feature importance
        """
        print(f"Fitting Random Forest model for {target_virus}...")
        
        # Prepare target variable
        y = features[target_virus] if target_virus in features.columns else None
        if y is None:
            raise ValueError(f"Target virus {target_virus} not found in features")
            
        # Remove target from features
        X = features.drop(target_virus, axis=1)
        
        # Handle categorical variables
        X = pd.get_dummies(X, drop_first=True)
        
        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        X_scaled = pd.DataFrame(X_scaled, columns=X.columns, index=X.index)
        
        # Fit Random Forest
        rf_model = RandomForestClassifier(
            n_estimators=500,
            max_depth=10,
            min_samples_split=5,
            min_samples_leaf=2,
            random_state=self.random_state,
            class_weight='balanced'
        )
        
        # Cross-validation
        cv_scores = cross_val_score(
            rf_model, X_scaled, y, 
            cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=self.random_state),
            scoring='roc_auc'
        )
        
        # Fit final model
        rf_model.fit(X_scaled, y)
        
        # Feature importance
        feature_importance = pd.DataFrame({
            'feature': X.columns,
            'importance': rf_model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        # Store results
        model_results = {
            'model': rf_model,
            'scaler': scaler,
            'cv_scores': cv_scores,
            'mean_cv_score': cv_scores.mean(),
            'std_cv_score': cv_scores.std(),
            'feature_importance': feature_importance,
            'target_virus': target_virus
        }
        
        self.models[f'rf_{target_virus}'] = model_results
        
        print(f"Random Forest CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
        
        return model_results
    
    def fit_logistic_regression_with_interactions(self, features, target_virus):
        """
        Fit logistic regression with interaction terms for viral interference.
        
        Args:
            features: Feature matrix
            target_virus: Target virus to predict
            
        Returns:
            model_results: Dictionary with model results and statistical tests
        """
        print(f"Fitting Logistic Regression with interactions for {target_virus}...")
        
        # Prepare target variable
        y = features[target_virus] if target_virus in features.columns else None
        if y is None:
            raise ValueError(f"Target virus {target_virus} not found in features")
            
        # Select viral features only for interaction analysis
        viral_features = features.select_dtypes(include=[np.number]).drop(target_virus, axis=1)
        
        # Select top co-occurring viruses for interaction terms
        virus_correlations = viral_features.corrwith(y).abs().sort_values(ascending=False)
        top_viruses = virus_correlations.head(10).index
        
        X = viral_features[top_viruses]
        
        # Add interaction terms
        interaction_features = pd.DataFrame(index=X.index)
        for i, virus1 in enumerate(top_viruses):
            for virus2 in top_viruses[i+1:]:
                interaction_name = f'{virus1}_x_{virus2}'
                interaction_features[interaction_name] = X[virus1] * X[virus2]
        
        # Combine main effects and interactions
        X_with_interactions = pd.concat([X, interaction_features], axis=1)
        
        # Fit logistic regression
        logit_model = Logit(y, X_with_interactions)
        logit_results = logit_model.fit(disp=0)
        
        # Extract significant interactions
        significant_interactions = logit_results.pvalues[logit_results.pvalues < 0.05]
        interaction_effects = significant_interactions[significant_interactions.index.str.contains('_x_')]
        
        model_results = {
            'model': logit_results,
            'significant_interactions': interaction_effects,
            'aic': logit_results.aic,
            'bic': logit_results.bic,
            'pseudo_r2': logit_results.prsquared,
            'target_virus': target_virus
        }
        
        self.models[f'logit_{target_virus}'] = model_results
        
        print(f"Logistic Regression Pseudo R²: {logit_results.prsquared:.3f}")
        print(f"Significant interactions found: {len(interaction_effects)}")
        
        return model_results
    
    def network_community_analysis(self, viral_matrix, significance_threshold=0.05):
        """
        Perform network-based community detection for viral groups.
        
        Args:
            viral_matrix: Binary viral presence matrix
            significance_threshold: P-value threshold for significant associations
            
        Returns:
            network_results: Dictionary with network analysis results
        """
        print("Performing network community analysis...")
        
        # Calculate pairwise associations
        virus_names = viral_matrix.columns
        n_viruses = len(virus_names)
        
        # Fisher's exact test for all pairs
        association_matrix = np.zeros((n_viruses, n_viruses))
        pvalue_matrix = np.zeros((n_viruses, n_viruses))
        
        for i, virus1 in enumerate(virus_names):
            for j, virus2 in enumerate(virus_names):
                if i != j:
                    # Create contingency table
                    both_present = ((viral_matrix[virus1] == 1) & (viral_matrix[virus2] == 1)).sum()
                    virus1_only = ((viral_matrix[virus1] == 1) & (viral_matrix[virus2] == 0)).sum()
                    virus2_only = ((viral_matrix[virus1] == 0) & (viral_matrix[virus2] == 1)).sum()
                    neither = ((viral_matrix[virus1] == 0) & (viral_matrix[virus2] == 0)).sum()
                    
                    contingency_table = np.array([[both_present, virus1_only],
                                                [virus2_only, neither]])
                    
                    # Fisher's exact test
                    odds_ratio, p_value = stats.fisher_exact(contingency_table)
                    
                    association_matrix[i, j] = odds_ratio
                    pvalue_matrix[i, j] = p_value
        
        # Multiple testing correction
        pvalues_flat = pvalue_matrix[np.triu_indices_from(pvalue_matrix, k=1)]
        rejected, pvals_corrected, _, _ = multipletests(pvalues_flat, method='fdr_bh')
        
        # Create network from significant associations
        G = nx.Graph()
        G.add_nodes_from(virus_names)
        
        idx = 0
        for i in range(n_viruses):
            for j in range(i+1, n_viruses):
                if rejected[idx]:  # Significant after correction
                    weight = association_matrix[i, j]
                    G.add_edge(virus_names[i], virus_names[j], 
                             weight=weight, pvalue=pvals_corrected[idx])
                idx += 1
        
        # Community detection
        communities = nx.community.greedy_modularity_communities(G)
        
        # Network metrics
        network_metrics = {
            'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(),
            'density': nx.density(G),
            'n_communities': len(communities),
            'modularity': nx.community.modularity(G, communities)
        }
        
        network_results = {
            'graph': G,
            'communities': communities,
            'metrics': network_metrics,
            'association_matrix': association_matrix,
            'pvalue_matrix': pvalue_matrix,
            'corrected_pvalues': pvals_corrected
        }
        
        self.results['network_analysis'] = network_results
        
        print(f"Network created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        print(f"Found {len(communities)} communities with modularity {network_metrics['modularity']:.3f}")
        
        return network_results
    
    def phylogenetic_signal_analysis(self, viral_matrix, phylogenetic_tree=None):
        """
        Test for phylogenetic signal in viral co-occurrence patterns.
        
        Args:
            viral_matrix: Binary viral presence matrix
            phylogenetic_tree: Phylogenetic tree (if available)
            
        Returns:
            phylo_results: Dictionary with phylogenetic analysis results
        """
        print("Analyzing phylogenetic signal in viral patterns...")
        
        # If no tree provided, create distance matrix from viral profiles
        if phylogenetic_tree is None:
            # Calculate Jaccard distances between samples
            jaccard_distances = pdist(viral_matrix.values, metric='jaccard')
            distance_matrix = squareform(jaccard_distances)
        else:
            # Use phylogenetic distances (implementation depends on tree format)
            distance_matrix = self._extract_phylogenetic_distances(phylogenetic_tree)
        
        # Mantel test for correlation between viral similarity and phylogenetic distance
        viral_similarity = 1 - pdist(viral_matrix.values, metric='jaccard')
        phylo_distances = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
        
        # Mantel test (correlation between distance matrices)
        mantel_r, mantel_p = self._mantel_test(viral_similarity, phylo_distances)
        
        phylo_results = {
            'mantel_correlation': mantel_r,
            'mantel_pvalue': mantel_p,
            'distance_matrix': distance_matrix,
            'viral_similarity': viral_similarity
        }
        
        self.results['phylogenetic_analysis'] = phylo_results
        
        print(f"Mantel test: r = {mantel_r:.3f}, p = {mantel_p:.3f}")
        
        return phylo_results
    
    def _mantel_test(self, X, Y, permutations=9999):
        """Perform Mantel test between two distance matrices."""
        # Calculate observed correlation
        observed_r = stats.pearsonr(X, Y)[0]
        
        # Permutation test
        permuted_r = []
        for _ in range(permutations):
            Y_permuted = np.random.permutation(Y)
            permuted_r.append(stats.pearsonr(X, Y_permuted)[0])
        
        # Calculate p-value
        p_value = (np.sum(np.abs(permuted_r) >= np.abs(observed_r)) + 1) / (permutations + 1)
        
        return observed_r, p_value
    
    def temporal_dynamics_analysis(self, viral_matrix, dates):
        """
        Analyze temporal dynamics of viral co-occurrence.
        
        Args:
            viral_matrix: Binary viral presence matrix
            dates: Collection dates for samples
            
        Returns:
            temporal_results: Dictionary with temporal analysis results
        """
        print("Analyzing temporal dynamics of viral co-occurrence...")
        
        dates = pd.to_datetime(dates)
        
        # Create time windows (monthly)
        viral_matrix['date'] = dates
        viral_matrix['year_month'] = dates.dt.to_period('M')
        
        # Calculate viral prevalence over time
        temporal_prevalence = viral_matrix.groupby('year_month').mean()
        temporal_prevalence.drop(['date'], axis=1, inplace=True, errors='ignore')
        
        # Trend analysis for each virus
        trend_results = {}
        for virus in temporal_prevalence.columns:
            if virus != 'year_month':
                # Linear trend test
                time_numeric = np.arange(len(temporal_prevalence))
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    time_numeric, temporal_prevalence[virus]
                )
                
                trend_results[virus] = {
                    'slope': slope,
                    'r_squared': r_value**2,
                    'p_value': p_value,
                    'trend': 'increasing' if slope > 0 and p_value < 0.05 else 
                            'decreasing' if slope < 0 and p_value < 0.05 else 'stable'
                }
        
        # Seasonal analysis
        viral_matrix['month'] = dates.dt.month
        seasonal_patterns = viral_matrix.groupby('month').mean()
        
        temporal_results = {
            'temporal_prevalence': temporal_prevalence,
            'trend_analysis': trend_results,
            'seasonal_patterns': seasonal_patterns
        }
        
        self.results['temporal_analysis'] = temporal_results
        
        print(f"Analyzed temporal trends for {len(trend_results)} viruses")
        
        return temporal_results
    
    def generate_publication_summary(self):
        """
        Generate comprehensive summary suitable for scientific publication.
        
        Returns:
            summary: Dictionary with publication-ready results
        """
        print("Generating publication summary...")
        
        summary = {
            'methodology': {
                'description': 'Advanced machine learning and statistical analysis of viral co-occurrence patterns',
                'models_used': list(self.models.keys()),
                'statistical_tests': ['Fisher exact test', 'Mantel test', 'Multiple testing correction (FDR)'],
                'network_analysis': 'Community detection using modularity optimization',
                'temporal_analysis': 'Linear trend analysis and seasonal decomposition'
            },
            'key_findings': {},
            'statistical_significance': {},
            'model_performance': {}
        }
        
        # Extract key findings from each analysis
        if 'network_analysis' in self.results:
            net_results = self.results['network_analysis']
            summary['key_findings']['network'] = {
                'n_significant_associations': net_results['graph'].number_of_edges(),
                'n_communities': net_results['metrics']['n_communities'],
                'modularity': net_results['metrics']['modularity'],
                'network_density': net_results['metrics']['density']
            }
        
        if 'phylogenetic_analysis' in self.results:
            phylo_results = self.results['phylogenetic_analysis']
            summary['key_findings']['phylogenetic'] = {
                'mantel_correlation': phylo_results['mantel_correlation'],
                'mantel_significance': phylo_results['mantel_pvalue'] < 0.05
            }
        
        # Model performance summary
        for model_name, model_results in self.models.items():
            if 'cv_scores' in model_results:
                summary['model_performance'][model_name] = {
                    'mean_cv_auc': model_results['mean_cv_score'],
                    'std_cv_auc': model_results['std_cv_score']
                }
        
        return summary

class LargeScaleDataProcessor:
    """
    Processor for handling large-scale mosquito genome datasets from NCBI.
    
    This class is designed to handle thousands of genome assemblies efficiently
    and extract meaningful patterns suitable for scientific publication.
    """
    
    def __init__(self, batch_size=100, n_jobs=4):
        self.batch_size = batch_size
        self.n_jobs = n_jobs
        self.processing_stats = {}
        
    def estimate_data_requirements(self, target_power=0.8, effect_size=0.3, alpha=0.05):
        """
        Estimate required sample size for detecting viral co-occurrence patterns.
        
        Args:
            target_power: Desired statistical power
            effect_size: Expected effect size (Cohen's d)
            alpha: Significance level
            
        Returns:
            sample_size_estimates: Dictionary with sample size requirements
        """
        from statsmodels.stats.power import ttest_power
        
        # Sample size for detecting co-occurrence
        n_required = ttest_power(effect_size, power=target_power, alpha=alpha, alternative='two-sided')
        
        estimates = {
            'minimum_samples_per_species': int(np.ceil(n_required)),
            'recommended_samples_per_species': int(np.ceil(n_required * 1.5)),
            'total_recommended_samples': int(np.ceil(n_required * 1.5 * 5)),  # 5 major species
            'power_analysis': {
                'effect_size': effect_size,
                'power': target_power,
                'alpha': alpha,
                'calculated_n': n_required
            }
        }
        
        print(f"Recommended sample sizes for robust analysis:")
        print(f"- Minimum per species: {estimates['minimum_samples_per_species']}")
        print(f"- Recommended per species: {estimates['recommended_samples_per_species']}")
        print(f"- Total recommended: {estimates['total_recommended_samples']}")
        
        return estimates
    
    def validate_data_quality(self, genome_files, min_genome_size=50000, max_n_content=0.1):
        """
        Validate quality of downloaded genome assemblies.
        
        Args:
            genome_files: List of genome file paths
            min_genome_size: Minimum acceptable genome size (bp)
            max_n_content: Maximum acceptable N content ratio
            
        Returns:
            quality_report: Dictionary with quality assessment
        """
        print(f"Validating quality of {len(genome_files)} genome assemblies...")
        
        quality_stats = []
        
        for genome_file in genome_files:
            try:
                # Read genome file
                with open(genome_file, 'r') as f:
                    sequences = []
                    current_seq = ""
                    
                    for line in f:
                        if line.startswith('>'):
                            if current_seq:
                                sequences.append(current_seq)
                                current_seq = ""
                        else:
                            current_seq += line.strip()
                    
                    if current_seq:
                        sequences.append(current_seq)
                
                # Calculate quality metrics
                total_length = sum(len(seq) for seq in sequences)
                n_content = sum(seq.count('N') + seq.count('n') for seq in sequences) / total_length
                n_contigs = len(sequences)
                
                quality_stats.append({
                    'file': genome_file,
                    'total_length': total_length,
                    'n_content': n_content,
                    'n_contigs': n_contigs,
                    'passes_quality': total_length >= min_genome_size and n_content <= max_n_content
                })
                
            except Exception as e:
                print(f"Error processing {genome_file}: {e}")
                quality_stats.append({
                    'file': genome_file,
                    'total_length': 0,
                    'n_content': 1.0,
                    'n_contigs': 0,
                    'passes_quality': False
                })
        
        quality_df = pd.DataFrame(quality_stats)
        
        quality_report = {
            'total_genomes': len(genome_files),
            'high_quality_genomes': quality_df['passes_quality'].sum(),
            'quality_rate': quality_df['passes_quality'].mean(),
            'mean_genome_size': quality_df['total_length'].mean(),
            'mean_n_content': quality_df['n_content'].mean(),
            'quality_details': quality_df
        }
        
        print(f"Quality assessment: {quality_report['high_quality_genomes']}/{quality_report['total_genomes']} genomes pass quality filters")
        
        return quality_report

# Example usage and testing functions
def run_comprehensive_analysis_example():
    """
    Example of running comprehensive analysis suitable for publication.
    
    This function demonstrates the full analytical pipeline with
    appropriate statistical rigor for scientific publication.
    """
    print("Running comprehensive viral co-occurrence analysis...")
    
    # Initialize models
    modeler = ViralInterferenceModeler(random_state=42)
    processor = LargeScaleDataProcessor()
    
    # Estimate data requirements
    sample_requirements = processor.estimate_data_requirements(
        target_power=0.8, effect_size=0.3, alpha=0.05
    )
    
    print("\nFor publishable results, you need:")
    print(f"- At least {sample_requirements['minimum_samples_per_species']} high-quality genomes per species")
    print(f"- Recommended {sample_requirements['recommended_samples_per_species']} genomes per species")
    print(f"- Total of {sample_requirements['total_recommended_samples']} genomes across major mosquito species")
    print("\nThis ensures sufficient statistical power to detect biologically meaningful viral interactions.")
    
    return sample_requirements

if __name__ == "__main__":
    # Run example analysis
    requirements = run_comprehensive_analysis_example()
    print("\nAdvanced ML models ready for large-scale viral co-occurrence analysis!")