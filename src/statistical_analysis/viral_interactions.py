"""Statistical analysis of viral interactions and co-occurrence patterns."""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact, chi2_contingency
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import networkx as nx
from typing import Dict, List, Tuple, Optional, Union
import warnings

from ..utils import get_logger, log_processing_step, get_config_value


class ViralInteractionAnalyzer:
    """Analyze viral co-occurrence patterns and potential interference."""
    
    def __init__(self):
        """Initialize the viral interaction analyzer."""
        self.logger = get_logger('statistical_analysis')
        
        # Load configuration parameters
        self.presence_threshold = get_config_value('statistical_analysis.presence_threshold_rpm', 1.0)
        self.fisher_alpha = get_config_value('statistical_analysis.fisher_test_alpha', 0.05)
        self.correction_method = get_config_value('statistical_analysis.multiple_testing_correction', 'fdr_bh')
        self.network_threshold = get_config_value('statistical_analysis.network_threshold', 0.3)
        
        self.logger.info("Viral interaction analyzer initialized")
    
    def analyze_cooccurrence(self, 
                           abundance_matrix: pd.DataFrame,
                           metadata: pd.DataFrame = None) -> Dict:
        """Analyze viral co-occurrence patterns across samples.
        
        Args:
            abundance_matrix: DataFrame with samples as rows, viruses as columns
            metadata: Optional metadata DataFrame with sample information
            
        Returns:
            Dictionary containing co-occurrence analysis results
        """
        self.logger.info("Starting viral co-occurrence analysis")
        
        results = {
            'pairwise_tests': None,
            'network': None,
            'summary_stats': {},
            'exclusion_candidates': [],
            'cooccurrence_candidates': []
        }
        
        with log_processing_step("Converting to presence/absence", self.logger):
            presence_matrix = self._convert_to_presence_absence(abundance_matrix)
            results['summary_stats']['presence_matrix_shape'] = presence_matrix.shape
            results['summary_stats']['total_virus_detections'] = presence_matrix.sum().sum()
        
        with log_processing_step("Pairwise statistical tests", self.logger):
            pairwise_results = self._pairwise_fisher_tests(presence_matrix)
            results['pairwise_tests'] = pairwise_results
        
        with log_processing_step("Multiple testing correction", self.logger):
            corrected_results = self._correct_multiple_testing(pairwise_results)
            results['pairwise_tests'] = corrected_results
        
        with log_processing_step("Identifying candidates", self.logger):
            exclusion_candidates = self._identify_exclusion_candidates(corrected_results)
            cooccurrence_candidates = self._identify_cooccurrence_candidates(corrected_results)
            results['exclusion_candidates'] = exclusion_candidates
            results['cooccurrence_candidates'] = cooccurrence_candidates
        
        with log_processing_step("Building interaction network", self.logger):
            network = self._build_interaction_network(corrected_results)
            results['network'] = network
        
        # Add summary statistics
        results['summary_stats'].update({
            'num_virus_pairs_tested': len(pairwise_results),
            'num_significant_associations': len(corrected_results[corrected_results['p_adj'] < self.fisher_alpha]),
            'num_exclusion_candidates': len(exclusion_candidates),
            'num_cooccurrence_candidates': len(cooccurrence_candidates)
        })
        
        self.logger.info(f"Co-occurrence analysis completed. Found {len(exclusion_candidates)} exclusion candidates")
        
        return results
    
    def _convert_to_presence_absence(self, abundance_matrix: pd.DataFrame) -> pd.DataFrame:
        """Convert abundance matrix to binary presence/absence.
        
        Args:
            abundance_matrix: Abundance data
            
        Returns:
            Binary presence/absence matrix
        """
        return (abundance_matrix >= self.presence_threshold).astype(int)
    
    def _pairwise_fisher_tests(self, presence_matrix: pd.DataFrame) -> pd.DataFrame:
        """Perform pairwise Fisher's exact tests for all virus pairs.
        
        Args:
            presence_matrix: Binary presence/absence matrix
            
        Returns:
            DataFrame with test results for each virus pair
        """
        viruses = presence_matrix.columns.tolist()
        results = []
        
        for i, virus_a in enumerate(viruses):
            for j, virus_b in enumerate(viruses[i+1:], i+1):
                # Create contingency table
                a_present = presence_matrix[virus_a]
                b_present = presence_matrix[virus_b]
                
                # 2x2 contingency table:
                # [[both_present, a_only], [b_only, neither]]
                both_present = ((a_present == 1) & (b_present == 1)).sum()
                a_only = ((a_present == 1) & (b_present == 0)).sum()
                b_only = ((a_present == 0) & (b_present == 1)).sum()
                neither = ((a_present == 0) & (b_present == 0)).sum()
                
                contingency_table = np.array([[both_present, a_only], [b_only, neither]])
                
                # Perform Fisher's exact test
                try:
                    odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')
                    
                    # Calculate additional statistics
                    total_samples = presence_matrix.shape[0]
                    a_prevalence = a_present.sum() / total_samples
                    b_prevalence = b_present.sum() / total_samples
                    
                    # Expected co-occurrence under independence
                    expected_cooccurrence = a_prevalence * b_prevalence * total_samples
                    observed_cooccurrence = both_present
                    
                    # Association type
                    if odds_ratio < 1:
                        association_type = 'negative'  # Potential exclusion
                    elif odds_ratio > 1:
                        association_type = 'positive'  # Co-occurrence
                    else:
                        association_type = 'independent'
                    
                    results.append({
                        'virus_a': virus_a,
                        'virus_b': virus_b,
                        'both_present': both_present,
                        'a_only': a_only,
                        'b_only': b_only,
                        'neither': neither,
                        'odds_ratio': odds_ratio,
                        'p_value': p_value,
                        'a_prevalence': a_prevalence,
                        'b_prevalence': b_prevalence,
                        'observed_cooccurrence': observed_cooccurrence,
                        'expected_cooccurrence': expected_cooccurrence,
                        'association_type': association_type
                    })
                    
                except Exception as e:
                    self.logger.warning(f"Fisher test failed for {virus_a} vs {virus_b}: {e}")
                    continue
        
        return pd.DataFrame(results)
    
    def _correct_multiple_testing(self, pairwise_results: pd.DataFrame) -> pd.DataFrame:
        """Apply multiple testing correction to p-values.
        
        Args:
            pairwise_results: DataFrame with pairwise test results
            
        Returns:
            DataFrame with corrected p-values
        """
        if len(pairwise_results) == 0:
            return pairwise_results
        
        # Apply multiple testing correction
        rejected, p_adj, alpha_sidak, alpha_bonf = multipletests(
            pairwise_results['p_value'], 
            alpha=self.fisher_alpha, 
            method=self.correction_method
        )
        
        pairwise_results = pairwise_results.copy()
        pairwise_results['p_adj'] = p_adj
        pairwise_results['significant'] = rejected
        pairwise_results['alpha_bonferroni'] = alpha_bonf
        
        return pairwise_results
    
    def _identify_exclusion_candidates(self, pairwise_results: pd.DataFrame) -> List[Dict]:
        """Identify candidate virus pairs showing exclusion patterns.
        
        Args:
            pairwise_results: DataFrame with corrected test results
            
        Returns:
            List of exclusion candidate dictionaries
        """
        # Filter for significant negative associations (potential exclusions)
        exclusion_mask = (
            (pairwise_results['significant']) & 
            (pairwise_results['association_type'] == 'negative') &
            (pairwise_results['odds_ratio'] < 0.5)  # Strong negative association
        )
        
        exclusion_candidates = []
        
        for _, row in pairwise_results[exclusion_mask].iterrows():
            candidate = {
                'virus_pair': (row['virus_a'], row['virus_b']),
                'odds_ratio': row['odds_ratio'],
                'p_value': row['p_value'],
                'p_adj': row['p_adj'],
                'both_present': row['both_present'],
                'expected_cooccurrence': row['expected_cooccurrence'],
                'exclusion_strength': 1 - row['odds_ratio'],  # Higher = stronger exclusion
                'evidence_level': self._assess_evidence_level(row)
            }
            exclusion_candidates.append(candidate)
        
        # Sort by exclusion strength
        exclusion_candidates.sort(key=lambda x: x['exclusion_strength'], reverse=True)
        
        return exclusion_candidates
    
    def _identify_cooccurrence_candidates(self, pairwise_results: pd.DataFrame) -> List[Dict]:
        """Identify candidate virus pairs showing co-occurrence patterns.
        
        Args:
            pairwise_results: DataFrame with corrected test results
            
        Returns:
            List of co-occurrence candidate dictionaries
        """
        # Filter for significant positive associations
        cooccurrence_mask = (
            (pairwise_results['significant']) & 
            (pairwise_results['association_type'] == 'positive') &
            (pairwise_results['odds_ratio'] > 2.0)  # Strong positive association
        )
        
        cooccurrence_candidates = []
        
        for _, row in pairwise_results[cooccurrence_mask].iterrows():
            candidate = {
                'virus_pair': (row['virus_a'], row['virus_b']),
                'odds_ratio': row['odds_ratio'],
                'p_value': row['p_value'],
                'p_adj': row['p_adj'],
                'both_present': row['both_present'],
                'expected_cooccurrence': row['expected_cooccurrence'],
                'cooccurrence_strength': row['odds_ratio'],
                'evidence_level': self._assess_evidence_level(row)
            }
            cooccurrence_candidates.append(candidate)
        
        # Sort by co-occurrence strength
        cooccurrence_candidates.sort(key=lambda x: x['cooccurrence_strength'], reverse=True)
        
        return cooccurrence_candidates
    
    def _assess_evidence_level(self, row: pd.Series) -> str:
        """Assess the evidence level for a viral interaction.
        
        Args:
            row: Row from pairwise results DataFrame
            
        Returns:
            Evidence level string
        """
        p_adj = row['p_adj']
        odds_ratio = row['odds_ratio']
        both_present = row['both_present']
        
        # Strong evidence criteria
        if p_adj < 0.001 and both_present >= 5:
            if odds_ratio < 0.1 or odds_ratio > 10:
                return 'strong'
            else:
                return 'moderate'
        elif p_adj < 0.01 and both_present >= 3:
            return 'moderate'
        elif p_adj < 0.05:
            return 'weak'
        else:
            return 'insufficient'
    
    def _build_interaction_network(self, pairwise_results: pd.DataFrame) -> nx.Graph:
        """Build a network of viral interactions.
        
        Args:
            pairwise_results: DataFrame with pairwise test results
            
        Returns:
            NetworkX graph of viral interactions
        """
        G = nx.Graph()
        
        # Add all viruses as nodes
        all_viruses = set(pairwise_results['virus_a'].tolist() + pairwise_results['virus_b'].tolist())
        G.add_nodes_from(all_viruses)
        
        # Add edges for significant associations
        significant_results = pairwise_results[pairwise_results['significant']]
        
        for _, row in significant_results.iterrows():
            virus_a = row['virus_a']
            virus_b = row['virus_b']
            
            # Edge attributes
            edge_attrs = {
                'weight': abs(np.log(row['odds_ratio'])),  # Log odds ratio magnitude
                'odds_ratio': row['odds_ratio'],
                'p_value': row['p_value'],
                'p_adj': row['p_adj'],
                'association_type': row['association_type'],
                'evidence_level': self._assess_evidence_level(row)
            }
            
            G.add_edge(virus_a, virus_b, **edge_attrs)
        
        return G
    
    def analyze_confounders(self, 
                          abundance_matrix: pd.DataFrame,
                          metadata: pd.DataFrame,
                          target_virus: str,
                          confounders: List[str] = None) -> Dict:
        """Analyze viral interactions while controlling for confounding factors.
        
        Args:
            abundance_matrix: Viral abundance matrix
            metadata: Sample metadata
            target_virus: Target virus to analyze
            confounders: List of confounder variables in metadata
            
        Returns:
            Dictionary with confounder analysis results
        """
        self.logger.info(f"Analyzing confounders for {target_virus}")
        
        if confounders is None:
            confounders = ['mosquito_species', 'location', 'collection_date']
        
        # Convert to presence/absence
        presence_matrix = self._convert_to_presence_absence(abundance_matrix)
        
        # Prepare data for logistic regression
        y = presence_matrix[target_virus]
        
        # Features: other viruses + confounders
        other_viruses = [v for v in presence_matrix.columns if v != target_virus]
        X_viral = presence_matrix[other_viruses]
        
        # Add confounder variables
        X_confounders = pd.DataFrame()
        for conf in confounders:
            if conf in metadata.columns:
                # Handle categorical variables
                if metadata[conf].dtype == 'object':
                    dummies = pd.get_dummies(metadata[conf], prefix=conf)
                    X_confounders = pd.concat([X_confounders, dummies], axis=1)
                else:
                    X_confounders[conf] = metadata[conf]
        
        # Combine viral and confounder features
        X = pd.concat([X_viral, X_confounders], axis=1)
        
        # Remove samples with missing data
        complete_mask = ~(X.isna().any(axis=1) | y.isna())
        X_clean = X[complete_mask]
        y_clean = y[complete_mask]
        
        results = {
            'target_virus': target_virus,
            'n_samples': len(y_clean),
            'viral_associations': {},
            'confounder_effects': {},
            'model_performance': {}
        }
        
        if len(y_clean) < 10 or y_clean.sum() < 3:
            self.logger.warning(f"Insufficient data for {target_virus} confounder analysis")
            return results
        
        try:
            # Fit logistic regression model
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_clean)
            
            model = LogisticRegression(random_state=42, max_iter=1000)
            model.fit(X_scaled, y_clean)
            
            # Extract results
            feature_names = X_clean.columns.tolist()
            coefficients = model.coef_[0]
            
            # Separate viral and confounder effects
            for i, feature in enumerate(feature_names):
                coef = coefficients[i]
                odds_ratio = np.exp(coef)
                
                if feature in other_viruses:
                    results['viral_associations'][feature] = {
                        'coefficient': coef,
                        'odds_ratio': odds_ratio,
                        'effect_direction': 'positive' if coef > 0 else 'negative'
                    }
                else:
                    results['confounder_effects'][feature] = {
                        'coefficient': coef,
                        'odds_ratio': odds_ratio,
                        'effect_direction': 'positive' if coef > 0 else 'negative'
                    }
            
            # Model performance
            train_score = model.score(X_scaled, y_clean)
            results['model_performance'] = {
                'accuracy': train_score,
                'n_features': len(feature_names),
                'intercept': model.intercept_[0]
            }
            
        except Exception as e:
            self.logger.error(f"Error in confounder analysis for {target_virus}: {e}")
            results['error'] = str(e)
        
        return results
    
    def permutation_test(self, 
                        abundance_matrix: pd.DataFrame,
                        virus_a: str,
                        virus_b: str,
                        n_permutations: int = 1000) -> Dict:
        """Perform permutation test for virus pair association.
        
        Args:
            abundance_matrix: Viral abundance matrix
            virus_a: First virus
            virus_b: Second virus
            n_permutations: Number of permutations
            
        Returns:
            Dictionary with permutation test results
        """
        self.logger.info(f"Performing permutation test for {virus_a} vs {virus_b}")
        
        presence_matrix = self._convert_to_presence_absence(abundance_matrix)
        
        # Observed test statistic (odds ratio)
        a_present = presence_matrix[virus_a]
        b_present = presence_matrix[virus_b]
        
        both_present = ((a_present == 1) & (b_present == 1)).sum()
        a_only = ((a_present == 1) & (b_present == 0)).sum()
        b_only = ((a_present == 0) & (b_present == 1)).sum()
        neither = ((a_present == 0) & (b_present == 0)).sum()
        
        contingency_table = np.array([[both_present, a_only], [b_only, neither]])
        observed_odds_ratio, _ = fisher_exact(contingency_table)
        
        # Permutation test
        permuted_odds_ratios = []
        
        for _ in range(n_permutations):
            # Shuffle virus B labels
            b_permuted = np.random.permutation(b_present)
            
            # Calculate permuted contingency table
            both_perm = ((a_present == 1) & (b_permuted == 1)).sum()
            a_only_perm = ((a_present == 1) & (b_permuted == 0)).sum()
            b_only_perm = ((a_present == 0) & (b_permuted == 1)).sum()
            neither_perm = ((a_present == 0) & (b_permuted == 0)).sum()
            
            perm_table = np.array([[both_perm, a_only_perm], [b_only_perm, neither_perm]])
            
            try:
                perm_odds_ratio, _ = fisher_exact(perm_table)
                permuted_odds_ratios.append(perm_odds_ratio)
            except:
                continue
        
        # Calculate p-value
        if observed_odds_ratio < 1:  # Testing for exclusion
            p_value = np.mean([or_perm <= observed_odds_ratio for or_perm in permuted_odds_ratios])
        else:  # Testing for co-occurrence
            p_value = np.mean([or_perm >= observed_odds_ratio for or_perm in permuted_odds_ratios])
        
        results = {
            'virus_a': virus_a,
            'virus_b': virus_b,
            'observed_odds_ratio': observed_odds_ratio,
            'permutation_p_value': p_value,
            'n_permutations': len(permuted_odds_ratios),
            'permuted_odds_ratios': permuted_odds_ratios
        }
        
        return results
    
    def export_results(self, results: Dict, output_path: str) -> None:
        """Export analysis results to files.
        
        Args:
            results: Analysis results dictionary
            output_path: Output directory path
        """
        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Export pairwise test results
        if results['pairwise_tests'] is not None:
            results['pairwise_tests'].to_csv(
                output_dir / 'pairwise_viral_associations.csv', 
                index=False
            )
        
        # Export exclusion candidates
        if results['exclusion_candidates']:
            exclusion_df = pd.DataFrame(results['exclusion_candidates'])
            exclusion_df.to_csv(
                output_dir / 'exclusion_candidates.csv', 
                index=False
            )
        
        # Export co-occurrence candidates
        if results['cooccurrence_candidates']:
            cooccurrence_df = pd.DataFrame(results['cooccurrence_candidates'])
            cooccurrence_df.to_csv(
                output_dir / 'cooccurrence_candidates.csv', 
                index=False
            )
        
        # Export network
        if results['network'] is not None:
            nx.write_gml(results['network'], output_dir / 'viral_interaction_network.gml')
        
        self.logger.info(f"Results exported to {output_dir}")