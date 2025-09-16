#!/usr/bin/env python3
"""
Viral Co-occurrence Analysis Module

This module analyzes viral co-occurrence patterns in mosquito genomes to:
1. Identify viruses that frequently appear together
2. Detect rare viruses that appear in isolation
3. Analyze viral exclusion patterns
4. Generate co-occurrence networks and statistics

Based on the user's requirements for understanding viral interactions
and identifying rare vs. common viral combinations.
"""

import pandas as pd
import numpy as np
import logging
from typing import List, Dict, Set, Tuple, Optional
from pathlib import Path
from collections import defaultdict, Counter
from dataclasses import dataclass
import itertools
from scipy.stats import fisher_exact, chi2_contingency
from scipy.cluster.hierarchy import linkage, dendrogram
import json

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class ViralCooccurrence:
    """Represents a viral co-occurrence relationship."""
    virus_a: str
    virus_b: str
    cooccurrence_count: int
    total_samples: int
    virus_a_alone: int
    virus_b_alone: int
    neither_virus: int
    
    @property
    def cooccurrence_frequency(self) -> float:
        """Frequency of co-occurrence."""
        return self.cooccurrence_count / self.total_samples if self.total_samples > 0 else 0.0
    
    @property
    def jaccard_index(self) -> float:
        """Jaccard similarity coefficient."""
        union = self.cooccurrence_count + self.virus_a_alone + self.virus_b_alone
        return self.cooccurrence_count / union if union > 0 else 0.0
    
    @property
    def phi_coefficient(self) -> float:
        """Phi coefficient (correlation for binary variables)."""
        # 2x2 contingency table
        a = self.cooccurrence_count  # both present
        b = self.virus_a_alone       # only A present
        c = self.virus_b_alone       # only B present
        d = self.neither_virus       # neither present
        
        numerator = (a * d) - (b * c)
        denominator = np.sqrt((a + b) * (c + d) * (a + c) * (b + d))
        
        return numerator / denominator if denominator > 0 else 0.0
    
    def fisher_exact_test(self) -> Tuple[float, float]:
        """Fisher's exact test for association."""
        # 2x2 contingency table
        table = [[self.cooccurrence_count, self.virus_a_alone],
                 [self.virus_b_alone, self.neither_virus]]
        
        try:
            odds_ratio, p_value = fisher_exact(table)
            return odds_ratio, p_value
        except:
            return 1.0, 1.0

@dataclass
class RareVirusAnalysis:
    """Analysis results for rare viruses."""
    virus_name: str
    occurrence_count: int
    total_samples: int
    prevalence: float
    isolation_score: float  # How often it appears alone
    exclusion_partners: List[str]  # Viruses it rarely co-occurs with
    
    @property
    def is_rare(self) -> bool:
        """Check if virus is considered rare (< 5% prevalence)."""
        return self.prevalence < 0.05
    
    @property
    def is_isolated(self) -> bool:
        """Check if virus tends to appear in isolation."""
        return self.isolation_score > 0.7

class ViralCooccurrenceAnalyzer:
    """
    Analyzes viral co-occurrence patterns in mosquito genome data.
    
    This analyzer:
    1. Processes viral detection results from multiple genomes
    2. Builds co-occurrence matrices
    3. Identifies significant associations and exclusions
    4. Detects rare and isolated viruses
    5. Generates network representations
    """
    
    def __init__(self, output_dir: str = "cooccurrence_analysis"):
        """
        Initialize the co-occurrence analyzer.
        
        Args:
            output_dir: Directory for saving analysis results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Analysis results storage
        self.viral_presence_matrix = None
        self.cooccurrence_results = []
        self.rare_virus_analyses = []
        self.exclusion_patterns = {}
        
        logger.info(f"Initialized ViralCooccurrenceAnalyzer in {self.output_dir}")
    
    def load_viral_detection_results(self, results_file: Path) -> pd.DataFrame:
        """
        Load viral detection results from CSV file.
        
        Args:
            results_file: Path to viral detection summary CSV
            
        Returns:
            DataFrame with viral detection data
        """
        try:
            df = pd.read_csv(results_file)
            logger.info(f"Loaded {len(df)} viral detection results")
            return df
        except Exception as e:
            logger.error(f"Error loading results file {results_file}: {e}")
            return pd.DataFrame()
    
    def extract_viral_presence_data(self, detection_results: List) -> pd.DataFrame:
        """
        Extract viral presence/absence data from detection results.
        
        Args:
            detection_results: List of ViralDetectionResult objects
            
        Returns:
            DataFrame with samples as rows and viruses as columns (binary presence/absence)
        """
        # Collect all unique viruses and samples
        all_viruses = set()
        sample_virus_data = {}
        
        for result in detection_results:
            sample_id = f"{result.species}_{Path(result.genome_file).stem}"
            
            # Extract virus names from hits
            sample_viruses = set()
            for hit in result.viral_hits:
                # Clean virus name from subject ID
                virus_name = self._clean_virus_name(hit.subject_id)
                sample_viruses.add(virus_name)
                all_viruses.add(virus_name)
            
            sample_virus_data[sample_id] = sample_viruses
        
        # Create presence/absence matrix
        all_viruses = sorted(list(all_viruses))
        samples = sorted(list(sample_virus_data.keys()))
        
        presence_matrix = []
        for sample in samples:
            row = [1 if virus in sample_virus_data[sample] else 0 for virus in all_viruses]
            presence_matrix.append(row)
        
        df = pd.DataFrame(presence_matrix, index=samples, columns=all_viruses)
        
        logger.info(f"Created presence matrix: {len(samples)} samples × {len(all_viruses)} viruses")
        
        return df
    
    def _clean_virus_name(self, subject_id: str) -> str:
        """
        Clean and standardize virus names from BLAST subject IDs.
        
        Args:
            subject_id: Raw subject ID from BLAST results
            
        Returns:
            Cleaned virus name
        """
        # Remove common prefixes and suffixes
        name = subject_id.lower()
        
        # Extract meaningful virus name
        if 'virus' in name:
            # Try to extract virus name
            parts = name.split()
            virus_parts = []
            
            for i, part in enumerate(parts):
                virus_parts.append(part)
                if 'virus' in part:
                    break
            
            cleaned = ' '.join(virus_parts[:3])  # Take first 3 parts max
        else:
            # Take first few meaningful parts
            parts = name.replace('_', ' ').split()
            cleaned = ' '.join(parts[:2])
        
        # Capitalize properly
        return ' '.join(word.capitalize() for word in cleaned.split())
    
    def calculate_cooccurrence_matrix(self, presence_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate viral co-occurrence matrix.
        
        Args:
            presence_df: Viral presence/absence DataFrame
            
        Returns:
            Co-occurrence matrix (virus × virus)
        """
        # Calculate co-occurrence counts
        cooccurrence_matrix = presence_df.T.dot(presence_df)
        
        # Set diagonal to occurrence counts (not co-occurrence)
        np.fill_diagonal(cooccurrence_matrix.values, presence_df.sum())
        
        logger.info(f"Calculated co-occurrence matrix: {cooccurrence_matrix.shape}")
        
        return cooccurrence_matrix
    
    def analyze_pairwise_cooccurrence(self, presence_df: pd.DataFrame) -> List[ViralCooccurrence]:
        """
        Analyze pairwise viral co-occurrence relationships.
        
        Args:
            presence_df: Viral presence/absence DataFrame
            
        Returns:
            List of ViralCooccurrence objects
        """
        viruses = presence_df.columns.tolist()
        total_samples = len(presence_df)
        cooccurrences = []
        
        # Analyze all virus pairs
        for i, virus_a in enumerate(viruses):
            for j, virus_b in enumerate(viruses[i+1:], i+1):
                
                # Calculate contingency table
                both_present = ((presence_df[virus_a] == 1) & (presence_df[virus_b] == 1)).sum()
                a_only = ((presence_df[virus_a] == 1) & (presence_df[virus_b] == 0)).sum()
                b_only = ((presence_df[virus_a] == 0) & (presence_df[virus_b] == 1)).sum()
                neither = ((presence_df[virus_a] == 0) & (presence_df[virus_b] == 0)).sum()
                
                cooccurrence = ViralCooccurrence(
                    virus_a=virus_a,
                    virus_b=virus_b,
                    cooccurrence_count=both_present,
                    total_samples=total_samples,
                    virus_a_alone=a_only,
                    virus_b_alone=b_only,
                    neither_virus=neither
                )
                
                cooccurrences.append(cooccurrence)
        
        logger.info(f"Analyzed {len(cooccurrences)} pairwise co-occurrences")
        
        return cooccurrences
    
    def identify_significant_associations(self, cooccurrences: List[ViralCooccurrence], 
                                        p_threshold: float = 0.05) -> List[ViralCooccurrence]:
        """
        Identify statistically significant viral associations.
        
        Args:
            cooccurrences: List of ViralCooccurrence objects
            p_threshold: P-value threshold for significance
            
        Returns:
            List of significant associations
        """
        significant = []
        
        for cooc in cooccurrences:
            if cooc.cooccurrence_count > 0:  # Only test pairs that co-occur
                odds_ratio, p_value = cooc.fisher_exact_test()
                
                if p_value < p_threshold and odds_ratio > 1.0:
                    significant.append(cooc)
        
        # Sort by strength of association
        significant.sort(key=lambda x: x.phi_coefficient, reverse=True)
        
        logger.info(f"Found {len(significant)} significant viral associations")
        
        return significant
    
    def identify_viral_exclusions(self, cooccurrences: List[ViralCooccurrence], 
                                exclusion_threshold: float = 0.1) -> List[ViralCooccurrence]:
        """
        Identify viral exclusion patterns (viruses that rarely co-occur).
        
        Args:
            cooccurrences: List of ViralCooccurrence objects
            exclusion_threshold: Maximum co-occurrence frequency for exclusion
            
        Returns:
            List of viral exclusions
        """
        exclusions = []
        
        for cooc in cooccurrences:
            # Check if both viruses are present in dataset but rarely co-occur
            virus_a_present = cooc.cooccurrence_count + cooc.virus_a_alone > 0
            virus_b_present = cooc.cooccurrence_count + cooc.virus_b_alone > 0
            
            if (virus_a_present and virus_b_present and 
                cooc.cooccurrence_frequency < exclusion_threshold):
                
                # Test for significant negative association
                odds_ratio, p_value = cooc.fisher_exact_test()
                
                if p_value < 0.05 and odds_ratio < 1.0:
                    exclusions.append(cooc)
        
        # Sort by strength of exclusion
        exclusions.sort(key=lambda x: abs(x.phi_coefficient), reverse=True)
        
        logger.info(f"Found {len(exclusions)} viral exclusion patterns")
        
        return exclusions
    
    def analyze_rare_viruses(self, presence_df: pd.DataFrame, 
                           rare_threshold: float = 0.05) -> List[RareVirusAnalysis]:
        """
        Analyze rare virus patterns and isolation tendencies.
        
        Args:
            presence_df: Viral presence/absence DataFrame
            rare_threshold: Prevalence threshold for rare viruses
            
        Returns:
            List of RareVirusAnalysis objects
        """
        rare_analyses = []
        total_samples = len(presence_df)
        
        for virus in presence_df.columns:
            occurrence_count = presence_df[virus].sum()
            prevalence = occurrence_count / total_samples
            
            if prevalence < rare_threshold and occurrence_count > 0:
                # Calculate isolation score
                virus_samples = presence_df[presence_df[virus] == 1]
                
                if len(virus_samples) > 0:
                    # Count samples where this virus appears alone
                    isolation_count = 0
                    for idx in virus_samples.index:
                        other_viruses = virus_samples.loc[idx].drop(virus).sum()
                        if other_viruses == 0:
                            isolation_count += 1
                    
                    isolation_score = isolation_count / len(virus_samples)
                    
                    # Find exclusion partners
                    exclusion_partners = []
                    for other_virus in presence_df.columns:
                        if other_virus != virus:
                            cooccurrence = ((presence_df[virus] == 1) & 
                                          (presence_df[other_virus] == 1)).sum()
                            if cooccurrence == 0 and presence_df[other_virus].sum() > 0:
                                exclusion_partners.append(other_virus)
                    
                    rare_analysis = RareVirusAnalysis(
                        virus_name=virus,
                        occurrence_count=occurrence_count,
                        total_samples=total_samples,
                        prevalence=prevalence,
                        isolation_score=isolation_score,
                        exclusion_partners=exclusion_partners
                    )
                    
                    rare_analyses.append(rare_analysis)
        
        # Sort by rarity (lowest prevalence first)
        rare_analyses.sort(key=lambda x: x.prevalence)
        
        logger.info(f"Identified {len(rare_analyses)} rare viruses")
        
        return rare_analyses
    
    def generate_cooccurrence_network(self, significant_associations: List[ViralCooccurrence]) -> Dict:
        """
        Generate network representation of viral co-occurrences.
        
        Args:
            significant_associations: List of significant viral associations
            
        Returns:
            Network data structure (nodes and edges)
        """
        nodes = set()
        edges = []
        
        for assoc in significant_associations:
            nodes.add(assoc.virus_a)
            nodes.add(assoc.virus_b)
            
            edge = {
                'source': assoc.virus_a,
                'target': assoc.virus_b,
                'weight': assoc.phi_coefficient,
                'cooccurrence_count': assoc.cooccurrence_count,
                'jaccard_index': assoc.jaccard_index
            }
            edges.append(edge)
        
        network = {
            'nodes': [{'id': node, 'label': node} for node in nodes],
            'edges': edges
        }
        
        logger.info(f"Generated network: {len(nodes)} nodes, {len(edges)} edges")
        
        return network
    
    def run_complete_analysis(self, detection_results: List) -> Dict:
        """
        Run complete co-occurrence analysis pipeline.
        
        Args:
            detection_results: List of ViralDetectionResult objects
            
        Returns:
            Dictionary containing all analysis results
        """
        logger.info("Starting complete viral co-occurrence analysis")
        
        # Extract presence/absence data
        presence_df = self.extract_viral_presence_data(detection_results)
        self.viral_presence_matrix = presence_df
        
        if presence_df.empty:
            logger.warning("No viral presence data found")
            return {}
        
        # Calculate co-occurrence matrix
        cooccurrence_matrix = self.calculate_cooccurrence_matrix(presence_df)
        
        # Analyze pairwise co-occurrences
        all_cooccurrences = self.analyze_pairwise_cooccurrence(presence_df)
        self.cooccurrence_results = all_cooccurrences
        
        # Identify significant associations
        significant_associations = self.identify_significant_associations(all_cooccurrences)
        
        # Identify exclusion patterns
        exclusion_patterns = self.identify_viral_exclusions(all_cooccurrences)
        
        # Analyze rare viruses
        rare_viruses = self.analyze_rare_viruses(presence_df)
        self.rare_virus_analyses = rare_viruses
        
        # Generate network
        network = self.generate_cooccurrence_network(significant_associations)
        
        # Compile results
        results = {
            'presence_matrix': presence_df,
            'cooccurrence_matrix': cooccurrence_matrix,
            'significant_associations': significant_associations,
            'exclusion_patterns': exclusion_patterns,
            'rare_viruses': rare_viruses,
            'network': network,
            'summary_stats': self._generate_summary_stats(presence_df, all_cooccurrences, rare_viruses)
        }
        
        # Save results
        self._save_analysis_results(results)
        
        logger.info("Completed viral co-occurrence analysis")
        
        return results
    
    def _generate_summary_stats(self, presence_df: pd.DataFrame, 
                              cooccurrences: List[ViralCooccurrence],
                              rare_viruses: List[RareVirusAnalysis]) -> Dict:
        """
        Generate summary statistics for the analysis.
        
        Args:
            presence_df: Viral presence/absence DataFrame
            cooccurrences: List of all co-occurrence analyses
            rare_viruses: List of rare virus analyses
            
        Returns:
            Dictionary of summary statistics
        """
        total_samples = len(presence_df)
        total_viruses = len(presence_df.columns)
        
        # Virus prevalence statistics
        virus_counts = presence_df.sum()
        avg_prevalence = virus_counts.mean() / total_samples
        
        # Co-occurrence statistics
        positive_cooccurrences = [c for c in cooccurrences if c.cooccurrence_count > 0]
        avg_cooccurrence_strength = np.mean([c.phi_coefficient for c in positive_cooccurrences]) if positive_cooccurrences else 0
        
        # Sample diversity
        viruses_per_sample = presence_df.sum(axis=1)
        avg_viruses_per_sample = viruses_per_sample.mean()
        
        summary = {
            'total_samples': total_samples,
            'total_viruses': total_viruses,
            'average_viral_prevalence': avg_prevalence,
            'rare_viruses_count': len(rare_viruses),
            'positive_cooccurrences': len(positive_cooccurrences),
            'average_cooccurrence_strength': avg_cooccurrence_strength,
            'average_viruses_per_sample': avg_viruses_per_sample,
            'max_viruses_per_sample': viruses_per_sample.max(),
            'samples_with_viruses': (viruses_per_sample > 0).sum(),
            'viral_richness_distribution': viruses_per_sample.value_counts().to_dict()
        }
        
        return summary
    
    def _save_analysis_results(self, results: Dict) -> None:
        """
        Save analysis results to files.
        
        Args:
            results: Dictionary containing analysis results
        """
        # Save presence matrix
        presence_file = self.output_dir / "viral_presence_matrix.csv"
        results['presence_matrix'].to_csv(presence_file)
        
        # Save co-occurrence matrix
        cooccurrence_file = self.output_dir / "viral_cooccurrence_matrix.csv"
        results['cooccurrence_matrix'].to_csv(cooccurrence_file)
        
        # Save significant associations
        if results['significant_associations']:
            assoc_data = []
            for assoc in results['significant_associations']:
                odds_ratio, p_value = assoc.fisher_exact_test()
                assoc_data.append({
                    'Virus_A': assoc.virus_a,
                    'Virus_B': assoc.virus_b,
                    'Cooccurrence_Count': assoc.cooccurrence_count,
                    'Cooccurrence_Frequency': assoc.cooccurrence_frequency,
                    'Jaccard_Index': assoc.jaccard_index,
                    'Phi_Coefficient': assoc.phi_coefficient,
                    'Odds_Ratio': odds_ratio,
                    'P_Value': p_value
                })
            
            assoc_df = pd.DataFrame(assoc_data)
            assoc_file = self.output_dir / "significant_viral_associations.csv"
            assoc_df.to_csv(assoc_file, index=False)
        
        # Save rare virus analysis
        if results['rare_viruses']:
            rare_data = []
            for rare in results['rare_viruses']:
                rare_data.append({
                    'Virus_Name': rare.virus_name,
                    'Occurrence_Count': rare.occurrence_count,
                    'Prevalence': rare.prevalence,
                    'Isolation_Score': rare.isolation_score,
                    'Is_Rare': rare.is_rare,
                    'Is_Isolated': rare.is_isolated,
                    'Exclusion_Partners_Count': len(rare.exclusion_partners),
                    'Exclusion_Partners': '; '.join(rare.exclusion_partners)
                })
            
            rare_df = pd.DataFrame(rare_data)
            rare_file = self.output_dir / "rare_virus_analysis.csv"
            rare_df.to_csv(rare_file, index=False)
        
        # Save network data
        network_file = self.output_dir / "viral_cooccurrence_network.json"
        with open(network_file, 'w') as f:
            json.dump(results['network'], f, indent=2)
        
        # Save summary statistics
        summary_file = self.output_dir / "analysis_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(results['summary_stats'], f, indent=2)
        
        logger.info(f"Analysis results saved to {self.output_dir}")


def main():
    """
    Example usage of the ViralCooccurrenceAnalyzer.
    """
    # This would typically be called with actual ViralDetectionResult objects
    # from the viral detection pipeline
    
    analyzer = ViralCooccurrenceAnalyzer("viral_cooccurrence_results")
    
    print("Viral Co-occurrence Analyzer initialized")
    print("To use this analyzer:")
    print("1. Run viral detection pipeline on mosquito genomes")
    print("2. Pass ViralDetectionResult objects to run_complete_analysis()")
    print("3. Results will be saved to the output directory")
    
    # Example of what the analysis would show:
    print("\nExample Analysis Output:")
    print("- Viral presence/absence matrix")
    print("- Significant viral co-occurrences")
    print("- Rare and isolated viruses")
    print("- Viral exclusion patterns")
    print("- Co-occurrence network for visualization")


if __name__ == "__main__":
    main()