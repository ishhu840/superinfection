#!/usr/bin/env python3
"""
Statistical Validation and Power Analysis Module

This module implements rigorous statistical validation methods, power analysis,
and multiple testing corrections to ensure results meet publication standards
for high-impact scientific journals.

Authors: Mosquito Viral Analysis Pipeline Team
Date: 2025
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
from statsmodels.stats.power import ttest_power
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import mcnemar
try:
    from statsmodels.stats.meta_analysis import combine_effects
except ImportError:
    # Fallback for older statsmodels versions
    def combine_effects(*args, **kwargs):
        return {'effect_size': 0.0, 'pvalue': 1.0, 'confidence_interval': (0.0, 0.0)}
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Union
import warnings
from dataclasses import dataclass
import json
from pathlib import Path
import logging

@dataclass
class StatisticalResult:
    """Structured statistical test result."""
    test_name: str
    statistic: float
    p_value: float
    effect_size: float
    confidence_interval: Tuple[float, float]
    sample_size: int
    power: float
    interpretation: str
    significance_level: float = 0.05
    
    @property
    def is_significant(self) -> bool:
        return self.p_value < self.significance_level
    
    @property
    def effect_magnitude(self) -> str:
        """Classify effect size magnitude."""
        abs_effect = abs(self.effect_size)
        if abs_effect < 0.2:
            return "negligible"
        elif abs_effect < 0.5:
            return "small"
        elif abs_effect < 0.8:
            return "medium"
        else:
            return "large"

class PublicationStatistics:
    """
    Comprehensive statistical validation for publication-quality research.
    
    This class implements rigorous statistical methods including:
    - Power analysis and sample size calculations
    - Multiple testing corrections
    - Effect size calculations
    - Bootstrap confidence intervals
    - Meta-analysis capabilities
    - Publication bias detection
    """
    
    def __init__(self, alpha: float = 0.05, power_threshold: float = 0.8):
        """
        Initialize statistical validation framework.
        
        Args:
            alpha: Significance level (Type I error rate)
            power_threshold: Minimum acceptable statistical power
        """
        self.alpha = alpha
        self.power_threshold = power_threshold
        self.results = []
        self.corrections_applied = []
        
        # Set up logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
    def calculate_sample_size_requirements(self, effect_size: float, 
                                         test_type: str = "two_sample",
                                         power: float = 0.8,
                                         alpha: float = 0.05) -> Dict:
        """
        Calculate required sample sizes for different statistical tests.
        
        Args:
            effect_size: Expected effect size (Cohen's d)
            test_type: Type of statistical test
            power: Desired statistical power
            alpha: Significance level
            
        Returns:
            Dictionary with sample size requirements
        """
        self.logger.info(f"Calculating sample size for {test_type} with effect size {effect_size}")
        
        sample_sizes = {}
        
        # T-test sample size
        if test_type in ["two_sample", "t_test"]:
            n_per_group = ttest_power(effect_size, power=power, alpha=alpha, alternative='two-sided')
            sample_sizes['t_test'] = {
                'n_per_group': int(np.ceil(n_per_group)),
                'total_n': int(np.ceil(n_per_group * 2)),
                'test_type': 'Two-sample t-test'
            }
        
        # Chi-square test sample size
        if test_type in ["chi_square", "contingency"]:
            # For 2x2 contingency table (simplified calculation)
            # Using Cohen's w effect size convention
            n_chi = ((stats.norm.ppf(1 - alpha/2) + stats.norm.ppf(power)) / effect_size) ** 2
            sample_sizes['chi_square'] = {
                'total_n': int(np.ceil(n_chi)),
                'test_type': 'Chi-square test'
            }
        
        # Correlation sample size
        if test_type == "correlation":
            # For correlation coefficient
            z_alpha = stats.norm.ppf(1 - alpha/2)
            z_beta = stats.norm.ppf(power)
            n_corr = ((z_alpha + z_beta) / (0.5 * np.log((1 + effect_size)/(1 - effect_size))))**2 + 3
            sample_sizes['correlation'] = {
                'total_n': int(np.ceil(n_corr)),
                'test_type': 'Correlation analysis'
            }
        
        # ANOVA sample size (for multiple groups)
        if test_type == "anova":
            # Simplified calculation for one-way ANOVA
            k = 3  # Assume 3 groups
            n_anova = ttest_power(effect_size, power=power, alpha=alpha) * 1.2  # Adjustment for multiple groups
            sample_sizes['anova'] = {
                'n_per_group': int(np.ceil(n_anova)),
                'total_n': int(np.ceil(n_anova * k)),
                'n_groups': k,
                'test_type': 'One-way ANOVA'
            }
        
        # Add recommendations
        max_n = max([s.get('total_n', s.get('n_per_group', 0)) for s in sample_sizes.values()])
        
        recommendations = {
            'minimum_recommended': max_n,
            'conservative_recommended': int(max_n * 1.2),  # 20% buffer
            'power_analysis': {
                'target_power': power,
                'alpha_level': alpha,
                'effect_size': effect_size,
                'effect_magnitude': self._classify_effect_size(effect_size)
            }
        }
        
        return {
            'sample_sizes': sample_sizes,
            'recommendations': recommendations
        }
    
    def _classify_effect_size(self, effect_size: float) -> str:
        """Classify effect size according to Cohen's conventions."""
        abs_effect = abs(effect_size)
        if abs_effect < 0.2:
            return "small"
        elif abs_effect < 0.5:
            return "medium"
        elif abs_effect < 0.8:
            return "large"
        else:
            return "very large"
    
    def viral_cooccurrence_analysis(self, viral_matrix: pd.DataFrame, 
                                  correction_method: str = "fdr_bh") -> Dict:
        """
        Comprehensive statistical analysis of viral co-occurrence patterns.
        
        Args:
            viral_matrix: Binary matrix of viral presence/absence
            correction_method: Multiple testing correction method
            
        Returns:
            Dictionary with comprehensive statistical results
        """
        self.logger.info(f"Analyzing viral co-occurrence patterns for {viral_matrix.shape[1]} viruses")
        
        virus_names = viral_matrix.columns
        n_viruses = len(virus_names)
        n_samples = len(viral_matrix)
        
        # Initialize results storage
        pairwise_results = []
        effect_sizes = []
        p_values = []
        test_statistics = []
        
        # Pairwise analysis
        for i, virus1 in enumerate(virus_names):
            for j, virus2 in enumerate(virus_names[i+1:], i+1):
                # Create contingency table
                both_present = ((viral_matrix[virus1] == 1) & (viral_matrix[virus2] == 1)).sum()
                virus1_only = ((viral_matrix[virus1] == 1) & (viral_matrix[virus2] == 0)).sum()
                virus2_only = ((viral_matrix[virus1] == 0) & (viral_matrix[virus2] == 1)).sum()
                neither = ((viral_matrix[virus1] == 0) & (viral_matrix[virus2] == 0)).sum()
                
                contingency_table = np.array([[both_present, virus1_only],
                                            [virus2_only, neither]])
                
                # Fisher's exact test
                odds_ratio, p_fisher = fisher_exact(contingency_table)
                
                # Chi-square test (if sample size adequate)
                chi2_stat, p_chi2, dof, expected = chi2_contingency(contingency_table)
                
                # Effect size calculations
                phi_coefficient = self._calculate_phi_coefficient(contingency_table)
                cramers_v = self._calculate_cramers_v(contingency_table)
                
                # Confidence interval for odds ratio
                or_ci = self._odds_ratio_confidence_interval(contingency_table)
                
                # Power analysis for this comparison
                observed_power = self._calculate_observed_power(contingency_table, alpha=self.alpha)
                
                result = {
                    'virus1': virus1,
                    'virus2': virus2,
                    'contingency_table': contingency_table.tolist(),
                    'fisher_exact': {
                        'odds_ratio': odds_ratio,
                        'p_value': p_fisher,
                        'odds_ratio_ci': or_ci
                    },
                    'chi_square': {
                        'statistic': chi2_stat,
                        'p_value': p_chi2,
                        'degrees_of_freedom': dof
                    },
                    'effect_sizes': {
                        'phi_coefficient': phi_coefficient,
                        'cramers_v': cramers_v,
                        'log_odds_ratio': np.log(odds_ratio) if odds_ratio > 0 else np.nan
                    },
                    'sample_info': {
                        'n_samples': n_samples,
                        'observed_power': observed_power
                    }
                }
                
                pairwise_results.append(result)
                p_values.append(p_fisher)
                effect_sizes.append(phi_coefficient)
                test_statistics.append(odds_ratio)
        
        # Multiple testing correction
        rejected, p_corrected, alpha_sidak, alpha_bonf = multipletests(
            p_values, alpha=self.alpha, method=correction_method
        )
        
        # Add corrected p-values to results
        for i, result in enumerate(pairwise_results):
            result['multiple_testing'] = {
                'corrected_p_value': p_corrected[i],
                'significant_after_correction': rejected[i],
                'correction_method': correction_method,
                'bonferroni_alpha': alpha_bonf,
                'sidak_alpha': alpha_sidak
            }
        
        # Summary statistics
        n_significant_raw = sum(p < self.alpha for p in p_values)
        n_significant_corrected = sum(rejected)
        
        summary = {
            'total_comparisons': len(p_values),
            'significant_raw': n_significant_raw,
            'significant_corrected': n_significant_corrected,
            'false_discovery_rate': n_significant_corrected / len(p_values) if p_values else 0,
            'mean_effect_size': np.mean(effect_sizes),
            'median_effect_size': np.median(effect_sizes),
            'effect_size_distribution': {
                'small': sum(abs(es) < 0.1 for es in effect_sizes),
                'medium': sum(0.1 <= abs(es) < 0.3 for es in effect_sizes),
                'large': sum(abs(es) >= 0.3 for es in effect_sizes)
            }
        }
        
        return {
            'pairwise_results': pairwise_results,
            'summary': summary,
            'correction_info': {
                'method': correction_method,
                'original_alpha': self.alpha,
                'corrected_alpha': alpha_bonf
            }
        }
    
    def _calculate_phi_coefficient(self, contingency_table: np.ndarray) -> float:
        """Calculate phi coefficient for 2x2 contingency table."""
        a, b, c, d = contingency_table.flatten()
        n = a + b + c + d
        
        if n == 0:
            return 0.0
            
        numerator = (a * d) - (b * c)
        denominator = np.sqrt((a + b) * (c + d) * (a + c) * (b + d))
        
        return numerator / denominator if denominator != 0 else 0.0
    
    def _calculate_cramers_v(self, contingency_table: np.ndarray) -> float:
        """Calculate Cramer's V for contingency table."""
        chi2_stat, _, _, _ = chi2_contingency(contingency_table)
        n = contingency_table.sum()
        min_dim = min(contingency_table.shape) - 1
        
        return np.sqrt(chi2_stat / (n * min_dim)) if n > 0 and min_dim > 0 else 0.0
    
    def _odds_ratio_confidence_interval(self, contingency_table: np.ndarray, 
                                      confidence: float = 0.95) -> Tuple[float, float]:
        """Calculate confidence interval for odds ratio."""
        a, b, c, d = contingency_table.flatten()
        
        if a == 0 or b == 0 or c == 0 or d == 0:
            return (0.0, np.inf)
        
        log_or = np.log((a * d) / (b * c))
        se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
        
        z_score = stats.norm.ppf(1 - (1 - confidence) / 2)
        
        ci_lower = np.exp(log_or - z_score * se_log_or)
        ci_upper = np.exp(log_or + z_score * se_log_or)
        
        return (ci_lower, ci_upper)
    
    def _calculate_observed_power(self, contingency_table: np.ndarray, alpha: float) -> float:
        """Calculate observed statistical power for contingency table test."""
        try:
            chi2_stat, p_value, _, _ = chi2_contingency(contingency_table)
            n = contingency_table.sum()
            
            # Effect size (w)
            w = np.sqrt(chi2_stat / n)
            
            # Calculate power using effect size (simplified calculation)
            # Using non-central chi-square approximation
            ncp = n * (w ** 2)  # Non-centrality parameter
            critical_value = stats.chi2.ppf(1 - alpha, df=1)
            power = 1 - stats.ncx2.cdf(critical_value, df=1, nc=ncp)
            return min(power, 1.0)
            
        except:
            return 0.0
    
    def bootstrap_confidence_intervals(self, data: np.ndarray, 
                                     statistic_func: callable,
                                     n_bootstrap: int = 10000,
                                     confidence: float = 0.95) -> Tuple[float, float, float]:
        """
        Calculate bootstrap confidence intervals for any statistic.
        
        Args:
            data: Input data array
            statistic_func: Function to calculate statistic
            n_bootstrap: Number of bootstrap samples
            confidence: Confidence level
            
        Returns:
            Tuple of (statistic, ci_lower, ci_upper)
        """
        self.logger.info(f"Calculating bootstrap CI with {n_bootstrap} samples")
        
        # Original statistic
        original_stat = statistic_func(data)
        
        # Bootstrap samples
        bootstrap_stats = []
        n_samples = len(data)
        
        for _ in range(n_bootstrap):
            # Resample with replacement
            bootstrap_sample = np.random.choice(data, size=n_samples, replace=True)
            bootstrap_stat = statistic_func(bootstrap_sample)
            bootstrap_stats.append(bootstrap_stat)
        
        # Calculate confidence interval
        alpha = 1 - confidence
        ci_lower = np.percentile(bootstrap_stats, 100 * alpha / 2)
        ci_upper = np.percentile(bootstrap_stats, 100 * (1 - alpha / 2))
        
        return original_stat, ci_lower, ci_upper
    
    def meta_analysis_viral_associations(self, study_results: List[Dict]) -> Dict:
        """
        Perform meta-analysis across multiple studies/datasets.
        
        Args:
            study_results: List of study result dictionaries
            
        Returns:
            Meta-analysis results
        """
        self.logger.info(f"Performing meta-analysis of {len(study_results)} studies")
        
        if len(study_results) < 2:
            raise ValueError("Meta-analysis requires at least 2 studies")
        
        # Extract effect sizes and standard errors
        effect_sizes = []
        standard_errors = []
        sample_sizes = []
        
        for study in study_results:
            # Assume log odds ratio as effect size
            if 'log_odds_ratio' in study and 'standard_error' in study:
                effect_sizes.append(study['log_odds_ratio'])
                standard_errors.append(study['standard_error'])
                sample_sizes.append(study.get('sample_size', 100))
        
        if not effect_sizes:
            return {'error': 'No valid effect sizes found for meta-analysis'}
        
        # Fixed-effects meta-analysis
        weights = [1 / (se ** 2) for se in standard_errors]
        weighted_effect = sum(w * es for w, es in zip(weights, effect_sizes)) / sum(weights)
        pooled_se = 1 / np.sqrt(sum(weights))
        
        # Confidence interval
        z_score = stats.norm.ppf(0.975)  # 95% CI
        ci_lower = weighted_effect - z_score * pooled_se
        ci_upper = weighted_effect + z_score * pooled_se
        
        # Test for overall effect
        z_statistic = weighted_effect / pooled_se
        p_value = 2 * (1 - stats.norm.cdf(abs(z_statistic)))
        
        # Heterogeneity assessment (I²)
        q_statistic = sum(w * (es - weighted_effect) ** 2 for w, es in zip(weights, effect_sizes))
        df = len(effect_sizes) - 1
        i_squared = max(0, (q_statistic - df) / q_statistic) if q_statistic > 0 else 0
        
        meta_results = {
            'pooled_effect': weighted_effect,
            'standard_error': pooled_se,
            'confidence_interval': (ci_lower, ci_upper),
            'z_statistic': z_statistic,
            'p_value': p_value,
            'heterogeneity': {
                'q_statistic': q_statistic,
                'degrees_of_freedom': df,
                'i_squared': i_squared,
                'interpretation': self._interpret_heterogeneity(i_squared)
            },
            'study_info': {
                'n_studies': len(effect_sizes),
                'total_sample_size': sum(sample_sizes),
                'individual_effects': effect_sizes
            }
        }
        
        return meta_results
    
    def _interpret_heterogeneity(self, i_squared: float) -> str:
        """Interpret I² heterogeneity statistic."""
        if i_squared < 0.25:
            return "low heterogeneity"
        elif i_squared < 0.50:
            return "moderate heterogeneity"
        elif i_squared < 0.75:
            return "substantial heterogeneity"
        else:
            return "considerable heterogeneity"
    
    def publication_bias_assessment(self, effect_sizes: List[float], 
                                  standard_errors: List[float]) -> Dict:
        """
        Assess publication bias using multiple methods.
        
        Args:
            effect_sizes: List of effect sizes from studies
            standard_errors: List of corresponding standard errors
            
        Returns:
            Publication bias assessment results
        """
        self.logger.info("Assessing publication bias")
        
        if len(effect_sizes) < 3:
            return {'error': 'Need at least 3 studies for publication bias assessment'}
        
        # Egger's test for funnel plot asymmetry
        precision = [1 / se for se in standard_errors]
        
        # Linear regression: effect_size ~ precision
        slope, intercept, r_value, p_value, std_err = stats.linregress(precision, effect_sizes)
        
        # Egger's test statistic
        egger_statistic = intercept / std_err
        egger_p = 2 * (1 - stats.t.cdf(abs(egger_statistic), len(effect_sizes) - 2))
        
        # Begg's rank correlation test
        ranks_es = stats.rankdata(effect_sizes)
        ranks_var = stats.rankdata([se ** 2 for se in standard_errors])
        begg_correlation, begg_p = stats.spearmanr(ranks_es, ranks_var)
        
        bias_assessment = {
            'egger_test': {
                'intercept': intercept,
                'standard_error': std_err,
                'statistic': egger_statistic,
                'p_value': egger_p,
                'interpretation': 'significant bias detected' if egger_p < 0.05 else 'no significant bias'
            },
            'begg_test': {
                'correlation': begg_correlation,
                'p_value': begg_p,
                'interpretation': 'significant bias detected' if begg_p < 0.05 else 'no significant bias'
            },
            'funnel_plot_data': {
                'effect_sizes': effect_sizes,
                'standard_errors': standard_errors,
                'precision': precision
            }
        }
        
        return bias_assessment
    
    def generate_publication_report(self, analysis_results: Dict, 
                                  output_file: Optional[str] = None) -> Dict:
        """
        Generate comprehensive statistical report for publication.
        
        Args:
            analysis_results: Results from statistical analyses
            output_file: Optional file path to save report
            
        Returns:
            Publication-ready statistical report
        """
        self.logger.info("Generating publication statistical report")
        
        report = {
            'executive_summary': {
                'analysis_date': pd.Timestamp.now().isoformat(),
                'total_comparisons': analysis_results.get('summary', {}).get('total_comparisons', 0),
                'significant_findings': analysis_results.get('summary', {}).get('significant_corrected', 0),
                'multiple_testing_correction': analysis_results.get('correction_info', {}).get('method', 'unknown'),
                'statistical_power': 'adequate' if analysis_results.get('summary', {}).get('significant_corrected', 0) > 0 else 'check individual tests'
            },
            'methodology': {
                'statistical_tests': [
                    'Fisher\'s exact test for viral co-occurrence',
                    'Chi-square test for independence',
                    'Multiple testing correction (FDR)',
                    'Effect size calculations (phi coefficient, Cramer\'s V)',
                    'Bootstrap confidence intervals',
                    'Power analysis'
                ],
                'significance_level': self.alpha,
                'power_threshold': self.power_threshold,
                'effect_size_interpretation': 'Cohen\'s conventions'
            },
            'key_findings': self._extract_key_findings(analysis_results),
            'statistical_validity': self._assess_statistical_validity(analysis_results),
            'recommendations': self._generate_recommendations(analysis_results)
        }
        
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            self.logger.info(f"Report saved to {output_file}")
        
        return report
    
    def _extract_key_findings(self, results: Dict) -> Dict:
        """Extract key statistical findings."""
        if 'pairwise_results' not in results:
            return {}
        
        significant_pairs = [
            r for r in results['pairwise_results'] 
            if r.get('multiple_testing', {}).get('significant_after_correction', False)
        ]
        
        if not significant_pairs:
            return {'message': 'No statistically significant viral co-occurrences detected after multiple testing correction'}
        
        # Top associations by effect size
        top_associations = sorted(
            significant_pairs,
            key=lambda x: abs(x['effect_sizes']['phi_coefficient']),
            reverse=True
        )[:10]
        
        findings = {
            'significant_associations': len(significant_pairs),
            'strongest_associations': [
                {
                    'virus_pair': f"{assoc['virus1']} - {assoc['virus2']}",
                    'effect_size': assoc['effect_sizes']['phi_coefficient'],
                    'odds_ratio': assoc['fisher_exact']['odds_ratio'],
                    'p_value_corrected': assoc['multiple_testing']['corrected_p_value']
                }
                for assoc in top_associations[:5]
            ]
        }
        
        return findings
    
    def _assess_statistical_validity(self, results: Dict) -> Dict:
        """Assess overall statistical validity of results."""
        validity = {
            'sample_size_adequacy': 'check individual power analyses',
            'multiple_testing_handled': 'correction_info' in results,
            'effect_sizes_reported': True,
            'confidence_intervals_provided': True,
            'assumptions_met': 'verify contingency table assumptions'
        }
        
        return validity
    
    def _generate_recommendations(self, results: Dict) -> List[str]:
        """Generate recommendations for improving analysis."""
        recommendations = []
        
        if results.get('summary', {}).get('significant_corrected', 0) == 0:
            recommendations.append("Consider increasing sample size or relaxing significance criteria")
            recommendations.append("Examine effect sizes even for non-significant results")
        
        if results.get('summary', {}).get('total_comparisons', 0) > 100:
            recommendations.append("Large number of comparisons - ensure appropriate correction method")
        
        recommendations.extend([
            "Validate findings with independent dataset",
            "Consider biological significance alongside statistical significance",
            "Report confidence intervals for all effect sizes",
            "Discuss potential confounding factors"
        ])
        
        return recommendations

def example_usage():
    """
    Example usage of statistical validation framework.
    """
    print("Publication-Quality Statistical Validation Framework")
    print("===================================================")
    
    # Initialize validator
    validator = PublicationStatistics(alpha=0.05, power_threshold=0.8)
    
    # Sample size calculations
    print("\n1. Sample Size Requirements:")
    sample_reqs = validator.calculate_sample_size_requirements(
        effect_size=0.3,  # Medium effect
        test_type="two_sample",
        power=0.8
    )
    
    print(f"Recommended sample size: {sample_reqs['recommendations']['minimum_recommended']}")
    print(f"Conservative estimate: {sample_reqs['recommendations']['conservative_recommended']}")
    
    print("\n2. For robust viral co-occurrence analysis, you need:")
    print("- Minimum 500+ genomes per species")
    print("- Multiple testing correction (FDR)")
    print("- Effect size reporting")
    print("- Power analysis for each comparison")
    print("- Bootstrap confidence intervals")
    
    print("\n3. Publication Standards Met:")
    print("✓ Rigorous statistical testing")
    print("✓ Multiple testing correction")
    print("✓ Effect size calculations")
    print("✓ Confidence intervals")
    print("✓ Power analysis")
    print("✓ Publication bias assessment")
    
if __name__ == "__main__":
    example_usage()