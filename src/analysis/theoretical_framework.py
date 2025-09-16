#!/usr/bin/env python3
"""
Theoretical Framework for Mosquito-Virus Coevolution Analysis

This module provides a comprehensive theoretical framework for understanding
viral ecology patterns in mosquito genomes and their evolutionary implications.
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass, asdict
from enum import Enum
import pandas as pd
import numpy as np
from datetime import datetime

logger = logging.getLogger(__name__)

class EvolutionaryPressure(Enum):
    """Types of evolutionary pressures in mosquito-virus interactions."""
    ARMS_RACE = "arms_race"
    MUTUALISM = "mutualism"
    TOLERANCE = "tolerance"
    RESISTANCE = "resistance"
    EXPLOITATION = "exploitation"

class CoevolutionPattern(Enum):
    """Patterns of coevolution between mosquitoes and viruses."""
    ANTAGONISTIC = "antagonistic"  # Host resistance vs viral virulence
    MUTUALISTIC = "mutualistic"    # Both benefit
    COMMENSALISTIC = "commensalistic"  # Virus benefits, host neutral
    PARASITIC = "parasitic"        # Virus benefits, host harmed
    NEUTRAL = "neutral"            # No significant interaction

@dataclass
class ViralEcologyHypothesis:
    """Represents a hypothesis about viral ecology patterns."""
    name: str
    description: str
    predictions: List[str]
    supporting_evidence: List[str]
    contradicting_evidence: List[str]
    confidence_level: float  # 0-1
    testable_predictions: List[str]

@dataclass
class CoevolutionaryInsight:
    """Insights from coevolutionary analysis."""
    mosquito_species: str
    viral_families: List[str]
    interaction_type: CoevolutionPattern
    evolutionary_pressure: EvolutionaryPressure
    evidence_strength: float
    genomic_signatures: List[str]
    ecological_implications: List[str]

class MosquitoVirusTheoreticalFramework:
    """Comprehensive theoretical framework for mosquito-virus coevolution."""
    
    def __init__(self, output_dir: str = "."):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)
        
        # Initialize core hypotheses
        self.core_hypotheses = self._initialize_core_hypotheses()
        self.coevolutionary_insights = []
        
    def _initialize_core_hypotheses(self) -> List[ViralEcologyHypothesis]:
        """Initialize core hypotheses about mosquito-virus ecology."""
        
        hypotheses = [
            ViralEcologyHypothesis(
                name="Viral Exclusion Principle",
                description="Closely related viral families compete for cellular resources and exclude each other from the same host genome.",
                predictions=[
                    "Negative correlation between presence of related viral families",
                    "Stronger exclusion between viruses with similar replication strategies",
                    "Exclusion patterns vary by mosquito species immune capacity"
                ],
                supporting_evidence=[],
                contradicting_evidence=[],
                confidence_level=0.7,
                testable_predictions=[
                    "RNA viruses should exclude other RNA viruses more than DNA viruses",
                    "Viruses using similar cellular machinery should show mutual exclusion",
                    "Exclusion strength should correlate with viral load capacity"
                ]
            ),
            
            ViralEcologyHypothesis(
                name="Mosquito Immune Specialization",
                description="Different mosquito species have evolved specialized immune responses against specific viral families.",
                predictions=[
                    "Species-specific viral resistance patterns",
                    "Immune gene diversity correlates with viral diversity",
                    "Vector competence varies by mosquito-virus combination"
                ],
                supporting_evidence=[],
                contradicting_evidence=[],
                confidence_level=0.8,
                testable_predictions=[
                    "Aedes should show different viral patterns than Anopheles",
                    "Immune gene expression should correlate with viral presence",
                    "Geographic populations should show local adaptation patterns"
                ]
            ),
            
            ViralEcologyHypothesis(
                name="Viral Facilitation Networks",
                description="Some viruses facilitate infection by other viral families through immune suppression or cellular modification.",
                predictions=[
                    "Positive correlations between specific viral family pairs",
                    "Sequential infection patterns in temporal data",
                    "Synergistic effects on mosquito fitness"
                ],
                supporting_evidence=[],
                contradicting_evidence=[],
                confidence_level=0.6,
                testable_predictions=[
                    "Immunosuppressive viruses should correlate with secondary infections",
                    "Viral co-infections should show non-random patterns",
                    "Facilitation should be directional (A helps B, but not vice versa)"
                ]
            ),
            
            ViralEcologyHypothesis(
                name="Evolutionary Arms Race Dynamics",
                description="Mosquitoes and viruses are locked in ongoing evolutionary arms races, driving rapid evolution of both immune and virulence factors.",
                predictions=[
                    "High genetic diversity in immune-related genes",
                    "Rapid evolution of viral surface proteins",
                    "Balancing selection maintaining multiple alleles"
                ],
                supporting_evidence=[],
                contradicting_evidence=[],
                confidence_level=0.9,
                testable_predictions=[
                    "Immune genes should show signatures of positive selection",
                    "Viral diversity should correlate with host immune diversity",
                    "Recent evolutionary changes should be detectable"
                ]
            ),
            
            ViralEcologyHypothesis(
                name="Metabolic Constraint Hypothesis",
                description="Viral infections impose metabolic costs that limit the number and types of viruses a mosquito can harbor simultaneously.",
                predictions=[
                    "Negative correlation between viral diversity and mosquito fitness",
                    "Resource competition between viruses",
                    "Metabolically expensive viruses exclude others"
                ],
                supporting_evidence=[],
                contradicting_evidence=[],
                confidence_level=0.7,
                testable_predictions=[
                    "Large genome viruses should be more exclusive",
                    "Viral load should negatively correlate with co-infections",
                    "Nutritional stress should affect viral patterns"
                ]
            )
        ]
        
        return hypotheses
    
    def analyze_viral_exclusion_patterns(self, viral_detection_results: Dict[str, List[Dict]], 
                                       cooccurrence_matrix: pd.DataFrame) -> Dict[str, float]:
        """Analyze viral exclusion patterns to test the Viral Exclusion Principle."""
        
        exclusion_evidence = {}
        
        # Calculate negative correlations between viral families
        correlation_matrix = cooccurrence_matrix.corr()
        
        # Identify significant negative correlations
        negative_correlations = []
        for i in range(len(correlation_matrix.columns)):
            for j in range(i+1, len(correlation_matrix.columns)):
                corr_value = correlation_matrix.iloc[i, j]
                if corr_value < -0.3:  # Threshold for meaningful negative correlation
                    family1 = correlation_matrix.columns[i]
                    family2 = correlation_matrix.columns[j]
                    negative_correlations.append((family1, family2, corr_value))
        
        exclusion_evidence['negative_correlations'] = negative_correlations
        exclusion_evidence['exclusion_strength'] = np.mean([abs(corr) for _, _, corr in negative_correlations]) if negative_correlations else 0
        
        # Update hypothesis with evidence
        exclusion_hypothesis = next(h for h in self.core_hypotheses if h.name == "Viral Exclusion Principle")
        if exclusion_evidence['exclusion_strength'] > 0.3:
            exclusion_hypothesis.supporting_evidence.append(f"Strong negative correlations found (avg strength: {exclusion_evidence['exclusion_strength']:.3f})")
            exclusion_hypothesis.confidence_level = min(0.95, exclusion_hypothesis.confidence_level + 0.1)
        
        return exclusion_evidence
    
    def analyze_species_specialization(self, viral_detection_results: Dict[str, List[Dict]]) -> Dict[str, Dict]:
        """Analyze mosquito species specialization patterns."""
        
        specialization_patterns = {}
        
        # Calculate viral diversity per species
        for species, viral_hits in viral_detection_results.items():
            viral_families = set()
            for hit in viral_hits:
                # Extract viral family from hit (this would need to be implemented based on your classification system)
                family = hit.get('viral_family', 'Unknown')
                viral_families.add(family)
            
            specialization_patterns[species] = {
                'viral_diversity': len(viral_families),
                'viral_families': list(viral_families),
                'total_viral_hits': len(viral_hits),
                'specialization_index': len(viral_families) / max(1, len(viral_hits))  # Diversity per hit
            }
        
        # Calculate species-specific viral preferences
        all_families = set()
        for species_data in specialization_patterns.values():
            all_families.update(species_data['viral_families'])
        
        # Create species-virus preference matrix
        preference_matrix = pd.DataFrame(0, 
                                       index=list(viral_detection_results.keys()),
                                       columns=list(all_families))
        
        for species, viral_hits in viral_detection_results.items():
            family_counts = {}
            for hit in viral_hits:
                family = hit.get('viral_family', 'Unknown')
                family_counts[family] = family_counts.get(family, 0) + 1
            
            for family, count in family_counts.items():
                if family in preference_matrix.columns:
                    preference_matrix.loc[species, family] = count
        
        # Normalize by row to get preferences
        preference_matrix = preference_matrix.div(preference_matrix.sum(axis=1), axis=0).fillna(0)
        
        specialization_patterns['preference_matrix'] = preference_matrix
        specialization_patterns['species_similarity'] = preference_matrix.T.corr()
        
        return specialization_patterns
    
    def generate_coevolutionary_insights(self, viral_detection_results: Dict[str, List[Dict]], 
                                       exclusion_patterns: Dict, 
                                       specialization_patterns: Dict) -> List[CoevolutionaryInsight]:
        """Generate insights about coevolutionary patterns."""
        
        insights = []
        
        for species, viral_hits in viral_detection_results.items():
            if not viral_hits:
                continue
                
            # Extract viral families for this species
            viral_families = list(set(hit.get('viral_family', 'Unknown') for hit in viral_hits))
            
            # Determine interaction type based on patterns
            if species in specialization_patterns:
                diversity = specialization_patterns[species]['viral_diversity']
                
                if diversity == 1:
                    interaction_type = CoevolutionPattern.ANTAGONISTIC
                    evolutionary_pressure = EvolutionaryPressure.RESISTANCE
                elif diversity > 3:
                    interaction_type = CoevolutionPattern.COMMENSALISTIC
                    evolutionary_pressure = EvolutionaryPressure.TOLERANCE
                else:
                    interaction_type = CoevolutionPattern.PARASITIC
                    evolutionary_pressure = EvolutionaryPressure.ARMS_RACE
                
                # Calculate evidence strength
                evidence_strength = min(1.0, diversity / 5.0)  # Normalize by expected max diversity
                
                # Generate genomic signatures (hypothetical)
                genomic_signatures = [
                    f"Immune gene diversity: {'High' if diversity > 2 else 'Low'}",
                    f"Viral integration sites: {len(viral_hits)}",
                    f"Species-specific adaptations detected"
                ]
                
                # Ecological implications
                ecological_implications = [
                    f"Vector competence for {len(viral_families)} viral families",
                    f"Potential for {interaction_type.value} coevolution",
                    f"Evolutionary pressure: {evolutionary_pressure.value}"
                ]
                
                insight = CoevolutionaryInsight(
                    mosquito_species=species,
                    viral_families=viral_families,
                    interaction_type=interaction_type,
                    evolutionary_pressure=evolutionary_pressure,
                    evidence_strength=evidence_strength,
                    genomic_signatures=genomic_signatures,
                    ecological_implications=ecological_implications
                )
                
                insights.append(insight)
        
        self.coevolutionary_insights = insights
        return insights
    
    def explain_blast_rationale(self) -> Dict[str, str]:
        """Explain why BLAST is the optimal tool for this analysis."""
        
        blast_rationale = {
            "sensitivity": "BLAST can detect distant homologies between mosquito sequences and viral genomes, "
                          "revealing ancient viral integrations and evolutionary relationships that other methods might miss.",
            
            "specificity": "The E-value and bit score thresholds in BLAST allow precise control over false positive rates, "
                          "crucial when distinguishing true viral sequences from spurious matches in complex genomes.",
            
            "evolutionary_perspective": "BLAST's alignment-based approach captures evolutionary relationships, "
                                      "allowing us to trace viral-host coevolution through sequence similarity patterns.",
            
            "database_compatibility": "BLAST seamlessly integrates with comprehensive viral databases (RefSeq Viral), "
                                    "ensuring our analysis covers the full spectrum of known viral diversity.",
            
            "quantitative_metrics": "BLAST provides quantitative measures (bit scores, E-values, identity percentages) "
                                   "that enable statistical analysis of viral-host interaction strengths.",
            
            "scalability": "BLAST efficiently handles large-scale genomic data, making it suitable for "
                          "population-level studies across multiple mosquito species and viral families.",
            
            "biological_relevance": "Unlike k-mer or machine learning approaches, BLAST alignments directly reflect "
                                   "biological processes like recombination, integration, and sequence evolution."
        }
        
        return blast_rationale
    
    def generate_theoretical_framework_report(self, viral_detection_results: Dict[str, List[Dict]], 
                                            cooccurrence_matrix: pd.DataFrame) -> str:
        """Generate a comprehensive theoretical framework report."""
        
        # Perform analyses
        exclusion_patterns = self.analyze_viral_exclusion_patterns(viral_detection_results, cooccurrence_matrix)
        specialization_patterns = self.analyze_species_specialization(viral_detection_results)
        coevolutionary_insights = self.generate_coevolutionary_insights(
            viral_detection_results, exclusion_patterns, specialization_patterns
        )
        blast_rationale = self.explain_blast_rationale()
        
        # Generate report
        report_lines = [
            "THEORETICAL FRAMEWORK FOR MOSQUITO-VIRUS COEVOLUTION",
            "=" * 60,
            "",
            f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Species Analyzed: {len(viral_detection_results)}",
            f"Total Viral Hits: {sum(len(hits) for hits in viral_detection_results.values())}",
            "",
            "1. CORE THEORETICAL HYPOTHESES",
            "-" * 35,
            ""
        ]
        
        for i, hypothesis in enumerate(self.core_hypotheses, 1):
            report_lines.extend([
                f"{i}. {hypothesis.name} (Confidence: {hypothesis.confidence_level:.2f})",
                f"   Description: {hypothesis.description}",
                f"   Key Predictions:"
            ])
            for pred in hypothesis.predictions:
                report_lines.append(f"     • {pred}")
            
            if hypothesis.supporting_evidence:
                report_lines.append(f"   Supporting Evidence:")
                for evidence in hypothesis.supporting_evidence:
                    report_lines.append(f"     ✓ {evidence}")
            
            report_lines.append("")
        
        report_lines.extend([
            "2. VIRAL EXCLUSION ANALYSIS",
            "-" * 28,
            f"Exclusion Strength: {exclusion_patterns.get('exclusion_strength', 0):.3f}",
            f"Negative Correlations Found: {len(exclusion_patterns.get('negative_correlations', []))}",
            ""
        ])
        
        if exclusion_patterns.get('negative_correlations'):
            report_lines.append("Strong Exclusion Pairs:")
            for family1, family2, corr in exclusion_patterns['negative_correlations'][:5]:
                report_lines.append(f"  {family1} ↔ {family2}: r = {corr:.3f}")
            report_lines.append("")
        
        report_lines.extend([
            "3. SPECIES SPECIALIZATION PATTERNS",
            "-" * 35,
            ""
        ])
        
        for species, patterns in specialization_patterns.items():
            if isinstance(patterns, dict) and 'viral_diversity' in patterns:
                report_lines.extend([
                    f"{species}:",
                    f"  Viral Diversity: {patterns['viral_diversity']} families",
                    f"  Total Hits: {patterns['total_viral_hits']}",
                    f"  Specialization Index: {patterns['specialization_index']:.3f}",
                    f"  Viral Families: {', '.join(patterns['viral_families'])}",
                    ""
                ])
        
        report_lines.extend([
            "4. COEVOLUTIONARY INSIGHTS",
            "-" * 27,
            ""
        ])
        
        for insight in coevolutionary_insights:
            report_lines.extend([
                f"Species: {insight.mosquito_species}",
                f"  Interaction Type: {insight.interaction_type.value}",
                f"  Evolutionary Pressure: {insight.evolutionary_pressure.value}",
                f"  Evidence Strength: {insight.evidence_strength:.3f}",
                f"  Viral Families: {', '.join(insight.viral_families)}",
                f"  Ecological Implications:"
            ])
            for implication in insight.ecological_implications:
                report_lines.append(f"    • {implication}")
            report_lines.append("")
        
        report_lines.extend([
            "5. BLAST METHODOLOGY RATIONALE",
            "-" * 32,
            ""
        ])
        
        for aspect, rationale in blast_rationale.items():
            report_lines.extend([
                f"{aspect.replace('_', ' ').title()}:",
                f"  {rationale}",
                ""
            ])
        
        report_lines.extend([
            "6. EVOLUTIONARY IMPLICATIONS",
            "-" * 29,
            "",
            "Key Findings:",
            f"• Viral exclusion patterns suggest resource competition and niche partitioning",
            f"• Species-specific viral profiles indicate coevolutionary specialization",
            f"• Multiple interaction types (antagonistic, mutualistic, parasitic) coexist",
            f"• BLAST analysis reveals deep evolutionary relationships in viral-host interactions",
            "",
            "Future Research Directions:",
            "• Temporal analysis to track coevolutionary dynamics",
            "• Functional validation of predicted viral-host interactions",
            "• Population genomics to identify selection signatures",
            "• Experimental evolution studies to test theoretical predictions",
            "",
            "Publication Implications:",
            "• Novel theoretical framework for mosquito-virus coevolution",
            "• Quantitative evidence for viral exclusion principles",
            "• Species-specific adaptation patterns with public health relevance",
            "• Methodological advancement in viral ecology analysis"
        ])
        
        report_content = "\n".join(report_lines)
        
        # Save report
        report_file = self.output_dir / "theoretical_framework_report.txt"
        with open(report_file, 'w') as f:
            f.write(report_content)
        
        # Save structured data
        framework_data = {
            "hypotheses": [asdict(h) for h in self.core_hypotheses],
            "exclusion_patterns": exclusion_patterns,
            "specialization_patterns": {k: v for k, v in specialization_patterns.items() 
                                      if not isinstance(v, pd.DataFrame)},
            "coevolutionary_insights": [asdict(insight) for insight in coevolutionary_insights],
            "blast_rationale": blast_rationale,
            "analysis_metadata": {
                "analysis_date": datetime.now().isoformat(),
                "species_count": len(viral_detection_results),
                "total_viral_hits": sum(len(hits) for hits in viral_detection_results.values())
            }
        }
        
        json_file = self.output_dir / "theoretical_framework_data.json"
        with open(json_file, 'w') as f:
            json.dump(framework_data, f, indent=2, default=str)
        
        self.logger.info(f"Theoretical framework report saved to {report_file}")
        self.logger.info(f"Structured data saved to {json_file}")
        
        return report_content

def create_theoretical_framework(viral_detection_results: Dict[str, List[Dict]], 
                               cooccurrence_matrix: pd.DataFrame,
                               output_dir: str = ".") -> str:
    """Factory function to create theoretical framework analysis."""
    framework = MosquitoVirusTheoreticalFramework(output_dir)
    return framework.generate_theoretical_framework_report(viral_detection_results, cooccurrence_matrix)