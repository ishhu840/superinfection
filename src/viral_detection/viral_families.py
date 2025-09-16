#!/usr/bin/env python3
"""
Comprehensive viral family classification system for mosquito-borne viruses.
This module provides detailed taxonomic information and detection parameters
for major viral families found in mosquito vectors.
"""

import logging
from typing import Dict, List, Set, Optional
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)

class ViralGenomeType(Enum):
    """Viral genome types for classification."""
    DNA_SINGLE_STRAND = "ssDNA"
    DNA_DOUBLE_STRAND = "dsDNA"
    RNA_SINGLE_STRAND_POSITIVE = "ssRNA(+)"
    RNA_SINGLE_STRAND_NEGATIVE = "ssRNA(-)"
    RNA_DOUBLE_STRAND = "dsRNA"
    RNA_SINGLE_STRAND_AMBISENSE = "ssRNA(±)"

class TransmissionMode(Enum):
    """Modes of viral transmission."""
    VECTOR_BORNE = "vector-borne"
    VERTICAL = "vertical"
    HORIZONTAL = "horizontal"
    ENVIRONMENTAL = "environmental"
    SEXUAL = "sexual"

@dataclass
class ViralFamily:
    """Comprehensive viral family information."""
    name: str
    genome_type: ViralGenomeType
    genome_size_range: tuple  # (min, max) in nucleotides
    envelope: bool
    transmission_modes: List[TransmissionMode]
    mosquito_hosts: List[str]
    vertebrate_hosts: List[str]
    geographic_distribution: List[str]
    pathogenicity: str  # "high", "moderate", "low", "unknown"
    detection_markers: List[str]  # Key genes/proteins for detection
    ncbi_taxonomy_id: Optional[int]
    representative_species: List[str]
    exclusion_patterns: List[str]  # Viruses that rarely co-occur
    
# Comprehensive viral family database
VIRAL_FAMILIES = {
    "Flaviviridae": ViralFamily(
        name="Flaviviridae",
        genome_type=ViralGenomeType.RNA_SINGLE_STRAND_POSITIVE,
        genome_size_range=(9000, 12000),
        envelope=True,
        transmission_modes=[TransmissionMode.VECTOR_BORNE, TransmissionMode.VERTICAL],
        mosquito_hosts=["Aedes aegypti", "Aedes albopictus", "Culex quinquefasciatus"],
        vertebrate_hosts=["Humans", "Primates", "Birds", "Mammals"],
        geographic_distribution=["Tropical", "Subtropical", "Temperate"],
        pathogenicity="high",
        detection_markers=["NS1", "NS3", "NS5", "E protein", "C protein"],
        ncbi_taxonomy_id=11050,
        representative_species=[
            "Dengue virus", "Zika virus", "Yellow fever virus", 
            "West Nile virus", "Japanese encephalitis virus", "Chikungunya virus"
        ],
        exclusion_patterns=["Bunyaviridae", "Reoviridae"]
    ),
    
    "Togaviridae": ViralFamily(
        name="Togaviridae",
        genome_type=ViralGenomeType.RNA_SINGLE_STRAND_POSITIVE,
        genome_size_range=(11000, 12000),
        envelope=True,
        transmission_modes=[TransmissionMode.VECTOR_BORNE],
        mosquito_hosts=["Aedes aegypti", "Aedes albopictus"],
        vertebrate_hosts=["Humans", "Mammals", "Birds"],
        geographic_distribution=["Tropical", "Subtropical"],
        pathogenicity="high",
        detection_markers=["nsP1", "nsP2", "nsP3", "nsP4", "E1", "E2"],
        ncbi_taxonomy_id=11018,
        representative_species=["Chikungunya virus", "Eastern equine encephalitis virus"],
        exclusion_patterns=["Flaviviridae"]
    ),
    
    "Bunyaviridae": ViralFamily(
        name="Bunyaviridae",
        genome_type=ViralGenomeType.RNA_SINGLE_STRAND_NEGATIVE,
        genome_size_range=(11000, 19000),
        envelope=True,
        transmission_modes=[TransmissionMode.VECTOR_BORNE, TransmissionMode.VERTICAL],
        mosquito_hosts=["Aedes aegypti", "Culex quinquefasciatus", "Anopheles gambiae"],
        vertebrate_hosts=["Humans", "Mammals", "Birds"],
        geographic_distribution=["Global"],
        pathogenicity="moderate",
        detection_markers=["N protein", "L protein", "Gc", "Gn"],
        ncbi_taxonomy_id=11266,
        representative_species=[
            "Rift Valley fever virus", "La Crosse virus", "Crimean-Congo hemorrhagic fever virus"
        ],
        exclusion_patterns=["Flaviviridae"]
    ),
    
    "Reoviridae": ViralFamily(
        name="Reoviridae",
        genome_type=ViralGenomeType.RNA_DOUBLE_STRAND,
        genome_size_range=(18000, 27000),
        envelope=False,
        transmission_modes=[TransmissionMode.VECTOR_BORNE, TransmissionMode.ENVIRONMENTAL],
        mosquito_hosts=["Culex quinquefasciatus", "Aedes aegypti"],
        vertebrate_hosts=["Humans", "Mammals", "Birds"],
        geographic_distribution=["Global"],
        pathogenicity="low",
        detection_markers=["VP1", "VP2", "VP3", "VP4", "VP6", "VP7"],
        ncbi_taxonomy_id=10880,
        representative_species=["Bluetongue virus", "Colorado tick fever virus"],
        exclusion_patterns=["Flaviviridae"]
    ),
    
    "Rhabdoviridae": ViralFamily(
        name="Rhabdoviridae",
        genome_type=ViralGenomeType.RNA_SINGLE_STRAND_NEGATIVE,
        genome_size_range=(11000, 15000),
        envelope=True,
        transmission_modes=[TransmissionMode.VECTOR_BORNE, TransmissionMode.HORIZONTAL],
        mosquito_hosts=["Culex quinquefasciatus", "Aedes aegypti"],
        vertebrate_hosts=["Mammals", "Birds", "Fish"],
        geographic_distribution=["Global"],
        pathogenicity="moderate",
        detection_markers=["N protein", "P protein", "M protein", "G protein", "L protein"],
        ncbi_taxonomy_id=11270,
        representative_species=["Vesicular stomatitis virus", "Rabies virus"],
        exclusion_patterns=[]
    ),
    
    "Parvoviridae": ViralFamily(
        name="Parvoviridae",
        genome_type=ViralGenomeType.DNA_SINGLE_STRAND,
        genome_size_range=(4000, 6000),
        envelope=False,
        transmission_modes=[TransmissionMode.VERTICAL, TransmissionMode.HORIZONTAL],
        mosquito_hosts=["Aedes aegypti", "Anopheles gambiae"],
        vertebrate_hosts=["Mammals", "Birds"],
        geographic_distribution=["Global"],
        pathogenicity="low",
        detection_markers=["NS1", "VP1", "VP2"],
        ncbi_taxonomy_id=10780,
        representative_species=["Aedes aegypti densovirus", "Anopheles gambiae densovirus"],
        exclusion_patterns=[]
    ),
    
    "Iridoviridae": ViralFamily(
        name="Iridoviridae",
        genome_type=ViralGenomeType.DNA_DOUBLE_STRAND,
        genome_size_range=(140000, 220000),
        envelope=True,
        transmission_modes=[TransmissionMode.HORIZONTAL, TransmissionMode.ENVIRONMENTAL],
        mosquito_hosts=["Aedes aegypti", "Culex quinquefasciatus"],
        vertebrate_hosts=["Fish", "Amphibians", "Reptiles"],
        geographic_distribution=["Global"],
        pathogenicity="low",
        detection_markers=["MCP", "DNA polymerase", "ATPase"],
        ncbi_taxonomy_id=10487,
        representative_species=["Invertebrate iridescent virus"],
        exclusion_patterns=[]
    )
}

class ViralFamilyClassifier:
    """Classifier for viral families in mosquito genomes."""
    
    def __init__(self):
        self.families = VIRAL_FAMILIES
        self.logger = logging.getLogger(__name__)
        
    def get_family_info(self, family_name: str) -> Optional[ViralFamily]:
        """Get detailed information about a viral family."""
        return self.families.get(family_name)
    
    def get_families_by_host(self, mosquito_species: str) -> List[ViralFamily]:
        """Get viral families that infect a specific mosquito species."""
        matching_families = []
        for family in self.families.values():
            if mosquito_species in family.mosquito_hosts:
                matching_families.append(family)
        return matching_families
    
    def get_families_by_pathogenicity(self, level: str) -> List[ViralFamily]:
        """Get viral families by pathogenicity level."""
        return [family for family in self.families.values() 
                if family.pathogenicity == level]
    
    def get_exclusion_patterns(self) -> Dict[str, List[str]]:
        """Get viral exclusion patterns for co-occurrence analysis."""
        exclusions = {}
        for name, family in self.families.items():
            if family.exclusion_patterns:
                exclusions[name] = family.exclusion_patterns
        return exclusions
    
    def classify_viral_hit(self, blast_hit: Dict) -> Optional[str]:
        """Classify a BLAST hit to a viral family based on sequence similarity."""
        # This would be implemented with more sophisticated logic
        # For now, simple keyword matching
        description = blast_hit.get('description', '').lower()
        
        for family_name, family in self.families.items():
            for species in family.representative_species:
                if any(word in description for word in species.lower().split()):
                    return family_name
        
        return None
    
    def get_detection_strategy(self, family_name: str) -> Dict[str, any]:
        """Get optimal detection strategy for a viral family."""
        family = self.families.get(family_name)
        if not family:
            return {}
        
        strategy = {
            "blast_parameters": {
                "evalue": 1e-5 if family.pathogenicity == "high" else 1e-3,
                "word_size": 11 if family.genome_type in [ViralGenomeType.RNA_SINGLE_STRAND_POSITIVE, 
                                                          ViralGenomeType.RNA_SINGLE_STRAND_NEGATIVE] else 28,
                "max_target_seqs": 100
            },
            "target_genes": family.detection_markers,
            "genome_characteristics": {
                "type": family.genome_type.value,
                "size_range": family.genome_size_range,
                "envelope": family.envelope
            }
        }
        
        return strategy
    
    def analyze_viral_ecology(self, detected_families: List[str], 
                            mosquito_species: str) -> Dict[str, any]:
        """Analyze viral ecology patterns and potential interactions."""
        analysis = {
            "detected_families": detected_families,
            "mosquito_species": mosquito_species,
            "pathogenicity_profile": {},
            "transmission_modes": set(),
            "exclusion_violations": [],
            "co_occurrence_patterns": [],
            "ecological_insights": []
        }
        
        # Analyze pathogenicity
        for family_name in detected_families:
            family = self.families.get(family_name)
            if family:
                analysis["pathogenicity_profile"][family_name] = family.pathogenicity
                analysis["transmission_modes"].update(family.transmission_modes)
        
        # Check for exclusion violations
        for family_name in detected_families:
            family = self.families.get(family_name)
            if family:
                for excluded in family.exclusion_patterns:
                    if excluded in detected_families:
                        analysis["exclusion_violations"].append((family_name, excluded))
        
        # Generate ecological insights
        if len(detected_families) > 1:
            analysis["ecological_insights"].append(
                "Multiple viral families detected - indicates diverse viral ecology"
            )
        
        if analysis["exclusion_violations"]:
            analysis["ecological_insights"].append(
                "Exclusion pattern violations detected - unusual co-occurrence"
            )
        
        high_path_count = sum(1 for p in analysis["pathogenicity_profile"].values() 
                             if p == "high")
        if high_path_count > 0:
            analysis["ecological_insights"].append(
                f"High pathogenicity viruses detected ({high_path_count}) - public health concern"
            )
        
        return analysis

def get_viral_family_classifier() -> ViralFamilyClassifier:
    """Factory function to get viral family classifier instance."""
    return ViralFamilyClassifier()