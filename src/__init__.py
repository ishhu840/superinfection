"""Viral Interference Analysis System

A computational system for studying viral interference and superinfection exclusion
in mosquitoes using AI/ML approaches and interactive visualization.
"""

__version__ = "0.1.0"
__author__ = "Viral Interference Research Team"
__email__ = "contact@example.com"

# Import main modules
from . import data_processing
from . import statistical_analysis
from . import machine_learning
from . import visualization
from . import utils

__all__ = [
    "data_processing",
    "statistical_analysis", 
    "machine_learning",
    "visualization",
    "utils"
]