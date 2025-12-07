"""
scPPIN: Single-cell Protein-Protein Interaction Network Analysis

A Python package for detecting functional modules in protein-protein interaction
networks by integrating single-cell RNA sequencing data.
"""

__version__ = '0.3.0'
__author__ = 'Shahroze Abbas'

# Main API - Class only
from .analyzer import scPPIN

__all__ = [
    'scPPIN',  # Only export the class
]
