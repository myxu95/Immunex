"""Allostery Analysis Module

Dynamic correlation analysis for studying cooperative motion
and allosteric communication pathways within proteins.

This module provides tools to:
- Calculate residue-residue contact time series
- Compute Pearson correlation between contact dynamics
- Identify synchronized and anti-correlated motions
- Visualize allosteric networks

Main Classes:
    ContactCorrelationAnalyzer: Dynamic cross-correlation of residue contacts
"""

from .contact_correlation import ContactCorrelationAnalyzer

__all__ = ['ContactCorrelationAnalyzer']
