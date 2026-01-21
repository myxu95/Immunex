"""
Angle analysis module for TCR-pMHC complexes

This module provides tools for:
1. Static angle clustering (identify binding modes)
2. MD trajectory angle stability analysis
"""

from .angle_clustering import AngleClusterer, TrajectoryStabilityAnalyzer

__all__ = ['AngleClusterer', 'TrajectoryStabilityAnalyzer']
