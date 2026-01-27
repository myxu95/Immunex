"""
Angle Analysis Module for MD Simulations

This module provides comprehensive angle calculation tools for molecular dynamics simulations,
particularly for protein-protein complexes like TCR-pMHC.

Submodules:
-----------
- principal_axes: Principal axes calculation based on inertia tensor
- docking_angles: TCR-pMHC docking angle analysis (Twist/Tilt/Swing)
- dihedral: General dihedral angle calculations (backbone phi/psi, sidechain chi)
- vector_angles: Vector angle utilities (principal axis angles, COM angles)
- trajectory_angles: MD trajectory angle analysis and time evolution

Key Features:
-------------
- Pure Python implementation (MDAnalysis + NumPy)
- Support for both static structures and MD trajectories
- Flexible chain/residue selection
- Batch processing capabilities
- Integration with AfterMD analysis pipeline

Example Usage:
--------------
Static structure analysis:
>>> from aftermd.analysis.angles import DockingAngleAnalyzer
>>> analyzer = DockingAngleAnalyzer('structure.pdb')
>>> twist, tilt, swing = analyzer.calculate_docking_angles(
...     mhc_selection='segid A',
...     tcr_selection='segid D E'
... )

MD trajectory analysis:
>>> from aftermd.analysis.angles import TrajectoryAngleAnalyzer
>>> analyzer = TrajectoryAngleAnalyzer('md.tpr', 'md.xtc')
>>> times, angles = analyzer.calculate_angle_evolution(
...     angle_type='docking',
...     output_file='docking_angles.csv'
... )

Notes:
------
- This is a complete redesign (2026-01-23)
- Old C++ tools archived in development/archived_scripts/
- For legacy code see: development/archived_scripts/README_ANGLES_ARCHIVE.md
"""

# Import implemented modules
from .principal_axes import PrincipalAxesCalculator
from .vector_angles import (
    angle_between_vectors,
    project_vector_onto_plane,
    dihedral_angle,
    signed_angle_between_vectors
)
from .docking_angles import DockingAngleAnalyzer

__all__ = [
    'PrincipalAxesCalculator',
    'angle_between_vectors',
    'project_vector_onto_plane',
    'dihedral_angle',
    'signed_angle_between_vectors',
    'DockingAngleAnalyzer',
]

# Version info
__version__ = '2.0.0'
__redesign_date__ = '2026-01-23'
