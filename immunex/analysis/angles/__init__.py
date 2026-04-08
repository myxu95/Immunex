"""
Angle Analysis Module for MD Simulations

This module provides docking angle calculation for TCR-pMHC complexes.

Public API (3 classes only):
-----------------------------
- DockingAngleAnalyzer: Main analyzer class
- DockingAngleInput: Standardized input parameters
- DockingAngleResult: Standardized output results

All internal implementation details are hidden within analyzer.py module.

Key Features:
-------------
- Sequence-based MHC region detection (robust to non-standard PDB numbering)
- ANARCI integration for TCR disulfide identification
- Geometric angle definitions (Crossing + Incident)
- Pure Python implementation (MDAnalysis + NumPy)
- Support for both static structures and MD trajectories
- Conforms to Immunex 6 design principles

Example Usage:
--------------
Basic trajectory analysis:
>>> from immunex.analysis.angles import DockingAngleAnalyzer, DockingAngleInput
>>> analyzer = DockingAngleAnalyzer()
>>> result = analyzer.analyze(DockingAngleInput(
...     topology='md.tpr',
...     trajectory='md_pbc.xtc',
...     stride=10,
...     output_dir='./results'
... ))
>>> print(f"Crossing: {result.statistics['crossing_mean']:.2f}°")

With progress tracking:
>>> def progress_callback(progress, message):
...     print(f"[{progress*100:.0f}%] {message}")
>>> analyzer.set_progress_callback(progress_callback)
>>> result = analyzer.analyze(input_params)

Notes:
------
- Refactored 2026-03-18: Consolidated from 9 files to 3 files
- All internal classes are hidden (prefixed with _)
- Only 3 classes are exported in public API
- Old modules archived for reference
"""

# Public API - Only 3 classes exported
from .angle_data_structures import DockingAngleInput, DockingAngleResult
from .analyzer import DockingAngleAnalyzer

__all__ = [
    'DockingAngleAnalyzer',
    'DockingAngleInput',
    'DockingAngleResult',
]

# Version info
__version__ = '4.3.0'
__refactoring_date__ = '2026-03-18'
__specification__ = 'Consolidated module with hidden internals'
