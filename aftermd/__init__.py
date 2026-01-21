"""
AfterMD: GROMACS MD Analysis Toolkit

A comprehensive toolkit for analyzing molecular dynamics simulation results
from GROMACS using MDAnalysis and GROMACS built-in tools.
"""

__version__ = "0.1.0"
__author__ = "Research Team"

# Utils - Tool functions
from .utils import BatchProcessor, PlotManager, PathManager

# Analysis modules - Core functionality
from .analysis import (
    trajectory, structure, quality,
    RMSDCalculator, RDFCalculator, RadiusGyrationCalculator,
    DistanceCalculator, HydrogenBondAnalyzer, RMSFAnalyzer,
    BFactorAnalyzer, ContactMapCalculator,
    GeometryAnalyzer, AtomInfoExtractor,
    pHLATCRAnalyzer, pHLATCRHydrogenBondAnalyzer, ComplexAngleAnalyzer
)

# Quality analysis modules
from .analysis.quality import (
    MDCompletenessChecker, StructureValidator,
    BatchTracker, QualityReporter, EnergyQualityChecker
)

# Preprocessing modules
# from .preprocessing import PBCProcessor  # TODO: Implement PBCProcessor

# Batch processing - Simple interface
# from .batch_process import process_md_tasks, discover_md_tasks, check_task_status  # TODO: Fix import

# SLURM cluster support
# from .utils.slurm_generator import generate_slurm_scripts_for_md_tasks  # TODO: Implement slurm_generator

__all__ = [
    # Utils
    "BatchProcessor",
    "PlotManager",
    "PathManager",
    # Analysis submodules
    "trajectory",
    "structure",
    "quality",
    # Trajectory analysis
    "RMSDCalculator",
    "RDFCalculator",
    "RadiusGyrationCalculator",
    "DistanceCalculator",
    "HydrogenBondAnalyzer",
    "RMSFAnalyzer",
    # pHLA-TCR specific analysis
    "pHLATCRAnalyzer",
    "pHLATCRHydrogenBondAnalyzer",
    "ComplexAngleAnalyzer",
    # Structure analysis
    "BFactorAnalyzer",
    "ContactMapCalculator",
    "GeometryAnalyzer",
    "AtomInfoExtractor",
    # Quality analysis
    "MDCompletenessChecker",
    "StructureValidator",
    "BatchTracker",
    "QualityReporter",
    # Preprocessing - TODO: Add when implemented
    # "PBCProcessor",
    # Batch processing - TODO: Fix imports
    # "process_md_tasks",
    # "discover_md_tasks",
    # "check_task_status",
    # SLURM cluster support - TODO: Implement
    # "generate_slurm_scripts_for_md_tasks"
]