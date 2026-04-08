"""
Immunex: GROMACS MD Analysis Toolkit

A comprehensive toolkit for analyzing molecular dynamics simulation results
from GROMACS using MDAnalysis and GROMACS built-in tools.
"""

__version__ = "0.1.0"
__author__ = "Research Team"

# Utils - Tool functions
from .utils import PlotManager, PathManager

# Analysis modules - Core functionality
from .analysis import (
    trajectory, structure, quality,
    PBCProcessor,
    RMSDCalculator, RDFCalculator, RadiusGyrationCalculator,
    DistanceCalculator, HydrogenBondAnalyzer, RMSFAnalyzer,
    BFactorAnalyzer, ContactMapCalculator,
    GeometryAnalyzer, AtomInfoExtractor,
    pHLATCRAnalyzer, pHLATCRHydrogenBondAnalyzer, ComplexAngleAnalyzer
)

# Quality analysis modules
from .analysis.quality import (
    MDCompletenessChecker, StructureValidator,
    BatchTracker, QualityReporter, EnergyQualityChecker,
    PostPBCValidator
)

# Trajectory convergence analysis
from .analysis.trajectory import RMSDConvergenceAnalyzer

# Pipeline modules
from .pipeline import PBCRMSDPipeline, QualityAssessmentPipeline

# Batch processing - Simple interface
from .pipeline import process_md_tasks, discover_md_tasks, check_task_status

# SLURM cluster support
from .cluster.slurm_generator import generate_slurm_scripts_for_md_tasks

__all__ = [
    # Utils
    "PlotManager",
    "PathManager",
    # Analysis submodules
    "trajectory",
    "structure",
    "quality",
    # Trajectory analysis
    "RMSDCalculator",
    "PBCProcessor",
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
    "PostPBCValidator",
    "RMSDConvergenceAnalyzer",
    # Pipeline modules
    "PBCRMSDPipeline",
    "QualityAssessmentPipeline",
    "process_md_tasks",
    "discover_md_tasks",
    "check_task_status",
    "generate_slurm_scripts_for_md_tasks"
]
