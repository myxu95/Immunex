# Lazy import submodules to avoid MDAnalysis dependency
from . import quality

# Import individual modules for convenience
from .trajectory import (
    pHLATCRAnalyzer,
    pHLATCRHydrogenBondAnalyzer,
    ComplexAngleAnalyzer,
    RMSDCalculator,
    RDFCalculator,
    RadiusGyrationCalculator,
    DistanceCalculator,
    HydrogenBondAnalyzer,
    RMSFAnalyzer,
    extract_sequence_from_topology,
    find_subsequence_position,
    calculate_cdr3_rmsf,
    analyze_phla_tcr_rmsf
)

from .structure import (
    BFactorAnalyzer,
    ContactMapCalculator,
    GeometryAnalyzer,
    AtomInfoExtractor
)

# Import quality analysis (no MDAnalysis dependency)
from .quality import (
    MDCompletenessChecker,
    StructureValidator,
    BatchTracker,
    QualityReporter
)

# Import allostery analysis
from .allostery import (
    ContactCorrelationAnalyzer
)

# Import free energy analysis
from .free_energy import (
    FELCalculator,
    FELVisualizer
)

__all__ = [
    # Submodules
    "trajectory",
    "structure",
    "quality",
    "allostery",
    "free_energy",
    # Trajectory analysis modules (lazy loaded)
    "RMSDCalculator",
    "RDFCalculator",
    "RadiusGyrationCalculator",
    "DistanceCalculator",
    "HydrogenBondAnalyzer",
    "RMSFAnalyzer",
    "extract_sequence_from_topology",
    "find_subsequence_position",
    "calculate_cdr3_rmsf",
    "analyze_phla_tcr_rmsf",
    "pHLATCRAnalyzer",
    "pHLATCRHydrogenBondAnalyzer",
    "ComplexAngleAnalyzer",
    # Structure analysis modules
    "BFactorAnalyzer",
    "ContactMapCalculator",
    "GeometryAnalyzer",
    "AtomInfoExtractor",
    # Quality analysis modules
    "MDCompletenessChecker",
    "StructureValidator",
    "BatchTracker",
    "QualityReporter",
    # Allostery analysis modules
    "ContactCorrelationAnalyzer",
    # Free energy analysis modules
    "FELCalculator",
    "FELVisualizer"
]

def __getattr__(name):
    """Lazy import for trajectory analysis classes that require MDAnalysis."""
    if name == "trajectory":
        from . import trajectory
        return trajectory
    elif name in ["RMSDCalculator", "RDFCalculator", "RadiusGyrationCalculator",
                   "DistanceCalculator", "HydrogenBondAnalyzer"]:
        from . import trajectory
        return getattr(trajectory, name)
    elif name == "structure":
        from . import structure
        return structure
    elif name == "allostery":
        from . import allostery
        return allostery
    elif name == "ContactCorrelationAnalyzer":
        from . import allostery
        return getattr(allostery, name)
    elif name == "free_energy":
        from . import free_energy
        return free_energy
    elif name in ["FELCalculator", "FELVisualizer"]:
        from . import free_energy
        return getattr(free_energy, name)
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")