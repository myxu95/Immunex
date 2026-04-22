"""
Pipeline Nodes - Layer 2 of the architecture.

Pipeline nodes are lightweight wrappers around Core Modules (Layer 1).
They handle context management, input validation, and error handling.

Subpackages currently group:
- interactions: contact and typed-interaction nodes
- topology: chain mapping, CDR detection, index generation, and identity nodes
- trajectory: preprocessing, RMSD, RMSF, and derived quality nodes
- interface: buried surface area and interface state nodes
- allostery: normal mode and perturbation-response nodes
- comparison: cross-condition/system comparison nodes
- geometry: docking and geometric descriptor nodes
"""

from .trajectory import (
    PreprocessNode,
    RMSDNode,
    RMSDPlotNode,
    RMSDQualityNode,
    RMSFNode,
    ResidueRMSFNode,
)
from .topology import (
    ChainIdentificationNode,
    IndexGenerationNode,
    CDRDetectionNode,
    BiologicalIdentityNode,
)
from .interface import (
    BSAAnalysisNode,
    InterfaceClusteringNode,
)
from .allostery import NormalModeNode
from .comparison import SystemComparisonNode
from .geometry import DockingAngleNode
from .interactions import (
    RRCSNode,
    ContactFrequencyNode,
    ContactAnnotationNode,
    ContactHeatmapNode,
    InteractionOccupancyNode,
    HydrogenBondPairNode,
    HydrogenBondAnnotationNode,
    HydrogenBondHeatmapNode,
    SaltBridgePairNode,
    SaltBridgeAnnotationNode,
    SaltBridgeHeatmapNode,
    HydrophobicContactPairNode,
    HydrophobicAnnotationNode,
    HydrophobicHeatmapNode,
    PiStackingPairNode,
    PiStackingAnnotationNode,
    PiStackingHeatmapNode,
    CationPiPairNode,
    CationPiAnnotationNode,
    CationPiHeatmapNode,
)

__all__ = [
    'PreprocessNode',
    'RMSDNode',
    'RMSDPlotNode',
    'RMSDQualityNode',
    'ChainIdentificationNode',
    'IndexGenerationNode',
    'CDRDetectionNode',
    'RMSFNode',
    'ResidueRMSFNode',
    'BiologicalIdentityNode',
    'BSAAnalysisNode',
    'NormalModeNode',
    'RRCSNode',
    'InterfaceClusteringNode',
    'SystemComparisonNode',
    'DockingAngleNode',
    'ContactFrequencyNode',
    'ContactAnnotationNode',
    'ContactHeatmapNode',
    'InteractionOccupancyNode',
    'HydrogenBondPairNode',
    'HydrogenBondAnnotationNode',
    'HydrogenBondHeatmapNode',
    'SaltBridgePairNode',
    'SaltBridgeAnnotationNode',
    'SaltBridgeHeatmapNode',
    'HydrophobicContactPairNode',
    'HydrophobicAnnotationNode',
    'HydrophobicHeatmapNode',
    'PiStackingPairNode',
    'PiStackingAnnotationNode',
    'PiStackingHeatmapNode',
    'CationPiPairNode',
    'CationPiAnnotationNode',
    'CationPiHeatmapNode',
]
