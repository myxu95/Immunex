"""
Immunex Pipeline Module

This module provides high-level pipelines that coordinate functional modules
to perform end-to-end MD trajectory processing and quality assessment.

New Architecture (2026-03-16):
- Pipeline base class (Layer 3: Orchestration)
- Pipeline nodes (Layer 2: Business logic wrappers)
- Standard pipelines (Pre-configured templates)
- BatchExecutor (Layer 4: Batch processing)

Legacy Pipelines:
- PBCRMSDPipeline: PBC correction + RMSD quality assessment pipeline
- QualityAssessmentPipeline: Comprehensive quality assessment pipeline

Author: Immunex Development Team
Date: 2026-03-15
"""

# New architecture (recommended)
from .base_pipeline import Pipeline
from .standard_pipelines import (
    StandardTrajectoryPipeline,
    QuickRMSDPipeline,
    PreprocessOnlyPipeline,
    PreprocessQualityPipeline,
)
from .analysis_pipelines import (
    TCRRMSDPipeline,
    CDRRMSDPipeline,
    pHLARMSDPipeline,
    HLAAlphaRMSDPipeline,
    InterfaceRMSDPipeline,
    AnnotatedRMSFPipeline,
    BiologicalIdentityPipeline,
    BSAPipeline,
    NormalModePipeline,
    InterfaceClusteringPipeline,
    DockingAnglePipeline,
    ContactFrequencyPipeline,
    HydrogenBondInteractionPipeline,
    SaltBridgeInteractionPipeline,
    HydrophobicInteractionPipeline,
    PiStackingInteractionPipeline,
    CationPiInteractionPipeline,
    AllosteryAnalysisPipeline,
    ComprehensiveAnalysisPipeline,
    CDRRMSFPipeline
)
from .analysis_workflows import (
    AnalysisPipeline,
    create_analysis_pipeline,
)
from .batch_executor import BatchExecutor
from .batch_workflow import (
    process_md_tasks,
    discover_md_tasks,
    check_task_status,
)
from .nodes import (
    PreprocessNode,
    RMSDNode,
    RMSDPlotNode,
    RMSDQualityNode,
    ChainIdentificationNode,
    IndexGenerationNode,
    CDRDetectionNode,
    RMSFNode,
    ResidueRMSFNode,
    BiologicalIdentityNode,
    BSAAnalysisNode,
    NormalModeNode,
    InterfaceClusteringNode,
    DockingAngleNode,
    ContactFrequencyNode,
    ContactAnnotationNode,
    ContactHeatmapNode,
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

# Legacy pipelines (backward compatibility)
from .pbc_rmsd_pipeline import PBCRMSDPipeline
from .quality_assessment_pipeline import QualityAssessmentPipeline

__all__ = [
    # New architecture
    'Pipeline',
    'BatchExecutor',
    'process_md_tasks',
    'discover_md_tasks',
    'check_task_status',
    'StandardTrajectoryPipeline',
    'QuickRMSDPipeline',
    'PreprocessOnlyPipeline',
    'PreprocessQualityPipeline',
    # Pipeline nodes
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
    'InterfaceClusteringNode',
    'DockingAngleNode',
    'ContactFrequencyNode',
    'ContactAnnotationNode',
    'ContactHeatmapNode',
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
    # Analysis pipelines
    'TCRRMSDPipeline',
    'CDRRMSDPipeline',
    'pHLARMSDPipeline',
    'HLAAlphaRMSDPipeline',
    'InterfaceRMSDPipeline',
    'AnnotatedRMSFPipeline',
    'BiologicalIdentityPipeline',
    'BSAPipeline',
    'NormalModePipeline',
    'InterfaceClusteringPipeline',
    'DockingAnglePipeline',
    'ContactFrequencyPipeline',
    'HydrogenBondInteractionPipeline',
    'SaltBridgeInteractionPipeline',
    'HydrophobicInteractionPipeline',
    'PiStackingInteractionPipeline',
    'CationPiInteractionPipeline',
    'AllosteryAnalysisPipeline',
    'ComprehensiveAnalysisPipeline',
    'CDRRMSFPipeline',
    'AnalysisPipeline',
    'create_analysis_pipeline',
    # Legacy
    'PBCRMSDPipeline',
    'QualityAssessmentPipeline'
]
