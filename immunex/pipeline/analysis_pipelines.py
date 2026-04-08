"""
Standard Analysis Pipelines - Pre-configured pipeline templates for common analyses.

This module provides commonly used analysis pipeline configurations for
TCR-pMHC complex analysis, eliminating code duplication across batch scripts.

Usage:
    >>> from pathlib import Path
    >>> from immunex.pipeline import TCRRMSDPipeline, BatchExecutor
    >>> from immunex.core import discover_tasks
    >>>
    >>> report = discover_tasks(Path("/data/simulations"))
    >>> executor = BatchExecutor(max_workers=4)
    >>> pipeline = TCRRMSDPipeline()
    >>> results = executor.execute_pipeline(report, pipeline)

Author: Immunex Development Team
Date: 2026-03-16
"""

from .base_pipeline import Pipeline
from .nodes import (
    PreprocessNode,
    RMSDNode,
    ChainIdentificationNode,
    IndexGenerationNode,
    CDRDetectionNode,
    RMSFNode,
    ResidueRMSFNode,
    BiologicalIdentityNode,
    BSAAnalysisNode,
    NormalModeNode,
    InterfaceClusteringNode,
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


class TCRRMSDPipeline(Pipeline):
    """
    TCR RMSD analysis pipeline.

    This pipeline performs:
    1. Chain identification using ANARCI
    2. Index file generation (pHLA and TCR groups)
    3. RMSD calculation for TCR aligned to pHLA

    Example:
        >>> from immunex.core import PipelineContext
        >>> context = PipelineContext(
        ...     system_id="1ao7",
        ...     topology="data/1ao7/md.tpr",
        ...     trajectory_raw="data/1ao7/md_processed.xtc",
        ...     structure_pdb="data/1ao7/md_converted.pdb"
        ... )
        >>> pipeline = TCRRMSDPipeline()
        >>> result = pipeline.execute(context)
        >>> print(f"Mean RMSD: {result.results['rmsd']['mean_rmsd']:.3f} nm")
    """

    def __init__(self,
                 use_anarci: bool = True,
                 rmsd_selection: str = "protein and name CA"):
        """
        Initialize TCR RMSD pipeline.

        Args:
            use_anarci: Whether to use ANARCI for chain identification
            rmsd_selection: Atom selection for RMSD calculation
        """
        nodes = [
            ChainIdentificationNode(
                method="anarci" if use_anarci else "heuristic",
                fallback_to_heuristic=True
            ),
            IndexGenerationNode(
                group_definitions={
                    'pHLA': 'chainID {mhc_alpha} {b2m} {peptide}',
                    'TCR': 'chainID {tcr_alpha} {tcr_beta}'
                },
                output_name="tcr_rmsd.ndx"
            ),
            RMSDNode(
                selection=rmsd_selection,
                method="mdanalysis"
            )
        ]
        super().__init__(nodes=nodes)


class CDRRMSDPipeline(Pipeline):
    """
    CDR region RMSD analysis pipeline.

    This pipeline performs:
    1. Chain identification
    2. Index file generation (CDR regions)
    3. RMSD calculation for CDR regions

    Note: Requires CDR detection capability (future implementation)

    Example:
        >>> pipeline = CDRRMSDPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self, cdr_regions: list = None):
        """
        Initialize CDR RMSD pipeline.

        Args:
            cdr_regions: List of CDR regions to analyze (e.g., ['CDR3_alpha', 'CDR3_beta'])
        """
        if cdr_regions is None:
            cdr_regions = ['CDR1', 'CDR2', 'CDR3']

        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            # TODO: Add CDRDetectionNode when implemented
            # TODO: Add IndexGenerationNode with CDR-specific groups
            # TODO: Add RMSDNode for each CDR region
        ]
        super().__init__(nodes=nodes)
        self.cdr_regions = cdr_regions


class pHLARMSDPipeline(Pipeline):
    """
    pHLA RMSD analysis pipeline.

    This pipeline calculates RMSD for pHLA complex components.

    Example:
        >>> pipeline = pHLARMSDPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self):
        """Initialize pHLA RMSD pipeline."""
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            IndexGenerationNode(
                group_definitions={
                    'pHLA': 'chainID {mhc_alpha} {b2m} {peptide}',
                    'MHC_alpha': 'chainID {mhc_alpha}',
                    'B2M': 'chainID {b2m}',
                    'Peptide': 'chainID {peptide}'
                },
                output_name="phla_rmsd.ndx"
            ),
            RMSDNode(selection="protein and name CA", method="mdanalysis")
        ]
        super().__init__(nodes=nodes)


class HLAAlphaRMSDPipeline(Pipeline):
    """
    HLA alpha chain RMSD analysis pipeline.

    This pipeline calculates RMSD specifically for the HLA alpha chain.

    Example:
        >>> pipeline = HLAAlphaRMSDPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self):
        """Initialize HLA alpha RMSD pipeline."""
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            IndexGenerationNode(
                group_definitions={
                    'HLA_alpha': 'chainID {mhc_alpha}'
                },
                output_name="hla_alpha.ndx"
            ),
            RMSDNode(selection="protein and name CA", method="mdanalysis")
        ]
        super().__init__(nodes=nodes)


class InterfaceRMSDPipeline(Pipeline):
    """
    Interface RMSD analysis pipeline.

    This pipeline calculates RMSD at the TCR-pMHC interface.

    Example:
        >>> pipeline = InterfaceRMSDPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self, interface_cutoff: float = 6.0):
        """
        Initialize interface RMSD pipeline.

        Args:
            interface_cutoff: Distance cutoff for interface definition (Angstrom)
        """
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            # TODO: Add InterfaceDetectionNode
            # TODO: Add RMSD calculation for interface residues
        ]
        super().__init__(nodes=nodes)
        self.interface_cutoff = interface_cutoff


class AnnotatedRMSFPipeline(Pipeline):
    """
    区域化 residue-level RMSF 分析 pipeline。

    这条主线用于单体系柔性分析，输出：
    - residue_rmsf.csv
    - region_rmsf_summary.csv
    - TCR / pHLA / region 三张图
    """

    def __init__(
        self,
        stride: int = 1,
        time_unit: str = "ps",
        selection: str | None = None,
        auto_identify_chains: bool = True,
        auto_detect_cdr: bool = True,
    ):
        nodes = []
        if auto_identify_chains:
            nodes.append(ChainIdentificationNode(method="anarci", fallback_to_heuristic=True))
        if auto_detect_cdr:
            nodes.append(CDRDetectionNode())
        nodes.append(
            ResidueRMSFNode(
                selection=selection,
                stride=stride,
                time_unit=time_unit,
            )
        )
        super().__init__(nodes=nodes)


class DockingAnglePipeline(Pipeline):
    """
    Docking angle analysis pipeline.

    This pipeline calculates TCR-pMHC docking angles (Crossing + Incident).
    Supports automatic chain identification and flexible configuration.

    Example:
        >>> from immunex.pipeline import DockingAnglePipeline
        >>> from immunex.core import PipelineContext
        >>>
        >>> context = PipelineContext(
        ...     system_id="1ao7",
        ...     topology="md.tpr",
        ...     trajectory_raw="md_pbc.xtc"
        ... )
        >>> pipeline = DockingAnglePipeline(stride=10)
        >>> result = pipeline.execute(context)
        >>> angles = result.results['docking_angles']
        >>> print(f"Crossing: {angles['statistics']['crossing_mean']:.2f}°")
    """

    def __init__(
        self,
        stride: int = 1,
        auto_identify_chains: bool = True,
        print_each_frame: bool = False,
    ):
        """
        Initialize docking angle pipeline.

        Args:
            stride: Frame stride for angle calculation
            auto_identify_chains: Whether to use automatic chain identification
            print_each_frame: 是否逐帧打印角度
        """
        from .nodes import DockingAngleNode

        nodes = []

        # Add chain identification if auto mode
        if auto_identify_chains:
            nodes.append(
                ChainIdentificationNode(method="anarci", fallback_to_heuristic=True)
            )

        # Add docking angle node
        nodes.append(
            DockingAngleNode(
                stride=stride,
                auto_identify_chains=auto_identify_chains,
                print_each_frame=print_each_frame,
            )
        )

        super().__init__(nodes=nodes)


class HydrogenBondInteractionPipeline(Pipeline):
    """pHLA-TCR 氢键相互作用分析 pipeline。"""

    def __init__(
        self,
        stride: int = 1,
        distance_cutoff: float = 3.5,
        angle_cutoff: float = 150.0,
    ):
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            CDRDetectionNode(),
            HydrogenBondPairNode(
                stride=stride,
                distance_cutoff=distance_cutoff,
                angle_cutoff=angle_cutoff,
            ),
            HydrogenBondAnnotationNode(),
            InteractionOccupancyNode(
                source_result_key="hbond_annotation",
                report_file_key="hbond_report_file",
                output_subdir="interactions/hydrogen_bonds/occupancy",
                family_name="hbond",
                result_key="hbond_occupancy",
            ),
            HydrogenBondHeatmapNode(),
        ]
        super().__init__(nodes=nodes)


class SaltBridgeInteractionPipeline(Pipeline):
    """pHLA-TCR 盐桥相互作用分析 pipeline。"""

    def __init__(self, stride: int = 1, distance_cutoff: float = 4.0):
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            CDRDetectionNode(),
            SaltBridgePairNode(stride=stride, distance_cutoff=distance_cutoff),
            SaltBridgeAnnotationNode(),
            InteractionOccupancyNode(
                source_result_key="salt_bridge_annotation",
                report_file_key="salt_bridge_report_file",
                output_subdir="interactions/salt_bridges/occupancy",
                family_name="saltbridge",
                result_key="saltbridge_occupancy",
            ),
            SaltBridgeHeatmapNode(),
        ]
        super().__init__(nodes=nodes)


class HydrophobicInteractionPipeline(Pipeline):
    """pHLA-TCR 疏水接触分析 pipeline。"""

    def __init__(self, stride: int = 1, distance_cutoff: float = 4.5):
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            CDRDetectionNode(),
            HydrophobicContactPairNode(stride=stride, distance_cutoff=distance_cutoff),
            HydrophobicAnnotationNode(),
            InteractionOccupancyNode(
                source_result_key="hydrophobic_annotation",
                report_file_key="hydrophobic_report_file",
                output_subdir="interactions/hydrophobic_contacts/occupancy",
                family_name="hydrophobic",
                result_key="hydrophobic_occupancy",
            ),
            HydrophobicHeatmapNode(),
        ]
        super().__init__(nodes=nodes)


class PiStackingInteractionPipeline(Pipeline):
    """pHLA-TCR pi-pi 相互作用分析 pipeline。"""

    def __init__(self, stride: int = 1, distance_cutoff: float = 6.5):
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            CDRDetectionNode(),
            PiStackingPairNode(stride=stride, distance_cutoff=distance_cutoff),
            PiStackingAnnotationNode(),
            InteractionOccupancyNode(
                source_result_key="pi_pi_annotation",
                report_file_key="pi_pi_report_file",
                output_subdir="interactions/pi_interactions/occupancy",
                family_name="pipi",
                result_key="pipi_occupancy",
            ),
            PiStackingHeatmapNode(),
        ]
        super().__init__(nodes=nodes)


class CationPiInteractionPipeline(Pipeline):
    """pHLA-TCR cation-pi 相互作用分析 pipeline。"""

    def __init__(self, stride: int = 1, distance_cutoff: float = 6.0, normal_angle_cutoff: float = 60.0):
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            CDRDetectionNode(),
            CationPiPairNode(
                stride=stride,
                distance_cutoff=distance_cutoff,
                normal_angle_cutoff=normal_angle_cutoff,
            ),
            CationPiAnnotationNode(),
            InteractionOccupancyNode(
                source_result_key="cation_pi_annotation",
                report_file_key="cation_pi_report_file",
                output_subdir="interactions/cation_pi_interactions/occupancy",
                family_name="cationpi",
                result_key="cationpi_occupancy",
            ),
            CationPiHeatmapNode(),
        ]
        super().__init__(nodes=nodes)


class ContactFrequencyPipeline(Pipeline):
    """
    Contact frequency analysis pipeline.

    This pipeline calculates residue-residue contact frequencies.

    Example:
        >>> pipeline = ContactFrequencyPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self, cutoff: float = 4.5, stride: int = 1, min_frequency: float = 0.0):
        """
        Initialize contact frequency pipeline.

        Args:
            cutoff: Distance cutoff for contacts (Angstrom)
        """
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            CDRDetectionNode(
                allow_fallback=True,
                include_cdrs=[1, 2, 3],
            ),
            ContactFrequencyNode(
                cutoff=cutoff,
                stride=stride,
                min_frequency=min_frequency,
            ),
            ContactAnnotationNode(),
            InteractionOccupancyNode(
                source_result_key="contact_annotation",
                report_file_key="contact_report_file",
                output_subdir="contacts/occupancy",
                family_name="contact",
                result_key="contact_occupancy",
            ),
            ContactHeatmapNode(),
        ]
        super().__init__(nodes=nodes)


class BSAPipeline(Pipeline):
    """pHLA-TCR buried surface area 分析 pipeline。"""

    def __init__(
        self,
        probe_radius: float = 1.4,
        stride: int = 1,
        time_unit: str = "ps",
        selection_a: str | None = None,
        selection_b: str | None = None,
        auto_identify_chains: bool = True,
    ):
        nodes = []
        if auto_identify_chains and (selection_a is None or selection_b is None):
            nodes.append(ChainIdentificationNode(method="anarci", fallback_to_heuristic=True))
        nodes.append(
            BSAAnalysisNode(
                selection_a=selection_a,
                selection_b=selection_b,
                probe_radius=probe_radius,
                stride=stride,
                time_unit=time_unit,
            )
        )
        super().__init__(nodes=nodes)


class BiologicalIdentityPipeline(Pipeline):
    """单体系基础生物学身份注释 pipeline。"""

    def __init__(
        self,
        cdr_metadata_path: str | None = None,
        auto_identify_chains: bool = True,
        auto_detect_cdr: bool = True,
    ):
        nodes = []
        if auto_identify_chains:
            nodes.append(ChainIdentificationNode(method="anarci", fallback_to_heuristic=True))
        if auto_detect_cdr:
            nodes.append(CDRDetectionNode())
        nodes.append(BiologicalIdentityNode(cdr_metadata_path=cdr_metadata_path))
        super().__init__(nodes=nodes)


class NormalModePipeline(Pipeline):
    """单体系 normal mode / PRS 分析 pipeline。"""

    def __init__(
        self,
        cutoff_angstrom: float = 10.0,
        n_low_modes: int = 10,
        prs_force_directions: int = 8,
        auto_identify_chains: bool = True,
        auto_detect_cdr: bool = True,
    ):
        nodes = []
        if auto_identify_chains:
            nodes.append(ChainIdentificationNode(method="anarci", fallback_to_heuristic=True))
        if auto_detect_cdr:
            nodes.append(CDRDetectionNode())
        nodes.append(
            NormalModeNode(
                cutoff_angstrom=cutoff_angstrom,
                n_low_modes=n_low_modes,
                prs_force_directions=prs_force_directions,
            )
        )
        super().__init__(nodes=nodes)


class InterfaceClusteringPipeline(Pipeline):
    """关键界面状态聚类 pipeline。"""

    def __init__(
        self,
        stride: int = 10,
        contact_cutoff_angstrom: float = 4.5,
        distance_cutoff: float = 0.35,
        linkage_method: str = "average",
        geometry_weight: float = 0.45,
        sidechain_weight: float = 0.25,
        interaction_weight: float = 0.30,
        auto_identify_chains: bool = True,
        auto_detect_cdr: bool = True,
    ):
        nodes = []
        if auto_identify_chains:
            nodes.append(ChainIdentificationNode(method="anarci", fallback_to_heuristic=True))
        if auto_detect_cdr:
            nodes.append(CDRDetectionNode())
        nodes.append(
            InterfaceClusteringNode(
                stride=stride,
                contact_cutoff_angstrom=contact_cutoff_angstrom,
                distance_cutoff=distance_cutoff,
                linkage_method=linkage_method,
                geometry_weight=geometry_weight,
                sidechain_weight=sidechain_weight,
                interaction_weight=interaction_weight,
            )
        )
        super().__init__(nodes=nodes)


class AllosteryAnalysisPipeline(Pipeline):
    """
    Allostery analysis pipeline.

    This pipeline performs contact correlation analysis for allosteric communication.

    Example:
        >>> pipeline = AllosteryAnalysisPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self):
        """Initialize allostery analysis pipeline."""
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            # TODO: Add AllosteryAnalysisNode when implemented
        ]
        super().__init__(nodes=nodes)


class ComprehensiveAnalysisPipeline(Pipeline):
    """
    Comprehensive analysis pipeline combining multiple analyses.

    This pipeline performs:
    1. Chain identification
    2. Index file generation
    3. RMSD analysis (TCR, pHLA)
    4. Docking angle analysis
    5. Contact frequency analysis

    Example:
        >>> pipeline = ComprehensiveAnalysisPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self):
        """Initialize comprehensive analysis pipeline."""
        nodes = [
            ChainIdentificationNode(method="anarci", fallback_to_heuristic=True),
            IndexGenerationNode(
                group_definitions={
                    'pHLA': 'chainID {mhc_alpha} {b2m} {peptide}',
                    'TCR': 'chainID {tcr_alpha} {tcr_beta}',
                    'MHC_alpha': 'chainID {mhc_alpha}',
                    'Peptide': 'chainID {peptide}'
                },
                output_name="comprehensive.ndx"
            ),
            RMSDNode(selection="protein and name CA", method="mdanalysis"),
            # TODO: Add DockingAngleNode
            # TODO: Add ContactFrequencyNode
        ]
        super().__init__(nodes=nodes)


class CDRRMSFPipeline(Pipeline):
    """
    CDR RMSF analysis pipeline.

    This pipeline combines ANARCI-based CDR detection with RMSF calculation
    to analyze flexibility of CDR regions (CDR1, CDR2, CDR3).

    This demonstrates the modular architecture:
    - Layer 1: CDRManager.detect_all_cdr_regions() + gmx rmsf
    - Layer 2: CDRDetectionNode + RMSFNode
    - Layer 3: CDRRMSFPipeline (this class)
    - Layer 4: CLI command (imn cdr-rmsf) + BatchExecutor

    Example:
        >>> from immunex.core import PipelineContext
        >>> from immunex.pipeline import CDRRMSFPipeline
        >>>
        >>> context = PipelineContext(
        ...     system_id="1ao7",
        ...     topology="data/1ao7/md.tpr",
        ...     trajectory_processed="data/1ao7/md_pbc.xtc"
        ... )
        >>>
        >>> pipeline = CDRRMSFPipeline(
        ...     chains={'TCR_alpha': 'chainID D', 'TCR_beta': 'chainID E'},
        ...     include_cdrs=[1, 2, 3]
        ... )
        >>> result = pipeline.execute(context)
        >>>
        >>> print(f"Detected {result.results['cdr_detection']['n_cdrs_detected']} CDR regions")
        >>> print(f"RMSF calculated for {result.results['cdr_rmsf']['n_cdrs']} regions")
    """

    def __init__(self,
                 chains: dict = None,
                 include_cdrs: list = None,
                 allow_fallback: bool = True,
                 numbering_scheme: str = "imgt",
                 residue_averaging: bool = True):
        """
        Initialize CDR RMSF pipeline.

        Args:
            chains: Dictionary mapping chain names to selection strings
                   e.g., {'TCR_alpha': 'chainID D', 'TCR_beta': 'chainID E'}
                   If None, uses default TCR chains (D, E)
            include_cdrs: List of CDR numbers to include (1, 2, 3)
                         Default: [1, 2, 3]
            allow_fallback: Allow regex fallback if ANARCI unavailable
            numbering_scheme: ANARCI numbering scheme (imgt/kabat/chothia)
            residue_averaging: Average RMSF by residue (-res flag)
        """
        if include_cdrs is None:
            include_cdrs = [1, 2, 3]

        nodes = [
            CDRDetectionNode(
                chains=chains,
                include_cdrs=include_cdrs,
                allow_fallback=allow_fallback,
                numbering_scheme=numbering_scheme
            ),
            RMSFNode(
                use_cdr_regions=True,
                residue_averaging=residue_averaging
            )
        ]
        super().__init__(nodes=nodes)
        self.include_cdrs = include_cdrs
