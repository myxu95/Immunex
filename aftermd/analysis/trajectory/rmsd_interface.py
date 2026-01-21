"""
RMSD Analysis Unified Interface

This module provides a high-level interface for RMSD calculations with flexible
component selection for alignment and analysis.

Supported components:
- HLA: MHC alpha and beta chains
- pHLA: MHC + peptide complex
- peptide: Antigenic peptide
- TCR: T-cell receptor (both chains)
- TCR_alpha: TCR alpha chain only
- TCR_beta: TCR beta chain only
- CDR3_alpha: CDR3 region of TCR alpha (requires sequence)
- CDR3_beta: CDR3 region of TCR beta (requires sequence)

Author: AfterMD Development Team
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import numpy as np
import pandas as pd

from ...utils.index_generator import IndexGenerator

logger = logging.getLogger(__name__)


class RMSDInterface:
    """
    Unified interface for RMSD calculations with flexible component selection.

    This class provides a streamlined API for calculating RMSD with different
    alignment and analysis components in pHLA-TCR complexes.

    Example:
        # Calculate TCR RMSD aligned to pHLA
        rmsd = RMSDInterface("md.tpr", "md_processed.xtc")
        result = rmsd.calculate(align="pHLA", calc="TCR")

        # Calculate CDR3 RMSD aligned to TCR
        result = rmsd.calculate(
            align="TCR",
            calc="CDR3_beta",
            cdr3_sequences={'beta': 'CASSLGQAYEQYF'}
        )
    """

    def __init__(self,
                 topology: str,
                 trajectory: str,
                 gmx_executable: str = "gmx",
                 auto_standardize: bool = True,
                 output_dir: Optional[str] = None):
        """
        Initialize RMSD Interface.

        Args:
            topology: Path to topology file (.tpr, .gro, .pdb)
            trajectory: Path to trajectory file (.xtc, .trr)
            gmx_executable: GROMACS executable command
            auto_standardize: Automatically standardize chain ordering
            output_dir: Output directory for results (default: same as trajectory)
        """
        self.topology = Path(topology)
        self.trajectory = Path(trajectory)
        self.gmx = gmx_executable

        if not self.topology.exists():
            raise FileNotFoundError(f"Topology not found: {self.topology}")
        if not self.trajectory.exists():
            raise FileNotFoundError(f"Trajectory not found: {self.trajectory}")

        self.output_dir = Path(output_dir) if output_dir else self.trajectory.parent
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.index_generator = IndexGenerator(
            topology=str(self.topology),
            gmx_executable=self.gmx,
            auto_standardize=auto_standardize
        )

        self.base_index_file = None
        self.cdr3_index_cache = {}
        self.fixed_components = ['HLA', 'HLA_alpha', 'HLA_beta', 'pHLA', 'peptide', 'TCR', 'TCR_alpha', 'TCR_beta']

        logger.info(f"RMSD Interface initialized")
        logger.info(f"  Topology: {self.topology.name}")
        logger.info(f"  Trajectory: {self.trajectory.name}")

        self._generate_base_index()

    def _generate_base_index(self) -> None:
        """
        Generate unified index file with all fixed components.

        Creates a single index file containing:
        - HLA (chains A, B)
        - pHLA (chains A, B, C)
        - peptide (chain C)
        - TCR (chains D, E)
        - TCR_alpha (chain D)
        - TCR_beta (chain E)

        CDR3 components are generated separately as they require sequences.
        """
        self.base_index_file = self.output_dir / "base_components.ndx"

        if self.base_index_file.exists():
            logger.info(f"Using existing base index: {self.base_index_file.name}")
            return

        logger.info("Generating unified base index for fixed components...")

        try:
            index_path = self.index_generator.generate_multi_component_index(
                components=self.fixed_components,
                output_file=str(self.base_index_file)
            )

            logger.info(f"Base index created: {self.base_index_file.name}")
            logger.info(f"  Components: {', '.join(self.fixed_components)}")

        except Exception as e:
            logger.error(f"Failed to generate base index: {e}")
            self.base_index_file = None
            raise

    def calculate(self,
                 align: str,
                 calc: str,
                 output_file: Optional[str] = None,
                 cdr3_sequences: Optional[Dict[str, str]] = None) -> Dict:
        """
        Calculate RMSD with specified alignment and calculation components.

        Args:
            align: Component for alignment (HLA/pHLA/peptide/TCR/etc.)
            calc: Component for RMSD calculation
            output_file: Output XVG file path (default: auto-generated)
            cdr3_sequences: CDR3 sequences dict {'alpha': seq, 'beta': seq}

        Returns:
            Dictionary with RMSD results:
                - 'time': Time array (ps)
                - 'rmsd': RMSD array (nm)
                - 'mean': Mean RMSD
                - 'std': Std RMSD
                - 'min': Min RMSD
                - 'max': Max RMSD
                - 'output_file': Path to XVG file

        Example:
            >>> result = rmsd.calculate(align="pHLA", calc="TCR")
            >>> print(f"Mean RMSD: {result['mean']:.3f} nm")
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"RMSD Calculation: Align={align}, Calc={calc}")
        logger.info(f"{'='*60}")

        align_index = self._get_or_create_index(align, cdr3_sequences)
        calc_index = self._get_or_create_index(calc, cdr3_sequences)

        if output_file is None:
            output_file = self.output_dir / f"rmsd_{align}_to_{calc}.xvg"
        else:
            output_file = Path(output_file)

        combined_index = self._merge_indices(align_index, calc_index)

        align_group, calc_group = self._get_group_names(align, calc)

        logger.info(f"Running gmx rms:")
        logger.info(f"  Align group: {align_group}")
        logger.info(f"  Calc group: {calc_group}")
        logger.info(f"  Output: {output_file.name}")

        self._run_gmx_rms(
            combined_index,
            output_file,
            align_group,
            calc_group
        )

        result = self._parse_rmsd_file(output_file)
        result['output_file'] = str(output_file)
        result['align_component'] = align
        result['calc_component'] = calc

        logger.info(f"\nRMSD Statistics:")
        logger.info(f"  Mean: {result['mean']:.4f} nm")
        logger.info(f"  Std:  {result['std']:.4f} nm")
        logger.info(f"  Min:  {result['min']:.4f} nm")
        logger.info(f"  Max:  {result['max']:.4f} nm")
        logger.info(f"{'='*60}\n")

        return result

    def _get_or_create_index(self,
                            component: str,
                            cdr3_sequences: Optional[Dict[str, str]] = None) -> Path:
        """
        Get index file for the specified component.

        For fixed components (HLA, pHLA, peptide, TCR, etc.):
            Returns the unified base index file.

        For CDR3 components:
            Generates separate index file (requires sequence).

        Args:
            component: Component name
            cdr3_sequences: CDR3 sequences if needed

        Returns:
            Path to index file containing the component
        """
        if component in self.fixed_components:
            logger.info(f"Using base index for {component}")
            return self.base_index_file

        if component.startswith("CDR3"):
            if component in self.cdr3_index_cache:
                logger.info(f"Using cached CDR3 index for {component}")
                return self.cdr3_index_cache[component]

            chain = 'alpha' if 'alpha' in component else 'beta'
            if not cdr3_sequences or chain not in cdr3_sequences:
                raise ValueError(f"CDR3 sequence for {chain} chain is required")

            output_file = self.output_dir / f"cdr3_{chain}_index.ndx"

            logger.info(f"Generating CDR3 index for {chain} chain...")

            index_path = self.index_generator.generate_cdr3_index(
                chain=chain,
                cdr3_sequence=cdr3_sequences[chain],
                output_file=str(output_file),
                ca_only=True
            )

            if index_path is None:
                raise RuntimeError(f"Failed to generate CDR3 index for {chain}")

            index_path = Path(index_path)
            self.cdr3_index_cache[component] = index_path

            logger.info(f"CDR3 index created: {index_path.name}")

            return index_path

        else:
            raise ValueError(f"Unknown component: {component}")

    def _merge_indices(self, index1: Path, index2: Path) -> Path:
        """
        Merge two index files into one if necessary.

        If both components are in the base index, no merge is needed.
        If one is a CDR3 index, merge it with the base index.

        Args:
            index1: First index file
            index2: Second index file

        Returns:
            Path to index file containing both components
        """
        if index1 == index2:
            return index1

        if index1 == self.base_index_file and index2 == self.base_index_file:
            return self.base_index_file

        merged_file = self.output_dir / "combined_rmsd_index.ndx"

        logger.info(f"Merging index files for gmx rms...")

        with open(merged_file, 'w') as outf:
            with open(index1, 'r') as inf:
                outf.write(inf.read())
                outf.write('\n')
            with open(index2, 'r') as inf:
                outf.write(inf.read())

        return merged_file

    def _get_group_names(self, align: str, calc: str) -> Tuple[str, str]:
        """
        Get group names from component names.

        Args:
            align: Alignment component
            calc: Calculation component

        Returns:
            Tuple of (align_group_name, calc_group_name)
        """
        def get_group_name(component: str) -> str:
            if component.startswith("CDR3"):
                chain = 'alpha' if 'alpha' in component else 'beta'
                return f"CDR3_TCR_{chain}_CA"
            else:
                return component

        return get_group_name(align), get_group_name(calc)

    def _run_gmx_rms(self,
                    index_file: Path,
                    output_file: Path,
                    align_group: str,
                    calc_group: str) -> None:
        """
        Run gmx rms command.

        Args:
            index_file: Index file path
            output_file: Output XVG file path
            align_group: Group name for alignment
            calc_group: Group name for RMSD calculation
        """
        cmd = [
            self.gmx, 'rms',
            '-s', str(self.topology),
            '-f', str(self.trajectory),
            '-n', str(index_file),
            '-o', str(output_file),
            '-tu', 'ns'
        ]

        stdin_input = f"{align_group}\n{calc_group}\n"

        try:
            process = subprocess.run(
                cmd,
                input=stdin_input,
                text=True,
                capture_output=True,
                timeout=600,
                check=True
            )

            logger.debug(f"gmx rms completed successfully")

        except subprocess.CalledProcessError as e:
            logger.error(f"gmx rms failed: {e.stderr}")
            raise RuntimeError(f"RMSD calculation failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("RMSD calculation timed out")

    def _parse_rmsd_file(self, rmsd_file: Path) -> Dict:
        """
        Parse RMSD XVG file.

        Args:
            rmsd_file: Path to XVG file

        Returns:
            Dictionary with time, rmsd, and statistics
        """
        data = []

        with open(rmsd_file, 'r') as f:
            for line in f:
                if line.startswith(('#', '@')):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    time = float(parts[0])
                    rmsd = float(parts[1])
                    data.append([time, rmsd])

        if not data:
            raise ValueError(f"No data found in {rmsd_file}")

        data = np.array(data)
        time_data = data[:, 0]
        rmsd_data = data[:, 1]

        return {
            'time': time_data,
            'rmsd': rmsd_data,
            'mean': float(np.mean(rmsd_data)),
            'std': float(np.std(rmsd_data)),
            'min': float(np.min(rmsd_data)),
            'max': float(np.max(rmsd_data)),
            'n_frames': len(rmsd_data)
        }

    def batch_calculate(self,
                       calculations: List[Dict],
                       output_csv: Optional[str] = None) -> pd.DataFrame:
        """
        Perform multiple RMSD calculations in batch.

        Args:
            calculations: List of calculation configs, each with:
                - 'align': Alignment component
                - 'calc': Calculation component
                - 'output_file': Optional output file name
                - 'cdr3_sequences': Optional CDR3 sequences
            output_csv: Optional CSV file to save summary

        Returns:
            DataFrame with RMSD statistics for all calculations

        Example:
            >>> calcs = [
            ...     {'align': 'pHLA', 'calc': 'TCR'},
            ...     {'align': 'pHLA', 'calc': 'peptide'},
            ...     {'align': 'TCR', 'calc': 'CDR3_beta',
            ...      'cdr3_sequences': {'beta': 'CASSLGQAYEQYF'}}
            ... ]
            >>> df = rmsd.batch_calculate(calcs)
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"Batch RMSD Calculation: {len(calculations)} tasks")
        logger.info(f"{'='*60}\n")

        results = []

        for i, calc_config in enumerate(calculations, 1):
            logger.info(f"Task {i}/{len(calculations)}: "
                       f"{calc_config['align']} -> {calc_config['calc']}")

            try:
                result = self.calculate(
                    align=calc_config['align'],
                    calc=calc_config['calc'],
                    output_file=calc_config.get('output_file'),
                    cdr3_sequences=calc_config.get('cdr3_sequences')
                )

                results.append({
                    'align_component': result['align_component'],
                    'calc_component': result['calc_component'],
                    'rmsd_mean_nm': result['mean'],
                    'rmsd_std_nm': result['std'],
                    'rmsd_min_nm': result['min'],
                    'rmsd_max_nm': result['max'],
                    'n_frames': result['n_frames'],
                    'output_file': result['output_file']
                })

            except Exception as e:
                logger.error(f"Failed: {e}")
                results.append({
                    'align_component': calc_config['align'],
                    'calc_component': calc_config['calc'],
                    'error': str(e)
                })

        df = pd.DataFrame(results)

        if output_csv:
            df.to_csv(output_csv, index=False)
            logger.info(f"\nBatch results saved: {output_csv}")

        logger.info(f"\n{'='*60}")
        logger.info(f"Batch Complete: {len(results)} results")
        logger.info(f"{'='*60}\n")

        return df

    @staticmethod
    def list_available_components() -> List[str]:
        """
        List all available component names.

        Returns:
            List of component names
        """
        return [
            'HLA', 'pHLA', 'peptide',
            'TCR', 'TCR_alpha', 'TCR_beta',
            'CDR3_alpha', 'CDR3_beta'
        ]
