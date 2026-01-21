"""
Buried SASA Calculator Module

Lightweight tool for calculating buried solvent-accessible surface area (SASA)
at protein-protein interfaces.
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from typing import Optional, Tuple
import logging
import os

logger = logging.getLogger(__name__)

# Try to import FreeSASA module (replacement for MDAnalysis SASA)
try:
    import freesasa
    SASA_AVAILABLE = True
except ImportError:
    SASA_AVAILABLE = False
    logger.warning("FreeSASA module not available. Install with: pip install freesasa")


class BuriedSASACalculator:
    """
    Calculate buried SASA time series at protein-protein interfaces.

    Buried SASA quantifies the surface area that becomes inaccessible to solvent
    upon complex formation:

        Buried SASA = (SASA_A + SASA_B - SASA_AB) / 2

    where:
        - SASA_A: Surface area of molecule A in isolation
        - SASA_B: Surface area of molecule B in isolation
        - SASA_AB: Surface area of the complex

    This is a lightweight tool designed for FEL analysis and other CV-based studies.

    Attributes:
        universe (MDAnalysis.Universe): Trajectory data
        topology (str): Path to topology file
        trajectory (str): Path to trajectory file

    Example:
        >>> calc = BuriedSASACalculator("md.tpr", "md.xtc")
        >>> calc.calculate_buried_sasa(
        ...     selection_a="segname PROD PROE",      # TCR
        ...     selection_b="segname PROA PROB PROC", # pMHC
        ...     output_file="interface_buried_sasa.csv"
        ... )

    Note:
        Requires FreeSASA library. Install with:
        pip install freesasa
    """

    def __init__(self, topology: str, trajectory: str):
        """
        Initialize BuriedSASACalculator.

        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)

        Raises:
            ValueError: If files cannot be loaded or SASA module unavailable
        """
        if not SASA_AVAILABLE:
            raise ValueError(
                "FreeSASA module not available. "
                "Install with: pip install freesasa"
            )

        self.topology = topology
        self.trajectory = trajectory

        try:
            self.universe = mda.Universe(topology, trajectory)
            logger.info(
                f"Buried SASA Calculator initialized: "
                f"{len(self.universe.trajectory)} frames, "
                f"{self.universe.atoms.n_atoms} atoms"
            )
        except Exception as e:
            logger.error(f"Failed to load trajectory: {e}")
            raise ValueError(f"Cannot initialize calculator: {e}")

    def calculate_buried_sasa(
        self,
        selection_a: str,
        selection_b: str,
        probe_radius: float = 1.4,
        output_file: Optional[str] = None,
        time_unit: str = 'ps',
        stride: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate buried SASA time series between two selections.

        Args:
            selection_a: MDAnalysis selection for molecule A
            selection_b: MDAnalysis selection for molecule B
            probe_radius: Radius of solvent probe in Angstroms (default: 1.4 for water)
            output_file: Path to save CSV output (optional)
            time_unit: Time unit in output ('ps' or 'ns')
            stride: Use every nth frame (default: 1, use all frames)

        Returns:
            Tuple of (times, buried_sasa_values):
                - times: Array of time values
                - buried_sasa_values: Array of buried SASA (Angstrom^2) per frame

        Raises:
            ValueError: If selections are invalid or return no atoms
        """
        logger.info(f"Calculating buried SASA between selections:")
        logger.info(f"  Selection A: {selection_a}")
        logger.info(f"  Selection B: {selection_b}")
        logger.info(f"  Probe radius: {probe_radius} A")
        logger.info(f"  Stride: {stride}")

        # Validate selections
        try:
            atoms_a = self.universe.select_atoms(selection_a)
            atoms_b = self.universe.select_atoms(selection_b)
            atoms_ab = self.universe.select_atoms(f"({selection_a}) or ({selection_b})")
        except Exception as e:
            logger.error(f"Selection failed: {e}")
            raise ValueError(f"Invalid selection: {e}")

        if len(atoms_a) == 0:
            raise ValueError(f"Selection A returned no atoms: {selection_a}")
        if len(atoms_b) == 0:
            raise ValueError(f"Selection B returned no atoms: {selection_b}")

        logger.info(f"  Molecule A: {len(atoms_a)} atoms")
        logger.info(f"  Molecule B: {len(atoms_b)} atoms")
        logger.info(f"  Complex: {len(atoms_ab)} atoms")

        # Calculate SASA for each component
        logger.info("\n[1/3] Calculating SASA for molecule A...")
        sasa_a = self._calculate_sasa(atoms_a, probe_radius, stride)

        logger.info("[2/3] Calculating SASA for molecule B...")
        sasa_b = self._calculate_sasa(atoms_b, probe_radius, stride)

        logger.info("[3/3] Calculating SASA for complex AB...")
        sasa_ab = self._calculate_sasa(atoms_ab, probe_radius, stride)

        # Calculate buried SASA
        # Formula: Buried SASA = (SASA_A + SASA_B - SASA_AB) / 2
        buried_sasa = (sasa_a + sasa_b - sasa_ab) / 2.0

        # Extract times
        n_frames = len(self.universe.trajectory[::stride])
        times = np.zeros(n_frames)
        for i, ts in enumerate(self.universe.trajectory[::stride]):
            if time_unit == 'ns':
                times[i] = ts.time / 1000.0
            else:  # ps
                times[i] = ts.time

        # Statistics
        mean_bsa = np.mean(buried_sasa)
        std_bsa = np.std(buried_sasa)
        min_bsa = np.min(buried_sasa)
        max_bsa = np.max(buried_sasa)

        logger.info(f"\nBuried SASA statistics:")
        logger.info(f"  Mean: {mean_bsa:.2f} ± {std_bsa:.2f} A^2")
        logger.info(f"  Range: [{min_bsa:.2f}, {max_bsa:.2f}] A^2")

        # Save to file if requested
        if output_file:
            self._save_results(
                times, buried_sasa,
                output_file, time_unit,
                selection_a, selection_b, probe_radius,
                sasa_a, sasa_b, sasa_ab
            )

        return times, buried_sasa

    def _calculate_sasa(
        self,
        atoms: mda.AtomGroup,
        probe_radius: float,
        stride: int
    ) -> np.ndarray:
        """
        Calculate SASA for an atom group across trajectory using FreeSASA.

        Args:
            atoms: MDAnalysis AtomGroup
            probe_radius: Probe radius in Angstroms
            stride: Frame stride

        Returns:
            Array of total SASA values per frame
        """
        import tempfile
        import os

        n_frames = len(self.universe.trajectory[::stride])
        sasa_values = np.zeros(n_frames)

        # Configure FreeSASA parameters
        parameters = freesasa.Parameters()
        parameters.setProbeRadius(probe_radius)

        # Process each frame
        for i, ts in enumerate(self.universe.trajectory[::stride]):
            if (i + 1) % 100 == 0 or i == 0:
                logger.info(f"    Frame {i + 1}/{n_frames}")

            # Write temporary PDB file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
                tmp_pdb_path = tmp_pdb.name
                tmp_pdb.write(self._atoms_to_pdb_string(atoms))

            # Calculate SASA using FreeSASA
            try:
                structure = freesasa.Structure(tmp_pdb_path)
                result = freesasa.calc(structure, parameters)
                sasa_values[i] = result.totalArea()
            except Exception as e:
                logger.warning(f"    Frame {i}: SASA calculation failed: {e}")
                sasa_values[i] = 0.0
            finally:
                # Clean up temporary file
                if os.path.exists(tmp_pdb_path):
                    os.unlink(tmp_pdb_path)

        return sasa_values

    def _atoms_to_pdb_string(self, atoms: mda.AtomGroup) -> str:
        """
        Convert MDAnalysis AtomGroup to PDB format string.

        Args:
            atoms: MDAnalysis AtomGroup

        Returns:
            PDB format string
        """
        pdb_lines = []

        for i, atom in enumerate(atoms, 1):
            # PDB ATOM record format
            # ATOM serial name altLoc resName chainID resSeq iCode x y z occupancy tempFactor element charge
            line = (
                f"ATOM  {i:5d} {atom.name:4s} "
                f"{atom.resname:3s} {atom.segid[0] if atom.segid else 'A':1s}"
                f"{atom.resid:4d}    "
                f"{atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}"
                f"  1.00  0.00           {atom.element if hasattr(atom, 'element') else atom.name[0]:2s}\n"
            )
            pdb_lines.append(line)

        pdb_lines.append("END\n")

        return ''.join(pdb_lines)

    def calculate_interface_ratio(
        self,
        selection_a: str,
        selection_b: str,
        probe_radius: float = 1.4,
        output_file: Optional[str] = None,
        time_unit: str = 'ps',
        stride: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate interface ratio (fraction of surface involved in binding).

        Interface ratio = Buried SASA / (SASA_A + SASA_B) * 2

        This provides a normalized metric between 0 and 1.

        Args:
            selection_a: MDAnalysis selection for molecule A
            selection_b: MDAnalysis selection for molecule B
            probe_radius: Radius of solvent probe in Angstroms
            output_file: Path to save CSV output (optional)
            time_unit: Time unit in output ('ps' or 'ns')
            stride: Use every nth frame

        Returns:
            Tuple of (times, interface_ratios)
        """
        # Get buried SASA components
        logger.info("Calculating interface ratio...")

        atoms_a = self.universe.select_atoms(selection_a)
        atoms_b = self.universe.select_atoms(selection_b)
        atoms_ab = self.universe.select_atoms(f"({selection_a}) or ({selection_b})")

        sasa_a = self._calculate_sasa(atoms_a, probe_radius, stride)
        sasa_b = self._calculate_sasa(atoms_b, probe_radius, stride)
        sasa_ab = self._calculate_sasa(atoms_ab, probe_radius, stride)

        buried_sasa = (sasa_a + sasa_b - sasa_ab) / 2.0

        # Calculate ratio
        interface_ratio = (2.0 * buried_sasa) / (sasa_a + sasa_b)

        # Extract times
        n_frames = len(self.universe.trajectory[::stride])
        times = np.zeros(n_frames)
        for i, ts in enumerate(self.universe.trajectory[::stride]):
            if time_unit == 'ns':
                times[i] = ts.time / 1000.0
            else:
                times[i] = ts.time

        logger.info(f"Interface ratio range: [{interface_ratio.min():.3f}, "
                   f"{interface_ratio.max():.3f}]")

        # Save if requested
        if output_file:
            self._save_results(
                times, interface_ratio,
                output_file, time_unit,
                selection_a, selection_b, probe_radius,
                sasa_a, sasa_b, sasa_ab,
                metric='interface_ratio'
            )

        return times, interface_ratio

    def _save_results(
        self,
        times: np.ndarray,
        values: np.ndarray,
        output_file: str,
        time_unit: str,
        selection_a: str,
        selection_b: str,
        probe_radius: float,
        sasa_a: np.ndarray,
        sasa_b: np.ndarray,
        sasa_ab: np.ndarray,
        metric: str = 'buried_sasa'
    ):
        """Save results to CSV file."""
        # Create output directory if needed
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # Create DataFrame with detailed columns
        df = pd.DataFrame({
            f'time_{time_unit}': times,
            metric: values,
            'sasa_a': sasa_a,
            'sasa_b': sasa_b,
            'sasa_ab': sasa_ab
        })

        # Save CSV
        df.to_csv(output_file, index=False)
        logger.info(f"\nResults saved to: {output_file}")

        # Save metadata
        metadata_file = output_file.replace('.csv', '_metadata.txt')
        with open(metadata_file, 'w') as f:
            f.write(f"Buried SASA Calculation Metadata\n")
            f.write(f"==================================\n")
            f.write(f"Topology: {self.topology}\n")
            f.write(f"Trajectory: {self.trajectory}\n")
            f.write(f"Selection A: {selection_a}\n")
            f.write(f"Selection B: {selection_b}\n")
            f.write(f"Probe radius: {probe_radius} A\n")
            f.write(f"Metric: {metric}\n")
            f.write(f"Frames: {len(times)}\n")
            f.write(f"Time range: {times[0]:.2f} - {times[-1]:.2f} {time_unit}\n")
            f.write(f"\nMean SASA_A: {np.mean(sasa_a):.2f} A^2\n")
            f.write(f"Mean SASA_B: {np.mean(sasa_b):.2f} A^2\n")
            f.write(f"Mean SASA_AB: {np.mean(sasa_ab):.2f} A^2\n")
            f.write(f"Mean {metric}: {np.mean(values):.2f}\n")
            f.write(f"Std {metric}: {np.std(values):.2f}\n")

        logger.info(f"Metadata saved to: {metadata_file}")
