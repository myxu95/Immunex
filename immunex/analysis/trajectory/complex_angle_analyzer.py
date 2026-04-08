#!/usr/bin/env python3
"""
Complex Angle Analyzer for pHLA-TCR Complexes (Fixed Version)

Comprehensive docking geometry analysis for TCR-pMHC MD trajectories:
- Crossing Angle: TCR long axis vs MHC groove axis
- Tilt Angle: TCR projection on MHC surface vs groove axis (VMD definition)
- Docking Angle: TCR rotation around peptide axis

**VERSION HISTORY**:
- v2.0.0 (2025-12-22): Added frame-to-frame direction alignment to eliminate
  180° sign flips caused by SVD direction ambiguity
- v2.1.0 (2026-01-15): CRITICAL FIX - Corrected incident_angle calculation
  to implement VMD's tilt angle definition (TCR projection vs groove axis)
  instead of incorrect TCR vs normal angle

Based on Rudolph & Wilson coordinate system framework.

Author: Immunex Development Team
Version: 2.1.0
"""

import MDAnalysis as mda
import numpy as np
import pandas as pd
import json
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple
from datetime import datetime
import warnings

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ComplexAngleAnalyzer:
    """
    Comprehensive TCR-pMHC docking angle analyzer with direction alignment.

    Calculates geometric parameters describing TCR orientation on MHC surface:
    1. Crossing Angle: TCR Vα/Vβ axis vs MHC α1/α2 helix axis (0-90°)
    2. Tilt Angle: TCR projection on MHC surface vs groove axis (0-90°, VMD definition)
    3. Docking Angle: TCR rotation around peptide axis (-180° to +180°)

    **KEY IMPROVEMENTS**:
    - v2.0.0: Frame-to-frame direction alignment for all axes
              Eliminates 180° sign flips from SVD direction ambiguity
              Ensures continuous angle trajectories
    - v2.1.0: CRITICAL FIX - Corrected tilt angle calculation
              Now uses VMD's projection method (TCR proj vs groove axis)
              Previous version incorrectly used TCR vs normal (3D angle)

    Chain Structure (pHLA-TCR standard):
    - Chain A: HLA-α (MHC heavy chain)
    - Chain B: HLA-β (β2-microglobulin)
    - Chain C: Peptide (antigen)
    - Chain D: TCR-α
    - Chain E: TCR-β
    """

    def __init__(self, task_name: str, trajectory: str, topology: str,
                 output_dir: Optional[str] = None):
        """
        Initialize complex angle analyzer.

        Args:
            task_name: Task identifier (e.g., "1ao7_run1")
            trajectory: Path to XTC trajectory file
            topology: Path to PDB/TPR topology file
            output_dir: Output directory (default: ./angle_results/{task_name})
        """
        self.task_name = task_name
        self.trajectory = Path(trajectory)
        self.topology = Path(topology)

        if output_dir is None:
            self.output_dir = Path("./angle_results") / task_name
        else:
            self.output_dir = Path(output_dir)

        # Create subdirectories
        self.subdirs = {
            'angles': self.output_dir / 'angles',
            'figures': self.output_dir / 'figures'
        }

        for subdir in self.subdirs.values():
            subdir.mkdir(parents=True, exist_ok=True)

        # Load universe
        logger.info(f"Loading trajectory: {self.trajectory}")
        logger.info(f"Loading topology: {self.topology}")

        self.u = mda.Universe(str(self.topology), str(self.trajectory))

        logger.info(f"Universe loaded successfully")
        logger.info(f"  Atoms: {len(self.u.atoms)}")
        logger.info(f"  Frames: {len(self.u.trajectory)}")
        logger.info(f"  Time range: {self.u.trajectory[0].time:.2f} - {self.u.trajectory[-1].time:.2f} ps")

        # Define chain selections
        self.chain_selections = {
            'HLA_alpha': 'chainID A and protein',
            'HLA_beta': 'chainID B and protein',
            'Peptide': 'chainID C and protein',
            'TCR_alpha': 'chainID D and protein',
            'TCR_beta': 'chainID E and protein'
        }

        # Validate chains
        self._validate_chains()

        # Store results
        self.results = {}

    def _validate_chains(self):
        """Validate that all expected chains are present."""
        logger.info("Validating chain structure...")

        for chain_name, selection in self.chain_selections.items():
            atoms = self.u.select_atoms(selection)
            if len(atoms) == 0:
                logger.warning(f"Chain not found: {chain_name} ({selection})")
            else:
                logger.info(f"  {chain_name}: {len(atoms)} atoms")

    def _extract_mhc_axis(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract MHC binding groove axis and surface normal using SVD.

        Method: Singular Value Decomposition on MHC α1 and α2 helix Cα atoms.
        - Largest eigenvalue direction → groove long axis
        - Smallest eigenvalue direction → surface normal (perpendicular to helices)

        Returns:
            Tuple of (mhc_axis, mhc_normal) as normalized vectors

        Reference:
            Rudolph & Wilson coordinate system (guideline.md Section 3)
        """
        # Select MHC α1 and α2 helix residues
        # Typical ranges for HLA-A: α1 (1-90), α2 (140-180)
        # Using broad selection to cover different HLA alleles
        mhc_helices = self.u.select_atoms(
            "(chainID A and protein and name CA and (resid 1:90 or resid 140:180))"
        )

        if len(mhc_helices) == 0:
            raise ValueError("No MHC helix Cα atoms found. Check chain A residue numbering.")

        # Get positions (already centered for current frame)
        coords = mhc_helices.positions

        # Center coordinates
        coords_centered = coords - coords.mean(axis=0)

        # SVD decomposition
        U, S, Vt = np.linalg.svd(coords_centered, full_matrices=False)

        # Extract principal axes
        mhc_axis = Vt[0]  # Largest variance = groove long axis
        mhc_normal = Vt[2]  # Smallest variance = surface normal

        # Normalize (should already be normalized from SVD, but ensure)
        mhc_axis = mhc_axis / np.linalg.norm(mhc_axis)
        mhc_normal = mhc_normal / np.linalg.norm(mhc_normal)

        return mhc_axis, mhc_normal

    def _extract_tcr_axis(self) -> np.ndarray:
        """
        Extract TCR variable domain axis.

        Method 1 (Primary): Vector connecting Vα and Vβ disulfide bond centers
        Method 2 (Fallback): SVD on combined Vα/Vβ Cα atoms

        Returns:
            Normalized TCR axis vector (pointing from Vα to Vβ)

        Notes:
            - Conserved disulfide bonds typically in CDR2 region
            - Vα: around residue 23 (IMGT numbering ~CDR1/FR2)
            - Vβ: around residue 23-57 (IMGT CDR2 region)
        """
        # Attempt Method 1: Disulfide bond centers
        try:
            # Select all Cys Cα atoms in TCR variable domains
            valpha_cys = self.u.select_atoms("chainID D and resname CYS and name CA")
            vbeta_cys = self.u.select_atoms("chainID E and resname CYS and name CA")

            if len(valpha_cys) >= 1 and len(vbeta_cys) >= 1:
                # Use center of geometry of all Cys (simple and robust)
                center_valpha = valpha_cys.center_of_geometry()
                center_vbeta = vbeta_cys.center_of_geometry()

                tcr_axis = center_vbeta - center_valpha
                tcr_axis = tcr_axis / np.linalg.norm(tcr_axis)

                return tcr_axis

        except Exception as e:
            logger.warning(f"Disulfide method failed: {e}, falling back to SVD method")

        # Fallback Method 2: SVD on Vα/Vβ Cα atoms
        logger.info("Using SVD method for TCR axis extraction")

        # Select all Vα and Vβ Cα atoms (variable domains only)
        # Assuming variable domains are first ~110 residues of each chain
        tcr_v_atoms = self.u.select_atoms(
            "((chainID D and name CA) or (chainID E and name CA)) and protein"
        )

        if len(tcr_v_atoms) == 0:
            raise ValueError("No TCR Cα atoms found")

        coords = tcr_v_atoms.positions
        coords_centered = coords - coords.mean(axis=0)

        U, S, Vt = np.linalg.svd(coords_centered, full_matrices=False)
        tcr_axis = Vt[0]  # Principal component

        tcr_axis = tcr_axis / np.linalg.norm(tcr_axis)

        return tcr_axis

    def _extract_peptide_axis(self) -> np.ndarray:
        """
        Extract peptide N→C terminal axis.

        Returns:
            Normalized vector from N-terminus to C-terminus
        """
        peptide = self.u.select_atoms("chainID C and protein and name CA")

        if len(peptide) < 2:
            raise ValueError("Peptide has fewer than 2 residues")

        # First and last Cα atoms
        n_term_pos = peptide[0].position
        c_term_pos = peptide[-1].position

        peptide_axis = c_term_pos - n_term_pos
        peptide_axis = peptide_axis / np.linalg.norm(peptide_axis)

        return peptide_axis

    def calculate_crossing_angle(self) -> pd.DataFrame:
        """
        Calculate TCR-MHC crossing angle over trajectory.

        Crossing Angle: Acute angle between TCR long axis and MHC groove axis.
        - Uses frame-to-frame direction alignment
        - Ensures 0-90° range (axial angle)
        - Typical canonical TCR: 30-70°

        **FIX**: Added direction alignment to prevent 180° flips from SVD ambiguity

        Returns:
            DataFrame with columns: time_ns, angle_deg

        Outputs:
            CSV: angles/crossing_angle.csv
        """
        logger.info("=" * 60)
        logger.info("Calculating Crossing Angle (TCR vs MHC axis)...")
        logger.info("=" * 60)

        times = []
        angles = []

        # Store previous frame directions for alignment
        prev_mhc_axis = None
        prev_tcr_axis = None

        for ts in self.u.trajectory:
            try:
                # Extract axes for current frame
                mhc_axis, _ = self._extract_mhc_axis()
                tcr_axis = self._extract_tcr_axis()

                # ✨ NEW: Direction alignment with previous frame
                if prev_mhc_axis is not None:
                    # Align MHC axis
                    if np.dot(mhc_axis, prev_mhc_axis) < 0:
                        mhc_axis = -mhc_axis
                    # Align TCR axis
                    if np.dot(tcr_axis, prev_tcr_axis) < 0:
                        tcr_axis = -tcr_axis

                # Save current directions for next frame
                prev_mhc_axis = mhc_axis.copy()
                prev_tcr_axis = tcr_axis.copy()

                # Calculate angle
                cos_angle = np.dot(tcr_axis, mhc_axis)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)

                angle_rad = np.arccos(cos_angle)
                angle_deg = np.degrees(angle_rad)

                # Ensure 0-90° range (axial angle)
                if angle_deg > 90.0:
                    angle_deg = 180.0 - angle_deg

                times.append(ts.time / 1000.0)
                angles.append(angle_deg)

            except Exception as e:
                logger.warning(f"Frame {ts.frame} failed: {e}")
                times.append(ts.time / 1000.0)
                angles.append(np.nan)

        # Create DataFrame
        df = pd.DataFrame({
            'time_ns': times,
            'angle_deg': angles
        })

        df = df.dropna()

        if len(df) == 0:
            logger.error("No valid crossing angles calculated!")
            return df

        # Save CSV
        output_file = self.subdirs['angles'] / 'crossing_angle.csv'
        df.to_csv(output_file, index=False)

        logger.info(f"Crossing angle saved: {output_file}")
        logger.info(f"  Valid frames: {len(df)}")
        logger.info(f"  Mean angle: {df['angle_deg'].mean():.2f}°")
        logger.info(f"  Std: {df['angle_deg'].std():.2f}°")
        logger.info(f"  Range: [{df['angle_deg'].min():.2f}, {df['angle_deg'].max():.2f}]°")

        # Store statistics
        self.results['crossing_angle'] = {
            'mean': float(df['angle_deg'].mean()),
            'std': float(df['angle_deg'].std()),
            'min': float(df['angle_deg'].min()),
            'max': float(df['angle_deg'].max()),
            'median': float(df['angle_deg'].median())
        }

        return df

    def calculate_incident_angle(self) -> pd.DataFrame:
        """
        Calculate TCR tilt angle (VMD definition - CORRECTED v2.1).

        Tilt Angle: Angle between TCR axis PROJECTION on MHC surface plane
                    and MHC groove axis (2D angle on surface plane).

        **CRITICAL FIX (v2.1.0 - 2026-01-15)**:
        - Previous implementation incorrectly calculated TCR vs MHC normal (3D angle)
        - Corrected to VMD's tilt angle definition: TCR projection vs groove axis
        - Method:
          1. Project TCR axis onto MHC surface plane (remove normal component)
          2. Calculate angle between projection and MHC groove axis

        Physical Interpretation:
        - Low angle (0-20°): TCR nearly aligned with groove axis
        - Medium angle (20-50°): Typical canonical TCR-pMHC docking
        - High angle (50-90°): TCR perpendicular to groove axis

        Typical values: 30-40° for canonical TCR-pMHC complexes

        Reference:
            VMD script: mkvmd_tiltangle.sh (lines 64-70)
            Rudolph & Wilson coordinate system

        **FIX**: Added direction alignment to prevent 180° flips

        Returns:
            DataFrame with columns: time_ns, angle_deg

        Outputs:
            CSV: angles/incident_angle.csv (name kept for backward compatibility)
        """
        logger.info("=" * 60)
        logger.info("Calculating Tilt Angle (VMD definition - CORRECTED)...")
        logger.info("=" * 60)

        times = []
        angles = []

        # Store previous frame directions for alignment
        prev_mhc_axis = None
        prev_mhc_normal = None
        prev_tcr_axis = None

        for ts in self.u.trajectory:
            try:
                # Extract axes
                mhc_axis, mhc_normal = self._extract_mhc_axis()
                tcr_axis = self._extract_tcr_axis()

                # ✨ Direction alignment with previous frame
                if prev_mhc_axis is not None:
                    if np.dot(mhc_axis, prev_mhc_axis) < 0:
                        mhc_axis = -mhc_axis
                    if np.dot(mhc_normal, prev_mhc_normal) < 0:
                        mhc_normal = -mhc_normal
                    if np.dot(tcr_axis, prev_tcr_axis) < 0:
                        tcr_axis = -tcr_axis

                prev_mhc_axis = mhc_axis.copy()
                prev_mhc_normal = mhc_normal.copy()
                prev_tcr_axis = tcr_axis.copy()

                # ========================================
                # CORRECTED VMD TILT ANGLE CALCULATION
                # ========================================
                # Step 1: Project TCR axis onto MHC surface plane
                # (Remove the component parallel to MHC normal)
                tcr_proj_on_plane = tcr_axis - np.dot(tcr_axis, mhc_normal) * mhc_normal
                tcr_proj_norm = np.linalg.norm(tcr_proj_on_plane)

                if tcr_proj_norm < 1e-6:
                    # TCR is perpendicular to MHC surface (very rare)
                    logger.warning(f"Frame {ts.frame}: TCR perpendicular to MHC surface")
                    times.append(ts.time / 1000.0)
                    angles.append(np.nan)
                    continue

                # Normalize the projection
                tcr_proj_normalized = tcr_proj_on_plane / tcr_proj_norm

                # Step 2: Calculate angle between projection and MHC groove axis
                cos_tilt = np.dot(tcr_proj_normalized, mhc_axis)
                cos_tilt = np.clip(cos_tilt, -1.0, 1.0)

                angle_rad = np.arccos(np.abs(cos_tilt))  # Use abs to ensure 0-90°
                angle_deg = np.degrees(angle_rad)

                # Ensure 0-90° range (axial angle)
                if angle_deg > 90.0:
                    angle_deg = 180.0 - angle_deg

                times.append(ts.time / 1000.0)
                angles.append(angle_deg)

            except Exception as e:
                logger.warning(f"Frame {ts.frame} failed: {e}")
                times.append(ts.time / 1000.0)
                angles.append(np.nan)

        # Create DataFrame
        df = pd.DataFrame({
            'time_ns': times,
            'angle_deg': angles
        })

        df = df.dropna()

        if len(df) == 0:
            logger.error("No valid tilt angles calculated!")
            return df

        # Save CSV (filename kept as incident_angle.csv for backward compatibility)
        output_file = self.subdirs['angles'] / 'incident_angle.csv'
        df.to_csv(output_file, index=False)

        logger.info(f"Tilt angle (VMD definition) saved: {output_file}")
        logger.info(f"  Valid frames: {len(df)}")
        logger.info(f"  Mean angle: {df['angle_deg'].mean():.2f}°")
        logger.info(f"  Std: {df['angle_deg'].std():.2f}°")
        logger.info(f"  Range: [{df['angle_deg'].min():.2f}, {df['angle_deg'].max():.2f}]°")

        # Store statistics
        self.results['incident_angle'] = {
            'mean': float(df['angle_deg'].mean()),
            'std': float(df['angle_deg'].std()),
            'min': float(df['angle_deg'].min()),
            'max': float(df['angle_deg'].max()),
            'median': float(df['angle_deg'].median())
        }

        return df

    def calculate_docking_angle(self) -> pd.DataFrame:
        """
        Calculate TCR docking angle (rotation around peptide axis).

        Docking Angle: Signed angle describing TCR rotation relative to peptide axis
        when projected onto MHC surface plane.

        **CRITICAL FIX**: Added direction alignment for mhc_normal to eliminate
        180° sign flips that were causing ~45% of tasks to show spurious angle jumps.

        Returns:
            DataFrame with columns: time_ns, angle_deg

        Outputs:
            CSV: angles/docking_angle.csv
        """
        logger.info("=" * 60)
        logger.info("Calculating Docking Angle (TCR rotation)...")
        logger.info("=" * 60)

        times = []
        angles = []

        # ✨ CRITICAL: Store previous mhc_normal for direction alignment
        prev_mhc_normal = None

        for ts in self.u.trajectory:
            try:
                # Extract axes
                mhc_axis, mhc_normal = self._extract_mhc_axis()
                tcr_axis = self._extract_tcr_axis()
                peptide_axis = self._extract_peptide_axis()

                # ✨ NEW: Align mhc_normal with previous frame
                if prev_mhc_normal is not None:
                    if np.dot(mhc_normal, prev_mhc_normal) < 0:
                        # Flip if direction reversed (SVD ambiguity)
                        mhc_normal = -mhc_normal

                # Save for next iteration
                prev_mhc_normal = mhc_normal.copy()

                # Project TCR axis onto MHC surface plane
                tcr_proj = tcr_axis - np.dot(tcr_axis, mhc_normal) * mhc_normal
                tcr_proj_norm = np.linalg.norm(tcr_proj)

                if tcr_proj_norm < 1e-6:
                    logger.warning(f"Frame {ts.frame}: TCR perpendicular to MHC")
                    times.append(ts.time / 1000.0)
                    angles.append(np.nan)
                    continue

                tcr_proj = tcr_proj / tcr_proj_norm

                # Project peptide axis onto MHC plane
                pep_proj = peptide_axis - np.dot(peptide_axis, mhc_normal) * mhc_normal
                pep_proj_norm = np.linalg.norm(pep_proj)

                if pep_proj_norm < 1e-6:
                    logger.warning(f"Frame {ts.frame}: Peptide perpendicular to MHC")
                    times.append(ts.time / 1000.0)
                    angles.append(np.nan)
                    continue

                pep_proj = pep_proj / pep_proj_norm

                # Calculate signed angle using atan2
                cos_angle = np.dot(tcr_proj, pep_proj)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)

                # Cross product for sign (now stable due to aligned mhc_normal)
                cross = np.cross(pep_proj, tcr_proj)
                sign = np.sign(np.dot(cross, mhc_normal))

                angle_rad = np.arccos(cos_angle)
                angle_deg = np.degrees(angle_rad) * sign

                times.append(ts.time / 1000.0)
                angles.append(angle_deg)

            except Exception as e:
                logger.warning(f"Frame {ts.frame} failed: {e}")
                times.append(ts.time / 1000.0)
                angles.append(np.nan)

        # Create DataFrame
        df = pd.DataFrame({
            'time_ns': times,
            'angle_deg': angles
        })

        df = df.dropna()

        if len(df) == 0:
            logger.error("No valid docking angles calculated!")
            return df

        # Save CSV
        output_file = self.subdirs['angles'] / 'docking_angle.csv'
        df.to_csv(output_file, index=False)

        logger.info(f"Docking angle saved: {output_file}")
        logger.info(f"  Valid frames: {len(df)}")
        logger.info(f"  Mean angle: {df['angle_deg'].mean():.2f}°")
        logger.info(f"  Std: {df['angle_deg'].std():.2f}°")
        logger.info(f"  Range: [{df['angle_deg'].min():.2f}, {df['angle_deg'].max():.2f}]°")

        # Store statistics
        self.results['docking_angle'] = {
            'mean': float(df['angle_deg'].mean()),
            'std': float(df['angle_deg'].std()),
            'min': float(df['angle_deg'].min()),
            'max': float(df['angle_deg'].max()),
            'median': float(df['angle_deg'].median())
        }

        return df

    def _save_summary(self):
        """Save analysis summary to JSON file."""
        summary = {
            'task_name': self.task_name,
            'trajectory': str(self.trajectory),
            'topology': str(self.topology),
            'n_atoms': len(self.u.atoms),
            'n_frames': len(self.u.trajectory),
            'time_range_ps': [float(self.u.trajectory[0].time),
                             float(self.u.trajectory[-1].time)],
            'timestamp': datetime.now().isoformat(),
            'version': '2.1.0 (VMD tilt angle corrected)',
            'results': self.results
        }

        summary_file = self.subdirs['angles'] / 'angle_stats.json'
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        logger.info(f"\nAnalysis summary saved: {summary_file}")

    def run_full_analysis(self) -> Dict:
        """
        Execute complete angle analysis workflow.

        Returns:
            Dictionary with all analysis results
        """
        logger.info("\n" + "=" * 70)
        logger.info(f"COMPLEX ANGLE ANALYSIS: {self.task_name}")
        logger.info("VERSION: 2.1.0 (VMD Tilt Angle Corrected)")
        logger.info("=" * 70)
        logger.info(f"Trajectory: {self.trajectory}")
        logger.info(f"Topology: {self.topology}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info("=" * 70 + "\n")

        start_time = datetime.now()

        try:
            # Phase 1: Crossing angle
            logger.info("\n>>> Phase 1/3: Crossing Angle")
            self.calculate_crossing_angle()

            # Phase 2: Tilt angle (corrected VMD definition)
            logger.info("\n>>> Phase 2/3: Tilt Angle (VMD definition - CORRECTED)")
            self.calculate_incident_angle()

            # Phase 3: Docking angle
            logger.info("\n>>> Phase 3/3: Docking Angle")
            self.calculate_docking_angle()

            # Save summary
            self._save_summary()

            elapsed = (datetime.now() - start_time).total_seconds()

            logger.info("\n" + "=" * 70)
            logger.info("ANALYSIS COMPLETE")
            logger.info("=" * 70)
            logger.info(f"Total time: {elapsed:.2f} seconds")
            logger.info(f"Output directory: {self.output_dir}")
            logger.info("=" * 70 + "\n")

            return self.results

        except Exception as e:
            logger.error(f"\nAnalysis failed: {e}")
            import traceback
            traceback.print_exc()
            raise


def main():
    """Main function for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Complex Angle Analyzer for pHLA-TCR Complexes (Fixed v2.1)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single task
  python complex_angle_analyzer.py \\
      1ao7_run1 \\
      /path/to/1ao7_processed.xtc \\
      /path/to/1ao7_standardized.pdb \\
      ./angle_results/1ao7_run1

  # Analyze with auto output directory
  python complex_angle_analyzer.py \\
      1ao7_run1 \\
      /path/to/trajectory.xtc \\
      /path/to/topology.pdb

Version: 2.1.0 (VMD Tilt Angle Corrected)
Changes:
  - v2.0.0: Fixed 180° sign flips via frame-to-frame direction alignment
  - v2.1.0: CRITICAL - Corrected tilt angle to use VMD projection method
        """
    )

    parser.add_argument('task_name', type=str,
                       help='Task name (e.g., 1ao7_run1)')
    parser.add_argument('trajectory', type=str,
                       help='Path to XTC trajectory file')
    parser.add_argument('topology', type=str,
                       help='Path to PDB/TPR topology file')
    parser.add_argument('output_dir', type=str, nargs='?', default=None,
                       help='Output directory (default: ./angle_results/{task_name})')

    args = parser.parse_args()

    # Create analyzer
    analyzer = ComplexAngleAnalyzer(
        task_name=args.task_name,
        trajectory=args.trajectory,
        topology=args.topology,
        output_dir=args.output_dir
    )

    # Run analysis
    try:
        results = analyzer.run_full_analysis()
        logger.info("SUCCESS!")
        return 0
    except Exception as e:
        logger.error(f"FAILED: {e}")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
