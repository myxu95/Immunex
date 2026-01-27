import numpy as np
import pandas as pd
import MDAnalysis as mda
from typing import Tuple, Optional
import logging

from .principal_axes import PrincipalAxesCalculator
from .vector_angles import (
    angle_between_vectors,
    project_vector_onto_plane,
    signed_angle_between_vectors
)

logger = logging.getLogger(__name__)


class DockingAngleAnalyzer:
    """
    TCR-pMHC docking angle analyzer.

    Calculates Twist, Tilt, and Swing angles based on structural invariants.
    Uses PCA to construct MHC reference coordinate system and disulfide bond
    anchors to define TCR orientation.

    Key features:
    - PCA-based reference frame eliminates overall motion
    - Disulfide bond anchors avoid CDR loop flexibility noise
    - Frame-by-frame quantification of binding conformations

    Angles:
    - Twist: Angle between TCR disulfide bond projection and MHC principal axis
    - Tilt: Angle between TCR V domain axis projection and MHC principal axis
    - Swing: Lateral deviation angle of TCR relative to peptide
    """

    def __init__(self, topology: str, trajectory: Optional[str] = None):
        """
        Initialize analyzer.

        Args:
            topology: Topology file (TPR/PDB/GRO)
            trajectory: Trajectory file (XTC/TRR), optional for single structure
        """
        if trajectory:
            self.universe = mda.Universe(topology, trajectory)
        else:
            self.universe = mda.Universe(topology)

        self.pca_calc = PrincipalAxesCalculator()

        logger.info(f"Loaded structure with {len(self.universe.atoms)} atoms")
        if trajectory:
            logger.info(f"Trajectory: {len(self.universe.trajectory)} frames")

    def _get_tcr_disulfide_vector(
        self,
        tcr_alpha_cys_selection: str,
        tcr_beta_cys_selection: str
    ) -> np.ndarray:
        """
        Calculate TCR disulfide bond vector.

        Uses conserved disulfide bonds in V domains:
        - V_alpha: typically Cys22-Cys92
        - V_beta: typically Cys23-Cys89

        Args:
            tcr_alpha_cys_selection: Selection for alpha chain cysteine pair
            tcr_beta_cys_selection: Selection for beta chain cysteine pair

        Returns:
            tcr_vector: Average disulfide bond vector
        """
        alpha_atoms = self.universe.select_atoms(tcr_alpha_cys_selection)
        beta_atoms = self.universe.select_atoms(tcr_beta_cys_selection)

        if len(alpha_atoms) != 2:
            raise ValueError(
                f"TCR alpha cysteine selection must return 2 atoms, got {len(alpha_atoms)}"
            )
        if len(beta_atoms) != 2:
            raise ValueError(
                f"TCR beta cysteine selection must return 2 atoms, got {len(beta_atoms)}"
            )

        vec_alpha = alpha_atoms.positions[1] - alpha_atoms.positions[0]
        vec_beta = beta_atoms.positions[1] - beta_atoms.positions[0]

        tcr_vector = (vec_alpha + vec_beta) / 2.0

        return tcr_vector

    def calculate_twist_angle(
        self,
        mhc_selection: str,
        tcr_alpha_cys_selection: str,
        tcr_beta_cys_selection: str
    ) -> float:
        """
        Calculate Twist angle.

        Definition: Angle between TCR disulfide bond vector projection on MHC plane
        and MHC principal axis.

        Algorithm:
        1. Calculate MHC principal axis and plane normal
        2. Calculate TCR disulfide bond vector
        3. Project TCR vector onto MHC plane
        4. Calculate angle between projection and MHC axis

        Args:
            mhc_selection: MHC binding groove selection (e.g., alpha1/alpha2 helices)
            tcr_alpha_cys_selection: TCR alpha chain cysteine pair
            tcr_beta_cys_selection: TCR beta chain cysteine pair

        Returns:
            twist_angle: Twist angle in degrees
        """
        mhc_atoms = self.universe.select_atoms(mhc_selection)
        mhc_axes, _ = self.pca_calc.calculate_principal_axes_from_atoms(mhc_atoms)

        mhc_axis = mhc_axes[0]
        mhc_normal = mhc_axes[2]

        tcr_disulfide = self._get_tcr_disulfide_vector(
            tcr_alpha_cys_selection,
            tcr_beta_cys_selection
        )

        tcr_projected = project_vector_onto_plane(tcr_disulfide, mhc_normal)
        tcr_projected = tcr_projected / np.linalg.norm(tcr_projected)

        twist_angle = angle_between_vectors(tcr_projected, mhc_axis)

        return twist_angle

    def calculate_tilt_angle(
        self,
        mhc_selection: str,
        tcr_v_selection: str
    ) -> float:
        """
        Calculate Tilt angle.

        Definition: Angle between TCR V domain principal axis projection on MHC plane
        and MHC principal axis.

        Algorithm:
        1. Calculate MHC coordinate system
        2. Calculate TCR V domain principal axis
        3. Project TCR axis onto MHC plane
        4. Calculate angle between projection and MHC axis

        Args:
            mhc_selection: MHC binding groove selection
            tcr_v_selection: TCR V domain selection (V_alpha + V_beta)

        Returns:
            tilt_angle: Tilt angle in degrees
        """
        mhc_atoms = self.universe.select_atoms(mhc_selection)
        mhc_axes, _ = self.pca_calc.calculate_principal_axes_from_atoms(mhc_atoms)

        mhc_axis = mhc_axes[0]
        mhc_normal = mhc_axes[2]

        tcr_atoms = self.universe.select_atoms(tcr_v_selection)
        tcr_axes, _ = self.pca_calc.calculate_principal_axes_from_atoms(tcr_atoms)

        tcr_axis = tcr_axes[0]

        tcr_projected = project_vector_onto_plane(tcr_axis, mhc_normal)
        tcr_projected = tcr_projected / np.linalg.norm(tcr_projected)

        tilt_angle = angle_between_vectors(tcr_projected, mhc_axis)

        return tilt_angle

    def calculate_swing_angle(
        self,
        mhc_selection: str,
        tcr_v_selection: str,
        peptide_selection: str
    ) -> float:
        """
        Calculate Swing angle.

        Definition: Lateral deviation angle of TCR relative to peptide.

        Algorithm:
        1. Calculate MHC principal axis
        2. Calculate TCR-peptide COM connection vector
        3. Project connection onto plane perpendicular to MHC axis
        4. Calculate lateral deviation angle

        Args:
            mhc_selection: MHC binding groove selection
            tcr_v_selection: TCR V domain selection
            peptide_selection: Peptide selection

        Returns:
            swing_angle: Swing angle in degrees (-180 to 180)
        """
        mhc_atoms = self.universe.select_atoms(mhc_selection)
        mhc_axes, _ = self.pca_calc.calculate_principal_axes_from_atoms(mhc_atoms)

        mhc_axis = mhc_axes[0]

        tcr_atoms = self.universe.select_atoms(tcr_v_selection)
        peptide_atoms = self.universe.select_atoms(peptide_selection)

        tcr_com = tcr_atoms.center_of_mass()
        peptide_com = peptide_atoms.center_of_mass()

        connection_vector = tcr_com - peptide_com

        projected = project_vector_onto_plane(connection_vector, mhc_axis)
        projected = projected / np.linalg.norm(projected)

        reference = mhc_axes[1]
        reference_projected = project_vector_onto_plane(reference, mhc_axis)
        reference_projected = reference_projected / np.linalg.norm(reference_projected)

        swing_angle = signed_angle_between_vectors(
            reference_projected,
            projected,
            mhc_axis
        )

        return swing_angle

    def calculate_docking_angles(
        self,
        mhc_selection: str,
        tcr_alpha_cys_selection: str,
        tcr_beta_cys_selection: str,
        tcr_v_selection: str,
        peptide_selection: str
    ) -> Tuple[float, float, float]:
        """
        Calculate docking angles for current frame.

        Args:
            mhc_selection: MHC binding groove selection
            tcr_alpha_cys_selection: TCR alpha chain cysteine pair
            tcr_beta_cys_selection: TCR beta chain cysteine pair
            tcr_v_selection: TCR V domain selection
            peptide_selection: Peptide selection

        Returns:
            twist, tilt, swing: Three angles in degrees
        """
        twist = self.calculate_twist_angle(
            mhc_selection,
            tcr_alpha_cys_selection,
            tcr_beta_cys_selection
        )

        tilt = self.calculate_tilt_angle(
            mhc_selection,
            tcr_v_selection
        )

        swing = self.calculate_swing_angle(
            mhc_selection,
            tcr_v_selection,
            peptide_selection
        )

        return twist, tilt, swing

    def calculate_docking_angles_trajectory(
        self,
        mhc_selection: str,
        tcr_alpha_cys_selection: str,
        tcr_beta_cys_selection: str,
        tcr_v_selection: str,
        peptide_selection: str,
        stride: int = 1,
        output_file: Optional[str] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate docking angles for trajectory.

        Args:
            mhc_selection: MHC binding groove selection
            tcr_alpha_cys_selection: TCR alpha chain cysteine pair
            tcr_beta_cys_selection: TCR beta chain cysteine pair
            tcr_v_selection: TCR V domain selection
            peptide_selection: Peptide selection
            stride: Frame stride (default: every frame)
            output_file: Output file path (CSV or XVG format)

        Returns:
            times: Time array (ps)
            twist: Twist angle array (degrees)
            tilt: Tilt angle array (degrees)
            swing: Swing angle array (degrees)
        """
        times = []
        twist_angles = []
        tilt_angles = []
        swing_angles = []

        n_frames = len(self.universe.trajectory[::stride])
        logger.info(f"Calculating docking angles for {n_frames} frames...")

        for i, ts in enumerate(self.universe.trajectory[::stride], 1):
            twist, tilt, swing = self.calculate_docking_angles(
                mhc_selection,
                tcr_alpha_cys_selection,
                tcr_beta_cys_selection,
                tcr_v_selection,
                peptide_selection
            )

            times.append(ts.time)
            twist_angles.append(twist)
            tilt_angles.append(tilt)
            swing_angles.append(swing)

            if i % 100 == 0:
                logger.info(f"Processed {i}/{n_frames} frames")

        times = np.array(times)
        twist = np.array(twist_angles)
        tilt = np.array(tilt_angles)
        swing = np.array(swing_angles)

        if output_file:
            self._save_angles(times, twist, tilt, swing, output_file)

        logger.info(f"Completed: {len(times)} frames processed")
        logger.info(f"Twist: {np.mean(twist):.2f} ± {np.std(twist):.2f} degrees")
        logger.info(f"Tilt: {np.mean(tilt):.2f} ± {np.std(tilt):.2f} degrees")
        logger.info(f"Swing: {np.mean(swing):.2f} ± {np.std(swing):.2f} degrees")

        return times, twist, tilt, swing

    def _save_angles(
        self,
        times: np.ndarray,
        twist: np.ndarray,
        tilt: np.ndarray,
        swing: np.ndarray,
        output_file: str
    ):
        """
        Save angle data to file.

        Args:
            times: Time array
            twist: Twist angle array
            tilt: Tilt angle array
            swing: Swing angle array
            output_file: Output file path (.csv or .xvg)
        """
        if output_file.endswith('.csv'):
            df = pd.DataFrame({
                'Time_ps': times,
                'Twist_deg': twist,
                'Tilt_deg': tilt,
                'Swing_deg': swing
            })
            df.to_csv(output_file, index=False)
            logger.info(f"Saved angles to {output_file}")

        elif output_file.endswith('.xvg'):
            with open(output_file, 'w') as f:
                f.write("# TCR-pMHC Docking Angles\n")
                f.write("# Generated by AfterMD DockingAngleAnalyzer\n")
                f.write("@ title \"Docking Angles vs Time\"\n")
                f.write("@ xaxis  label \"Time (ps)\"\n")
                f.write("@ yaxis  label \"Angle (degrees)\"\n")
                f.write("@ TYPE xy\n")
                f.write("@ s0 legend \"Twist\"\n")
                f.write("@ s1 legend \"Tilt\"\n")
                f.write("@ s2 legend \"Swing\"\n")
                for t, tw, ti, sw in zip(times, twist, tilt, swing):
                    f.write(f"{t:12.3f} {tw:8.3f} {ti:8.3f} {sw:8.3f}\n")
            logger.info(f"Saved angles to {output_file}")

        else:
            logger.warning(f"Unknown file format: {output_file}. Supported: .csv, .xvg")
