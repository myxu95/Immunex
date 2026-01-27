import numpy as np
import MDAnalysis as mda
from typing import Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class PrincipalAxesCalculator:
    """
    Principal axes calculator based on PCA.

    Used to dynamically construct molecular reference coordinate systems,
    eliminating overall rotation and translation in MD trajectories.

    The algorithm:
    1. Calculate center of mass (COM)
    2. Center coordinates: centered = positions - COM
    3. Construct inertia tensor (3x3 symmetric matrix)
    4. Eigenvalue decomposition: eigenvals, eigenvecs = np.linalg.eigh(I)
    5. Sort by moments of inertia: principal axis = eigenvector with largest eigenvalue
    """

    def __init__(self):
        """Initialize calculator."""
        pass

    def calculate_inertia_tensor(
        self,
        positions: np.ndarray,
        masses: np.ndarray
    ) -> np.ndarray:
        """
        Calculate inertia tensor.

        The inertia tensor is defined as:
        I_xx = sum(m * (y^2 + z^2))
        I_yy = sum(m * (x^2 + z^2))
        I_zz = sum(m * (x^2 + y^2))
        I_xy = -sum(m * x * y)
        I_xz = -sum(m * x * z)
        I_yz = -sum(m * y * z)

        Args:
            positions: Atom coordinates (N, 3)
            masses: Atom masses (N,)

        Returns:
            inertia_tensor: 3x3 symmetric matrix
        """
        com = np.average(positions, weights=masses, axis=0)
        centered = positions - com

        x = centered[:, 0]
        y = centered[:, 1]
        z = centered[:, 2]

        Ixx = np.sum(masses * (y**2 + z**2))
        Iyy = np.sum(masses * (x**2 + z**2))
        Izz = np.sum(masses * (x**2 + y**2))
        Ixy = -np.sum(masses * x * y)
        Ixz = -np.sum(masses * x * z)
        Iyz = -np.sum(masses * y * z)

        I = np.array([
            [Ixx, Ixy, Ixz],
            [Ixy, Iyy, Iyz],
            [Ixz, Iyz, Izz]
        ])

        return I

    def calculate_principal_axes(
        self,
        positions: np.ndarray,
        masses: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate principal axes via PCA.

        Args:
            positions: Atom coordinates (N, 3)
            masses: Atom masses (N,), defaults to unit mass

        Returns:
            axes: Principal axis directions (3, 3), each row is a principal axis
            moments: Principal moments of inertia (3,), in descending order
        """
        if masses is None:
            masses = np.ones(len(positions))

        I = self.calculate_inertia_tensor(positions, masses)

        eigenvals, eigenvecs = np.linalg.eigh(I)

        idx = np.argsort(eigenvals)[::-1]
        moments = eigenvals[idx]
        axes = eigenvecs[:, idx].T

        return axes, moments

    def calculate_principal_axes_from_atoms(
        self,
        atoms: mda.AtomGroup
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate principal axes from MDAnalysis atom group.

        Args:
            atoms: MDAnalysis AtomGroup

        Returns:
            axes: Principal axis directions (3, 3)
            moments: Principal moments of inertia (3,)
        """
        positions = atoms.positions
        masses = atoms.masses
        return self.calculate_principal_axes(positions, masses)

    def calculate_principal_axes_trajectory(
        self,
        atoms: mda.AtomGroup,
        universe: mda.Universe,
        stride: int = 1
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate principal axes for trajectory frames.

        Args:
            atoms: MDAnalysis AtomGroup (selection)
            universe: MDAnalysis Universe
            stride: Frame stride

        Returns:
            times: Time array (n_frames,)
            axes_trajectory: Principal axes for each frame (n_frames, 3, 3)
            moments_trajectory: Principal moments for each frame (n_frames, 3)
        """
        times = []
        axes_list = []
        moments_list = []

        logger.info(f"Calculating principal axes for {len(universe.trajectory)} frames...")

        for ts in universe.trajectory[::stride]:
            axes, moments = self.calculate_principal_axes_from_atoms(atoms)
            times.append(ts.time)
            axes_list.append(axes)
            moments_list.append(moments)

        times = np.array(times)
        axes_trajectory = np.array(axes_list)
        moments_trajectory = np.array(moments_list)

        logger.info(f"Completed: {len(times)} frames processed")

        return times, axes_trajectory, moments_trajectory
