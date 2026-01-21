"""
Contact Number Calculator Module

Lightweight tool for calculating contact number time series between two selections.
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from typing import Optional, Tuple
import logging
import os

logger = logging.getLogger(__name__)


class ContactNumberCalculator:
    """
    Calculate contact number time series between two atom selections.

    This is a lightweight tool designed for FEL analysis and other CV-based
    studies. It computes the number of atom pairs within a distance cutoff
    across the trajectory.

    Attributes:
        universe (MDAnalysis.Universe): Trajectory data
        topology (str): Path to topology file
        trajectory (str): Path to trajectory file

    Example:
        >>> calc = ContactNumberCalculator("md.tpr", "md.xtc")
        >>> calc.calculate_contact_number(
        ...     selection1="segname PROE and resid 103:115",  # CDR3beta
        ...     selection2="segname PROC",                     # peptide
        ...     cutoff=4.5,
        ...     output_file="cdr3beta_peptide_contacts.csv"
        ... )
    """

    def __init__(self, topology: str, trajectory: str):
        """
        Initialize ContactNumberCalculator.

        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)

        Raises:
            ValueError: If files cannot be loaded
        """
        self.topology = topology
        self.trajectory = trajectory

        try:
            self.universe = mda.Universe(topology, trajectory)
            logger.info(
                f"Contact Number Calculator initialized: "
                f"{len(self.universe.trajectory)} frames, "
                f"{self.universe.atoms.n_atoms} atoms"
            )
        except Exception as e:
            logger.error(f"Failed to load trajectory: {e}")
            raise ValueError(f"Cannot initialize calculator: {e}")

    def calculate_contact_number(
        self,
        selection1: str,
        selection2: str,
        cutoff: float = 4.5,
        heavy_atoms_only: bool = True,
        output_file: Optional[str] = None,
        time_unit: str = 'ps',
        use_pbc: bool = False
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate contact number time series between two selections.

        A contact is defined as any atom pair (one from each selection) within
        the cutoff distance. The total number of such contacts is computed for
        each frame.

        Args:
            selection1: MDAnalysis selection string for first group
            selection2: MDAnalysis selection string for second group
            cutoff: Distance cutoff for contacts (Angstroms)
            heavy_atoms_only: Only consider non-hydrogen atoms
            output_file: Path to save CSV output (optional)
            time_unit: Time unit in output ('ps' or 'ns')
            use_pbc: Use periodic boundary conditions for distance calculation.
                     Set to False (default) for PBC-processed trajectories.
                     Set to True only for raw trajectories with periodic boundaries.

        Returns:
            Tuple of (times, contact_numbers):
                - times: Array of time values
                - contact_numbers: Array of contact counts per frame

        Raises:
            ValueError: If selections are invalid or return no atoms
        """
        logger.info(f"Calculating contact number between selections:")
        logger.info(f"  Selection 1: {selection1}")
        logger.info(f"  Selection 2: {selection2}")
        logger.info(f"  Cutoff: {cutoff} A")

        # Apply selections
        try:
            if heavy_atoms_only:
                sel1 = f"({selection1}) and not name H*"
                sel2 = f"({selection2}) and not name H*"
            else:
                sel1 = selection1
                sel2 = selection2

            group1 = self.universe.select_atoms(sel1)
            group2 = self.universe.select_atoms(sel2)
        except Exception as e:
            logger.error(f"Selection failed: {e}")
            raise ValueError(f"Invalid selection: {e}")

        if len(group1) == 0:
            raise ValueError(f"Selection 1 returned no atoms: {selection1}")
        if len(group2) == 0:
            raise ValueError(f"Selection 2 returned no atoms: {selection2}")

        logger.info(f"  Group 1: {len(group1)} atoms")
        logger.info(f"  Group 2: {len(group2)} atoms")
        logger.info(f"  Maximum possible contacts: {len(group1) * len(group2)}")

        # Calculate contact numbers across trajectory
        n_frames = len(self.universe.trajectory)
        times = np.zeros(n_frames)
        contact_numbers = np.zeros(n_frames, dtype=int)

        logger.info(f"Processing {n_frames} frames...")

        for frame_idx, ts in enumerate(self.universe.trajectory):
            # Progress logging
            if (frame_idx + 1) % 100 == 0 or frame_idx == 0:
                logger.info(f"  Frame {frame_idx + 1}/{n_frames}")

            # Calculate distance matrix between groups
            # For PBC-processed trajectories, use box=None
            # For raw trajectories with periodic boundaries, use box=ts.dimensions
            box = ts.dimensions if use_pbc else None
            dist_matrix = distance_array(
                group1.positions,
                group2.positions,
                box=box
            )

            # Count contacts within cutoff
            n_contacts = np.sum(dist_matrix < cutoff)
            contact_numbers[frame_idx] = n_contacts

            # Store time
            if time_unit == 'ns':
                times[frame_idx] = ts.time / 1000.0
            else:  # ps
                times[frame_idx] = ts.time

        # Statistics
        mean_contacts = np.mean(contact_numbers)
        std_contacts = np.std(contact_numbers)
        min_contacts = np.min(contact_numbers)
        max_contacts = np.max(contact_numbers)

        logger.info(f"\nContact number statistics:")
        logger.info(f"  Mean: {mean_contacts:.2f} ± {std_contacts:.2f}")
        logger.info(f"  Range: [{min_contacts}, {max_contacts}]")

        # Save to file if requested
        if output_file:
            self._save_results(
                times, contact_numbers,
                output_file, time_unit,
                selection1, selection2, cutoff
            )

        return times, contact_numbers

    def calculate_contact_fraction(
        self,
        selection1: str,
        selection2: str,
        cutoff: float = 4.5,
        heavy_atoms_only: bool = True,
        output_file: Optional[str] = None,
        time_unit: str = 'ps'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate contact fraction (normalized contact number) time series.

        Contact fraction = contact_number / max_possible_contacts
        where max_possible_contacts = len(group1) * len(group2)

        This provides a normalized metric between 0 and 1.

        Args:
            selection1: MDAnalysis selection string for first group
            selection2: MDAnalysis selection string for second group
            cutoff: Distance cutoff for contacts (Angstroms)
            heavy_atoms_only: Only consider non-hydrogen atoms
            output_file: Path to save CSV output (optional)
            time_unit: Time unit in output ('ps' or 'ns')

        Returns:
            Tuple of (times, contact_fractions)
        """
        # Get contact numbers
        times, contact_numbers = self.calculate_contact_number(
            selection1, selection2, cutoff,
            heavy_atoms_only, output_file=None, time_unit=time_unit
        )

        # Calculate max possible contacts
        if heavy_atoms_only:
            sel1 = f"({selection1}) and not name H*"
            sel2 = f"({selection2}) and not name H*"
        else:
            sel1 = selection1
            sel2 = selection2

        group1 = self.universe.select_atoms(sel1)
        group2 = self.universe.select_atoms(sel2)
        max_contacts = len(group1) * len(group2)

        # Normalize
        contact_fractions = contact_numbers / max_contacts

        logger.info(f"Contact fraction range: [{contact_fractions.min():.3f}, "
                   f"{contact_fractions.max():.3f}]")

        # Save if requested
        if output_file:
            self._save_results(
                times, contact_fractions,
                output_file, time_unit,
                selection1, selection2, cutoff,
                metric='contact_fraction'
            )

        return times, contact_fractions

    def calculate_contact_frequency(
        self,
        selection1: str,
        selection2: str,
        cutoff: float = 4.5,
        heavy_atoms_only: bool = True,
        output_file: Optional[str] = None,
        time_unit: str = 'ps'
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Calculate contact frequency: fraction of frames with at least one contact.

        Contact frequency = (number of frames with contact) / (total frames)
        A frame is considered to have contact if ANY atom pair from the two
        selections is within the cutoff distance.

        Args:
            selection1: MDAnalysis selection string for first group
            selection2: MDAnalysis selection string for second group
            cutoff: Distance cutoff for contacts (Angstroms)
            heavy_atoms_only: Only consider non-hydrogen atoms
            output_file: Path to save CSV output (optional)
            time_unit: Time unit in output ('ps' or 'ns')

        Returns:
            Tuple of (times, contact_status, frequency):
                - times: Array of time values
                - contact_status: Binary array (1=contact, 0=no contact) per frame
                - frequency: Overall contact frequency (0.0 to 1.0)
        """
        logger.info(f"Calculating contact frequency between selections:")
        logger.info(f"  Selection 1: {selection1}")
        logger.info(f"  Selection 2: {selection2}")
        logger.info(f"  Cutoff: {cutoff} A")

        # Apply selections
        try:
            if heavy_atoms_only:
                sel1 = f"({selection1}) and not name H*"
                sel2 = f"({selection2}) and not name H*"
            else:
                sel1 = selection1
                sel2 = selection2

            group1 = self.universe.select_atoms(sel1)
            group2 = self.universe.select_atoms(sel2)
        except Exception as e:
            logger.error(f"Selection failed: {e}")
            raise ValueError(f"Invalid selection: {e}")

        if len(group1) == 0:
            raise ValueError(f"Selection 1 returned no atoms: {selection1}")
        if len(group2) == 0:
            raise ValueError(f"Selection 2 returned no atoms: {selection2}")

        logger.info(f"  Group 1: {len(group1)} atoms")
        logger.info(f"  Group 2: {len(group2)} atoms")

        # Calculate contact status across trajectory
        n_frames = len(self.universe.trajectory)
        times = np.zeros(n_frames)
        contact_status = np.zeros(n_frames, dtype=int)

        logger.info(f"Processing {n_frames} frames...")

        for frame_idx, ts in enumerate(self.universe.trajectory):
            if (frame_idx + 1) % 100 == 0 or frame_idx == 0:
                logger.info(f"  Frame {frame_idx + 1}/{n_frames}")

            # Calculate distance matrix between groups
            # For PBC-processed trajectories, use box=None
            # For raw trajectories with periodic boundaries, use box=ts.dimensions
            box = ts.dimensions if use_pbc else None
            dist_matrix = distance_array(
                group1.positions,
                group2.positions,
                box=box
            )

            # Check if ANY contact exists (at least one pair within cutoff)
            has_contact = np.any(dist_matrix < cutoff)
            contact_status[frame_idx] = 1 if has_contact else 0

            # Store time
            if time_unit == 'ns':
                times[frame_idx] = ts.time / 1000.0
            else:
                times[frame_idx] = ts.time

        # Calculate overall frequency
        frequency = np.sum(contact_status) / n_frames

        logger.info(f"\nContact frequency results:")
        logger.info(f"  Frames with contact: {np.sum(contact_status)}/{n_frames}")
        logger.info(f"  Contact frequency: {frequency:.4f} ({frequency*100:.2f}%)")

        # Save to file if requested
        if output_file:
            self._save_contact_frequency_results(
                times, contact_status, frequency,
                output_file, time_unit,
                selection1, selection2, cutoff
            )

        return times, contact_status, frequency

    def _save_contact_frequency_results(
        self,
        times: np.ndarray,
        contact_status: np.ndarray,
        frequency: float,
        output_file: str,
        time_unit: str,
        selection1: str,
        selection2: str,
        cutoff: float
    ):
        """Save contact frequency results to CSV file."""
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # Create DataFrame
        df = pd.DataFrame({
            f'Time_{time_unit}': times,
            'Contact_Status': contact_status
        })

        # Save CSV
        df.to_csv(output_file, index=False)
        logger.info(f"\nResults saved to: {output_file}")

        # Save summary metadata
        metadata_file = output_file.replace('.csv', '_summary.txt')
        with open(metadata_file, 'w') as f:
            f.write(f"Contact Frequency Analysis Summary\n")
            f.write(f"===================================\n\n")
            f.write(f"Topology: {self.topology}\n")
            f.write(f"Trajectory: {self.trajectory}\n\n")
            f.write(f"Selection 1: {selection1}\n")
            f.write(f"Selection 2: {selection2}\n")
            f.write(f"Cutoff: {cutoff} A\n\n")
            f.write(f"Total frames: {len(times)}\n")
            f.write(f"Frames with contact: {np.sum(contact_status)}\n")
            f.write(f"Frames without contact: {len(times) - np.sum(contact_status)}\n\n")
            f.write(f"Contact Frequency: {frequency:.4f} ({frequency*100:.2f}%)\n\n")
            f.write(f"Time range: {times[0]:.2f} - {times[-1]:.2f} {time_unit}\n")

        logger.info(f"Summary saved to: {metadata_file}")

    def _save_results(
        self,
        times: np.ndarray,
        values: np.ndarray,
        output_file: str,
        time_unit: str,
        selection1: str,
        selection2: str,
        cutoff: float,
        metric: str = 'contact_number'
    ):
        """Save results to CSV file."""
        # Create output directory if needed
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # Create DataFrame
        df = pd.DataFrame({
            f'time_{time_unit}': times,
            metric: values
        })

        # Save CSV
        df.to_csv(output_file, index=False)
        logger.info(f"\nResults saved to: {output_file}")

        # Also save metadata as comment in header
        metadata_file = output_file.replace('.csv', '_metadata.txt')
        with open(metadata_file, 'w') as f:
            f.write(f"Contact Number Calculation Metadata\n")
            f.write(f"=====================================\n")
            f.write(f"Topology: {self.topology}\n")
            f.write(f"Trajectory: {self.trajectory}\n")
            f.write(f"Selection 1: {selection1}\n")
            f.write(f"Selection 2: {selection2}\n")
            f.write(f"Cutoff: {cutoff} A\n")
            f.write(f"Metric: {metric}\n")
            f.write(f"Frames: {len(times)}\n")
            f.write(f"Time range: {times[0]:.2f} - {times[-1]:.2f} {time_unit}\n")
            f.write(f"Mean value: {np.mean(values):.2f}\n")
            f.write(f"Std value: {np.std(values):.2f}\n")

        logger.info(f"Metadata saved to: {metadata_file}")
