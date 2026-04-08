import numpy as np
import pandas as pd
import MDAnalysis as mda
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class DistanceCalculator:
    """Distance calculation module for MD trajectories."""
    
    def __init__(self, topology: str, trajectory: str):
        """
        Initialize DistanceCalculator.
        
        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)
        """
        self.topology = topology
        self.trajectory = trajectory
        
        try:
            self.universe = mda.Universe(topology, trajectory)
            logger.info(f"Distance Calculator initialized with {len(self.universe.trajectory)} frames")
        except Exception as e:
            logger.error(f"Failed to load trajectory for distance calculation: {e}")
            raise
    
    def calculate_center_of_mass_distance(self,
                                        selection1: str,
                                        selection2: str,
                                        output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate distance between centers of mass of two selections over time.
        
        Args:
            selection1: First atom selection
            selection2: Second atom selection
            output_file: Optional output file
            
        Returns:
            Tuple of (time, distance) arrays
        """
        group1 = self.universe.select_atoms(selection1)
        group2 = self.universe.select_atoms(selection2)
        
        logger.info(f"Calculating COM distance between {len(group1)} and {len(group2)} atoms")
        
        times = []
        distances = []
        
        for ts in self.universe.trajectory:
            com1 = group1.center_of_mass()
            com2 = group2.center_of_mass()
            distance = np.linalg.norm(com1 - com2)
            distances.append(distance)
            times.append(ts.time)
        
        times = np.array(times)
        distances = np.array(distances)
        
        if output_file:
            self._save_data(times, distances, output_file, "COM Distance")
        
        logger.info(f"COM distance calculation completed. Mean distance: {np.mean(distances):.3f} Å")
        return times, distances
    
    def calculate_minimum_distance(self,
                                 selection1: str,
                                 selection2: str,
                                 output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate minimum distance between two selections over time.
        
        Args:
            selection1: First atom selection
            selection2: Second atom selection
            output_file: Optional output file
            
        Returns:
            Tuple of (time, minimum_distance) arrays
        """
        group1 = self.universe.select_atoms(selection1)
        group2 = self.universe.select_atoms(selection2)
        
        times = []
        min_distances = []
        
        for ts in self.universe.trajectory:
            pos1 = group1.positions
            pos2 = group2.positions
            
            # Calculate all pairwise distances
            diff = pos1[:, np.newaxis, :] - pos2[np.newaxis, :, :]
            distances = np.sqrt(np.sum(diff**2, axis=2))
            min_dist = np.min(distances)
            
            min_distances.append(min_dist)
            times.append(ts.time)
        
        times = np.array(times)
        min_distances = np.array(min_distances)
        
        if output_file:
            self._save_data(times, min_distances, output_file, "Minimum Distance")
        
        logger.info(f"Minimum distance calculation completed. Mean min distance: {np.mean(min_distances):.3f} Å")
        return times, min_distances
    
    def calculate_atom_distance(self,
                              atom1_selection: str,
                              atom2_selection: str,
                              output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate distance between two specific atoms over time.
        
        Args:
            atom1_selection: First atom selection (should select single atom)
            atom2_selection: Second atom selection (should select single atom)
            output_file: Optional output file
            
        Returns:
            Tuple of (time, distance) arrays
        """
        atom1 = self.universe.select_atoms(atom1_selection)
        atom2 = self.universe.select_atoms(atom2_selection)
        
        if len(atom1) != 1 or len(atom2) != 1:
            raise ValueError("Atom selections must each select exactly one atom")
        
        times = []
        distances = []
        
        for ts in self.universe.trajectory:
            pos1 = atom1.positions[0]
            pos2 = atom2.positions[0]
            distance = np.linalg.norm(pos1 - pos2)
            distances.append(distance)
            times.append(ts.time)
        
        times = np.array(times)
        distances = np.array(distances)
        
        if output_file:
            self._save_data(times, distances, output_file, "Atom Distance")
        
        logger.info(f"Atom distance calculation completed")
        return times, distances
    
    def _save_data(self, times: np.ndarray, distances: np.ndarray, filename: str, distance_type: str):
        """Save distance data to file."""
        if filename.endswith('.csv'):
            df = pd.DataFrame({'Time (ps)': times, f'{distance_type} (Angstrom)': distances})
            df.to_csv(filename, index=False)
        elif filename.endswith('.xvg'):
            self._write_xvg_file(filename, times, distances, distance_type)
        else:
            np.savetxt(filename, np.column_stack([times, distances]), 
                      header=f"Time (ps)\t{distance_type} (Angstrom)")
    
    def _write_xvg_file(self, filename: str, times: np.ndarray, distances: np.ndarray, distance_type: str):
        """Write distance data to XVG format file."""
        with open(filename, 'w') as f:
            f.write(f"# {distance_type} calculation by Immunex\n")
            f.write(f"@ title \"{distance_type} vs Time\"\n")
            f.write("@ xaxis  label \"Time (ps)\"\n")
            f.write("@ yaxis  label \"Distance (Angstrom)\"\n")
            f.write("@ TYPE xy\n")
            
            for t, d in zip(times, distances):
                f.write(f"{t:12.6f} {d:12.6f}\n")