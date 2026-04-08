import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from typing import Optional, Dict
import logging

logger = logging.getLogger(__name__)


class HydrogenBondAnalyzer:
    """Hydrogen bond analysis module for MD trajectories."""
    
    def __init__(self, topology: str, trajectory: str):
        """
        Initialize HydrogenBondAnalyzer.
        
        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)
        """
        self.topology = topology
        self.trajectory = trajectory
        
        try:
            self.universe = mda.Universe(topology, trajectory)
            logger.info(f"H-bond Analyzer initialized with {len(self.universe.trajectory)} frames")
        except Exception as e:
            logger.error(f"Failed to load trajectory for H-bond analysis: {e}")
            raise
    
    def calculate(self,
                 donors_sel: str = "protein",
                 acceptors_sel: str = "protein", 
                 distance_cutoff: float = 3.0,
                 angle_cutoff: float = 150.0,
                 output_file: Optional[str] = None) -> Dict:
        """
        Calculate hydrogen bonds analysis.
        
        Args:
            donors_sel: Donor atom selection
            acceptors_sel: Acceptor atom selection
            distance_cutoff: Maximum H-A distance in Angstroms
            angle_cutoff: Minimum D-H-A angle in degrees
            output_file: Optional output file for time series
            
        Returns:
            Dictionary with hydrogen bond analysis results
        """
        logger.info(f"Analyzing H-bonds: donors='{donors_sel}', acceptors='{acceptors_sel}'")
        logger.info(f"Cutoffs: distance={distance_cutoff}Å, angle={angle_cutoff}°")
        
        hbonds = HydrogenBondAnalysis(
            universe=self.universe,
            donors_sel=donors_sel,
            acceptors_sel=acceptors_sel,
            d_a_cutoff=distance_cutoff,
            d_h_a_angle_cutoff=angle_cutoff
        )
        
        hbonds.run()
        
        if output_file:
            times = np.array([ts.time for ts in self.universe.trajectory])
            self._save_time_series(times, hbonds.count_by_time, output_file)
        
        results = {
            "count_by_time": hbonds.count_by_time,
            "count_by_ids": hbonds.count_by_ids,
            "hbonds": hbonds.hbonds,
            "mean_count": np.mean(hbonds.count_by_time),
            "std_count": np.std(hbonds.count_by_time),
            "max_count": np.max(hbonds.count_by_time),
            "min_count": np.min(hbonds.count_by_time)
        }
        
        logger.info(f"H-bond analysis completed. Mean count: {results['mean_count']:.2f}")
        return results
    
    def calculate_intermolecular(self,
                               mol1_sel: str,
                               mol2_sel: str,
                               distance_cutoff: float = 3.0,
                               angle_cutoff: float = 150.0,
                               output_file: Optional[str] = None) -> Dict:
        """
        Calculate intermolecular hydrogen bonds between two molecular groups.
        
        Args:
            mol1_sel: First molecule selection
            mol2_sel: Second molecule selection
            distance_cutoff: Maximum H-A distance in Angstroms
            angle_cutoff: Minimum D-H-A angle in degrees
            output_file: Optional output file
            
        Returns:
            Dictionary with intermolecular H-bond results
        """
        # Find donors and acceptors in each molecule
        mol1_donors = f"({mol1_sel}) and (name N* or name O*)"
        mol1_acceptors = f"({mol1_sel}) and (name N* or name O*)"
        mol2_donors = f"({mol2_sel}) and (name N* or name O*)"
        mol2_acceptors = f"({mol2_sel}) and (name N* or name O*)"
        
        # Calculate mol1 donors -> mol2 acceptors
        hbonds_1to2 = HydrogenBondAnalysis(
            universe=self.universe,
            donors_sel=mol1_donors,
            acceptors_sel=mol2_acceptors,
            d_a_cutoff=distance_cutoff,
            d_h_a_angle_cutoff=angle_cutoff
        )
        hbonds_1to2.run()
        
        # Calculate mol2 donors -> mol1 acceptors
        hbonds_2to1 = HydrogenBondAnalysis(
            universe=self.universe,
            donors_sel=mol2_donors,
            acceptors_sel=mol1_acceptors,
            d_a_cutoff=distance_cutoff,
            d_h_a_angle_cutoff=angle_cutoff
        )
        hbonds_2to1.run()
        
        # Combine results
        total_count = hbonds_1to2.count_by_time + hbonds_2to1.count_by_time
        
        if output_file:
            times = np.array([ts.time for ts in self.universe.trajectory])
            self._save_intermolecular_data(times, hbonds_1to2.count_by_time, 
                                         hbonds_2to1.count_by_time, total_count, output_file)
        
        results = {
            "mol1_to_mol2_count": hbonds_1to2.count_by_time,
            "mol2_to_mol1_count": hbonds_2to1.count_by_time,
            "total_count": total_count,
            "mean_total": np.mean(total_count),
            "std_total": np.std(total_count)
        }
        
        logger.info(f"Intermolecular H-bond analysis completed. Mean total: {results['mean_total']:.2f}")
        return results
    
    def _save_time_series(self, times: np.ndarray, counts: np.ndarray, filename: str):
        """Save H-bond time series data."""
        if filename.endswith('.csv'):
            df = pd.DataFrame({'Time (ps)': times, 'H-bond Count': counts})
            df.to_csv(filename, index=False)
        elif filename.endswith('.xvg'):
            self._write_xvg_file(filename, times, counts, "H-bond Count")
        else:
            np.savetxt(filename, np.column_stack([times, counts]), 
                      header="Time (ps)\tH-bond Count")
    
    def _save_intermolecular_data(self, times: np.ndarray, count_1to2: np.ndarray,
                                count_2to1: np.ndarray, total: np.ndarray, filename: str):
        """Save intermolecular H-bond data."""
        if filename.endswith('.csv'):
            df = pd.DataFrame({
                'Time (ps)': times,
                'Mol1_to_Mol2': count_1to2,
                'Mol2_to_Mol1': count_2to1,
                'Total': total
            })
            df.to_csv(filename, index=False)
        else:
            data = np.column_stack([times, count_1to2, count_2to1, total])
            np.savetxt(filename, data, 
                      header="Time (ps)\tMol1_to_Mol2\tMol2_to_Mol1\tTotal")
    
    def _write_xvg_file(self, filename: str, times: np.ndarray, counts: np.ndarray, ylabel: str):
        """Write H-bond data to XVG format file."""
        with open(filename, 'w') as f:
            f.write("# Hydrogen bond analysis by Immunex\n")
            f.write("@ title \"Hydrogen Bonds vs Time\"\n")
            f.write("@ xaxis  label \"Time (ps)\"\n")
            f.write(f"@ yaxis  label \"{ylabel}\"\n")
            f.write("@ TYPE xy\n")
            
            for t, count in zip(times, counts):
                f.write(f"{t:12.6f} {count:12.0f}\n")