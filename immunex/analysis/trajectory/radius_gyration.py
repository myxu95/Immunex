import numpy as np
import pandas as pd
import MDAnalysis as mda
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class RadiusGyrationCalculator:
    """Radius of gyration calculation module."""
    
    def __init__(self, topology: str, trajectory: str):
        """
        Initialize RadiusGyrationCalculator.
        
        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)
        """
        self.topology = topology
        self.trajectory = trajectory
        
        try:
            self.universe = mda.Universe(topology, trajectory)
            logger.info(f"Radius of Gyration Calculator initialized with {len(self.universe.trajectory)} frames")
        except Exception as e:
            logger.error(f"Failed to load trajectory for Rg: {e}")
            raise
    
    def calculate(self,
                 selection: str = "protein",
                 output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate radius of gyration over trajectory.
        
        Args:
            selection: Atom selection
            output_file: Optional output file
            
        Returns:
            Tuple of (time, radius_of_gyration) arrays
        """
        atoms = self.universe.select_atoms(selection)
        logger.info(f"Calculating Rg for {len(atoms)} atoms with selection: {selection}")
        
        times = []
        rg_values = []
        
        for ts in self.universe.trajectory:
            rg = atoms.radius_of_gyration()
            rg_values.append(rg)
            times.append(ts.time)
        
        times = np.array(times)
        rg_values = np.array(rg_values)
        
        if output_file:
            self._save_data(times, rg_values, output_file)
        
        logger.info(f"Radius of gyration calculation completed. Mean Rg: {np.mean(rg_values):.3f} nm")
        return times, rg_values
    
    def calculate_components(self,
                           selection: str = "protein",
                           output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate radius of gyration components (Rgx, Rgy, Rgz).
        
        Args:
            selection: Atom selection
            output_file: Optional output file
            
        Returns:
            Tuple of (time, Rgx, Rgy, Rgz) arrays
        """
        atoms = self.universe.select_atoms(selection)
        
        times = []
        rgx_values = []
        rgy_values = []
        rgz_values = []
        
        for ts in self.universe.trajectory:
            positions = atoms.positions
            com = atoms.center_of_mass()
            rel_pos = positions - com
            
            # Calculate component-wise radius of gyration
            rgx = np.sqrt(np.mean(rel_pos[:, 0]**2))
            rgy = np.sqrt(np.mean(rel_pos[:, 1]**2))
            rgz = np.sqrt(np.mean(rel_pos[:, 2]**2))
            
            rgx_values.append(rgx)
            rgy_values.append(rgy)
            rgz_values.append(rgz)
            times.append(ts.time)
        
        times = np.array(times)
        rgx_values = np.array(rgx_values)
        rgy_values = np.array(rgy_values)
        rgz_values = np.array(rgz_values)
        
        if output_file:
            self._save_components_data(times, rgx_values, rgy_values, rgz_values, output_file)
        
        logger.info(f"Radius of gyration components calculation completed")
        return times, rgx_values, rgy_values, rgz_values
    
    def _save_data(self, times: np.ndarray, rg_values: np.ndarray, filename: str):
        """Save Rg data to file."""
        if filename.endswith('.csv'):
            df = pd.DataFrame({'Time (ps)': times, 'Rg (nm)': rg_values})
            df.to_csv(filename, index=False)
        elif filename.endswith('.xvg'):
            self._write_xvg_file(filename, times, rg_values)
        else:
            np.savetxt(filename, np.column_stack([times, rg_values]), 
                      header="Time (ps)\tRg (nm)")
    
    def _save_components_data(self, times: np.ndarray, rgx: np.ndarray, 
                            rgy: np.ndarray, rgz: np.ndarray, filename: str):
        """Save Rg components data to file."""
        if filename.endswith('.csv'):
            df = pd.DataFrame({
                'Time (ps)': times, 
                'Rgx (nm)': rgx, 
                'Rgy (nm)': rgy, 
                'Rgz (nm)': rgz
            })
            df.to_csv(filename, index=False)
        else:
            data = np.column_stack([times, rgx, rgy, rgz])
            np.savetxt(filename, data, header="Time (ps)\tRgx (nm)\tRgy (nm)\tRgz (nm)")
    
    def _write_xvg_file(self, filename: str, times: np.ndarray, rg_values: np.ndarray):
        """Write Rg data to XVG format file."""
        with open(filename, 'w') as f:
            f.write("# Radius of gyration calculation by Immunex\n")
            f.write("@ title \"Radius of Gyration vs Time\"\n")
            f.write("@ xaxis  label \"Time (ps)\"\n")
            f.write("@ yaxis  label \"Rg (nm)\"\n")
            f.write("@ TYPE xy\n")
            
            for t, rg in zip(times, rg_values):
                f.write(f"{t:12.6f} {rg:12.6f}\n")