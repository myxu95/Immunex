import numpy as np
import pandas as pd
import subprocess
from typing import Optional, Tuple
import logging
from ...utils.group_selector import GroupSelector

# Optional MDAnalysis import
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    mda = None
    rms = None
    align = None

logger = logging.getLogger(__name__)


class RMSDCalculator:
    """RMSD calculation module for MD trajectories."""
    
    def __init__(self, topology: str, trajectory: str, gmx_executable: str = "gmx"):
        """
        Initialize RMSDCalculator.

        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)
            gmx_executable: GROMACS executable command
        """
        self.topology = topology
        self.trajectory = trajectory
        self.gmx = gmx_executable

        # Initialize group selector
        self.group_selector = GroupSelector(topology, gmx_executable)

        # Load MDAnalysis universe if available
        self.universe = None
        if HAS_MDANALYSIS:
            try:
                self.universe = mda.Universe(topology, trajectory)
                logger.info(f"RMSD Calculator initialized with {len(self.universe.trajectory)} frames (MDAnalysis available)")
            except Exception as e:
                logger.warning(f"Failed to load trajectory with MDAnalysis: {e}")
                logger.info("MDAnalysis methods will not be available, but GROMACS methods will work")
        else:
            logger.info("RMSD Calculator initialized (MDAnalysis not available, GROMACS methods only)")
    
    def calculate_mdanalysis(self,
                           selection: str = "protein and name CA",
                           reference_frame: int = 0,
                           output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate RMSD using MDAnalysis.

        Args:
            selection: Atom selection for RMSD calculation
            reference_frame: Reference frame index
            output_file: Optional output file for RMSD data

        Returns:
            Tuple of (time, rmsd) arrays

        Raises:
            RuntimeError: If MDAnalysis is not available
        """
        if not HAS_MDANALYSIS or self.universe is None:
            raise RuntimeError("MDAnalysis is not available or Universe failed to load. Use calculate_gromacs() instead.")

        atoms = self.universe.select_atoms(selection)

        # Set reference frame and store reference coordinates
        self.universe.trajectory[reference_frame]
        reference_coords = atoms.positions.copy()

        rmsd_values = []
        times = []

        # Calculate RMSD for each frame
        for ts in self.universe.trajectory:
            # Calculate RMSD with optimal alignment (without modifying coordinates)
            rmsd_val = rms.rmsd(atoms.positions, reference_coords, superposition=True)
            rmsd_values.append(rmsd_val)
            times.append(ts.time)

        times = np.array(times)
        rmsd_values = np.array(rmsd_values)

        if output_file:
            self._save_data(times, rmsd_values, output_file)

        logger.info(f"MDAnalysis RMSD completed. Mean RMSD: {np.mean(rmsd_values):.3f} nm")
        return times, rmsd_values
    
    def calculate_gromacs(self,
                         selection_group: Optional[str] = None,
                         output_file: str = "rmsd.xvg",
                         rmsd_type: str = "backbone") -> str:
        """
        Calculate RMSD using GROMACS gmx rms command with intelligent group selection.
        
        Args:
            selection_group: Specific group name/number (if None, auto-select based on rmsd_type)
            output_file: Output XVG file path
            rmsd_type: Type of RMSD calculation ('backbone', 'protein', 'calpha', 'heavy')
            
        Returns:
            Path to output XVG file
        """
        # Generate input string for fit and calculation groups
        if selection_group is None:
            # Use intelligent two-step group selection
            stdin_input = self.group_selector.get_rmsd_input_string(rmsd_type)
            fit_group, calc_group = self.group_selector.get_rmsd_groups(rmsd_type)
            logger.info(f"Auto-selected RMSD groups - Fit: {fit_group}, Calc: {calc_group}")
        else:
            # Use the same group for both fit and calculation
            stdin_input = f"{selection_group}\n{selection_group}\n"
            logger.info(f"Using custom group for both fit and calc: {selection_group}")
        
        cmd = [
            "rms",
            "-f", self.trajectory,
            "-s", self.topology,
            "-o", output_file
        ]
        
        try:
            result = subprocess.run(
                [self.gmx] + cmd,
                input=stdin_input,
                text=True,
                capture_output=True,
                check=True
            )
            logger.info(f"GROMACS RMSD calculation completed: {output_file}")
            return output_file
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS RMSD calculation failed: {e}")
            logger.error(f"stderr: {e.stderr}")
            raise
    
    def calculate_gromacs_custom_groups(self,
                                       reference_group: str,
                                       analysis_group: str, 
                                       output_file: str = "rmsd.xvg") -> str:
        """
        Calculate RMSD using GROMACS with custom reference and analysis groups.
        
        Args:
            reference_group: Reference group for fitting
            analysis_group: Analysis group for RMSD calculation
            output_file: Output XVG file path
            
        Returns:
            Path to output XVG file
        """
        stdin_input = f"{reference_group}\n{analysis_group}\n"
        
        cmd = [
            "rms",
            "-f", self.trajectory,
            "-s", self.topology,
            "-o", output_file
        ]
        
        try:
            result = subprocess.run(
                [self.gmx] + cmd,
                input=stdin_input,
                text=True,
                capture_output=True,
                check=True
            )
            logger.info(f"GROMACS RMSD calculation completed with custom groups: {output_file}")
            return output_file
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS RMSD calculation failed: {e}")
            logger.error(f"stderr: {e.stderr}")
            raise
    
    def _save_data(self, times: np.ndarray, rmsd_values: np.ndarray, filename: str):
        """Save RMSD data to file."""
        if filename.endswith('.csv'):
            df = pd.DataFrame({'Time (ps)': times, 'RMSD (nm)': rmsd_values})
            df.to_csv(filename, index=False)
        elif filename.endswith('.xvg'):
            self._write_xvg_file(filename, times, rmsd_values)
        else:
            np.savetxt(filename, np.column_stack([times, rmsd_values]), 
                      header="Time (ps)\tRMSD (nm)")
    
    def _write_xvg_file(self, filename: str, times: np.ndarray, rmsd_values: np.ndarray):
        """Write data to XVG format file."""
        with open(filename, 'w') as f:
            f.write("# RMSD calculation by AfterMD\n")
            f.write("@ title \"RMSD vs Time\"\n")
            f.write("@ xaxis  label \"Time (ps)\"\n")
            f.write("@ yaxis  label \"RMSD (nm)\"\n")
            f.write("@ TYPE xy\n")
            
            for t, rmsd in zip(times, rmsd_values):
                f.write(f"{t:12.6f} {rmsd:12.6f}\n")