import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class RDFCalculator:
    """Radial Distribution Function calculation module."""
    
    def __init__(self, topology: str, trajectory: str):
        """
        Initialize RDFCalculator.
        
        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)
        """
        self.topology = topology
        self.trajectory = trajectory
        
        try:
            self.universe = mda.Universe(topology, trajectory)
            logger.info(f"RDF Calculator initialized with {len(self.universe.trajectory)} frames")
        except Exception as e:
            logger.error(f"Failed to load trajectory for RDF: {e}")
            raise
    
    def calculate(self,
                 selection1: str,
                 selection2: str,
                 nbins: int = 75,
                 range_max: float = 15.0,
                 output_file: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate radial distribution function.
        
        Args:
            selection1: First atom selection
            selection2: Second atom selection
            nbins: Number of bins
            range_max: Maximum distance in Angstroms
            output_file: Optional output file
            
        Returns:
            Tuple of (bins, rdf) arrays
        """
        g1 = self.universe.select_atoms(selection1)
        g2 = self.universe.select_atoms(selection2)
        
        logger.info(f"Calculating RDF between {len(g1)} and {len(g2)} atoms")
        
        rdf_analysis = InterRDF(g1, g2, nbins=nbins, range=(0.0, range_max))
        rdf_analysis.run()
        
        if output_file:
            self._save_data(rdf_analysis.bins, rdf_analysis.rdf, 
                          output_file, selection1, selection2)
        
        logger.info(f"RDF calculation completed for {selection1} - {selection2}")
        return rdf_analysis.bins, rdf_analysis.rdf
    
    def _save_data(self, bins: np.ndarray, rdf: np.ndarray, filename: str,
                  sel1: str, sel2: str):
        """Save RDF data to file."""
        if filename.endswith('.csv'):
            df = pd.DataFrame({'Distance (Angstrom)': bins, 'g(r)': rdf})
            df.to_csv(filename, index=False)
        elif filename.endswith('.xvg'):
            self._write_xvg_file(filename, bins, rdf, sel1, sel2)
        else:
            np.savetxt(filename, np.column_stack([bins, rdf]), 
                      header="Distance (Angstrom)\tg(r)")
    
    def _write_xvg_file(self, filename: str, bins: np.ndarray, rdf: np.ndarray,
                       sel1: str, sel2: str):
        """Write RDF data to XVG format file."""
        with open(filename, 'w') as f:
            f.write(f"# RDF calculation by Immunex: {sel1} - {sel2}\n")
            f.write("@ title \"Radial Distribution Function\"\n")
            f.write("@ xaxis  label \"Distance (Angstrom)\"\n")
            f.write("@ yaxis  label \"g(r)\"\n")
            f.write("@ TYPE xy\n")
            
            for r, g in zip(bins, rdf):
                f.write(f"{r:12.6f} {g:12.6f}\n")