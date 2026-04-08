#!/usr/bin/env python3
"""
pHLA-TCR Complex Trajectory Analysis Module

Comprehensive analysis toolkit for pHLA-TCR protein complex MD trajectories.
Focuses on interface stability, binding interactions, and dynamic correlations.

Author: Immunex Development Team
Date: 2025-12-07
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, distances, contacts
from typing import Dict, List, Tuple, Optional
import logging
from pathlib import Path
import json
from datetime import datetime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class pHLATCRAnalyzer:
    """
    Comprehensive analysis pipeline for pHLA-TCR protein complexes.

    Chain structure:
    - Chain A: HLA-α (MHC heavy chain, ~275 residues)
    - Chain B: HLA-β (β2-microglobulin, ~100 residues)
    - Chain C: Peptide (antigen peptide, ~9 residues)
    - Chain D: TCR-α (~200 residues)
    - Chain E: TCR-β (~250 residues)
    """

    def __init__(self,
                 task_name: str,
                 trajectory: str,
                 topology: str,
                 output_dir: Optional[str] = None):
        """
        Initialize pHLA-TCR analyzer.

        Args:
            task_name: Task identifier (e.g., "1ao7_run1")
            trajectory: Path to trajectory file (.xtc)
            topology: Path to topology file (.tpr or .pdb)
            output_dir: Output directory for results
        """
        self.task_name = task_name
        self.trajectory = Path(trajectory)
        self.topology = Path(topology)

        # Setup output directory
        if output_dir is None:
            output_dir = Path(f"./analysis_results/{task_name}")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        self.subdirs = {
            'stability': self.output_dir / 'stability',
            'interface': self.output_dir / 'interface',
            'hydrogen_bonds': self.output_dir / 'hydrogen_bonds',
            'contacts': self.output_dir / 'contacts',
            'dynamics': self.output_dir / 'dynamics',
            'figures': self.output_dir / 'figures'
        }
        for subdir in self.subdirs.values():
            subdir.mkdir(exist_ok=True)

        # Load trajectory
        try:
            self.u = mda.Universe(str(topology), str(trajectory))
            logger.info(f"Loaded trajectory: {len(self.u.trajectory)} frames")
            logger.info(f"Total atoms: {self.u.atoms.n_atoms}")
        except Exception as e:
            logger.error(f"Failed to load trajectory: {e}")
            raise

        # Define chain selections
        self.chain_selections = {
            'HLA_alpha': 'chainID A and protein',
            'HLA_beta': 'chainID B and protein',
            'Peptide': 'chainID C and protein',
            'TCR_alpha': 'chainID D and protein',
            'TCR_beta': 'chainID E and protein'
        }

        # Validate chains exist
        self._validate_chains()

        # Store analysis results
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

    def analyze_global_rmsd(self,
                           selection: str = "protein and name CA",
                           reference_frame: int = 0) -> pd.DataFrame:
        """
        Calculate global RMSD for the entire complex.

        Args:
            selection: Atom selection for RMSD
            reference_frame: Reference frame index

        Returns:
            DataFrame with time and RMSD values
        """
        logger.info("Calculating global RMSD...")

        # Select atoms
        atoms = self.u.select_atoms(selection)

        # Set reference
        reference = self.u.trajectory[reference_frame]
        ref_positions = atoms.positions.copy()

        # Calculate RMSD for each frame
        rmsd_values = []
        times = []

        for ts in self.u.trajectory:
            # Align to reference
            rms.rmsd(atoms.positions, ref_positions, superposition=True)
            rmsd_val = rms.rmsd(atoms.positions, ref_positions)

            rmsd_values.append(rmsd_val)
            times.append(ts.time)

        # Create DataFrame
        df = pd.DataFrame({
            'time_ps': times,
            'time_ns': np.array(times) / 1000.0,
            'rmsd_angstrom': rmsd_values
        })

        # Save results
        output_file = self.subdirs['stability'] / 'global_rmsd.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"Global RMSD saved to {output_file}")

        # Store in results
        self.results['global_rmsd'] = {
            'mean': np.mean(rmsd_values),
            'std': np.std(rmsd_values),
            'min': np.min(rmsd_values),
            'max': np.max(rmsd_values)
        }

        return df

    def analyze_chain_rmsd(self, reference_frame: int = 0) -> Dict[str, pd.DataFrame]:
        """
        Calculate RMSD for each individual chain.

        Args:
            reference_frame: Reference frame index

        Returns:
            Dictionary of DataFrames for each chain
        """
        logger.info("Calculating per-chain RMSD...")

        chain_results = {}

        for chain_name, selection in self.chain_selections.items():
            logger.info(f"  Processing {chain_name}...")

            # Select CA atoms for this chain
            ca_selection = f"{selection} and name CA"
            atoms = self.u.select_atoms(ca_selection)

            if len(atoms) == 0:
                logger.warning(f"    No atoms found for {chain_name}")
                continue

            # Set reference
            reference = self.u.trajectory[reference_frame]
            ref_positions = atoms.positions.copy()

            # Calculate RMSD
            rmsd_values = []
            times = []

            for ts in self.u.trajectory:
                rmsd_val = rms.rmsd(atoms.positions, ref_positions, superposition=True)
                rmsd_values.append(rmsd_val)
                times.append(ts.time)

            # Create DataFrame
            df = pd.DataFrame({
                'time_ps': times,
                'time_ns': np.array(times) / 1000.0,
                'rmsd_angstrom': rmsd_values
            })

            chain_results[chain_name] = df

        # Combine all chains into one DataFrame
        combined_df = pd.DataFrame({'time_ns': chain_results[list(chain_results.keys())[0]]['time_ns']})

        for chain_name, df in chain_results.items():
            combined_df[f'{chain_name}_rmsd'] = df['rmsd_angstrom']

        # Save combined results
        output_file = self.subdirs['stability'] / 'chain_rmsd.csv'
        combined_df.to_csv(output_file, index=False)
        logger.info(f"Chain RMSD saved to {output_file}")

        # Store statistics
        self.results['chain_rmsd'] = {}
        for chain_name, df in chain_results.items():
            self.results['chain_rmsd'][chain_name] = {
                'mean': np.mean(df['rmsd_angstrom']),
                'std': np.std(df['rmsd_angstrom']),
                'min': np.min(df['rmsd_angstrom']),
                'max': np.max(df['rmsd_angstrom'])
            }

        return chain_results

    def analyze_interface_distances(self) -> pd.DataFrame:
        """
        Calculate key interface distances over time.

        Key interfaces:
        1. Peptide to HLA binding groove (COM distance)
        2. TCR to Peptide (minimum distance)
        3. TCR to HLA (COM distance)

        Returns:
            DataFrame with interface distances
        """
        logger.info("Calculating interface distances...")

        # Define interface pairs
        interfaces = {
            'Peptide_HLA_COM': {
                'group1': 'chainID C',
                'group2': 'chainID A or chainID B'
            },
            'TCR_Peptide_min': {
                'group1': 'chainID D or chainID E',
                'group2': 'chainID C'
            },
            'TCR_HLA_COM': {
                'group1': 'chainID D or chainID E',
                'group2': 'chainID A or chainID B'
            }
        }

        results_data = {'time_ns': []}

        for ts in self.u.trajectory:
            results_data['time_ns'].append(ts.time / 1000.0)

        for interface_name, selections in interfaces.items():
            logger.info(f"  Processing {interface_name}...")

            group1 = self.u.select_atoms(selections['group1'])
            group2 = self.u.select_atoms(selections['group2'])

            if len(group1) == 0 or len(group2) == 0:
                logger.warning(f"    Empty selection for {interface_name}")
                continue

            distances_list = []

            for ts in self.u.trajectory:
                if '_min' in interface_name:
                    # Minimum distance between any atoms
                    dist_array = distances.distance_array(
                        group1.positions,
                        group2.positions
                    )
                    min_dist = np.min(dist_array)
                    distances_list.append(min_dist)
                else:
                    # Center of mass distance
                    com1 = group1.center_of_mass()
                    com2 = group2.center_of_mass()
                    com_dist = np.linalg.norm(com1 - com2)
                    distances_list.append(com_dist)

            results_data[interface_name] = distances_list

        # Create DataFrame
        df = pd.DataFrame(results_data)

        # Save results
        output_file = self.subdirs['interface'] / 'interface_distances.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"Interface distances saved to {output_file}")

        # Store statistics
        self.results['interface_distances'] = {}
        for col in df.columns:
            if col != 'time_ns':
                self.results['interface_distances'][col] = {
                    'mean': np.mean(df[col]),
                    'std': np.std(df[col]),
                    'min': np.min(df[col]),
                    'max': np.max(df[col])
                }

        return df

    def analyze_contacts(self, cutoff: float = 4.5) -> Dict:
        """
        Calculate contact frequency maps for key interfaces.

        Args:
            cutoff: Distance cutoff for defining contacts (Angstroms)

        Returns:
            Dictionary with contact analysis results
        """
        logger.info(f"Calculating contacts (cutoff={cutoff} Å)...")

        # Define interface pairs for contact analysis
        interface_pairs = {
            'Peptide_HLA': {
                'group1': 'chainID C',
                'group2': 'chainID A or chainID B'
            },
            'TCR_Peptide': {
                'group1': 'chainID D or chainID E',
                'group2': 'chainID C'
            },
            'TCR_HLA': {
                'group1': 'chainID D or chainID E',
                'group2': 'chainID A or chainID B'
            }
        }

        contact_results = {}

        for interface_name, selections in interface_pairs.items():
            logger.info(f"  Processing {interface_name}...")

            group1 = self.u.select_atoms(selections['group1'])
            group2 = self.u.select_atoms(selections['group2'])

            if len(group1) == 0 or len(group2) == 0:
                logger.warning(f"    Empty selection for {interface_name}")
                continue

            # Track contacts over time
            n_frames = len(self.u.trajectory)
            contact_counts = []
            contact_fractions = []
            times = []

            for ts in self.u.trajectory:
                # Calculate distance matrix
                dist_matrix = distances.distance_array(
                    group1.positions,
                    group2.positions
                )

                # Count contacts
                n_contacts = np.sum(dist_matrix < cutoff)
                contact_counts.append(n_contacts)

                # Calculate contact fraction
                total_pairs = len(group1) * len(group2)
                contact_fraction = n_contacts / total_pairs if total_pairs > 0 else 0
                contact_fractions.append(contact_fraction)

                times.append(ts.time / 1000.0)

            # Save time series
            df = pd.DataFrame({
                'time_ns': times,
                'contact_count': contact_counts,
                'contact_fraction': contact_fractions
            })

            output_file = self.subdirs['contacts'] / f'{interface_name}_contacts.csv'
            df.to_csv(output_file, index=False)

            contact_results[interface_name] = {
                'mean_contacts': np.mean(contact_counts),
                'std_contacts': np.std(contact_counts),
                'mean_fraction': np.mean(contact_fractions),
                'dataframe': df
            }

        logger.info(f"Contact analysis saved to {self.subdirs['contacts']}")

        # Store in results
        self.results['contacts'] = {
            k: {key: val for key, val in v.items() if key != 'dataframe'}
            for k, v in contact_results.items()
        }

        return contact_results

    def analyze_radius_gyration(self) -> pd.DataFrame:
        """
        Calculate radius of gyration for the complex and individual chains.

        Returns:
            DataFrame with Rg values over time
        """
        logger.info("Calculating radius of gyration...")

        results_data = {'time_ns': []}

        # Whole complex Rg
        protein = self.u.select_atoms("protein")
        rg_complex = []

        for ts in self.u.trajectory:
            results_data['time_ns'].append(ts.time / 1000.0)
            rg_complex.append(protein.radius_of_gyration())

        results_data['Rg_complex'] = rg_complex

        # Per-chain Rg
        for chain_name, selection in self.chain_selections.items():
            atoms = self.u.select_atoms(selection)

            if len(atoms) == 0:
                continue

            rg_values = []
            for ts in self.u.trajectory:
                rg_values.append(atoms.radius_of_gyration())

            results_data[f'Rg_{chain_name}'] = rg_values

        # Create DataFrame
        df = pd.DataFrame(results_data)

        # Save results
        output_file = self.subdirs['stability'] / 'radius_gyration.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"Radius of gyration saved to {output_file}")

        # Store statistics
        self.results['radius_gyration'] = {}
        for col in df.columns:
            if col != 'time_ns':
                self.results['radius_gyration'][col] = {
                    'mean': np.mean(df[col]),
                    'std': np.std(df[col])
                }

        return df

    def analyze_global_rmsf(self, selection: str = "protein and name CA") -> pd.DataFrame:
        """
        Calculate global RMSF (Root Mean Square Fluctuation) for the entire protein complex.

        RMSF quantifies per-atom flexibility by measuring deviation from mean position
        across the trajectory. Higher values indicate more flexible regions.

        Biological Significance:
        - Low RMSF (<1.0 Å): Rigid regions, structural core, stable binding
        - Medium RMSF (1.0-2.5 Å): Moderate flexibility, functional regions
        - High RMSF (>2.5 Å): Highly flexible regions, loops, entropic freedom

        Args:
            selection: Atom selection for RMSF calculation
                      Default: "protein and name CA" (1 CA per residue for residue-level analysis)

        Returns:
            DataFrame with columns:
            - residue_id: Residue index (1-based PDB numbering)
            - residue_name: Three-letter residue code (e.g., 'ALA', 'GLY')
            - atom_name: Atom name (e.g., 'CA')
            - rmsf_angstrom: RMSF value in Angstroms

        Outputs:
            CSV: dynamics/global_rmsf.csv
            Statistics: self.results['global_rmsf']

        Note:
            RMSF calculation assumes trajectory is pre-aligned. Our trajectories
            are already aligned via PBC correction (gmx trjconv -fit rot+trans).
        """
        logger.info("Calculating global RMSF...")

        # Check trajectory length
        if len(self.u.trajectory) < 2:
            logger.error("RMSF requires at least 2 frames in trajectory")
            return pd.DataFrame()

        # Select atoms
        atoms = self.u.select_atoms(selection)
        n_atoms = len(atoms)

        if n_atoms == 0:
            logger.warning(f"Empty selection: {selection}")
            return pd.DataFrame()

        logger.info(f"  Selected {n_atoms} atoms for RMSF calculation")

        # Import RMSF analyzer
        from MDAnalysis.analysis.rms import RMSF

        # Run RMSF analysis
        # Note: RMSF assumes trajectory is already aligned (which ours is)
        rmsf_analyzer = RMSF(atoms).run()
        rmsf_values = rmsf_analyzer.results.rmsf  # in Angstroms

        # Build DataFrame with residue information
        data = {
            'residue_id': atoms.resids,         # 1-based residue numbering
            'residue_name': atoms.resnames,     # e.g., 'ALA', 'GLY'
            'atom_name': atoms.names,           # e.g., 'CA'
            'rmsf_angstrom': rmsf_values
        }

        df = pd.DataFrame(data)

        # Save results
        output_file = self.subdirs['dynamics'] / 'global_rmsf.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"  Global RMSF saved to {output_file}")

        # Calculate and store statistics
        self.results['global_rmsf'] = {
            'mean': float(np.mean(rmsf_values)),
            'std': float(np.std(rmsf_values)),
            'min': float(np.min(rmsf_values)),
            'max': float(np.max(rmsf_values)),
            'median': float(np.median(rmsf_values))
        }

        logger.info(f"  Mean RMSF: {self.results['global_rmsf']['mean']:.2f} Å")
        logger.info(f"  Max RMSF: {self.results['global_rmsf']['max']:.2f} Å")

        # Validation: warn if unusually high RMSF (may indicate misalignment)
        if self.results['global_rmsf']['mean'] > 5.0:
            logger.warning(
                f"Unusually high mean RMSF detected ({self.results['global_rmsf']['mean']:.2f} Å > 5.0 Å). "
                "This may indicate the trajectory is not properly aligned. "
                "Ensure 'gmx trjconv -fit rot+trans' was applied during PBC correction."
            )

        return df

    def analyze_chain_rmsf(self) -> Dict[str, pd.DataFrame]:
        """
        Calculate RMSF for each individual chain.

        Provides per-chain flexibility analysis to identify chain-specific
        dynamic behavior and compare rigidity across different protein components.

        Biological Insights:
        - Peptide: Low RMSF indicates stable binding (lock-in mechanism)
        - HLA chains: Moderate RMSF expected in loops, low in helices
        - TCR chains: CDR loops show highest RMSF (entropic clamp effect)

        Returns:
            Dictionary mapping chain names to DataFrames with RMSF data
            Keys: 'Peptide', 'HLA_alpha', 'HLA_beta', 'TCR_alpha', 'TCR_beta'

        Outputs:
            CSV: dynamics/chain_rmsf.csv (combined data for all chains)
            Statistics: self.results['chain_rmsf'][chain_name]

        Note:
            Only CA atoms are used for residue-level analysis (1 atom per residue)
        """
        logger.info("Calculating per-chain RMSF...")

        # Check trajectory length
        if len(self.u.trajectory) < 2:
            logger.error("RMSF requires at least 2 frames in trajectory")
            return {}

        # Import RMSF analyzer
        from MDAnalysis.analysis.rms import RMSF

        chain_results = {}

        # Process each chain
        for chain_name, selection in self.chain_selections.items():
            logger.info(f"  Processing {chain_name}...")

            # Select CA atoms for this chain
            ca_selection = f"{selection} and name CA"
            atoms = self.u.select_atoms(ca_selection)

            if len(atoms) == 0:
                logger.warning(f"    No CA atoms found for {chain_name}")
                continue

            # Run RMSF analysis
            rmsf_analyzer = RMSF(atoms).run()
            rmsf_values = rmsf_analyzer.results.rmsf

            # Build DataFrame for this chain
            df = pd.DataFrame({
                'residue_id': atoms.resids,
                'residue_name': atoms.resnames,
                'chain_id': [chain_name] * len(atoms),
                'rmsf_angstrom': rmsf_values
            })

            chain_results[chain_name] = df

            mean_rmsf = float(np.mean(rmsf_values))
            logger.info(f"    {len(atoms)} residues, Mean RMSF: {mean_rmsf:.2f} Å")

        # Combine all chains into one CSV
        if chain_results:
            combined_df = pd.concat([df for df in chain_results.values()],
                                    ignore_index=True)

            # Save combined results
            output_file = self.subdirs['dynamics'] / 'chain_rmsf.csv'
            combined_df.to_csv(output_file, index=False)
            logger.info(f"  Chain RMSF saved to {output_file}")

            # Store statistics per chain
            self.results['chain_rmsf'] = {}
            for chain_name, df in chain_results.items():
                self.results['chain_rmsf'][chain_name] = {
                    'mean': float(np.mean(df['rmsf_angstrom'])),
                    'std': float(np.std(df['rmsf_angstrom'])),
                    'min': float(np.min(df['rmsf_angstrom'])),
                    'max': float(np.max(df['rmsf_angstrom'])),
                    'median': float(np.median(df['rmsf_angstrom']))
                }

        return chain_results

    def _extract_sequence_from_topology(self, chain_selection: str) -> str:
        """
        Extract amino acid sequence from MD topology.

        Args:
            chain_selection: Chain selection string (e.g., "chainID D and protein")

        Returns:
            Single-letter amino acid sequence string
        """
        atoms = self.u.select_atoms(f"{chain_selection} and name CA")

        # Three-letter to single-letter conversion
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        sequence = ''.join([aa_map.get(res, 'X') for res in atoms.resnames])

        return sequence

    def _find_subsequence_position(self, full_seq: str, subseq: str) -> tuple:
        """
        Find subsequence position in full sequence.

        Args:
            full_seq: Full sequence
            subseq: Subsequence (CDR3)

        Returns:
            (start_index, end_index) or None if not found
        """
        pos = full_seq.find(subseq)

        if pos == -1:
            return None

        return (pos, pos + len(subseq))

    def analyze_cdr3_rmsf(self, cdr3_sequences: dict) -> dict:
        """
        Calculate RMSF for CDR3 regions.

        Args:
            cdr3_sequences: CDR3 sequence dictionary
                {
                    'cdr3_alpha': str,  # TCR-α CDR3 sequence
                    'cdr3_beta': str    # TCR-β CDR3 sequence
                }

        Returns:
            CDR3 RMSF results dictionary

        Outputs:
            CSV: dynamics/cdr3_alpha_rmsf.csv
            CSV: dynamics/cdr3_beta_rmsf.csv
            Statistics: self.results['cdr3_rmsf']
        """
        logger.info("Calculating CDR3-specific RMSF...")

        # Check trajectory
        if len(self.u.trajectory) < 2:
            logger.error("RMSF requires at least 2 frames")
            return {}

        # Import RMSF analyzer
        from MDAnalysis.analysis.rms import RMSF

        cdr3_results = {}
        self.results['cdr3_rmsf'] = {}

        # Process CDR3α
        if 'cdr3_alpha' in cdr3_sequences and cdr3_sequences['cdr3_alpha']:
            cdr3_alpha_seq = cdr3_sequences['cdr3_alpha']
            logger.info(f"  Processing CDR3α: {cdr3_alpha_seq}")

            # Extract TCR-α full sequence
            tcra_full_seq = self._extract_sequence_from_topology(
                self.chain_selections['TCR_alpha']
            )

            # Find CDR3 position
            pos = self._find_subsequence_position(tcra_full_seq, cdr3_alpha_seq)

            if pos is None:
                logger.warning(f"    CDR3α sequence not found in TCR-α chain")
                logger.warning(f"    TCR-α sequence: {tcra_full_seq[:50]}...")
            else:
                start_idx, end_idx = pos
                logger.info(f"    CDR3α position: sequence index {start_idx}-{end_idx}")

                # Get corresponding residues
                tcra_ca = self.u.select_atoms(f"{self.chain_selections['TCR_alpha']} and name CA")
                cdr3_residues = tcra_ca[start_idx:end_idx]

                logger.info(f"    Residue range: {cdr3_residues.resids[0]}-{cdr3_residues.resids[-1]}")

                # Calculate RMSF
                rmsf_analyzer = RMSF(cdr3_residues).run()
                rmsf_values = rmsf_analyzer.results.rmsf

                # Build DataFrame
                df = pd.DataFrame({
                    'residue_id': cdr3_residues.resids,
                    'residue_name': cdr3_residues.resnames,
                    'residue_seq_index': range(len(cdr3_residues)),
                    'rmsf_angstrom': rmsf_values
                })

                # Save CSV
                output_file = self.subdirs['dynamics'] / 'cdr3_alpha_rmsf.csv'
                df.to_csv(output_file, index=False)
                logger.info(f"    CDR3α RMSF saved to {output_file}")

                # Store statistics
                self.results['cdr3_rmsf']['CDR3_alpha'] = {
                    'mean': float(np.mean(rmsf_values)),
                    'std': float(np.std(rmsf_values)),
                    'min': float(np.min(rmsf_values)),
                    'max': float(np.max(rmsf_values)),
                    'median': float(np.median(rmsf_values)),
                    'sequence': cdr3_alpha_seq,
                    'residue_range': (int(cdr3_residues.resids[0]),
                                     int(cdr3_residues.resids[-1]))
                }

                cdr3_results['CDR3_alpha'] = df

                logger.info(f"    Mean RMSF: {self.results['cdr3_rmsf']['CDR3_alpha']['mean']:.2f} Å")

        # Process CDR3β
        if 'cdr3_beta' in cdr3_sequences and cdr3_sequences['cdr3_beta']:
            cdr3_beta_seq = cdr3_sequences['cdr3_beta']
            logger.info(f"  Processing CDR3β: {cdr3_beta_seq}")

            tcrb_full_seq = self._extract_sequence_from_topology(
                self.chain_selections['TCR_beta']
            )

            pos = self._find_subsequence_position(tcrb_full_seq, cdr3_beta_seq)

            if pos is None:
                logger.warning(f"    CDR3β sequence not found in TCR-β chain")
                logger.warning(f"    TCR-β sequence: {tcrb_full_seq[:50]}...")
            else:
                start_idx, end_idx = pos
                logger.info(f"    CDR3β position: sequence index {start_idx}-{end_idx}")

                tcrb_ca = self.u.select_atoms(f"{self.chain_selections['TCR_beta']} and name CA")
                cdr3_residues = tcrb_ca[start_idx:end_idx]

                logger.info(f"    Residue range: {cdr3_residues.resids[0]}-{cdr3_residues.resids[-1]}")

                rmsf_analyzer = RMSF(cdr3_residues).run()
                rmsf_values = rmsf_analyzer.results.rmsf

                df = pd.DataFrame({
                    'residue_id': cdr3_residues.resids,
                    'residue_name': cdr3_residues.resnames,
                    'residue_seq_index': range(len(cdr3_residues)),
                    'rmsf_angstrom': rmsf_values
                })

                output_file = self.subdirs['dynamics'] / 'cdr3_beta_rmsf.csv'
                df.to_csv(output_file, index=False)
                logger.info(f"    CDR3β RMSF saved to {output_file}")

                self.results['cdr3_rmsf']['CDR3_beta'] = {
                    'mean': float(np.mean(rmsf_values)),
                    'std': float(np.std(rmsf_values)),
                    'min': float(np.min(rmsf_values)),
                    'max': float(np.max(rmsf_values)),
                    'median': float(np.median(rmsf_values)),
                    'sequence': cdr3_beta_seq,
                    'residue_range': (int(cdr3_residues.resids[0]),
                                     int(cdr3_residues.resids[-1]))
                }

                cdr3_results['CDR3_beta'] = df

                logger.info(f"    Mean RMSF: {self.results['cdr3_rmsf']['CDR3_beta']['mean']:.2f} Å")

        return cdr3_results

    def run_full_analysis(self) -> Dict:
        """
        Execute complete analysis workflow.

        Returns:
            Dictionary with all analysis results
        """
        logger.info(f"Starting full analysis for {self.task_name}")
        logger.info(f"Output directory: {self.output_dir}")

        start_time = datetime.now()

        try:
            # 1. Structural stability
            logger.info("\n=== Phase 1: Structural Stability ===")
            self.analyze_global_rmsd()
            self.analyze_chain_rmsd()
            self.analyze_radius_gyration()

            # 1b. Flexibility analysis
            logger.info("\n=== Phase 1b: Flexibility Analysis ===")
            self.analyze_global_rmsf()
            self.analyze_chain_rmsf()

            # 2. Interface analysis
            logger.info("\n=== Phase 2: Interface Analysis ===")
            self.analyze_interface_distances()
            self.analyze_contacts()

            # 3. Save summary
            self._save_summary()

            elapsed = (datetime.now() - start_time).total_seconds()
            logger.info(f"\nAnalysis completed in {elapsed:.2f} seconds")
            logger.info(f"Results saved to: {self.output_dir}")

            return self.results

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise

    def _save_summary(self):
        """Save analysis summary to JSON file."""
        summary = {
            'task_name': self.task_name,
            'trajectory': str(self.trajectory),
            'topology': str(self.topology),
            'n_frames': len(self.u.trajectory),
            'n_atoms': self.u.atoms.n_atoms,
            'timestamp': datetime.now().isoformat(),
            'results': self.results
        }

        output_file = self.output_dir / 'analysis_summary.json'
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)

        logger.info(f"Summary saved to {output_file}")


def main():
    """Example usage."""
    import sys

    if len(sys.argv) < 4:
        print("Usage: python phla_tcr_analyzer.py <task_name> <trajectory> <topology> [output_dir]")
        sys.exit(1)

    task_name = sys.argv[1]
    trajectory = sys.argv[2]
    topology = sys.argv[3]
    output_dir = sys.argv[4] if len(sys.argv) > 4 else None

    # Create analyzer
    analyzer = pHLATCRAnalyzer(
        task_name=task_name,
        trajectory=trajectory,
        topology=topology,
        output_dir=output_dir
    )

    # Run analysis
    results = analyzer.run_full_analysis()

    print("\n=== Analysis Summary ===")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
