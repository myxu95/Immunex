#!/usr/bin/env python3
"""
CDR Loop Selector - Reusable utility for identifying TCR CDR regions

This module provides functionality to locate CDR loops in MD trajectories
using exact sequence matching from a reference CSV file.
"""
import pandas as pd
import numpy as np
import MDAnalysis as mda
from typing import Dict, Optional, Tuple, List
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class CDRSelector:
    """
    Select CDR regions from MD trajectories using exact sequence matching

    This class provides methods to identify CDR loops by matching exact sequences
    from a reference CSV file against all chains in a PDB structure.

    Attributes:
        cdr_reference: Dictionary containing CDR sequences for each PDB
        aa_code: Amino acid three-letter to one-letter code mapping

    Example:
        >>> selector = CDRSelector("tcr_cdrs_output.csv")
        >>> cdr_info = selector.get_cdr_regions(
        ...     pdb_file="md_converted.pdb",
        ...     pdb_id="1AO7"
        ... )
        >>> alpha_cdr3 = cdr_info['alpha']['cdr3']
        >>> print(f"CDR3 chain: {alpha_cdr3['chain_id']}")
        >>> print(f"CDR3 residues: {alpha_cdr3['residue_range']}")
    """

    def __init__(self, cdr_csv_file: str):
        """
        Initialize CDR selector with reference data

        Args:
            cdr_csv_file: Path to CSV file containing CDR sequences
                         Expected columns: PDB_ID, alpha_cdr1, alpha_cdr2, alpha_cdr3,
                                         beta_cdr1, beta_cdr2, beta_cdr3
        """
        self.cdr_reference = self._load_cdr_reference(cdr_csv_file)

        # Amino acid code mapping
        self.aa_code = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }

    def _load_cdr_reference(self, csv_file: str) -> Dict:
        """Load CDR sequences from CSV file"""
        df = pd.read_csv(csv_file)

        cdr_dict = {}
        for _, row in df.iterrows():
            pdb_id = row['PDB_ID'].upper()

            cdr_dict[pdb_id] = {
                'alpha': {
                    'sequence': row['alpha_variable_domain'],
                    'cdr1': row['alpha_cdr1'] if pd.notna(row['alpha_cdr1']) else None,
                    'cdr2': row['alpha_cdr2'] if pd.notna(row['alpha_cdr2']) else None,
                    'cdr3': row['alpha_cdr3'] if pd.notna(row['alpha_cdr3']) else None,
                },
                'beta': {
                    'sequence': row['beta_variable_domain'],
                    'cdr1': row['beta_cdr1'] if pd.notna(row['beta_cdr1']) else None,
                    'cdr2': row['beta_cdr2'] if pd.notna(row['beta_cdr2']) else None,
                    'cdr3': row['beta_cdr3'] if pd.notna(row['beta_cdr3']) else None,
                }
            }

        return cdr_dict

    def get_chain_sequences(self, pdb_file: str) -> Dict[str, Dict]:
        """
        Extract sequences from all chains in PDB file

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary mapping chain IDs to their sequences and residues
        """
        u = mda.Universe(pdb_file)

        chain_sequences = {}
        for chain_id in set(u.atoms.chainIDs):
            if not chain_id:
                continue

            chain_atoms = u.select_atoms(f'chainID {chain_id} and name CA')
            if len(chain_atoms) > 0:
                seq = ''.join([self.aa_code.get(r.resname, 'X') for r in chain_atoms.residues])
                chain_sequences[chain_id] = {
                    'sequence': seq,
                    'residues': chain_atoms.residues
                }

        return chain_sequences

    def find_sequence_in_chain(self, chain_sequence: str, target_sequence: str) -> Optional[Tuple[int, int]]:
        """
        Find target sequence within a chain sequence

        Args:
            chain_sequence: Full sequence of the chain
            target_sequence: CDR sequence to find

        Returns:
            Tuple of (start_residue, end_residue) in 1-based numbering, or None if not found
        """
        if not target_sequence or pd.isna(target_sequence):
            return None

        pos = chain_sequence.find(target_sequence)
        if pos == -1:
            return None

        # Convert to 1-based residue numbering
        start = pos + 1
        end = pos + len(target_sequence)

        return (start, end)

    def get_cdr_regions(self, pdb_file: str, pdb_id: str) -> Dict:
        """
        Get all CDR regions for a specific PDB

        Args:
            pdb_file: Path to PDB file
            pdb_id: PDB identifier (e.g., "1AO7")

        Returns:
            Dictionary with CDR information for alpha and beta chains:
            {
                'alpha': {
                    'cdr1': {'chain_id': 'D', 'sequence': 'DRGSQS',
                            'residue_range': [26, 31], 'ca_indices': [...]},
                    'cdr2': {...},
                    'cdr3': {...}
                },
                'beta': {...}
            }
        """
        pdb_id = pdb_id.upper()

        if pdb_id not in self.cdr_reference:
            raise ValueError(f"PDB {pdb_id} not found in reference data")

        cdr_ref = self.cdr_reference[pdb_id]
        u = mda.Universe(pdb_file)
        chain_sequences = self.get_chain_sequences(pdb_file)

        cdr_regions = {}

        for chain_type in ['alpha', 'beta']:
            chain_data = {}

            # Search for each CDR in all chains
            for cdr_num in [1, 2, 3]:
                cdr_key = f'cdr{cdr_num}'
                cdr_seq = cdr_ref[chain_type][cdr_key]

                if not cdr_seq or pd.isna(cdr_seq):
                    continue

                # Search all chains for this CDR sequence
                found = False
                for chain_id, chain_info in chain_sequences.items():
                    cdr_range = self.find_sequence_in_chain(chain_info['sequence'], cdr_seq)

                    if cdr_range:
                        start, end = cdr_range

                        # Get CA atom indices
                        cdr_atoms = u.select_atoms(
                            f'chainID {chain_id} and name CA and resid {start}:{end}'
                        )

                        chain_data[cdr_key] = {
                            'chain_id': chain_id,
                            'sequence': cdr_seq,
                            'residue_range': [int(start), int(end)],
                            'residue_count': end - start + 1,
                            'ca_count': len(cdr_atoms),
                            'ca_indices': [int(a.index + 1) for a in cdr_atoms]
                        }

                        logger.info(f"  {chain_type.upper()} CDR{cdr_num}: {cdr_seq}")
                        logger.info(f"    Chain {chain_id}, Residues: {start}-{end} ({len(cdr_atoms)} CA atoms)")
                        found = True
                        break

                if not found:
                    logger.warning(f"  {chain_type.upper()} CDR{cdr_num}: Sequence '{cdr_seq}' not found")

            cdr_regions[chain_type] = chain_data

        return cdr_regions

    def create_selection_string(self, cdr_info: Dict, chain_type: str, cdr_number: int,
                               atom_name: str = 'CA') -> Optional[str]:
        """
        Create MDAnalysis selection string for a specific CDR

        Args:
            cdr_info: CDR information dictionary from get_cdr_regions()
            chain_type: 'alpha' or 'beta'
            cdr_number: 1, 2, or 3
            atom_name: Atom name to select (default: 'CA')

        Returns:
            Selection string like "chainID D and name CA and resid 26:31"
            or None if CDR not found

        Example:
            >>> cdr_info = selector.get_cdr_regions("md.pdb", "1AO7")
            >>> sel_str = selector.create_selection_string(cdr_info, 'alpha', 3)
            >>> atoms = u.select_atoms(sel_str)
        """
        cdr_key = f'cdr{cdr_number}'

        if chain_type not in cdr_info or cdr_key not in cdr_info[chain_type]:
            return None

        cdr_data = cdr_info[chain_type][cdr_key]
        chain_id = cdr_data['chain_id']
        start, end = cdr_data['residue_range']

        return f"chainID {chain_id} and name {atom_name} and resid {start}:{end}"

    def get_all_cdr_atoms(self, universe: mda.Universe, cdr_info: Dict,
                         atom_name: str = 'CA') -> Dict:
        """
        Get MDAnalysis AtomGroups for all CDRs

        Args:
            universe: MDAnalysis Universe object
            cdr_info: CDR information from get_cdr_regions()
            atom_name: Atom name to select (default: 'CA')

        Returns:
            Dictionary of AtomGroups:
            {
                'alpha_cdr1': AtomGroup,
                'alpha_cdr2': AtomGroup,
                'alpha_cdr3': AtomGroup,
                'beta_cdr1': AtomGroup,
                'beta_cdr2': AtomGroup,
                'beta_cdr3': AtomGroup
            }
        """
        cdr_atoms = {}

        for chain_type in ['alpha', 'beta']:
            for cdr_num in [1, 2, 3]:
                key = f'{chain_type}_cdr{cdr_num}'
                sel_str = self.create_selection_string(cdr_info, chain_type, cdr_num, atom_name)

                if sel_str:
                    cdr_atoms[key] = universe.select_atoms(sel_str)

        return cdr_atoms
