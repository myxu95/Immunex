#!/usr/bin/env python3
"""
Robust Chain Detection for MD Systems

A more robust approach to chain detection and group selection for PBC processing,
with multiple fallback strategies and better error handling.
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import subprocess

logger = logging.getLogger(__name__)


class RobustChainDetector:
    """Robust chain detector with multiple fallback strategies."""

    def __init__(self, topology_file: str, gro_file: Optional[str] = None, gmx_executable: str = "gmx"):
        self.topology_file = topology_file
        self.gro_file = gro_file
        self.gmx = gmx_executable

        # Protein residue types (expanded list)
        self.protein_residues = {
            # Standard amino acids
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            # Modified residues
            'ACE', 'NME', 'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN',
            # Alternative names
            'HSD', 'HSE', 'HSP', 'CSS', 'CYM'
        }

    def get_center_group_with_fallback(self) -> Tuple[str, str]:
        """
        Get appropriate center group using multiple strategies with fallbacks.

        Returns:
            Tuple of (group_selection, group_description)
        """
        strategies = [
            self._try_shortest_chain_from_tpr,
            self._try_shortest_chain_from_gro,
            self._try_manual_chain_selection,
            self._try_backbone_group,
            self._try_calpha_group
        ]

        for i, strategy in enumerate(strategies, 1):
            try:
                logger.info(f"尝试策略 {i}: {strategy.__name__}")
                result = strategy()
                if result:
                    group_sel, description = result
                    logger.info(f"策略 {i} 成功: {description}")
                    return group_sel, description

            except Exception as e:
                logger.warning(f"策略 {i} 失败: {e}")
                continue

        # Ultimate fallback - but with warning
        logger.error("所有策略失败，使用最后的备用方案")
        return "4", "Backbone (Last Resort - May Not Be Optimal)"

    def _try_shortest_chain_from_tpr(self) -> Optional[Tuple[str, str]]:
        """Strategy 1: Try to create shortest chain index from TPR file."""
        try:
            # Try gmx make_ndx with different approaches
            approaches = [
                self._gmx_make_ndx_standard,
                self._gmx_make_ndx_quiet,
                self._gmx_make_ndx_minimal
            ]

            for approach in approaches:
                try:
                    chains = approach()
                    if chains:
                        shortest_chain = min(chains, key=lambda x: len(chains[x]))
                        # Create index file with shortest chain
                        index_file = self._create_shortest_chain_index(chains[shortest_chain])
                        if index_file:
                            return "18", f"Shortest Chain from TPR ({len(chains[shortest_chain])} atoms)"
                except Exception as e:
                    logger.debug(f"TPR approach failed: {e}")
                    continue

            return None

        except Exception as e:
            logger.warning(f"TPR strategy failed: {e}")
            return None

    def _try_shortest_chain_from_gro(self) -> Optional[Tuple[str, str]]:
        """Strategy 2: Parse chains from GRO file with multiple parsing methods."""
        if not self.gro_file or not Path(self.gro_file).exists():
            return None

        try:
            parsing_methods = [
                self._parse_gro_strict,
                self._parse_gro_flexible,
                self._parse_gro_residue_based
            ]

            for method in parsing_methods:
                try:
                    chains = method()
                    if chains and len(chains) > 1:  # Only use if multiple chains detected
                        shortest_chain_atoms = min(chains.values(), key=len)
                        index_file = self._create_shortest_chain_index(shortest_chain_atoms)
                        if index_file:
                            return "18", f"Shortest Chain from GRO ({len(shortest_chain_atoms)} atoms)"
                except Exception as e:
                    logger.debug(f"GRO parsing method failed: {e}")
                    continue

            return None

        except Exception as e:
            logger.warning(f"GRO strategy failed: {e}")
            return None

    def _try_manual_chain_selection(self) -> Optional[Tuple[str, str]]:
        """Strategy 3: Try to manually identify reasonable groups."""
        try:
            # Get available groups
            groups = self._get_available_groups()
            if not groups:
                return None

            # Look for chain-specific groups (common GROMACS patterns)
            chain_patterns = [
                r'Chain[_\s]*[A-Z]',
                r'Protein[_\s]*Chain[_\s]*[A-Z]',
                r'[A-Z][_\s]*chain',
                r'chain[_\s]*[A-Z]'
            ]

            chain_groups = []
            for group_id, group_info in groups.items():
                group_name = group_info['name']
                for pattern in chain_patterns:
                    if re.search(pattern, group_name, re.IGNORECASE):
                        chain_groups.append((group_id, group_info))
                        break

            if chain_groups:
                # Select smallest chain group
                smallest_chain = min(chain_groups, key=lambda x: x[1]['atoms'])
                return str(smallest_chain[0]), f"Manual Chain Selection: {smallest_chain[1]['name']}"

            return None

        except Exception as e:
            logger.warning(f"Manual selection strategy failed: {e}")
            return None

    def _try_backbone_group(self) -> Optional[Tuple[str, str]]:
        """Strategy 4: Use backbone group if available."""
        try:
            groups = self._get_available_groups()

            # Look for backbone-like groups
            backbone_patterns = ['Backbone', 'MainChain', 'Protein-H', 'C-alpha']

            for pattern in backbone_patterns:
                for group_id, group_info in groups.items():
                    if pattern.lower() in group_info['name'].lower():
                        return str(group_id), f"Backbone Group: {group_info['name']}"

            return None

        except Exception as e:
            logger.warning(f"Backbone strategy failed: {e}")
            return None

    def _try_calpha_group(self) -> Optional[Tuple[str, str]]:
        """Strategy 5: Use C-alpha group."""
        try:
            groups = self._get_available_groups()

            for group_id, group_info in groups.items():
                if 'alpha' in group_info['name'].lower() or 'ca' in group_info['name'].lower():
                    return str(group_id), f"C-alpha Group: {group_info['name']}"

            return None

        except Exception as e:
            logger.warning(f"C-alpha strategy failed: {e}")
            return None

    def _get_available_groups(self) -> Dict[int, Dict[str, Union[str, int]]]:
        """Get available groups from topology."""
        try:
            result = subprocess.run(
                [self.gmx, "make_ndx", "-f", self.topology_file, "-o", "/dev/null"],
                input="q\n",
                text=True,
                capture_output=True,
                timeout=30
            )

            # Parse groups from output
            groups = {}
            for line in result.stderr.split('\n'):
                # Look for group definitions
                match = re.search(r'^\s*(\d+)\s+(.+?):\s*(\d+)\s+atoms', line)
                if match:
                    group_id = int(match.group(1))
                    group_name = match.group(2).strip()
                    atom_count = int(match.group(3))

                    groups[group_id] = {
                        'name': group_name,
                        'atoms': atom_count
                    }

            return groups

        except Exception as e:
            logger.error(f"Failed to get available groups: {e}")
            return {}

    def _parse_gro_flexible(self) -> Dict[int, List[int]]:
        """Flexible GRO parsing with multiple format handling."""
        chains = {}
        current_chain = 0
        current_atoms = []
        last_resid = None

        try:
            with open(self.gro_file, 'r') as f:
                lines = f.readlines()

            atom_lines = lines[2:-1]  # Skip header, count, and box line

            for line in atom_lines:
                if len(line.strip()) < 30:  # Minimum reasonable line length
                    continue

                try:
                    # Multiple parsing approaches
                    atom_info = self._parse_gro_line_flexible(line)
                    if not atom_info:
                        continue

                    resid, resname, atom_id = atom_info

                    if resname in self.protein_residues:
                        # Chain detection logic with more flexibility
                        if last_resid is not None:
                            # Allow for larger gaps in residue numbering
                            if abs(resid - last_resid) > 5:  # More lenient gap detection
                                if len(current_atoms) >= 10:  # Lower threshold
                                    chains[current_chain] = current_atoms
                                    current_chain += 1
                                current_atoms = []

                        current_atoms.append(atom_id)
                        last_resid = resid

                except Exception as e:
                    logger.debug(f"Failed to parse line: {line.strip()}: {e}")
                    continue

            # Add final chain
            if len(current_atoms) >= 10:
                chains[current_chain] = current_atoms

            logger.info(f"Flexible parsing detected {len(chains)} chains")
            return chains

        except Exception as e:
            logger.error(f"Flexible GRO parsing failed: {e}")
            return {}

    def _parse_gro_line_flexible(self, line: str) -> Optional[Tuple[int, str, int]]:
        """Parse GRO line with multiple format attempts."""
        try:
            # Attempt 1: Standard format
            if len(line) >= 44:
                resid_resname = line[0:5].strip()
                atom_id = int(line[10:15].strip())

                # Try different residue parsing methods
                for pattern in [r'(\d+)([A-Z]+)', r'(\d+)\s+([A-Z]+)', r'(\d+)(.+)']:
                    match = re.match(pattern, resid_resname)
                    if match:
                        resid = int(match.group(1))
                        resname = match.group(2).strip()
                        return resid, resname, atom_id

            # Attempt 2: Space-separated format
            parts = line.split()
            if len(parts) >= 6:
                # Look for patterns like "123ALA" or "123 ALA"
                resid_part = parts[0]
                match = re.search(r'(\d+)([A-Z]+)', resid_part)
                if match:
                    resid = int(match.group(1))
                    resname = match.group(2)
                    atom_id = int(parts[2])  # Usually third field
                    return resid, resname, atom_id

            return None

        except Exception as e:
            logger.debug(f"Line parsing failed: {e}")
            return None

    # Placeholder methods for other parsing strategies
    def _gmx_make_ndx_standard(self):
        """Standard gmx make_ndx approach."""
        # Implementation would go here
        return None

    def _gmx_make_ndx_quiet(self):
        """Quiet gmx make_ndx approach."""
        # Implementation would go here
        return None

    def _gmx_make_ndx_minimal(self):
        """Minimal gmx make_ndx approach."""
        # Implementation would go here
        return None

    def _parse_gro_strict(self):
        """Strict GRO parsing."""
        # Implementation would go here
        return None

    def _parse_gro_residue_based(self):
        """Residue-based GRO parsing."""
        # Implementation would go here
        return None

    def _create_shortest_chain_index(self, atoms: List[int]) -> Optional[str]:
        """Create index file with shortest chain."""
        # Implementation would go here
        return None


def get_robust_center_group(topology_file: str,
                          gro_file: Optional[str] = None,
                          gmx_executable: str = "gmx") -> Tuple[str, str]:
    """
    Get center group using robust detection with multiple fallback strategies.

    Returns:
        Tuple of (group_selection, description)
    """
    detector = RobustChainDetector(topology_file, gro_file, gmx_executable)
    return detector.get_center_group_with_fallback()