#!/usr/bin/env python3
"""
Unit tests for TopologyChainIdentifier

Tests the topology-based chain identification functionality without
requiring actual MD data files.

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
import unittest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import subprocess

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.topology.topology_chain_identifier import (
    TopologyChainInfo,
    TopologyChainMapping,
    TopologyChainIdentifier
)


class TestTopologyChainInfo(unittest.TestCase):
    """Test TopologyChainInfo dataclass."""

    def test_chain_info_creation(self):
        """Test TopologyChainInfo creation."""
        chain = TopologyChainInfo(
            original_group_id=18,
            original_group_name='ch0_Protein',
            atom_count=1200,
            residue_count=150,
            assigned_component='HLA_alpha',
            confidence=1.0
        )

        self.assertEqual(chain.original_group_id, 18)
        self.assertEqual(chain.original_group_name, 'ch0_Protein')
        self.assertEqual(chain.atom_count, 1200)
        self.assertEqual(chain.residue_count, 150)
        self.assertEqual(chain.assigned_component, 'HLA_alpha')
        self.assertEqual(chain.confidence, 1.0)


class TestTopologyChainMapping(unittest.TestCase):
    """Test TopologyChainMapping dataclass."""

    def test_mapping_success(self):
        """Test successful mapping result."""
        chain1 = TopologyChainInfo(18, 'ch0_Protein', 1200, 150, 'HLA_alpha', 1.0)
        chain2 = TopologyChainInfo(19, 'ch1_Protein', 800, 100, 'beta2m', 1.0)

        mapping = TopologyChainMapping(
            success=True,
            chains={18: chain1, 19: chain2},
            component_map={'HLA_alpha': 18, 'beta2m': 19}
        )

        self.assertTrue(mapping.success)
        self.assertEqual(len(mapping.chains), 2)
        self.assertEqual(mapping.component_map['HLA_alpha'], 18)
        self.assertIsNone(mapping.error_message)

    def test_mapping_failure(self):
        """Test failed mapping result."""
        mapping = TopologyChainMapping(
            success=False,
            chains={},
            component_map={},
            error_message="No chains found"
        )

        self.assertFalse(mapping.success)
        self.assertEqual(len(mapping.chains), 0)
        self.assertEqual(mapping.error_message, "No chains found")


class TestTopologyChainIdentifier(unittest.TestCase):
    """Test TopologyChainIdentifier core functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.identifier = TopologyChainIdentifier()

    def test_initialization(self):
        """Test identifier initialization."""
        self.assertEqual(self.identifier.gmx, "gmx")

        custom_identifier = TopologyChainIdentifier(gmx_executable="gmx_mpi")
        self.assertEqual(custom_identifier.gmx, "gmx_mpi")

    def test_length_ranges(self):
        """Test length range definitions."""
        ranges = TopologyChainIdentifier.LENGTH_RANGES

        self.assertIn('peptide', ranges)
        self.assertIn('beta2m', ranges)
        self.assertIn('TCR_alpha', ranges)
        self.assertIn('TCR_beta', ranges)
        self.assertIn('HLA_alpha', ranges)

        # Check peptide range
        self.assertEqual(ranges['peptide']['min_res'], 5)
        self.assertEqual(ranges['peptide']['max_res'], 22)
        self.assertEqual(ranges['peptide']['typical'], 9)

    def test_parse_splitch_output(self):
        """Test parsing of gmx make_ndx splitch output."""
        stderr_output = """
Splitting Protein into chains (using column 21 of the pdb file).
There are 5 chains and 0 residues with unknown chain ID.

  0 System              : 50000 atoms
  1 Protein             :  5000 atoms
 18 ch0_Protein         :  2200 atoms
 19 ch1_Protein         :   800 atoms
 20 ch2_Protein         :    70 atoms
 21 ch3_Protein         :  1440 atoms
 22 ch4_Protein         :  1920 atoms
"""

        chain_groups = self.identifier._parse_splitch_output(stderr_output)

        self.assertEqual(len(chain_groups), 5)

        # Check first chain
        self.assertEqual(chain_groups[0]['group_id'], 18)
        self.assertEqual(chain_groups[0]['name'], 'ch0_Protein')
        self.assertEqual(chain_groups[0]['atom_count'], 2200)

        # Check last chain
        self.assertEqual(chain_groups[4]['group_id'], 22)
        self.assertEqual(chain_groups[4]['name'], 'ch4_Protein')
        self.assertEqual(chain_groups[4]['atom_count'], 1920)

    def test_extract_chain_info(self):
        """Test chain information extraction."""
        chain_groups = [
            {'group_id': 18, 'name': 'ch0_Protein', 'atom_count': 2200},
            {'group_id': 19, 'name': 'ch1_Protein', 'atom_count': 800},
            {'group_id': 20, 'name': 'ch2_Protein', 'atom_count': 70}
        ]

        chains_info = self.identifier._extract_chain_info(chain_groups, "test.tpr")

        self.assertEqual(len(chains_info), 3)

        # Check residue count estimation (atom_count / 8)
        self.assertEqual(chains_info[0].residue_count, 2200 // 8)  # 275
        self.assertEqual(chains_info[1].residue_count, 800 // 8)   # 100
        self.assertEqual(chains_info[2].residue_count, 70 // 8)    # 8

        # Check other fields
        self.assertEqual(chains_info[0].original_group_id, 18)
        self.assertEqual(chains_info[0].original_group_name, 'ch0_Protein')
        self.assertEqual(chains_info[0].atom_count, 2200)

    def test_assign_components_by_length(self):
        """Test component assignment based on length."""
        # Create chains with typical lengths
        chains_info = [
            TopologyChainInfo(18, 'ch0_Protein', 2200, 275, '', 0.0),  # HLA_alpha (longest)
            TopologyChainInfo(19, 'ch1_Protein', 800, 100, '', 0.0),   # beta2m (2nd shortest)
            TopologyChainInfo(20, 'ch2_Protein', 70, 9, '', 0.0),      # peptide (shortest)
            TopologyChainInfo(21, 'ch3_Protein', 1440, 180, '', 0.0),  # TCR_alpha (middle)
            TopologyChainInfo(22, 'ch4_Protein', 1920, 240, '', 0.0)   # TCR_beta (2nd longest)
        ]

        assigned = self.identifier._assign_components_by_length(chains_info)

        # Check assignments (sorted by length)
        self.assertEqual(len(assigned), 5)

        # Find assignments by group_id
        group_20 = assigned[20]  # Shortest
        self.assertEqual(group_20.assigned_component, 'peptide')
        self.assertEqual(group_20.confidence, 1.0)

        group_19 = assigned[19]  # 2nd shortest
        self.assertEqual(group_19.assigned_component, 'beta2m')

        group_21 = assigned[21]  # Middle
        self.assertEqual(group_21.assigned_component, 'TCR_alpha')

        group_22 = assigned[22]  # 2nd longest
        self.assertEqual(group_22.assigned_component, 'TCR_beta')

        group_18 = assigned[18]  # Longest
        self.assertEqual(group_18.assigned_component, 'HLA_alpha')

    def test_assign_components_out_of_range(self):
        """Test component assignment with out-of-range lengths."""
        # Create chains with unusual lengths
        chains_info = [
            TopologyChainInfo(18, 'ch0_Protein', 2400, 300, '', 0.0),  # Too long for HLA
            TopologyChainInfo(19, 'ch1_Protein', 800, 100, '', 0.0),
            TopologyChainInfo(20, 'ch2_Protein', 200, 25, '', 0.0),    # Too long for peptide
            TopologyChainInfo(21, 'ch3_Protein', 1440, 180, '', 0.0),
            TopologyChainInfo(22, 'ch4_Protein', 1920, 240, '', 0.0)
        ]

        assigned = self.identifier._assign_components_by_length(chains_info)

        # Check that out-of-range chains have lower confidence
        group_20 = assigned[20]  # Out of peptide range
        self.assertEqual(group_20.assigned_component, 'peptide')
        self.assertEqual(group_20.confidence, 0.5)

        group_18 = assigned[18]  # Out of HLA_alpha range
        self.assertEqual(group_18.assigned_component, 'HLA_alpha')
        self.assertEqual(group_18.confidence, 0.5)

    def test_get_component_selection_string(self):
        """Test getting GROMACS selection string for components."""
        # Create mock mapping
        chain1 = TopologyChainInfo(18, 'ch0_Protein', 2200, 275, 'HLA_alpha', 1.0)
        chain2 = TopologyChainInfo(19, 'ch1_Protein', 800, 100, 'beta2m', 1.0)
        chain3 = TopologyChainInfo(20, 'ch2_Protein', 70, 9, 'peptide', 1.0)
        chain4 = TopologyChainInfo(21, 'ch3_Protein', 1440, 180, 'TCR_alpha', 1.0)
        chain5 = TopologyChainInfo(22, 'ch4_Protein', 1920, 240, 'TCR_beta', 1.0)

        mapping = TopologyChainMapping(
            success=True,
            chains={18: chain1, 19: chain2, 20: chain3, 21: chain4, 22: chain5},
            component_map={
                'HLA_alpha': 18,
                'beta2m': 19,
                'peptide': 20,
                'TCR_alpha': 21,
                'TCR_beta': 22
            }
        )

        # Test individual components
        self.assertEqual(
            self.identifier.get_component_selection_string('peptide', mapping),
            '20'
        )
        self.assertEqual(
            self.identifier.get_component_selection_string('HLA_alpha', mapping),
            '18'
        )

        # Test combined components
        phla_selection = self.identifier.get_component_selection_string('pHLA', mapping)
        self.assertIn('18', phla_selection)  # HLA_alpha
        self.assertIn('19', phla_selection)  # beta2m
        self.assertIn('20', phla_selection)  # peptide
        self.assertIn('|', phla_selection)

        tcr_selection = self.identifier.get_component_selection_string('TCR', mapping)
        self.assertIn('21', tcr_selection)  # TCR_alpha
        self.assertIn('22', tcr_selection)  # TCR_beta
        self.assertIn('|', tcr_selection)

        hla_selection = self.identifier.get_component_selection_string('HLA', mapping)
        self.assertIn('18', hla_selection)  # HLA_alpha
        self.assertIn('19', hla_selection)  # beta2m

    def test_get_component_selection_missing(self):
        """Test getting selection string for missing component."""
        mapping = TopologyChainMapping(
            success=True,
            chains={},
            component_map={'peptide': 20}
        )

        # Missing component
        result = self.identifier.get_component_selection_string('HLA_alpha', mapping)
        self.assertIsNone(result)

        # Failed mapping
        failed_mapping = TopologyChainMapping(
            success=False,
            chains={},
            component_map={},
            error_message="Test error"
        )
        result = self.identifier.get_component_selection_string('peptide', failed_mapping)
        self.assertIsNone(result)


def run_tests():
    """Run all tests."""
    print("=" * 70)
    print("Running TopologyChainIdentifier Unit Tests")
    print("=" * 70)

    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestTopologyChainInfo))
    suite.addTests(loader.loadTestsFromTestCase(TestTopologyChainMapping))
    suite.addTests(loader.loadTestsFromTestCase(TestTopologyChainIdentifier))

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    print("\n" + "=" * 70)
    if result.wasSuccessful():
        print("All tests passed!")
    else:
        print(f"Tests failed: {len(result.failures)} failures, {len(result.errors)} errors")
    print("=" * 70)

    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
