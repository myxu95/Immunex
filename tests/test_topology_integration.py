#!/usr/bin/env python3
"""
Test topology-based chain identification integration with IndexManager

This test verifies that the TopologyChainIdentifier is properly integrated
into the IndexManager and that the new methods work correctly.

Author: Immunex Development Team
Date: 2026-03-17
"""

import sys
import unittest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PipelineContext
from immunex.analysis import IndexManager
from immunex.analysis import TopologyChainIdentifier, TopologyChainInfo, TopologyChainMapping


class TestTopologyIntegration(unittest.TestCase):
    """Test integration of TopologyChainIdentifier with IndexManager."""

    def setUp(self):
        """Set up test context."""
        self.context = PipelineContext(
            system_id="test_system",
            topology="test.tpr",
            trajectory_raw="test.xtc",
            output_dir="/tmp/test_output"
        )

    def test_indexmanager_has_topology_methods(self):
        """Test that IndexManager has topology identification methods."""
        index_mgr = IndexManager(self.context)

        # Check methods exist
        self.assertTrue(hasattr(index_mgr, 'identify_chains_from_topology'))
        self.assertTrue(callable(index_mgr.identify_chains_from_topology))

        self.assertTrue(hasattr(index_mgr, 'ensure_base_index_from_topology'))
        self.assertTrue(callable(index_mgr.ensure_base_index_from_topology))

    def test_topology_chain_identifier_imported(self):
        """Test that TopologyChainIdentifier can be imported."""
        from immunex.analysis import TopologyChainIdentifier, TopologyChainInfo, TopologyChainMapping

        # Verify classes exist
        self.assertIsNotNone(TopologyChainIdentifier)
        self.assertIsNotNone(TopologyChainInfo)
        self.assertIsNotNone(TopologyChainMapping)

    @patch('immunex.analysis.topology.topology_chain_identifier.TopologyChainIdentifier.identify_chains_from_topology')
    def test_identify_chains_from_topology_success(self, mock_identify):
        """Test successful topology-based chain identification."""
        # Mock successful result
        mock_result = Mock(spec=TopologyChainMapping)
        mock_result.success = True
        mock_result.component_map = {
            'peptide': 18,
            'beta2m': 19,
            'TCR_alpha': 20,
            'TCR_beta': 21,
            'HLA_alpha': 22
        }
        mock_result.chains = {
            18: Mock(assigned_component='peptide', residue_count=9, atom_count=72),
            19: Mock(assigned_component='beta2m', residue_count=99, atom_count=792),
            20: Mock(assigned_component='TCR_alpha', residue_count=180, atom_count=1440),
            21: Mock(assigned_component='TCR_beta', residue_count=240, atom_count=1920),
            22: Mock(assigned_component='HLA_alpha', residue_count=275, atom_count=2200)
        }
        mock_result.error_message = None

        mock_identify.return_value = mock_result

        # Create IndexManager and identify chains
        index_mgr = IndexManager(self.context)
        result = index_mgr.identify_chains_from_topology()

        # Verify result
        self.assertTrue(result.success)
        self.assertEqual(len(result.component_map), 5)
        self.assertIn('peptide', result.component_map)
        self.assertIn('HLA_alpha', result.component_map)

        # Verify metadata updated
        self.assertIn('topology_chain_mapping', self.context.metadata)

        # Verify chain standardization status
        self.assertTrue(index_mgr._chain_standardized)
        validation_result = index_mgr._chain_validation_result
        self.assertEqual(validation_result['method'], 'topology_based')
        self.assertTrue(validation_result['standardized'])
        self.assertTrue(validation_result['valid'])

    @patch('immunex.analysis.topology.topology_chain_identifier.TopologyChainIdentifier.identify_chains_from_topology')
    def test_identify_chains_from_topology_failure(self, mock_identify):
        """Test failed topology-based chain identification."""
        # Mock failed result
        mock_result = Mock(spec=TopologyChainMapping)
        mock_result.success = False
        mock_result.component_map = {}
        mock_result.chains = {}
        mock_result.error_message = "Failed to split chains: gmx command not found"

        mock_identify.return_value = mock_result

        # Create IndexManager and identify chains
        index_mgr = IndexManager(self.context)
        result = index_mgr.identify_chains_from_topology()

        # Verify result
        self.assertFalse(result.success)
        self.assertIsNotNone(result.error_message)

        # Verify chain standardization status
        self.assertFalse(index_mgr._chain_standardized)
        validation_result = index_mgr._chain_validation_result
        self.assertEqual(validation_result['method'], 'topology_based')
        self.assertFalse(validation_result['standardized'])

    def test_topology_chain_info_dataclass(self):
        """Test TopologyChainInfo dataclass."""
        chain_info = TopologyChainInfo(
            original_group_id=18,
            original_group_name='ch0_Protein',
            atom_count=72,
            residue_count=9,
            assigned_component='peptide',
            confidence=1.0
        )

        self.assertEqual(chain_info.original_group_id, 18)
        self.assertEqual(chain_info.assigned_component, 'peptide')
        self.assertEqual(chain_info.residue_count, 9)
        self.assertEqual(chain_info.confidence, 1.0)

    def test_topology_chain_mapping_dataclass(self):
        """Test TopologyChainMapping dataclass."""
        mapping = TopologyChainMapping(
            success=True,
            chains={
                18: Mock(assigned_component='peptide', residue_count=9)
            },
            component_map={'peptide': 18},
            error_message=None
        )

        self.assertTrue(mapping.success)
        self.assertEqual(len(mapping.chains), 1)
        self.assertIn('peptide', mapping.component_map)
        self.assertIsNone(mapping.error_message)

    def test_length_ranges_defined(self):
        """Test that TopologyChainIdentifier has correct length ranges."""
        identifier = TopologyChainIdentifier()

        # Verify LENGTH_RANGES attribute exists
        self.assertTrue(hasattr(identifier, 'LENGTH_RANGES'))

        length_ranges = identifier.LENGTH_RANGES

        # Verify all components defined
        self.assertIn('peptide', length_ranges)
        self.assertIn('beta2m', length_ranges)
        self.assertIn('TCR_alpha', length_ranges)
        self.assertIn('TCR_beta', length_ranges)
        self.assertIn('HLA_alpha', length_ranges)

        # Verify range structure
        peptide_range = length_ranges['peptide']
        self.assertIn('min_res', peptide_range)
        self.assertIn('max_res', peptide_range)
        self.assertIn('typical', peptide_range)

        # Verify logical ordering (shortest to longest)
        self.assertLess(length_ranges['peptide']['typical'],
                       length_ranges['beta2m']['typical'])
        self.assertLess(length_ranges['TCR_alpha']['typical'],
                       length_ranges['TCR_beta']['typical'])
        self.assertLess(length_ranges['TCR_beta']['typical'],
                       length_ranges['HLA_alpha']['typical'])


class TestTopologyIdentifierMethods(unittest.TestCase):
    """Test specific methods of TopologyChainIdentifier."""

    def setUp(self):
        """Set up identifier."""
        self.identifier = TopologyChainIdentifier()

    def test_identifier_initialization(self):
        """Test identifier initialization."""
        identifier = TopologyChainIdentifier(gmx_executable="gmx")
        self.assertEqual(identifier.gmx, "gmx")

        identifier_custom = TopologyChainIdentifier(gmx_executable="gmx_mpi")
        self.assertEqual(identifier_custom.gmx, "gmx_mpi")

    def test_get_component_selection_string_individual(self):
        """Test getting selection string for individual components."""
        # Mock mapping
        mock_mapping = Mock(spec=TopologyChainMapping)
        mock_mapping.success = True
        mock_mapping.component_map = {
            'peptide': 18,
            'HLA_alpha': 22,
            'TCR_alpha': 20
        }

        # Test individual component
        selection = self.identifier.get_component_selection_string('peptide', mock_mapping)
        self.assertEqual(selection, '18')

        selection = self.identifier.get_component_selection_string('HLA_alpha', mock_mapping)
        self.assertEqual(selection, '22')

    def test_get_component_selection_string_combined(self):
        """Test getting selection string for combined components."""
        # Mock mapping
        mock_mapping = Mock(spec=TopologyChainMapping)
        mock_mapping.success = True
        mock_mapping.component_map = {
            'HLA_alpha': 22,
            'beta2m': 19,
            'peptide': 18,
            'TCR_alpha': 20,
            'TCR_beta': 21
        }

        # Test pHLA (HLA_alpha + beta2m + peptide)
        selection = self.identifier.get_component_selection_string('pHLA', mock_mapping)
        # Should be "22 | 19 | 18" (order from component_map)
        self.assertIn('22', selection)
        self.assertIn('19', selection)
        self.assertIn('18', selection)
        self.assertIn('|', selection)

        # Test TCR (TCR_alpha + TCR_beta)
        selection = self.identifier.get_component_selection_string('TCR', mock_mapping)
        self.assertIn('20', selection)
        self.assertIn('21', selection)

        # Test HLA (HLA_alpha + beta2m)
        selection = self.identifier.get_component_selection_string('HLA', mock_mapping)
        self.assertIn('22', selection)
        self.assertIn('19', selection)

    def test_get_component_selection_string_not_found(self):
        """Test handling of non-existent component."""
        mock_mapping = Mock(spec=TopologyChainMapping)
        mock_mapping.success = True
        mock_mapping.component_map = {'peptide': 18}

        selection = self.identifier.get_component_selection_string('nonexistent', mock_mapping)
        self.assertIsNone(selection)

    def test_get_component_selection_string_failed_mapping(self):
        """Test handling of failed mapping."""
        mock_mapping = Mock(spec=TopologyChainMapping)
        mock_mapping.success = False

        selection = self.identifier.get_component_selection_string('peptide', mock_mapping)
        self.assertIsNone(selection)


def run_tests():
    """Run all tests."""
    print("=" * 70)
    print("Running Topology Integration Tests")
    print("=" * 70)

    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestTopologyIntegration))
    suite.addTests(loader.loadTestsFromTestCase(TestTopologyIdentifierMethods))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    print("\n" + "=" * 70)
    if result.wasSuccessful():
        print("All tests passed!")
        print("Topology-based chain identification is properly integrated.")
    else:
        print(f"Tests failed: {len(result.failures)} failures, {len(result.errors)} errors")
    print("=" * 70)

    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
