#!/usr/bin/env python3
"""
Unit tests for IndexManager

This tests the core functionality of the unified IndexManager without
requiring actual MD data files.

Author: Immunex Development Team
Date: 2026-03-17
"""

import sys
import unittest
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PipelineContext
from immunex.analysis import IndexManager, StandardComponent, GroupInfo


class TestStandardComponent(unittest.TestCase):
    """Test StandardComponent Enum."""

    def test_component_chains(self):
        """Test that components have correct chain definitions."""
        # Static components
        self.assertEqual(StandardComponent.HLA.chains, ['A', 'B'])
        self.assertEqual(StandardComponent.PHLA.chains, ['A', 'B', 'C'])
        self.assertEqual(StandardComponent.PEPTIDE.chains, ['C'])
        self.assertEqual(StandardComponent.TCR.chains, ['D', 'E'])
        self.assertEqual(StandardComponent.TCR_ALPHA.chains, ['D'])
        self.assertEqual(StandardComponent.TCR_BETA.chains, ['E'])

        # Dynamic components - CDR loops
        self.assertEqual(StandardComponent.CDR1_ALPHA.chains, ['D'])
        self.assertEqual(StandardComponent.CDR2_ALPHA.chains, ['D'])
        self.assertEqual(StandardComponent.CDR3_ALPHA.chains, ['D'])
        self.assertEqual(StandardComponent.CDR1_BETA.chains, ['E'])
        self.assertEqual(StandardComponent.CDR2_BETA.chains, ['E'])
        self.assertEqual(StandardComponent.CDR3_BETA.chains, ['E'])

        # Dynamic components - Interface
        self.assertEqual(StandardComponent.INTERFACE.chains, ['A', 'B', 'C', 'D', 'E'])

    def test_component_names(self):
        """Test component name formatting."""
        # Static components
        self.assertEqual(StandardComponent.PHLA.component_name, 'pHLA')
        self.assertEqual(StandardComponent.HLA.component_name, 'HLA')
        self.assertEqual(StandardComponent.HLA_ALPHA.component_name, 'HLA_alpha')
        self.assertEqual(StandardComponent.TCR_BETA.component_name, 'TCR_beta')

        # Dynamic components
        self.assertEqual(StandardComponent.CDR3_ALPHA.component_name, 'CDR3_alpha')
        self.assertEqual(StandardComponent.CDR1_BETA.component_name, 'CDR1_beta')
        self.assertEqual(StandardComponent.INTERFACE.component_name, 'INTERFACE')

    def test_enum_iteration(self):
        """Test that we can iterate over all components."""
        components = list(StandardComponent)
        # 8 static + 6 CDR + 1 Interface = 15 total
        self.assertEqual(len(components), 15)

        component_names = [c.component_name for c in components]
        # Static components
        self.assertIn('pHLA', component_names)
        self.assertIn('TCR', component_names)
        # Dynamic components
        self.assertIn('CDR3_alpha', component_names)
        self.assertIn('CDR3_beta', component_names)
        self.assertIn('INTERFACE', component_names)

    def test_component_dynamic_properties(self):
        """Test dynamic component detection."""
        # Static components should not be dynamic
        self.assertFalse(StandardComponent.HLA.is_dynamic)
        self.assertFalse(StandardComponent.PHLA.is_dynamic)
        self.assertFalse(StandardComponent.TCR.is_dynamic)

        # CDR components should be dynamic and require sequence
        self.assertTrue(StandardComponent.CDR3_ALPHA.is_dynamic)
        self.assertTrue(StandardComponent.CDR3_ALPHA.requires_sequence)
        self.assertFalse(StandardComponent.CDR3_ALPHA.requires_distance_calculation)

        self.assertTrue(StandardComponent.CDR1_BETA.is_dynamic)
        self.assertTrue(StandardComponent.CDR1_BETA.requires_sequence)

        # Interface should be dynamic and require distance calculation
        self.assertTrue(StandardComponent.INTERFACE.is_dynamic)
        self.assertFalse(StandardComponent.INTERFACE.requires_sequence)
        self.assertTrue(StandardComponent.INTERFACE.requires_distance_calculation)


class TestGroupInfo(unittest.TestCase):
    """Test GroupInfo dataclass."""

    def test_group_info_creation(self):
        """Test GroupInfo creation."""
        group = GroupInfo(
            group_id=20,
            group_name='pHLA',
            atom_count=1500,
            component=StandardComponent.PHLA
        )

        self.assertEqual(group.group_id, 20)
        self.assertEqual(group.group_name, 'pHLA')
        self.assertEqual(group.atom_count, 1500)
        self.assertEqual(group.component, StandardComponent.PHLA)

    def test_group_info_without_component(self):
        """Test GroupInfo without component association."""
        group = GroupInfo(
            group_id=25,
            group_name='CDR3_TCR_beta',
            atom_count=15
        )

        self.assertEqual(group.group_id, 25)
        self.assertIsNone(group.component)


class TestIndexManager(unittest.TestCase):
    """Test IndexManager core functionality."""

    def setUp(self):
        """Set up test context."""
        self.context = PipelineContext(
            system_id="test_system",
            topology="test.tpr",
            trajectory_raw="test.xtc",
            structure_pdb="test.pdb",
            output_dir="/tmp/test_output"
        )

    def test_index_manager_initialization(self):
        """Test IndexManager initialization."""
        index_mgr = IndexManager(self.context)

        self.assertEqual(index_mgr.context.system_id, "test_system")
        self.assertIsNone(index_mgr.base_index_file)
        self.assertEqual(len(index_mgr.group_registry), 0)
        self.assertEqual(len(index_mgr.dynamic_index_files), 0)

    def test_lazy_initialization_via_context(self):
        """Test lazy initialization through PipelineContext."""
        # IndexManager should be created lazily
        self.assertIsNone(self.context._index_manager)

        # Access through property
        index_mgr = self.context.index_manager

        # Should be initialized now
        self.assertIsNotNone(index_mgr)
        self.assertIsInstance(index_mgr, IndexManager)

        # Second access should return same instance
        index_mgr2 = self.context.index_manager
        self.assertIs(index_mgr, index_mgr2)

    def test_parse_ndx_content(self):
        """Test parsing of .ndx file content."""
        # Create a test index manager
        index_mgr = self.context.index_manager

        # Simulate parsing ndx content
        test_ndx_content = """[ System ]
   1    2    3    4    5

[ Protein ]
   1    2    3

[ pHLA ]
  10   20   30   40   50

[ TCR ]
  60   70   80
"""
        # Create temporary ndx file
        test_ndx_file = Path("/tmp/test_index.ndx")
        with open(test_ndx_file, 'w') as f:
            f.write(test_ndx_content)

        # Parse the file
        index_mgr._parse_index_file(test_ndx_file)

        # Verify groups were registered
        self.assertIn('System', index_mgr.group_registry)
        self.assertIn('Protein', index_mgr.group_registry)
        self.assertIn('pHLA', index_mgr.group_registry)
        self.assertIn('TCR', index_mgr.group_registry)

        # Verify group IDs (sequential)
        self.assertEqual(index_mgr.group_registry['System'].group_id, 0)
        self.assertEqual(index_mgr.group_registry['Protein'].group_id, 1)
        self.assertEqual(index_mgr.group_registry['pHLA'].group_id, 2)
        self.assertEqual(index_mgr.group_registry['TCR'].group_id, 3)

        # Verify atom counts
        self.assertEqual(index_mgr.group_registry['System'].atom_count, 5)
        self.assertEqual(index_mgr.group_registry['pHLA'].atom_count, 5)
        self.assertEqual(index_mgr.group_registry['TCR'].atom_count, 3)

        # Clean up
        test_ndx_file.unlink()

    def test_group_registry_lookup(self):
        """Test group registry lookup after manual registration."""
        index_mgr = self.context.index_manager

        # Manually register some groups
        index_mgr._register_group(20, 'pHLA', 1500)
        index_mgr._register_group(21, 'TCR', 800)
        index_mgr._register_group(22, 'peptide', 100)

        # Test get_group_id
        self.assertEqual(index_mgr.group_registry['pHLA'].group_id, 20)
        self.assertEqual(index_mgr.group_registry['TCR'].group_id, 21)
        self.assertEqual(index_mgr.group_registry['peptide'].group_id, 22)

        # Test get_group_info
        phla_info = index_mgr.group_registry['pHLA']
        self.assertEqual(phla_info.group_name, 'pHLA')
        self.assertEqual(phla_info.atom_count, 1500)

    def test_repr(self):
        """Test string representation."""
        index_mgr = self.context.index_manager
        repr_str = repr(index_mgr)

        self.assertIn('IndexManager', repr_str)
        self.assertIn('test_system', repr_str)

    def test_chain_validation_methods_exist(self):
        """Test that chain validation methods exist."""
        index_mgr = self.context.index_manager

        # Test methods exist and are callable
        self.assertTrue(hasattr(index_mgr, 'get_chain_validation_report'))
        self.assertTrue(callable(index_mgr.get_chain_validation_report))

        self.assertTrue(hasattr(index_mgr, 'is_chain_standardized'))
        self.assertTrue(callable(index_mgr.is_chain_standardized))

        # Test initial state
        self.assertIsNone(index_mgr._chain_standardized)
        self.assertIsNone(index_mgr._chain_validation_result)


def run_tests():
    """Run all tests."""
    print("=" * 70)
    print("Running IndexManager Unit Tests")
    print("=" * 70)

    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestStandardComponent))
    suite.addTests(loader.loadTestsFromTestCase(TestGroupInfo))
    suite.addTests(loader.loadTestsFromTestCase(TestIndexManager))

    # Run tests
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
