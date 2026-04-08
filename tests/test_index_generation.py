#!/usr/bin/env python3
"""
Unit tests for standardized IndexGenerator

Tests the index generation module following 6 design principles.

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
import unittest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import tempfile

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.topology import (
    IndexGenerator,
    IndexGenerationInput,
    IndexGenerationResult,
    IndexGenerationMethod,
    ComponentDefinition,
    ComponentIndexInfo,
    InvalidInputError,
    TopologyFileError,
    GROMACSCommandError,
    SideEffectTracker
)


class TestErrorDefinitions(unittest.TestCase):
    """Test error definitions (Principle 4: Clear Errors)."""

    def test_invalid_input_error(self):
        """Test InvalidInputError creation."""
        error = InvalidInputError('topology', '/invalid/path', 'File does not exist')
        
        self.assertEqual(error.param_name, 'topology')
        self.assertEqual(error.value, '/invalid/path')
        self.assertEqual(error.reason, 'File does not exist')
        self.assertIn('topology', str(error))
        self.assertIn('File does not exist', str(error))

    def test_gromacs_command_error(self):
        """Test GROMACSCommandError with suggestion."""
        error = GROMACSCommandError(
            'gmx make_ndx',
            'Chain A not found',
            'Check chain IDs in PDB file'
        )
        
        self.assertEqual(error.command, 'gmx make_ndx')
        self.assertEqual(error.stderr, 'Chain A not found')
        self.assertIn('Suggestion', str(error))


class TestInputDefinitions(unittest.TestCase):
    """Test input definitions (Principle 1: Clear Inputs)."""

    def setUp(self):
        """Create temporary topology file for testing."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_topology = Path(self.temp_dir) / "test.pdb"
        self.temp_topology.touch()

    def test_component_definition(self):
        """Test ComponentDefinition creation."""
        comp = ComponentDefinition(
            name="pHLA",
            chains=['A', 'B', 'C'],
            ca_only=False
        )
        
        self.assertEqual(comp.name, "pHLA")
        self.assertEqual(comp.chains, ['A', 'B', 'C'])
        self.assertFalse(comp.ca_only)

    def test_input_validation_success(self):
        """Test successful input validation."""
        input_params = IndexGenerationInput(
            topology=str(self.temp_topology),
            method=IndexGenerationMethod.PDB_BASED,
            components=[
                ComponentDefinition(name="pHLA", chains=['A', 'B', 'C'])
            ],
            output_file=str(Path(self.temp_dir) / "test.ndx")
        )
        
        # Should not raise
        input_params.validate()

    def test_input_validation_missing_topology(self):
        """Test validation fails for missing topology."""
        input_params = IndexGenerationInput(
            topology="/nonexistent/file.pdb",
            method=IndexGenerationMethod.PDB_BASED,
            components=[ComponentDefinition(name="pHLA", chains=['A'])],
            output_file="test.ndx"
        )
        
        with self.assertRaises(InvalidInputError) as ctx:
            input_params.validate()
        
        self.assertEqual(ctx.exception.param_name, 'topology')
        self.assertIn('does not exist', ctx.exception.reason)

    def test_input_validation_invalid_extension(self):
        """Test validation fails for invalid topology extension."""
        temp_invalid = Path(self.temp_dir) / "test.txt"
        temp_invalid.touch()
        
        input_params = IndexGenerationInput(
            topology=str(temp_invalid),
            method=IndexGenerationMethod.PDB_BASED,
            components=[ComponentDefinition(name="pHLA", chains=['A'])],
            output_file="test.ndx"
        )
        
        with self.assertRaises(InvalidInputError) as ctx:
            input_params.validate()
        
        self.assertEqual(ctx.exception.param_name, 'topology')
        self.assertIn('Invalid extension', ctx.exception.reason)

    def test_input_validation_empty_components(self):
        """Test validation fails for empty components."""
        input_params = IndexGenerationInput(
            topology=str(self.temp_topology),
            method=IndexGenerationMethod.PDB_BASED,
            components=[],
            output_file="test.ndx"
        )
        
        with self.assertRaises(InvalidInputError) as ctx:
            input_params.validate()
        
        self.assertEqual(ctx.exception.param_name, 'components')

    def test_input_validation_pdb_based_requires_chains(self):
        """Test PDB-based method requires chains."""
        input_params = IndexGenerationInput(
            topology=str(self.temp_topology),
            method=IndexGenerationMethod.PDB_BASED,
            components=[ComponentDefinition(name="pHLA")],  # Missing chains
            output_file="test.ndx"
        )
        
        with self.assertRaises(InvalidInputError):
            input_params.validate()


class TestOutputDefinitions(unittest.TestCase):
    """Test output definitions (Principle 2: Clear Outputs)."""

    def test_component_index_info(self):
        """Test ComponentIndexInfo creation."""
        info = ComponentIndexInfo(
            name="pHLA",
            group_id=18,
            atom_count=1500,
            residue_count=200,
            chains=['A', 'B', 'C']
        )
        
        self.assertEqual(info.name, "pHLA")
        self.assertEqual(info.group_id, 18)
        self.assertEqual(info.atom_count, 1500)

    def test_generation_result_success(self):
        """Test successful result creation."""
        result = IndexGenerationResult(
            success=True,
            output_file="test.ndx",
            components=[
                ComponentIndexInfo(name="pHLA", group_id=18, atom_count=1500)
            ],
            processing_stats={'n_components': 1, 'processing_time_sec': 1.2}
        )
        
        self.assertTrue(result.success)
        self.assertEqual(result.output_file, "test.ndx")
        self.assertEqual(len(result.components), 1)
        self.assertIsNone(result.error_message)

    def test_generation_result_failure(self):
        """Test failed result creation."""
        result = IndexGenerationResult(
            success=False,
            error_message="Topology file not found"
        )
        
        self.assertFalse(result.success)
        self.assertIsNone(result.output_file)
        self.assertEqual(result.error_message, "Topology file not found")

    def test_result_to_dict(self):
        """Test result serialization to dict."""
        result = IndexGenerationResult(
            success=True,
            output_file="test.ndx",
            components=[
                ComponentIndexInfo(name="pHLA", group_id=18, atom_count=1500)
            ]
        )
        
        result_dict = result.to_dict()
        
        self.assertIsInstance(result_dict, dict)
        self.assertTrue(result_dict['success'])
        self.assertEqual(result_dict['output_file'], "test.ndx")
        self.assertEqual(len(result_dict['components']), 1)


class TestSideEffectTracker(unittest.TestCase):
    """Test side effect tracking (Principle 3: Clear Side Effects)."""

    def test_track_file_creation(self):
        """Test file creation tracking."""
        tracker = SideEffectTracker()
        tracker.track_file_creation("/tmp/test.ndx")
        
        self.assertEqual(len(tracker.files_created), 1)
        self.assertEqual(tracker.files_created[0], "/tmp/test.ndx")

    def test_track_command(self):
        """Test command execution tracking."""
        tracker = SideEffectTracker()
        tracker.track_command(['gmx', 'make_ndx'], 0, "")
        
        self.assertEqual(len(tracker.commands_executed), 1)
        self.assertIn('gmx make_ndx', tracker.commands_executed[0]['command'])
        self.assertEqual(tracker.commands_executed[0]['returncode'], 0)

    def test_get_summary(self):
        """Test side effect summary."""
        tracker = SideEffectTracker()
        tracker.track_file_creation("/tmp/test1.ndx")
        tracker.track_file_creation("/tmp/test2.ndx")
        tracker.track_file_deletion("/tmp/temp.ndx")
        tracker.track_command(['gmx', 'make_ndx'], 0)
        
        summary = tracker.get_summary()
        
        self.assertEqual(summary['files_created'], 2)
        self.assertEqual(summary['files_deleted'], 1)
        self.assertEqual(summary['commands_executed'], 1)


class TestIndexGenerator(unittest.TestCase):
    """Test IndexGenerator core functionality (Principles 5 & 6: Testable & Schedulable)."""

    def setUp(self):
        """Set up test fixtures."""
        self.generator = IndexGenerator()
        self.temp_dir = tempfile.mkdtemp()
        self.temp_topology = Path(self.temp_dir) / "test.pdb"
        
        # Create minimal PDB file
        with open(self.temp_topology, 'w') as f:
            f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
            f.write("END\n")

    def test_initialization(self):
        """Test generator initialization."""
        self.assertIsNotNone(self.generator)
        self.assertIsNone(self.generator.progress_callback)
        self.assertFalse(self.generator.is_cancelled())

    def test_progress_callback(self):
        """Test progress callback (Principle 6: Schedulable)."""
        progress_calls = []
        
        def callback(progress, message):
            progress_calls.append((progress, message))
        
        self.generator.set_progress_callback(callback)
        self.generator._report_progress(0.5, "Test message")
        
        self.assertEqual(len(progress_calls), 1)
        self.assertEqual(progress_calls[0][0], 0.5)
        self.assertEqual(progress_calls[0][1], "Test message")

    def test_cancellation(self):
        """Test operation cancellation (Principle 6: Schedulable)."""
        self.assertFalse(self.generator.is_cancelled())
        
        self.generator.cancel()
        
        self.assertTrue(self.generator.is_cancelled())

    def test_generate_invalid_input(self):
        """Test generation with invalid input."""
        input_params = IndexGenerationInput(
            topology="/nonexistent/file.pdb",
            method=IndexGenerationMethod.PDB_BASED,
            components=[ComponentDefinition(name="pHLA", chains=['A'])],
            output_file="test.ndx"
        )
        
        result = self.generator.generate(input_params)
        
        self.assertFalse(result.success)
        self.assertIsNotNone(result.error_message)
        self.assertIn('does not exist', result.error_message)

    def test_generate_cancelled(self):
        """Test generation when cancelled."""
        input_params = IndexGenerationInput(
            topology=str(self.temp_topology),
            method=IndexGenerationMethod.PDB_BASED,
            components=[ComponentDefinition(name="pHLA", chains=['A'])],
            output_file=str(Path(self.temp_dir) / "test.ndx")
        )
        
        # Cancel before generation
        self.generator.cancel()
        
        result = self.generator.generate(input_params)
        
        self.assertFalse(result.success)
        self.assertIn('cancelled', result.error_message.lower())


def run_tests():
    """Run all tests."""
    print("=" * 70)
    print("Running Standardized IndexGenerator Unit Tests")
    print("=" * 70)

    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestErrorDefinitions))
    suite.addTests(loader.loadTestsFromTestCase(TestInputDefinitions))
    suite.addTests(loader.loadTestsFromTestCase(TestOutputDefinitions))
    suite.addTests(loader.loadTestsFromTestCase(TestSideEffectTracker))
    suite.addTests(loader.loadTestsFromTestCase(TestIndexGenerator))

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
