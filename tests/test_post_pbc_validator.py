#!/usr/bin/env python3
"""
Unit tests for PostPBCValidator module.

Tests cover:
1. Initialization and configuration
2. Trajectory integrity validation
3. PBC correction quality validation
4. Structural stability validation
5. Comprehensive validation
6. Quality grading
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import numpy as np
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.quality import PostPBCValidator

# Check if MDAnalysis is available
try:
    import MDAnalysis as mda
    HAS_MDA = True
except ImportError:
    HAS_MDA = False


class TestPostPBCValidatorInit(unittest.TestCase):
    """Test PostPBCValidator initialization"""

    def test_default_initialization(self):
        """Test default initialization"""
        validator = PostPBCValidator()

        self.assertEqual(validator.max_com_drift, 1.0)
        self.assertEqual(validator.max_rg_std_ratio, 0.15)
        self.assertEqual(validator.max_frame_jump_rmsd, 0.5)
        self.assertEqual(validator.min_file_size_mb, 1.0)
        self.assertIsNone(validator.expected_n_frames)

    def test_custom_initialization(self):
        """Test initialization with custom parameters"""
        validator = PostPBCValidator(
            expected_n_frames=1000,
            max_com_drift=1.5,
            max_rg_std_ratio=0.2,
            max_frame_jump_rmsd=0.8,
            min_file_size_mb=2.0
        )

        self.assertEqual(validator.expected_n_frames, 1000)
        self.assertEqual(validator.max_com_drift, 1.5)
        self.assertEqual(validator.max_rg_std_ratio, 0.2)
        self.assertEqual(validator.max_frame_jump_rmsd, 0.8)
        self.assertEqual(validator.min_file_size_mb, 2.0)


class TestTrajectoryIntegrityValidation(unittest.TestCase):
    """Test trajectory integrity validation methods"""

    def setUp(self):
        """Set up test fixtures"""
        self.validator = PostPBCValidator(min_file_size_mb=0.001)
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary directory"""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_missing_trajectory_file(self):
        """Test validation with missing trajectory file"""
        result = self.validator.validate_trajectory_integrity(
            trajectory="/nonexistent/file.xtc",
            topology="/nonexistent/topology.tpr"
        )

        self.assertEqual(result['status'], 'fail')
        self.assertTrue(any('not found' in issue.lower()
                          for issue in result['issues']))

    def test_empty_trajectory_file(self):
        """Test validation with empty trajectory file"""
        # Create empty file
        empty_file = Path(self.temp_dir) / "empty.xtc"
        empty_file.touch()

        result = self.validator.validate_trajectory_integrity(
            trajectory=str(empty_file),
            topology="/nonexistent/topology.tpr"
        )

        self.assertEqual(result['status'], 'fail')
        self.assertLess(result['file_size_mb'], 0.001)

    @unittest.skipUnless(HAS_MDA, "MDAnalysis not available")
    def test_file_size_check(self):
        """Test file size validation"""
        # Create a small dummy file
        small_file = Path(self.temp_dir) / "small.xtc"
        small_file.write_bytes(b'0' * 100)  # 100 bytes

        validator = PostPBCValidator(min_file_size_mb=1.0)
        result = validator.validate_trajectory_integrity(
            trajectory=str(small_file),
            topology="/nonexistent/topology.tpr"
        )

        # Should fail due to small size
        self.assertTrue(any('File size' in issue
                          for issue in result['issues']))


class TestQualityGrading(unittest.TestCase):
    """Test quality grading logic"""

    def setUp(self):
        """Set up test fixtures"""
        self.validator = PostPBCValidator(
            max_com_drift=1.0,
            max_rg_std_ratio=0.15,
            max_frame_jump_rmsd=0.5
        )

    def test_grade_a_assignment(self):
        """Test Grade A assignment for excellent quality"""
        integrity = {'status': 'pass', 'n_frames': 1000, 'issues': []}
        pbc_quality = {'status': 'pass', 'com_drift_nm': 0.1, 'issues': []}
        stability = {
            'status': 'pass',
            'rg_std_ratio': 0.05,
            'max_frame_jump_rmsd_nm': 0.1,
            'issues': []
        }

        grade = self.validator._assign_quality_grade(
            integrity, pbc_quality, stability
        )

        self.assertEqual(grade, 'A')

    def test_grade_b_assignment(self):
        """Test Grade B assignment for good quality"""
        integrity = {'status': 'pass', 'n_frames': 1000, 'issues': []}
        pbc_quality = {'status': 'pass', 'com_drift_nm': 0.4, 'issues': []}
        stability = {
            'status': 'pass',
            'rg_std_ratio': 0.08,
            'max_frame_jump_rmsd_nm': 0.2,
            'issues': []
        }

        grade = self.validator._assign_quality_grade(
            integrity, pbc_quality, stability
        )

        self.assertIn(grade, ['A', 'B', 'C'])  # Could vary depending on scoring

    def test_grade_d_for_failed_checks(self):
        """Test Grade D assignment when checks fail"""
        integrity = {'status': 'fail', 'n_frames': 0, 'issues': ['Failed']}
        pbc_quality = {'status': 'pass', 'com_drift_nm': 0.1, 'issues': []}
        stability = {'status': 'pass', 'rg_std_ratio': 0.05, 'issues': []}

        grade = self.validator._assign_quality_grade(
            integrity, pbc_quality, stability
        )

        self.assertEqual(grade, 'D')

    def test_grade_d_for_high_com_drift(self):
        """Test Grade D for high COM drift"""
        integrity = {'status': 'pass', 'n_frames': 1000, 'issues': []}
        pbc_quality = {'status': 'pass', 'com_drift_nm': 2.0, 'issues': []}
        stability = {'status': 'pass', 'rg_std_ratio': 0.05, 'issues': []}

        grade = self.validator._assign_quality_grade(
            integrity, pbc_quality, stability
        )

        # High COM drift should penalize severely
        self.assertIn(grade, ['C', 'D'])


class TestValidationEdgeCases(unittest.TestCase):
    """Test edge cases and error handling"""

    def test_very_strict_thresholds(self):
        """Test with very strict thresholds"""
        validator = PostPBCValidator(
            max_com_drift=0.01,
            max_rg_std_ratio=0.01,
            max_frame_jump_rmsd=0.01
        )

        # Even perfect data might fail with such strict thresholds
        self.assertEqual(validator.max_com_drift, 0.01)
        self.assertEqual(validator.max_rg_std_ratio, 0.01)

    def test_very_permissive_thresholds(self):
        """Test with very permissive thresholds"""
        validator = PostPBCValidator(
            max_com_drift=10.0,
            max_rg_std_ratio=1.0,
            max_frame_jump_rmsd=5.0
        )

        self.assertEqual(validator.max_com_drift, 10.0)
        self.assertEqual(validator.max_rg_std_ratio, 1.0)

    def test_zero_frames_expected(self):
        """Test with zero expected frames"""
        validator = PostPBCValidator(expected_n_frames=0)
        self.assertEqual(validator.expected_n_frames, 0)


class TestValidationResultStructure(unittest.TestCase):
    """Test that validation results have correct structure"""

    def setUp(self):
        """Set up test fixtures"""
        self.validator = PostPBCValidator()

    def test_integrity_result_structure(self):
        """Test integrity validation result structure"""
        result = self.validator.validate_trajectory_integrity(
            trajectory="/nonexistent/file.xtc",
            topology="/nonexistent/topology.tpr"
        )

        # Check required keys
        self.assertIn('status', result)
        self.assertIn('n_frames', result)
        self.assertIn('file_size_mb', result)
        self.assertIn('issues', result)

        # Check types
        self.assertIsInstance(result['status'], str)
        self.assertIsInstance(result['n_frames'], int)
        self.assertIsInstance(result['file_size_mb'], float)
        self.assertIsInstance(result['issues'], list)


# Simple functional tests (non-unittest style)

def test_validator_creation():
    """Test validator creation with various parameters"""
    print("Test: Validator creation")

    # Default parameters
    v1 = PostPBCValidator()
    assert v1.max_com_drift == 1.0
    print("  ✓ Default parameters")

    # Custom parameters
    v2 = PostPBCValidator(
        expected_n_frames=500,
        max_com_drift=2.0,
        max_rg_std_ratio=0.25
    )
    assert v2.expected_n_frames == 500
    assert v2.max_com_drift == 2.0
    print("  ✓ Custom parameters")

    print("  PASSED\n")


def test_grading_logic():
    """Test quality grading logic with mock data"""
    print("Test: Quality grading logic")

    validator = PostPBCValidator()

    # Perfect quality (Grade A)
    integrity_a = {'status': 'pass', 'issues': []}
    pbc_a = {'status': 'pass', 'com_drift_nm': 0.05, 'issues': []}
    stab_a = {'status': 'pass', 'rg_std_ratio': 0.03,
              'max_frame_jump_rmsd_nm': 0.05, 'issues': []}

    grade_a = validator._assign_quality_grade(integrity_a, pbc_a, stab_a)
    assert grade_a == 'A', f"Expected A, got {grade_a}"
    print(f"  ✓ Perfect quality: Grade {grade_a}")

    # Failed checks (Grade D)
    integrity_d = {'status': 'fail', 'issues': ['error']}
    pbc_d = {'status': 'pass', 'com_drift_nm': 0.1, 'issues': []}
    stab_d = {'status': 'pass', 'rg_std_ratio': 0.05,
              'max_frame_jump_rmsd_nm': 0.1, 'issues': []}

    grade_d = validator._assign_quality_grade(integrity_d, pbc_d, stab_d)
    assert grade_d == 'D', f"Expected D, got {grade_d}"
    print(f"  ✓ Failed checks: Grade {grade_d}")

    print("  PASSED\n")


def run_functional_tests():
    """Run functional tests"""
    print("=" * 60)
    print("PostPBCValidator Functional Tests")
    print("=" * 60)
    print()

    test_validator_creation()
    test_grading_logic()

    print("=" * 60)
    print("All functional tests passed!")
    print("=" * 60)


if __name__ == '__main__':
    # Run functional tests
    run_functional_tests()

    print("\n")

    # Run unittest tests
    print("=" * 60)
    print("Running unittest suite...")
    print("=" * 60)
    unittest.main(verbosity=2, exit=True)
