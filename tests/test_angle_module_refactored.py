#!/usr/bin/env python3
"""
Test suite for refactored angle analysis module (without pytest).

Tests the new simplified API with only 3 public classes:
- DockingAngleAnalyzer
- DockingAngleInput
- DockingAngleResult

Author: Immunex Development Team
Date: 2026-03-19
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_module_imports():
    """Test that all 3 public classes can be imported."""
    print("\n" + "=" * 70)
    print("Test 1: Module Imports")
    print("=" * 70)

    try:
        from immunex.analysis.angles import (
            DockingAngleAnalyzer,
            DockingAngleInput,
            DockingAngleResult
        )
        print("✓ All 3 public classes imported successfully")
        print(f"  - DockingAngleAnalyzer: {DockingAngleAnalyzer}")
        print(f"  - DockingAngleInput: {DockingAngleInput}")
        print(f"  - DockingAngleResult: {DockingAngleResult}")
        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        return False


def test_public_api_count():
    """Test that only 3 classes are exported."""
    print("\n" + "=" * 70)
    print("Test 2: Public API Count")
    print("=" * 70)

    try:
        from immunex.analysis import angles

        exported = angles.__all__
        print(f"Exported classes: {exported}")

        if len(exported) == 3:
            print(f"✓ Correct number of exported classes: {len(exported)}")
            return True
        else:
            print(f"✗ Expected 3 classes, got {len(exported)}")
            return False
    except Exception as e:
        print(f"✗ Test failed: {e}")
        return False


def test_internal_classes_hidden():
    """Test that internal classes are not accessible."""
    print("\n" + "=" * 70)
    print("Test 3: Internal Classes Hidden")
    print("=" * 70)

    internal_classes = [
        '_GeometryUtils',
        '_SequenceAligner',
        '_MHCGrooveDetector',
        '_TCRAxisCalculator'
    ]

    all_hidden = True
    for class_name in internal_classes:
        try:
            exec(f"from immunex.analysis.angles import {class_name}")
            print(f"✗ {class_name} should not be accessible")
            all_hidden = False
        except ImportError:
            print(f"✓ {class_name} correctly hidden")

    return all_hidden


def test_old_api_removed():
    """Test that old API classes are no longer accessible."""
    print("\n" + "=" * 70)
    print("Test 4: Old API Removed")
    print("=" * 70)

    old_classes = [
        'PlaneFitter',
        'MHCSequenceAligner',
        'DockingAnglePrimaryAnalyzer',
        'ConservedCysteineDetector'
    ]

    all_removed = True
    for class_name in old_classes:
        try:
            exec(f"from immunex.analysis.angles import {class_name}")
            print(f"✗ {class_name} should be removed")
            all_removed = False
        except ImportError:
            print(f"✓ {class_name} correctly removed")

    return all_removed


def test_docking_angle_input_creation():
    """Test DockingAngleInput creation."""
    print("\n" + "=" * 70)
    print("Test 5: DockingAngleInput Creation")
    print("=" * 70)

    try:
        from immunex.analysis.angles import DockingAngleInput

        # Minimal parameters
        input1 = DockingAngleInput(topology="test.pdb")
        print(f"✓ Minimal input created: topology={input1.topology}")

        # Full parameters
        input2 = DockingAngleInput(
            topology="md.tpr",
            trajectory="md.xtc",
            mhc_selection="chainID A",
            tcr_alpha_selection="chainID D",
            tcr_beta_selection="chainID E",
            stride=10,
            output_dir="./results"
        )
        print(f"✓ Full input created with stride={input2.stride}")

        # Test to_dict
        input_dict = input1.to_dict()
        print(f"✓ Input.to_dict() works: {list(input_dict.keys())[:3]}...")

        return True
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_docking_angle_result_creation():
    """Test DockingAngleResult creation."""
    print("\n" + "=" * 70)
    print("Test 6: DockingAngleResult Creation")
    print("=" * 70)

    try:
        from immunex.analysis.angles import DockingAngleResult
        import numpy as np

        # Success result with trajectory data
        result1 = DockingAngleResult(
            success=True,
            times=np.array([0, 100, 200]),
            crossing_angles=np.array([45.0, 46.0, 44.5]),
            incident_angles=np.array([60.0, 61.0, 59.5]),
            statistics={'crossing_mean': 45.17}
        )
        print(f"✓ Success result created: {len(result1.times)} frames")

        # Failure result
        result2 = DockingAngleResult(
            success=False,
            error_message="Test error"
        )
        print(f"✓ Failure result created: {result2.error_message}")

        # Single frame result
        result3 = DockingAngleResult(
            success=True,
            crossing_angle=45.5,
            incident_angle=60.2
        )
        print(f"✓ Single frame result: crossing={result3.crossing_angle}°")

        # Test to_dict and get_summary
        result_dict = result1.to_dict()
        summary = result1.get_summary()
        print(f"✓ Result.to_dict() and get_summary() work")

        return True
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_docking_angle_analyzer_creation():
    """Test DockingAngleAnalyzer creation and methods."""
    print("\n" + "=" * 70)
    print("Test 7: DockingAngleAnalyzer Basic Functionality")
    print("=" * 70)

    try:
        from immunex.analysis.angles import DockingAngleAnalyzer

        # Create analyzer
        analyzer = DockingAngleAnalyzer()
        print(f"✓ Analyzer created: {analyzer}")

        # Test methods exist
        assert hasattr(analyzer, 'analyze')
        assert hasattr(analyzer, 'set_progress_callback')
        assert hasattr(analyzer, 'cancel')
        assert hasattr(analyzer, 'is_cancelled')
        print("✓ All required methods exist")

        # Test progress callback
        callback_called = []
        def test_callback(progress, message):
            callback_called.append((progress, message))

        analyzer.set_progress_callback(test_callback)
        print("✓ Progress callback set successfully")

        # Test cancellation
        assert analyzer.is_cancelled() is False
        analyzer.cancel()
        assert analyzer.is_cancelled() is True
        print("✓ Cancellation mechanism works")

        return True
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_analyzer_with_invalid_input():
    """Test analyzer with invalid input (should handle gracefully)."""
    print("\n" + "=" * 70)
    print("Test 8: Analyzer Error Handling")
    print("=" * 70)

    try:
        from immunex.analysis.angles import (
            DockingAngleAnalyzer,
            DockingAngleInput
        )

        analyzer = DockingAngleAnalyzer()

        # Test with nonexistent file
        input_params = DockingAngleInput(topology="nonexistent_file.pdb")
        result = analyzer.analyze(input_params)

        if not result.success:
            print(f"✓ Analyzer correctly handles invalid input")
            print(f"  Error message: {result.error_message[:60]}...")
            return True
        else:
            print("✗ Analyzer should return failure for invalid input")
            return False

    except Exception as e:
        print(f"✗ Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_analyzer_with_cancellation():
    """Test analyzer respects cancellation."""
    print("\n" + "=" * 70)
    print("Test 9: Analyzer Cancellation Handling")
    print("=" * 70)

    try:
        from immunex.analysis.angles import (
            DockingAngleAnalyzer,
            DockingAngleInput
        )

        analyzer = DockingAngleAnalyzer()
        analyzer.cancel()  # Cancel before starting

        input_params = DockingAngleInput(topology="test.pdb")
        result = analyzer.analyze(input_params)

        if not result.success and 'cancel' in result.error_message.lower():
            print(f"✓ Analyzer correctly handles cancellation")
            print(f"  Error message: {result.error_message}")
            return True
        else:
            print("✗ Analyzer should return cancellation error")
            return False

    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_typical_usage_pattern():
    """Test the typical recommended usage pattern."""
    print("\n" + "=" * 70)
    print("Test 10: Typical Usage Pattern")
    print("=" * 70)

    try:
        from immunex.analysis.angles import (
            DockingAngleAnalyzer,
            DockingAngleInput
        )

        # This is the recommended usage pattern
        print("Testing pattern:")
        print("  analyzer = DockingAngleAnalyzer()")
        print("  input_params = DockingAngleInput(...)")
        print("  result = analyzer.analyze(input_params)")

        analyzer = DockingAngleAnalyzer()

        # Track progress
        progress_updates = []
        def callback(progress, message):
            progress_updates.append((progress, message))

        analyzer.set_progress_callback(callback)

        input_params = DockingAngleInput(
            topology="test.pdb",
            stride=10
        )

        result = analyzer.analyze(input_params)

        # Check result structure
        assert hasattr(result, 'success')
        assert hasattr(result, 'error_message')
        assert isinstance(result.success, bool)

        print(f"✓ Usage pattern works correctly")
        print(f"  - Result.success: {result.success}")
        print(f"  - Progress updates received: {len(progress_updates)}")
        if progress_updates:
            print(f"  - First progress: {progress_updates[0][1]}")

        return True

    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_module_version():
    """Test module version info."""
    print("\n" + "=" * 70)
    print("Test 11: Module Version Info")
    print("=" * 70)

    try:
        from immunex.analysis import angles

        print(f"Module version: {angles.__version__}")
        print(f"Refactoring date: {angles.__refactoring_date__}")

        assert angles.__version__ == '4.0.0'
        assert angles.__refactoring_date__ == '2026-03-18'

        print("✓ Module version info correct")
        return True

    except Exception as e:
        print(f"✗ Test failed: {e}")
        return False


def main():
    """Run all tests."""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 15 + "Refactored Angle Module Test Suite" + " " * 19 + "║")
    print("╚" + "=" * 68 + "╝")

    tests = [
        ("Module Imports", test_module_imports),
        ("Public API Count", test_public_api_count),
        ("Internal Classes Hidden", test_internal_classes_hidden),
        ("Old API Removed", test_old_api_removed),
        ("DockingAngleInput Creation", test_docking_angle_input_creation),
        ("DockingAngleResult Creation", test_docking_angle_result_creation),
        ("DockingAngleAnalyzer Creation", test_docking_angle_analyzer_creation),
        ("Analyzer Error Handling", test_analyzer_with_invalid_input),
        ("Analyzer Cancellation", test_analyzer_with_cancellation),
        ("Typical Usage Pattern", test_typical_usage_pattern),
        ("Module Version", test_module_version),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n✗ {test_name} failed with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))

    # Print summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status:8s} - {test_name}")

    print("=" * 70)
    print(f"Results: {passed}/{total} tests passed ({passed/total*100:.0f}%)")
    print("=" * 70)

    if passed == total:
        print("\n🎉 All tests passed! The refactored module is ready to use.")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please review the errors above.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
