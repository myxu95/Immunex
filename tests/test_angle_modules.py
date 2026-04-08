#!/usr/bin/env python3
"""
Unit tests for angle analysis modules.

Tests:
1. Principal axes calculation (orthogonality, normalization)
2. Vector angle utilities (angle range, projection)
3. Docking angle calculations (basic validation)
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.angles import (
    PrincipalAxesCalculator,
    angle_between_vectors,
    project_vector_onto_plane,
    signed_angle_between_vectors
)


def test_principal_axes_orthogonality():
    """Test that principal axes are orthogonal."""
    print("Test 1: Principal axes orthogonality")

    positions = np.random.rand(50, 3) * 10
    masses = np.ones(50)

    calc = PrincipalAxesCalculator()
    axes, moments = calc.calculate_principal_axes(positions, masses)

    dot_product_matrix = np.dot(axes, axes.T)
    identity = np.eye(3)

    is_orthogonal = np.allclose(dot_product_matrix, identity, atol=1e-6)

    print(f"  Axes shape: {axes.shape}")
    print(f"  Moments: {moments}")
    print(f"  Dot product matrix:\n{dot_product_matrix}")
    print(f"  Orthogonal: {is_orthogonal}")

    assert is_orthogonal, "Principal axes are not orthogonal"
    print("  PASSED\n")


def test_principal_axes_sorted():
    """Test that moments are sorted in descending order."""
    print("Test 2: Principal moments sorting")

    positions = np.random.rand(100, 3) * 20
    masses = np.random.rand(100) + 0.5

    calc = PrincipalAxesCalculator()
    axes, moments = calc.calculate_principal_axes(positions, masses)

    is_sorted = np.all(moments[:-1] >= moments[1:])

    print(f"  Moments: {moments}")
    print(f"  Sorted (descending): {is_sorted}")

    assert is_sorted, "Principal moments are not in descending order"
    print("  PASSED\n")


def test_angle_between_vectors():
    """Test angle calculation between vectors."""
    print("Test 3: Angle between vectors")

    v1 = np.array([1, 0, 0])
    v2 = np.array([0, 1, 0])
    angle = angle_between_vectors(v1, v2)

    print(f"  v1: {v1}, v2: {v2}")
    print(f"  Angle: {angle:.2f} degrees (expected: 90)")

    assert np.isclose(angle, 90.0, atol=0.1), f"Expected 90 degrees, got {angle}"

    v3 = np.array([1, 1, 0]) / np.sqrt(2)
    angle2 = angle_between_vectors(v1, v3)

    print(f"  v1: {v1}, v3: {v3}")
    print(f"  Angle: {angle2:.2f} degrees (expected: 45)")

    assert np.isclose(angle2, 45.0, atol=0.1), f"Expected 45 degrees, got {angle2}"
    print("  PASSED\n")


def test_project_vector_onto_plane():
    """Test vector projection onto plane."""
    print("Test 4: Vector projection onto plane")

    vector = np.array([1, 1, 1])
    normal = np.array([0, 0, 1])

    projected = project_vector_onto_plane(vector, normal)

    print(f"  Original vector: {vector}")
    print(f"  Plane normal: {normal}")
    print(f"  Projected: {projected}")
    print(f"  Expected: [1, 1, 0]")

    expected = np.array([1, 1, 0])
    assert np.allclose(projected, expected), f"Expected {expected}, got {projected}"

    dot_product = np.dot(projected, normal)
    print(f"  Dot product with normal: {dot_product:.6f} (should be ~0)")

    assert np.isclose(dot_product, 0.0, atol=1e-6), "Projected vector not perpendicular to normal"
    print("  PASSED\n")


def test_signed_angle():
    """Test signed angle calculation."""
    print("Test 5: Signed angle between vectors")

    v1 = np.array([1, 0, 0])
    v2 = np.array([0, 1, 0])
    normal = np.array([0, 0, 1])

    angle = signed_angle_between_vectors(v1, v2, normal)

    print(f"  v1: {v1}, v2: {v2}")
    print(f"  Normal: {normal}")
    print(f"  Signed angle: {angle:.2f} degrees (expected: 90)")

    assert np.isclose(angle, 90.0, atol=0.1), f"Expected 90 degrees, got {angle}"

    v3 = np.array([0, -1, 0])
    angle2 = signed_angle_between_vectors(v1, v3, normal)

    print(f"  v1: {v1}, v3: {v3}")
    print(f"  Signed angle: {angle2:.2f} degrees (expected: -90)")

    assert np.isclose(angle2, -90.0, atol=0.1), f"Expected -90 degrees, got {angle2}"
    print("  PASSED\n")


def test_inertia_tensor_symmetry():
    """Test that inertia tensor is symmetric."""
    print("Test 6: Inertia tensor symmetry")

    positions = np.random.rand(30, 3) * 15
    masses = np.random.rand(30) + 1.0

    calc = PrincipalAxesCalculator()
    I = calc.calculate_inertia_tensor(positions, masses)

    is_symmetric = np.allclose(I, I.T, atol=1e-10)

    print(f"  Inertia tensor:\n{I}")
    print(f"  Transpose:\n{I.T}")
    print(f"  Symmetric: {is_symmetric}")

    assert is_symmetric, "Inertia tensor is not symmetric"
    print("  PASSED\n")


def test_angle_range():
    """Test that angles are in valid range."""
    print("Test 7: Angle range validation")

    n_tests = 100
    valid_angles = []

    for _ in range(n_tests):
        v1 = np.random.rand(3) - 0.5
        v2 = np.random.rand(3) - 0.5

        if np.linalg.norm(v1) < 0.1 or np.linalg.norm(v2) < 0.1:
            continue

        angle = angle_between_vectors(v1, v2)
        valid_angles.append(0 <= angle <= 180)

    all_valid = all(valid_angles)

    print(f"  Tested {len(valid_angles)} random vector pairs")
    print(f"  All angles in [0, 180]: {all_valid}")

    assert all_valid, "Some angles are outside valid range [0, 180]"
    print("  PASSED\n")


def run_all_tests():
    """Run all unit tests."""
    print("=" * 60)
    print("Running Angle Analysis Module Unit Tests")
    print("=" * 60 + "\n")

    tests = [
        test_principal_axes_orthogonality,
        test_principal_axes_sorted,
        test_angle_between_vectors,
        test_project_vector_onto_plane,
        test_signed_angle,
        test_inertia_tensor_symmetry,
        test_angle_range,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  FAILED: {e}\n")
            failed += 1
        except Exception as e:
            print(f"  ERROR: {e}\n")
            failed += 1

    print("=" * 60)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
