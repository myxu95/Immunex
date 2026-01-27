#!/usr/bin/env python3
"""
Quick demonstration of angle analysis module functionality.

This script demonstrates the core features without requiring actual MD files.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.angles import (
    PrincipalAxesCalculator,
    angle_between_vectors,
    project_vector_onto_plane,
    signed_angle_between_vectors
)


def demo_principal_axes():
    """Demonstrate principal axes calculation."""
    print("=" * 70)
    print("Demo 1: Principal Axes Calculation (PCA)")
    print("=" * 70)

    positions = np.array([
        [1, 0, 0], [2, 0, 0], [3, 0, 0],
        [0, 1, 0], [0, 2, 0],
        [0, 0, 1]
    ], dtype=float)
    masses = np.array([12, 12, 12, 14, 14, 16], dtype=float)

    calc = PrincipalAxesCalculator()
    axes, moments = calc.calculate_principal_axes(positions, masses)

    print("\nInput: 6 atoms along x, y, z directions")
    print(f"Positions shape: {positions.shape}")
    print(f"Masses: {masses}")

    print("\nPrincipal Axes (each row is a principal axis):")
    print(axes)

    print(f"\nPrincipal Moments (descending order):")
    print(f"  I1 = {moments[0]:.2f}")
    print(f"  I2 = {moments[1]:.2f}")
    print(f"  I3 = {moments[2]:.2f}")

    print("\nOrthogonality check (should be identity matrix):")
    orthogonality = np.dot(axes, axes.T)
    print(orthogonality)
    print()


def demo_vector_angles():
    """Demonstrate vector angle calculations."""
    print("=" * 70)
    print("Demo 2: Vector Angle Calculations")
    print("=" * 70)

    v1 = np.array([1, 0, 0])
    v2 = np.array([1, 1, 0]) / np.sqrt(2)
    v3 = np.array([0, 1, 0])

    angle_12 = angle_between_vectors(v1, v2)
    angle_13 = angle_between_vectors(v1, v3)

    print(f"\nv1 = {v1}")
    print(f"v2 = {v2}")
    print(f"v3 = {v3}")

    print(f"\nAngle(v1, v2) = {angle_12:.2f} degrees (expected: 45)")
    print(f"Angle(v1, v3) = {angle_13:.2f} degrees (expected: 90)")
    print()


def demo_projection():
    """Demonstrate vector projection onto plane."""
    print("=" * 70)
    print("Demo 3: Vector Projection onto Plane")
    print("=" * 70)

    vector = np.array([1, 1, 1])
    plane_normal = np.array([0, 0, 1])

    projected = project_vector_onto_plane(vector, plane_normal)

    print(f"\nOriginal vector: {vector}")
    print(f"Plane normal: {plane_normal}")
    print(f"Projected vector: {projected}")
    print(f"Expected: [1, 1, 0]")

    perpendicularity = np.dot(projected, plane_normal)
    print(f"\nProjected · Normal = {perpendicularity:.10f} (should be ~0)")
    print()


def demo_signed_angle():
    """Demonstrate signed angle calculation."""
    print("=" * 70)
    print("Demo 4: Signed Angle in Plane")
    print("=" * 70)

    v1 = np.array([1, 0, 0])
    v2 = np.array([0, 1, 0])
    v3 = np.array([0, -1, 0])
    normal = np.array([0, 0, 1])

    angle_pos = signed_angle_between_vectors(v1, v2, normal)
    angle_neg = signed_angle_between_vectors(v1, v3, normal)

    print(f"\nReference: v1 = {v1}")
    print(f"Normal: {normal}")

    print(f"\nRotating from v1 to v2 = {v2}")
    print(f"  Signed angle = {angle_pos:.2f} degrees (expected: +90)")

    print(f"\nRotating from v1 to v3 = {v3}")
    print(f"  Signed angle = {angle_neg:.2f} degrees (expected: -90)")
    print()


def demo_tcr_pmhc_concept():
    """Demonstrate TCR-pMHC angle concept."""
    print("=" * 70)
    print("Demo 5: TCR-pMHC Docking Angle Concept")
    print("=" * 70)

    print("\nMHC Binding Groove:")
    mhc_groove = np.array([
        [0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0],
        [0, 1, 0], [1, 1, 0], [2, 1, 0], [3, 1, 0]
    ], dtype=float)

    calc = PrincipalAxesCalculator()
    masses = np.ones(len(mhc_groove))
    mhc_axes, _ = calc.calculate_principal_axes(mhc_groove, masses)

    print(f"  {len(mhc_groove)} atoms representing binding groove")
    print(f"  Principal axis (main direction): {mhc_axes[0]}")
    print(f"  Normal (perpendicular): {mhc_axes[2]}")

    print("\nTCR Disulfide Bond Vector:")
    tcr_vector = np.array([0.7, 0.7, 0.1])
    print(f"  Vector: {tcr_vector}")

    print("\nTwist Angle Calculation:")
    print("  1. Project TCR vector onto MHC plane")
    tcr_projected = project_vector_onto_plane(tcr_vector, mhc_axes[2])
    tcr_projected = tcr_projected / np.linalg.norm(tcr_projected)
    print(f"     Projected: {tcr_projected}")

    print("  2. Calculate angle with MHC axis")
    twist = angle_between_vectors(tcr_projected, mhc_axes[0])
    print(f"     Twist angle: {twist:.2f} degrees")

    print("\nThis represents the TCR binding orientation!")
    print()


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 70)
    print("AfterMD Angle Analysis Module - Quick Demonstration")
    print("=" * 70 + "\n")

    demo_principal_axes()
    demo_vector_angles()
    demo_projection()
    demo_signed_angle()
    demo_tcr_pmhc_concept()

    print("=" * 70)
    print("Demonstration complete!")
    print("=" * 70)
    print("\nFor real MD analysis:")
    print("  - See: examples/docking_angles_usage.py")
    print("  - Read: development/ANGLES_IMPLEMENTATION_SUMMARY.md")
    print("  - Tests: development/test_angle_modules.py")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
