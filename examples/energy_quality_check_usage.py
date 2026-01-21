#!/usr/bin/env python3
"""
Energy Quality Check Usage Example

This example demonstrates how to use the EnergyQualityChecker to assess
MD simulation quality based on energy file (md.edr) analysis.

Advantages:
- Fast screening (no PBC processing needed)
- Independent quality assessment
- Reliable thermodynamic indicators
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.quality import EnergyQualityChecker


def example_1_basic_usage():
    """Example 1: Basic energy quality check"""
    print("=" * 80)
    print("Example 1: Basic Energy Quality Check")
    print("=" * 80)

    # Initialize checker with default parameters
    checker = EnergyQualityChecker(
        target_temperature=300.0,
        temp_tolerance=5.0,
        gmx_executable="gmx"
    )

    # Path to energy file
    edr_file = "/path/to/your/md.edr"

    # Check temperature stability
    temp_result = checker.check_temperature_stability(edr_file)

    print(f"\nTemperature Quality Check:")
    print(f"  Grade: {temp_result.get('grade', 'N/A')}")
    print(f"  Assessment: {temp_result.get('assessment', 'N/A')}")

    if temp_result.get('status') == 'success':
        stats = temp_result['statistics']
        print(f"  Mean Temperature: {stats['mean']:.2f} K")
        print(f"  Std Deviation: {stats['std']:.2f} K")
        print(f"  Temperature Range: [{stats['min']:.2f}, {stats['max']:.2f}] K")


def example_2_comprehensive_check():
    """Example 2: Comprehensive energy quality check"""
    print("\n" + "=" * 80)
    print("Example 2: Comprehensive Energy Quality Check")
    print("=" * 80)

    # Initialize checker
    checker = EnergyQualityChecker(
        target_temperature=300.0,
        temp_tolerance=5.0,
        target_pressure=1.0,
        pressure_tolerance=50.0,
        energy_drift_threshold=2.0
    )

    # Path to energy file
    edr_file = "/path/to/your/md.edr"

    # Comprehensive check
    result = checker.comprehensive_energy_check(edr_file)

    print(f"\nOverall Energy Quality:")
    print(f"  Grade: {result.get('energy_grade', 'N/A')}")
    print(f"  Score: {result.get('score', 0):.1f}/100")
    print(f"  Status: {result.get('overall_status', 'N/A')}")

    # Print summary
    if 'summary' in result:
        print(f"\nSummary: {result['summary']}")

    # Print issues
    if result.get('issues'):
        print(f"\nIssues Detected:")
        for issue in result['issues']:
            print(f"  - {issue}")

    # Detailed results
    print(f"\nDetailed Results:")
    details = result.get('details', {})

    # Temperature
    if 'temperature' in details:
        temp = details['temperature']
        print(f"  Temperature: Grade {temp.get('grade', 'N/A')} - {temp.get('assessment', 'N/A')}")

    # Pressure
    if 'pressure' in details and details['pressure'].get('status') == 'success':
        pressure = details['pressure']
        print(f"  Pressure: Grade {pressure.get('grade', 'N/A')} - {pressure.get('assessment', 'N/A')}")

    # Energy conservation
    if 'energy_conservation' in details:
        energy = details['energy_conservation']
        print(f"  Energy Conservation: Grade {energy.get('grade', 'N/A')} - {energy.get('assessment', 'N/A')}")

    # Energy balance
    if 'energy_balance' in details:
        balance = details['energy_balance']
        print(f"  Energy Balance: Grade {balance.get('grade', 'N/A')} - {balance.get('assessment', 'N/A')}")


def example_3_batch_screening():
    """Example 3: Batch energy quality screening"""
    print("\n" + "=" * 80)
    print("Example 3: Batch Energy Quality Screening")
    print("=" * 80)

    # Initialize checker
    checker = EnergyQualityChecker(
        target_temperature=300.0,
        temp_tolerance=5.0,
        energy_drift_threshold=2.0
    )

    # Simulate multiple MD directories
    md_directories = [
        "/path/to/md1",
        "/path/to/md2",
        "/path/to/md3",
        "/path/to/md4",
        "/path/to/md5"
    ]

    # Results storage
    results = []
    qualified_count = 0

    print(f"\nScreening {len(md_directories)} MD simulations...")
    print("-" * 80)

    for i, md_dir in enumerate(md_directories, 1):
        edr_file = Path(md_dir) / "prod" / "md.edr"

        # Alternative: check root directory
        if not edr_file.exists():
            edr_file = Path(md_dir) / "md.edr"

        if not edr_file.exists():
            print(f"[{i:2d}] {Path(md_dir).name:<20} - SKIP (no edr file)")
            continue

        # Quick check
        result = checker.comprehensive_energy_check(str(edr_file))

        grade = result.get('energy_grade', 'D')
        score = result.get('score', 0)

        # Decision: qualify if grade is A, B, or C
        is_qualified = grade in ['A', 'B', 'C']

        if is_qualified:
            qualified_count += 1
            status_icon = "✓"
        else:
            status_icon = "✗"

        print(f"[{i:2d}] {Path(md_dir).name:<20} Grade: {grade} ({score:.1f}/100) {status_icon}")

        results.append({
            'md_dir': md_dir,
            'grade': grade,
            'score': score,
            'qualified': is_qualified,
            'result': result
        })

    # Summary
    print("-" * 80)
    print(f"\nScreening Summary:")
    print(f"  Total Checked: {len(results)}")
    print(f"  Qualified (A/B/C): {qualified_count}")
    print(f"  Failed (D): {len(results) - qualified_count}")
    print(f"  Qualification Rate: {qualified_count/len(results)*100:.1f}%")

    # Grade distribution
    grade_counts = {}
    for r in results:
        grade = r['grade']
        grade_counts[grade] = grade_counts.get(grade, 0) + 1

    print(f"\nGrade Distribution:")
    for grade in ['A', 'B', 'C', 'D']:
        count = grade_counts.get(grade, 0)
        print(f"  Grade {grade}: {count} ({count/len(results)*100:.1f}%)")


def example_4_individual_checks():
    """Example 4: Individual energy component checks"""
    print("\n" + "=" * 80)
    print("Example 4: Individual Energy Component Checks")
    print("=" * 80)

    checker = EnergyQualityChecker()
    edr_file = "/path/to/your/md.edr"

    # 1. Temperature check
    print("\n1. Temperature Stability Check:")
    temp_result = checker.check_temperature_stability(edr_file)
    if temp_result.get('status') == 'success':
        print(f"   Grade: {temp_result['grade']}")
        print(f"   Mean: {temp_result['statistics']['mean']:.2f} K")
        print(f"   Std: {temp_result['statistics']['std']:.2f} K")
        print(f"   Drift: {temp_result['drift']['drift_percentage']:.3f}%")

    # 2. Pressure check (NPT only)
    print("\n2. Pressure Stability Check:")
    pressure_result = checker.check_pressure_stability(edr_file)
    if pressure_result.get('status') == 'success':
        print(f"   Grade: {pressure_result['grade']}")
        print(f"   Mean: {pressure_result['statistics']['mean']:.1f} bar")
        print(f"   Std: {pressure_result['statistics']['std']:.1f} bar")
    elif pressure_result.get('status') == 'skipped':
        print(f"   Skipped: {pressure_result['message']}")

    # 3. Energy conservation check
    print("\n3. Total Energy Conservation Check:")
    energy_result = checker.check_total_energy_conservation(edr_file)
    if energy_result.get('status') == 'success':
        print(f"   Grade: {energy_result['grade']}")
        print(f"   Energy Drift: {energy_result['drift']['drift_percentage']:.3f}%")
        print(f"   Assessment: {energy_result['drift']['assessment']}")

    # 4. Potential/Kinetic balance check
    print("\n4. Potential/Kinetic Energy Balance Check:")
    balance_result = checker.check_potential_kinetic_balance(edr_file)
    if balance_result.get('status') == 'success':
        print(f"   Grade: {balance_result['grade']}")
        print(f"   Pot/Kin Ratio: {balance_result['ratio']['mean']:.2f}")
        print(f"   Ratio CV: {balance_result['ratio']['coefficient_of_variation']:.3f}")


def example_5_custom_thresholds():
    """Example 5: Using custom quality thresholds"""
    print("\n" + "=" * 80)
    print("Example 5: Custom Quality Thresholds")
    print("=" * 80)

    # Strict thresholds for high-quality simulations
    strict_checker = EnergyQualityChecker(
        target_temperature=310.0,        # Different target temperature
        temp_tolerance=2.0,              # Stricter tolerance
        target_pressure=1.0,
        pressure_tolerance=20.0,         # Stricter pressure tolerance
        energy_drift_threshold=0.5       # Stricter drift threshold
    )

    # Lenient thresholds for exploratory simulations
    lenient_checker = EnergyQualityChecker(
        target_temperature=300.0,
        temp_tolerance=10.0,             # More lenient
        target_pressure=1.0,
        pressure_tolerance=100.0,        # More lenient
        energy_drift_threshold=5.0       # More lenient
    )

    edr_file = "/path/to/your/md.edr"

    print("\nStrict Quality Check (Publication-grade):")
    strict_result = strict_checker.comprehensive_energy_check(edr_file)
    print(f"  Grade: {strict_result.get('energy_grade', 'N/A')}")
    print(f"  Score: {strict_result.get('score', 0):.1f}/100")

    print("\nLenient Quality Check (Exploratory):")
    lenient_result = lenient_checker.comprehensive_energy_check(edr_file)
    print(f"  Grade: {lenient_result.get('energy_grade', 'N/A')}")
    print(f"  Score: {lenient_result.get('score', 0):.1f}/100")


def main():
    """Run all examples"""
    print("\n")
    print("╔" + "═" * 78 + "╗")
    print("║" + " " * 20 + "Energy Quality Checker - Usage Examples" + " " * 19 + "║")
    print("╚" + "═" * 78 + "╝")

    # Run examples
    try:
        # example_1_basic_usage()
        # example_2_comprehensive_check()
        # example_3_batch_screening()
        # example_4_individual_checks()
        # example_5_custom_thresholds()

        print("\n" + "=" * 80)
        print("Note: Examples are commented out by default.")
        print("Uncomment the examples you want to run and update file paths.")
        print("=" * 80)

    except Exception as e:
        print(f"\nError: {e}")
        print("Make sure to update file paths to your actual MD simulation directories.")


if __name__ == "__main__":
    main()
