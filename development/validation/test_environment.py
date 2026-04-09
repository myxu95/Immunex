#!/usr/bin/env python3
"""
Test Immunex Environment

Verify all required dependencies are correctly installed.

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
from pathlib import Path


def test_imports():
    """Test all critical imports."""
    print("=" * 70)
    print("Testing Critical Dependencies")
    print("=" * 70)
    print()
    
    tests = []
    
    # Core scientific
    print("1. Core Scientific Libraries")
    print("-" * 70)
    tests.append(test_import("numpy", "NumPy"))
    tests.append(test_import("pandas", "Pandas"))
    tests.append(test_import("scipy", "SciPy"))
    tests.append(test_import("matplotlib", "Matplotlib"))
    print()
    
    # MD analysis
    print("2. MD Analysis Tools")
    print("-" * 70)
    tests.append(test_import("MDAnalysis", "MDAnalysis"))
    print()
    
    # Structure analysis
    print("3. Structure Analysis")
    print("-" * 70)
    tests.append(test_import("Bio", "BioPython"))
    tests.append(test_import("biotite", "Biotite"))
    print()
    
    # CDR detection (critical)
    print("4. CDR Detection (CRITICAL)")
    print("-" * 70)
    anarci_ok = test_import("anarci", "ANARCI")
    tests.append(anarci_ok)
    if anarci_ok:
        try:
            from anarci import anarci
            result = anarci([('test', 'EVQLVESGGGLVQPGGSLRLSCAASGFTFS')], 
                           scheme='imgt', output=False)
            if result:
                print("  ANARCI functional test: PASSED")
            else:
                print("  ANARCI functional test: WARNING (no result)")
        except Exception as e:
            print(f"  ANARCI functional test: FAILED ({e})")
    print()
    
    # Visualization
    print("5. Visualization")
    print("-" * 70)
    tests.append(test_import("seaborn", "Seaborn"))
    tests.append(test_import("plotly", "Plotly"))
    print()
    
    # Optional but recommended
    print("6. Optional Dependencies")
    print("-" * 70)
    test_import("freesasa", "FreeSASA", optional=True)
    test_import("nglview", "NGLView", optional=True)
    print()
    
    # Immunex itself
    print("7. Immunex Package")
    print("-" * 70)
    sys.path.insert(0, str(Path(__file__).parent.parent))
    immunex_ok = test_import("immunex", "Immunex")
    tests.append(immunex_ok)
    if immunex_ok:
        try:
            import immunex
            if hasattr(immunex, '__version__'):
                print(f"  Immunex version: {immunex.__version__}")
            else:
                print("  Immunex version: development")
        except:
            pass
    print()
    
    # Summary
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    
    total = len(tests)
    passed = sum(tests)
    failed = total - passed
    
    print(f"\nTotal tests: {total}")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    
    if failed == 0:
        print("\n✓ All critical dependencies installed successfully!")
        return 0
    else:
        print(f"\n✗ {failed} critical dependencies missing!")
        print("Please check the installation.")
        return 1


def test_import(module_name: str, display_name: str, optional: bool = False) -> bool:
    """
    Test if a module can be imported.
    
    Returns:
        bool: True if import successful
    """
    try:
        module = __import__(module_name)
        version = getattr(module, '__version__', 'unknown')
        status = "OK" if not optional else "OK (optional)"
        print(f"  {display_name:20s} - {status:15s} (v{version})")
        return True
    except ImportError as e:
        if optional:
            print(f"  {display_name:20s} - NOT INSTALLED (optional)")
            return True
        else:
            print(f"  {display_name:20s} - FAILED: {e}")
            return False


def test_gromacs():
    """Test GROMACS availability."""
    print("=" * 70)
    print("Testing GROMACS")
    print("=" * 70)
    print()
    
    import subprocess
    
    try:
        result = subprocess.run(['gmx', '--version'], 
                              capture_output=True, 
                              text=True, 
                              timeout=5)
        
        if result.returncode == 0:
            # Extract version
            for line in result.stdout.split('\n'):
                if 'GROMACS version' in line:
                    print(f"  {line.strip()}")
                    break
            print("  GROMACS: OK")
            return True
        else:
            print("  GROMACS: Command failed")
            return False
            
    except FileNotFoundError:
        print("  GROMACS: NOT FOUND in PATH")
        print("  Note: GROMACS is optional if using system installation")
        return False
    except Exception as e:
        print(f"  GROMACS: Error - {e}")
        return False


def main():
    """Run all tests."""
    print("\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + " " * 15 + "Immunex Environment Test" + " " * 29 + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "═" * 68 + "╝")
    print()
    
    # Test imports
    import_status = test_imports()
    
    # Test GROMACS
    print()
    test_gromacs()
    
    print()
    print("=" * 70)
    
    if import_status == 0:
        print("Environment is ready for use!")
        print()
        print("Next steps:")
        print("  1. Run index generation test:")
        print("     python development/validation/test_cdr_detection.py")
        print()
        print("  2. Generate comprehensive index:")
        print("     python development/validation/test_comprehensive_index.py")
        print()
        print("  3. Start using Immunex:")
        print("     python examples/index_generation_usage.py")
    else:
        print("Please fix the missing dependencies before proceeding.")
    
    print("=" * 70)
    print()
    
    return import_status


if __name__ == "__main__":
    sys.exit(main())
