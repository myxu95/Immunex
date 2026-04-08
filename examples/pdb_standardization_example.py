#!/usr/bin/env python3
"""
Example usage of PDB Chain Standardizer module

This script demonstrates various use cases of the PDBChainStandardizer class
for standardizing chain identifiers in protein complex PDB files.
"""

from pathlib import Path
from immunex.analysis import PDBChainStandardizer

def example_single_file():
    """Example: Standardize a single PDB file."""
    print("="*80)
    print("Example 1: Single File Standardization")
    print("="*80)

    standardizer = PDBChainStandardizer()

    # Process single PDB file
    result = standardizer.process_single(
        input_pdb="input/1ao7.pdb",
        output_pdb="output/1ao7_standardized.pdb",
        task_name="1ao7_run1"
    )

    print(f"\nResult:")
    print(f"  Status: {result.status}")
    print(f"  Chains found: {result.num_chains}")
    print(f"  Mapping: {result.chain_mapping}")
    print(f"  Processing time: {result.processing_time:.2f}s")


def example_batch_processing():
    """Example: Batch process multiple PDB files."""
    print("\n" + "="*80)
    print("Example 2: Batch Processing")
    print("="*80)

    standardizer = PDBChainStandardizer()

    # Prepare input/output pairs
    input_dir = Path("input")
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    pairs = []
    for pdb_file in input_dir.glob("*.pdb"):
        task_name = pdb_file.stem
        output_file = output_dir / f"{task_name}_standardized.pdb"
        pairs.append((pdb_file, output_file, task_name))

    # Batch process with 4 parallel workers
    results = standardizer.batch_process(pairs, n_processes=4)

    # Save report
    standardizer.save_report(results, "standardization_report.csv")

    # Print summary
    ok_count = sum(1 for r in results if r.status == "OK")
    print(f"\n✓ Successfully standardized: {ok_count}/{len(results)} files")


def example_custom_configuration():
    """Example: Custom standardizer configuration."""
    print("\n" + "="*80)
    print("Example 3: Custom Configuration")
    print("="*80)

    # Custom chain order for different complex type
    standardizer = PDBChainStandardizer(
        standard_order=['A', 'B', 'C'],
        expected_chain_count=3,
        max_residue_threshold=500
    )

    result = standardizer.process_single(
        input_pdb="input/trimer_complex.pdb",
        output_pdb="output/trimer_standardized.pdb",
        task_name="trimer"
    )

    print(f"Custom configuration result: {result.status}")


def example_validation():
    """Example: Validate standardized PDB files."""
    print("\n" + "="*80)
    print("Example 4: Validation")
    print("="*80)

    standardizer = PDBChainStandardizer()

    # Validate a PDB file
    is_valid, message = standardizer.validate_standardized_pdb(
        "output/1ao7_standardized.pdb"
    )

    print(f"Validation result: {is_valid}")
    print(f"Message: {message}")


def example_analyze_only():
    """Example: Analyze chains without standardization."""
    print("\n" + "="*80)
    print("Example 5: Chain Analysis Only")
    print("="*80)

    standardizer = PDBChainStandardizer()

    # Analyze chain composition
    chains = standardizer.analyze_chains("input/1ao7.pdb")

    print(f"\nFound {len(chains)} protein chains:")
    for i, chain in enumerate(chains, 1):
        print(f"  {i}. Chain {chain.chain_id}: "
              f"{chain.residue_count} residues, "
              f"{chain.atom_count} atoms")

    # Create mapping (without writing output)
    mapping, status = standardizer.create_mapping(chains)

    if status == "OK":
        print(f"\nProposed mapping:")
        for old, new in mapping.items():
            print(f"  {old} → {new}")


def example_workflow_with_checks():
    """Example: Complete workflow with error handling."""
    print("\n" + "="*80)
    print("Example 6: Complete Workflow with Error Handling")
    print("="*80)

    standardizer = PDBChainStandardizer()
    pdb_files = ["1ao7.pdb", "2vlk.pdb", "3kxf.pdb"]

    for pdb_name in pdb_files:
        input_pdb = Path("input") / pdb_name
        output_pdb = Path("output") / f"{input_pdb.stem}_std.pdb"

        try:
            # Process file
            result = standardizer.process_single(
                input_pdb=input_pdb,
                output_pdb=output_pdb,
                task_name=input_pdb.stem,
                skip_if_standard=True
            )

            # Handle different statuses
            if result.status == "OK":
                print(f"✓ {pdb_name}: Standardized successfully")

                # Validate output
                is_valid, msg = standardizer.validate_standardized_pdb(output_pdb)
                if is_valid:
                    print(f"  Validation: PASS")
                else:
                    print(f"  Validation: FAIL - {msg}")

            elif result.status == "ALREADY_STANDARD":
                print(f"✓ {pdb_name}: Already in standard order")

            elif result.status == "MULTICHAIN":
                print(f"⚠ {pdb_name}: Too many chains ({result.num_chains})")

            elif result.status == "INSUFFICIENT":
                print(f"⚠ {pdb_name}: Too few chains ({result.num_chains})")

            elif result.status == "NO_PDB":
                print(f"✗ {pdb_name}: File not found")

            else:
                print(f"✗ {pdb_name}: Error - {result.error_message}")

        except Exception as e:
            print(f"✗ {pdb_name}: Exception - {str(e)}")


def example_process_from_csv():
    """Example: Process tasks from CSV file list."""
    print("\n" + "="*80)
    print("Example 7: Process from CSV Task List")
    print("="*80)

    import csv

    standardizer = PDBChainStandardizer()

    # Read task list from CSV
    pairs = []
    with open("task_list.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            input_pdb = Path(row['input_pdb'])
            output_pdb = Path(row['output_dir']) / f"{input_pdb.stem}_std.pdb"
            pairs.append((input_pdb, output_pdb, row['task_name']))

    # Batch process
    results = standardizer.batch_process(pairs, n_processes=4)

    # Categorize results
    success = [r for r in results if r.status in ["OK", "ALREADY_STANDARD"]]
    warnings = [r for r in results if r.status in ["MULTICHAIN", "INSUFFICIENT"]]
    errors = [r for r in results if r.status in ["ERROR", "NO_PDB"]]

    print(f"\nProcessing summary:")
    print(f"  Success: {len(success)}")
    print(f"  Warnings: {len(warnings)}")
    print(f"  Errors: {len(errors)}")

    # Save detailed report
    standardizer.save_report(results, "detailed_report.csv")


if __name__ == "__main__":
    print("PDB Chain Standardizer - Usage Examples")
    print("========================================\n")

    # Run examples
    try:
        example_single_file()
    except Exception as e:
        print(f"Example 1 failed: {e}")

    try:
        example_batch_processing()
    except Exception as e:
        print(f"Example 2 failed: {e}")

    try:
        example_custom_configuration()
    except Exception as e:
        print(f"Example 3 failed: {e}")

    try:
        example_validation()
    except Exception as e:
        print(f"Example 4 failed: {e}")

    try:
        example_analyze_only()
    except Exception as e:
        print(f"Example 5 failed: {e}")

    try:
        example_workflow_with_checks()
    except Exception as e:
        print(f"Example 6 failed: {e}")

    try:
        example_process_from_csv()
    except Exception as e:
        print(f"Example 7 failed: {e}")

    print("\n" + "="*80)
    print("All examples completed!")
    print("="*80)
