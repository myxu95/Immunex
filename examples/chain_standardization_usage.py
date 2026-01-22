"""
Chain Standardization and Index Generation Usage Examples

This script demonstrates the usage of enhanced chain standardization features:
1. Dual PDB output (protein-only and full system)
2. Chain-based index generation (supports PDB and TPR input)
3. Integration with PBC processing workflow

Author: AfterMD Development Team
Date: 2026-01-21
"""

import os
from pathlib import Path
from aftermd.utils import (
    PDBChainStandardizer,
    ChainBasedIndexGenerator,
    generate_peptide_index_from_pdb,
    generate_peptide_index_from_tpr
)


# ============================================================================
# Example 1: Dual PDB Output (Protein-Only + Full System)
# ============================================================================

def example_dual_pdb_output():
    """
    Generate two versions of standardized PDB:
    - protein-only version (for structure analysis)
    - full system version (for MD trajectory analysis)
    """
    print("=" * 80)
    print("Example 1: Dual PDB Output")
    print("=" * 80)

    # Initialize standardizer
    standardizer = PDBChainStandardizer()

    # Input: Original PDB file (any chain arrangement)
    input_pdb = "input/7n1e_complex.pdb"
    output_prefix = "output/7n1e_std"

    # Process and generate dual outputs
    results = standardizer.process_single_dual_output(
        input_pdb=input_pdb,
        output_prefix=output_prefix,
        force=True
    )

    # Check results
    if results['protein']['status'] == 'success':
        print(f"\nProtein-only PDB generated:")
        print(f"  Output: {results['protein']['output_file']}")
        print(f"  Chains: {results['protein']['chain_mapping']}")
        print(f"    A: HLA-alpha (longest)")
        print(f"    B: beta2-microglobulin")
        print(f"    C: Peptide (shortest)")
        print(f"    D: TCR-alpha")
        print(f"    E: TCR-beta")

    if results['full']['status'] == 'success':
        print(f"\nFull system PDB generated:")
        print(f"  Output: {results['full']['output_file']}")
        print(f"  Contains: Protein + Water + Ions (protein chains standardized)")

    return results


# ============================================================================
# Example 2: Index Generation from PDB File
# ============================================================================

def example_index_from_pdb():
    """
    Generate GROMACS index file from standardized PDB file.
    Use explicit chain ID (Chain C = peptide).
    """
    print("\n" + "=" * 80)
    print("Example 2: Index Generation from PDB")
    print("=" * 80)

    # Use protein-only standardized PDB
    pdb_file = "output/7n1e_std_protein.pdb"
    output_dir = "output/indices"

    # Method 1: Using convenience function
    print("\nMethod 1: Convenience function")
    index_file = generate_peptide_index_from_pdb(
        pdb_file=pdb_file,
        output_dir=output_dir,
        peptide_chain_id='C'  # Chain C is always peptide after standardization
    )
    print(f"Generated index file: {index_file}")

    # Method 2: Using class (for more control)
    print("\nMethod 2: Using ChainBasedIndexGenerator class")
    generator = ChainBasedIndexGenerator(pdb_file)
    index_file = generator.generate_peptide_index(
        output_dir=output_dir,
        peptide_chain_id='C',
        index_name="peptide_chain.ndx"
    )
    print(f"Generated index file: {index_file}")

    return index_file


# ============================================================================
# Example 3: Index Generation from TPR File
# ============================================================================

def example_index_from_tpr():
    """
    Generate GROMACS index file from TPR file.
    TPR preserves chain IDs from standardized PDB.
    """
    print("\n" + "=" * 80)
    print("Example 3: Index Generation from TPR")
    print("=" * 80)

    # Use TPR file from GROMACS simulation
    tpr_file = "md_tasks/7n1e/md.tpr"
    output_dir = "md_tasks/7n1e/indices"

    # Method 1: Convenience function
    print("\nMethod 1: Convenience function")
    index_file = generate_peptide_index_from_tpr(
        tpr_file=tpr_file,
        output_dir=output_dir,
        peptide_chain_id='C'
    )
    print(f"Generated index file: {index_file}")
    print("Note: PDB was extracted from TPR temporarily (auto-cleanup)")

    # Method 2: Using class
    print("\nMethod 2: Using ChainBasedIndexGenerator class")
    generator = ChainBasedIndexGenerator(tpr_file)
    index_file = generator.generate_peptide_index(
        output_dir=output_dir,
        peptide_chain_id='C',
        index_name="peptide_chain.ndx"
    )
    print(f"Generated index file: {index_file}")

    return index_file


# ============================================================================
# Example 4: Complete Workflow Integration
# ============================================================================

def example_complete_workflow():
    """
    Complete workflow: PDB standardization -> Index generation -> PBC correction
    """
    print("\n" + "=" * 80)
    print("Example 4: Complete Workflow")
    print("=" * 80)

    # Step 1: Standardize PDB (dual output)
    print("\n[Step 1] PDB Standardization")
    print("-" * 80)

    standardizer = PDBChainStandardizer()
    results = standardizer.process_single_dual_output(
        input_pdb="input/7n1e_complex.pdb",
        output_prefix="output/7n1e_std",
        force=True
    )

    protein_pdb = results['protein']['output_file']
    full_pdb = results['full']['output_file']
    print(f"Standardized PDBs generated:")
    print(f"  Protein-only: {protein_pdb}")
    print(f"  Full system:  {full_pdb}")

    # Step 2: Prepare MD system (GROMACS pdb2gmx, editconf, solvate, genion)
    # User performs this step manually or via script
    print("\n[Step 2] MD System Preparation (User performs)")
    print("-" * 80)
    print("Commands:")
    print(f"  gmx pdb2gmx -f {full_pdb} -o md.gro -water tip3p")
    print("  gmx editconf -f md.gro -o md_box.gro -box 8 8 8")
    print("  gmx solvate -cp md_box.gro -o md_solv.gro")
    print("  gmx grompp -f ions.mdp -o ions.tpr")
    print("  gmx genion -s ions.tpr -o md_system.gro")
    print("  gmx grompp -f md.mdp -o md.tpr")
    print("  gmx mdrun -s md.tpr -deffnm md")

    # Step 3: Generate index for PBC correction
    print("\n[Step 3] Index Generation for PBC Correction")
    print("-" * 80)

    # Option A: From TPR (recommended if MD already run)
    tpr_file = "md_tasks/7n1e/md.tpr"
    if Path(tpr_file).exists():
        print(f"Using TPR file: {tpr_file}")
        index_file = generate_peptide_index_from_tpr(
            tpr_file=tpr_file,
            output_dir="md_tasks/7n1e",
            peptide_chain_id='C'
        )
    # Option B: From protein-only PDB (before MD)
    else:
        print(f"Using protein-only PDB: {protein_pdb}")
        index_file = generate_peptide_index_from_pdb(
            pdb_file=protein_pdb,
            output_dir="output/indices",
            peptide_chain_id='C'
        )

    print(f"Generated index file: {index_file}")

    # Step 4: PBC correction using generated index
    print("\n[Step 4] PBC Correction (User performs)")
    print("-" * 80)
    print("Commands:")
    print(f"  gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc \\")
    print(f"    -n {index_file} -center -pbc atom")
    print("  # Select: Peptide_Chain_C for centering, System for output")
    print(f"  gmx trjconv -s md.tpr -f md_center.xtc -o md_whole.xtc -pbc whole")
    print(f"  gmx trjconv -s md.tpr -f md_whole.xtc -o md_pbc.xtc -fit rot+trans")

    print("\n" + "=" * 80)
    print("Workflow completed successfully!")
    print("=" * 80)


# ============================================================================
# Example 5: Batch Processing Multiple PDBs
# ============================================================================

def example_batch_processing():
    """
    Batch process multiple PDB files: standardization + index generation
    """
    print("\n" + "=" * 80)
    print("Example 5: Batch Processing")
    print("=" * 80)

    pdb_files = [
        "input/7n1e_complex.pdb",
        "input/1ao7_complex.pdb",
        "input/3pl6_complex.pdb"
    ]

    standardizer = PDBChainStandardizer()

    for pdb_file in pdb_files:
        pdb_name = Path(pdb_file).stem
        print(f"\nProcessing: {pdb_name}")
        print("-" * 80)

        # Standardize
        output_prefix = f"output/{pdb_name}_std"
        results = standardizer.process_single_dual_output(
            input_pdb=pdb_file,
            output_prefix=output_prefix,
            force=True
        )

        if results['protein']['status'] == 'success':
            protein_pdb = results['protein']['output_file']
            print(f"  Standardized: {protein_pdb}")

            # Generate index
            index_file = generate_peptide_index_from_pdb(
                pdb_file=protein_pdb,
                output_dir=f"output/{pdb_name}_indices"
            )
            print(f"  Index file:   {index_file}")
        else:
            print(f"  Failed: {results['protein']['message']}")


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    print("\n")
    print("=" * 80)
    print("AfterMD Chain Standardization Usage Examples")
    print("=" * 80)

    # Create output directories
    os.makedirs("output/indices", exist_ok=True)

    # Run examples
    try:
        # Example 1: Dual PDB output
        example_dual_pdb_output()

        # Example 2: Index from PDB
        example_index_from_pdb()

        # Example 3: Index from TPR
        # example_index_from_tpr()  # Uncomment if TPR file exists

        # Example 4: Complete workflow
        example_complete_workflow()

        # Example 5: Batch processing
        # example_batch_processing()  # Uncomment for batch processing

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    print("\n")
    print("=" * 80)
    print("Examples completed")
    print("=" * 80)
