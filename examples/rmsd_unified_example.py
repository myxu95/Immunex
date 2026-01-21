#!/usr/bin/env python3
"""
Example: Using the Unified RMSD Interface

This example demonstrates the new hybrid index approach:
- Fixed components (HLA, pHLA, peptide, TCR) use unified base index
- CDR3 components are generated on-demand

Author: AfterMD Development Team
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.trajectory import RMSDInterface


def main():
    """Demonstrate RMSD calculation with unified index approach."""

    print("="*80)
    print("RMSD Calculation Example - Hybrid Index Approach")
    print("="*80)

    # Example files (replace with your actual files)
    topology = "md.tpr"
    trajectory = "md_processed.xtc"
    output_dir = "rmsd_example_output"

    # CDR3 sequences for your system
    cdr3_sequences = {
        'alpha': 'CAVRDDSNYQLIW',
        'beta': 'CASSLGQAYEQYF'
    }

    print("\n1. Initializing RMSD Interface...")
    print("   This will generate the unified base index with all fixed components")

    rmsd = RMSDInterface(
        topology=topology,
        trajectory=trajectory,
        output_dir=output_dir,
        auto_standardize=True
    )

    print("\n2. Calculate TCR RMSD aligned to pHLA (using base index)")
    result_tcr = rmsd.calculate(align="pHLA", calc="TCR")

    print(f"\nTCR RMSD Results:")
    print(f"  Mean: {result_tcr['mean']:.4f} nm")
    print(f"  Std:  {result_tcr['std']:.4f} nm")
    print(f"  Output: {result_tcr['output_file']}")

    print("\n3. Calculate peptide RMSD aligned to pHLA (using base index)")
    result_peptide = rmsd.calculate(align="pHLA", calc="peptide")

    print(f"\nPeptide RMSD Results:")
    print(f"  Mean: {result_peptide['mean']:.4f} nm")
    print(f"  Std:  {result_peptide['std']:.4f} nm")

    print("\n4. Calculate CDR3β RMSD aligned to TCR")
    print("   This will generate a separate cdr3_beta_index.ndx")

    result_cdr3 = rmsd.calculate(
        align="TCR",
        calc="CDR3_beta",
        cdr3_sequences=cdr3_sequences
    )

    print(f"\nCDR3β RMSD Results:")
    print(f"  Mean: {result_cdr3['mean']:.4f} nm")
    print(f"  Std:  {result_cdr3['std']:.4f} nm")

    print("\n5. Batch calculation with multiple components")

    calculations = [
        {'align': 'pHLA', 'calc': 'TCR'},
        {'align': 'pHLA', 'calc': 'peptide'},
        {'align': 'pHLA', 'calc': 'TCR_alpha'},
        {'align': 'pHLA', 'calc': 'TCR_beta'},
        {
            'align': 'TCR',
            'calc': 'CDR3_beta',
            'cdr3_sequences': cdr3_sequences
        }
    ]

    df = rmsd.batch_calculate(
        calculations=calculations,
        output_csv="rmsd_batch_summary.csv"
    )

    print("\nBatch Results:")
    print(df[['align_component', 'calc_component', 'rmsd_mean_nm', 'rmsd_std_nm']])

    print("\n" + "="*80)
    print("Example Complete!")
    print("="*80)

    print("\nGenerated Files:")
    print(f"  {output_dir}/")
    print(f"    ├── base_components.ndx         (6 fixed components)")
    print(f"    ├── cdr3_beta_index.ndx         (CDR3β only)")
    print(f"    ├── combined_rmsd_index.ndx     (merged for CDR3 analysis)")
    print(f"    ├── rmsd_pHLA_to_TCR.xvg")
    print(f"    ├── rmsd_pHLA_to_peptide.xvg")
    print(f"    ├── rmsd_TCR_to_CDR3_beta.xvg")
    print(f"    └── rmsd_batch_summary.csv")

    print("\nKey Points:")
    print("  • Base index contains: HLA, pHLA, peptide, TCR, TCR_alpha, TCR_beta")
    print("  • CDR3 indices generated only when needed")
    print("  • Minimal file clutter compared to individual index files")


if __name__ == "__main__":
    main()
