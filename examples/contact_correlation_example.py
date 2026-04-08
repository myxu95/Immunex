"""
Contact Correlation Analysis Example

This script demonstrates how to use the ContactCorrelationAnalyzer
to study allosteric communication and cooperative motion within proteins.

The analysis computes dynamic cross-correlation of residue contacts,
revealing which residue pairs move in a synchronized (positive correlation)
or anti-correlated (negative correlation) manner.
"""

import sys
from pathlib import Path

# Add project root to path if needed
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from immunex.analysis.allostery import ContactCorrelationAnalyzer


def basic_analysis_example():
    """Basic contact correlation analysis example"""

    print("=" * 70)
    print("EXAMPLE 1: Basic Contact Correlation Analysis")
    print("=" * 70)

    # Initialize analyzer with your trajectory files
    topology_file = "path/to/md.tpr"
    trajectory_file = "path/to/md_pbc.xtc"

    # Update these paths to your actual files
    print("\nNOTE: Update topology_file and trajectory_file paths in this script")
    print(f"  topology_file = '{topology_file}'")
    print(f"  trajectory_file = '{trajectory_file}'")

    # Uncomment below when you have actual files
    """
    analyzer = ContactCorrelationAnalyzer(topology_file, trajectory_file)

    # Run complete analysis with default parameters
    results = analyzer.run_full_analysis(
        selection="protein",           # Analyze entire protein
        output_dir="./allostery_results",
        cutoff=4.5,                   # Contact distance threshold (Angstrom)
        seq_dist_cutoff=3,            # Exclude |i-j| < 3 (sequence neighbors)
        min_frequency=0.1,            # Only contacts present in >= 10% of frames
        plot_heatmap=True             # Generate red-blue heatmap
    )

    # Print summary
    print(f"\\nAnalysis Results:")
    print(f"  Persistent contacts found: {results['n_contacts']}")
    print(f"  Trajectory frames analyzed: {results['n_frames']}")
    print(f"  Mean correlation: {results['mean_correlation']:.3f}")
    print(f"  High positive correlations (>0.7): {results['n_high_pos_correlation']}")
    print(f"  High negative correlations (<-0.7): {results['n_high_neg_correlation']}")

    print(f"\\nOutput files:")
    for key, filepath in results['files'].items():
        print(f"  {key}: {filepath}")
    """


def advanced_analysis_example():
    """Advanced analysis with custom parameters"""

    print("\\n")
    print("=" * 70)
    print("EXAMPLE 2: Advanced Analysis with Custom Parameters")
    print("=" * 70)

    topology_file = "path/to/md.tpr"
    trajectory_file = "path/to/md_pbc.xtc"

    print("\\nNOTE: Update file paths in this script before running")

    # Uncomment below when you have actual files
    """
    analyzer = ContactCorrelationAnalyzer(topology_file, trajectory_file)

    # Step-by-step analysis for more control

    # Step 1: Calculate contact matrix with custom parameters
    print("\\nStep 1: Calculating contact matrix...")
    contact_matrix, contact_labels, times = analyzer.calculate_contact_matrix(
        selection="protein and resid 1-100",  # Analyze specific residue range
        cutoff=5.0,                           # Larger cutoff (5.0 A)
        seq_dist_cutoff=5,                    # More stringent sequence filter
        min_frequency=0.2                     # Only very persistent contacts (20%)
    )

    print(f"  Found {len(contact_labels)} persistent contacts")
    print(f"  Contact matrix shape: {contact_matrix.shape}")

    # Step 2: Calculate correlation matrix
    print("\\nStep 2: Calculating correlation matrix...")
    correlation_matrix = analyzer.calculate_correlation_matrix(contact_matrix)

    # Step 3: Generate heatmap with custom size
    print("\\nStep 3: Generating heatmap...")
    analyzer.plot_correlation_heatmap(
        correlation_matrix=correlation_matrix,
        contact_labels=contact_labels,
        output_file="./custom_correlation_heatmap.png",
        figsize=(16, 14)  # Larger figure
    )

    # Step 4: Save results
    print("\\nStep 4: Saving results...")
    saved_files = analyzer.save_results(
        contact_matrix=contact_matrix,
        correlation_matrix=correlation_matrix,
        contact_labels=contact_labels,
        times=times,
        output_dir="./custom_allostery_results"
    )

    print("\\nAnalysis complete!")
    """


def analyze_specific_region():
    """Analyze allosteric communication in a specific protein region"""

    print("\\n")
    print("=" * 70)
    print("EXAMPLE 3: Domain-Specific Analysis (e.g., TCR-pMHC)")
    print("=" * 70)

    topology_file = "path/to/md.tpr"
    trajectory_file = "path/to/md_pbc.xtc"

    print("\\nNOTE: Update file paths in this script before running")

    # Uncomment below when you have actual files
    """
    analyzer = ContactCorrelationAnalyzer(topology_file, trajectory_file)

    # Example 1: Analyze TCR alpha chain
    print("\\nAnalyzing TCR alpha chain...")
    results_tcr_alpha = analyzer.run_full_analysis(
        selection="segname PROD",  # TCR alpha chain
        output_dir="./allostery_tcr_alpha",
        cutoff=4.5,
        plot_heatmap=True
    )

    # Example 2: Analyze peptide
    print("\\nAnalyzing peptide...")
    results_peptide = analyzer.run_full_analysis(
        selection="segname PROC",  # Peptide
        output_dir="./allostery_peptide",
        cutoff=4.5,
        plot_heatmap=True
    )

    # Example 3: Analyze entire complex
    print("\\nAnalyzing entire complex...")
    results_complex = analyzer.run_full_analysis(
        selection="protein",  # Entire complex
        output_dir="./allostery_complex",
        cutoff=4.5,
        min_frequency=0.15,  # More stringent for large system
        plot_heatmap=True
    )

    print("\\nComparison:")
    print(f"  TCR alpha contacts: {results_tcr_alpha['n_contacts']}")
    print(f"  Peptide contacts: {results_peptide['n_contacts']}")
    print(f"  Complex contacts: {results_complex['n_contacts']}")
    """


def interpret_results():
    """Guide for interpreting correlation results"""

    print("\\n")
    print("=" * 70)
    print("INTERPRETING RESULTS")
    print("=" * 70)

    print("""
The contact correlation analysis produces several outputs:

1. **Contact Labels (contact_labels.csv)**:
   - Lists all persistent residue-residue contacts
   - Includes residue IDs, names, and contact frequency
   - Sorted by frequency (most stable contacts first)

2. **Correlation Matrix (correlation_matrix.csv / .npz)**:
   - N_contacts × N_contacts symmetric matrix
   - Values range from -1 to +1:
     * +1: Perfect positive correlation (contacts form/break together)
     *  0: No correlation (independent motion)
     * -1: Perfect negative correlation (anti-correlated)

3. **Contact Time Series (contact_time_series.npz)**:
   - Binary matrix (N_contacts × T_frames)
   - 1 = contact present, 0 = contact absent
   - Useful for detailed temporal analysis

4. **Correlation Heatmap (correlation_heatmap.png)**:
   - Visual representation of correlation matrix
   - Red blocks: Positively correlated contact groups (cooperative motion)
   - Blue blocks: Negatively correlated groups (competitive motion)
   - White: No correlation (independent regions)

**Biological Interpretation**:

- **Positive correlation clusters**:
  Indicate regions that move together as functional units
  (e.g., domain breathing, hinge motions)

- **Negative correlation patterns**:
  Suggest compensatory motions or allosteric switches
  (e.g., one site opens while another closes)

- **High correlation between distant contacts**:
  Evidence of long-range allosteric communication pathways

**Key Metrics**:

- n_contacts: Total number of persistent contacts analyzed
- mean_correlation: Overall level of coordination
- n_high_pos_correlation (>0.7): Strongly cooperative pairs
- n_high_neg_correlation (<-0.7): Strongly anti-correlated pairs

**Recommended Parameters**:

For typical protein analysis:
- cutoff: 4.5 Å (standard residue contact)
- seq_dist_cutoff: 3 (exclude trivial sequence neighbors)
- min_frequency: 0.1 (10% persistence threshold)

For large complexes or noisy trajectories:
- Increase min_frequency (e.g., 0.2) to focus on stable contacts
- Increase cutoff (e.g., 5.0 Å) to capture more interactions

For detailed local analysis:
- Decrease cutoff (e.g., 4.0 Å) for tighter contacts only
- Decrease seq_dist_cutoff (e.g., 2) if interested in local clusters
""")


def main():
    """Run all examples"""

    print("\\n" + "=" * 70)
    print("Immunex Contact Correlation Analysis Examples")
    print("=" * 70)

    # Run examples
    basic_analysis_example()
    advanced_analysis_example()
    analyze_specific_region()
    interpret_results()

    print("\\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)
    print("\\nTo run analysis on your data:")
    print("1. Update topology_file and trajectory_file paths in this script")
    print("2. Uncomment the analysis code blocks")
    print("3. Run: python examples/contact_correlation_example.py")
    print("=" * 70 + "\\n")


if __name__ == "__main__":
    main()
