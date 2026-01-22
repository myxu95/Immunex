#!/usr/bin/env python3
"""
Step 1: Static Angle Clustering

Perform clustering analysis on TCR-pMHC static structures based on
Docking angle and Incident angle to identify distinct binding modes.
"""

import sys
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.angles import AngleClusterer


def main():
    parser = argparse.ArgumentParser(
        description='TCR-pMHC Static Angle Clustering Analysis (Step 1)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Auto-determine optimal number of clusters
  python step1_cluster_static_angles.py \\
      --input classI_complexes.csv \\
      --output angle_clustering_results

  # Use specific number of clusters
  python step1_cluster_static_angles.py \\
      --input classI_complexes.csv \\
      --output angle_clustering_results \\
      --n-clusters 4

  # Use hierarchical clustering
  python step1_cluster_static_angles.py \\
      --input classI_complexes.csv \\
      --output angle_clustering_results \\
      --algorithm hierarchical \\
      --n-clusters 4
        '''
    )

    parser.add_argument('--input', required=True,
                       help='Input CSV file with angle data (must have: PDB ID, Docking angle, Incident angle)')
    parser.add_argument('--output', required=True,
                       help='Output directory for results')
    parser.add_argument('--n-clusters', default='auto',
                       help='Number of clusters (integer or "auto" for automatic selection, default: auto)')
    parser.add_argument('--algorithm', choices=['kmeans', 'hierarchical'], default='kmeans',
                       help='Clustering algorithm (default: kmeans)')
    parser.add_argument('--k-range', default='2,11',
                       help='Range for automatic k selection (default: 2,11)')
    parser.add_argument('--random-seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')

    args = parser.parse_args()

    # Parse k-range
    k_min, k_max = map(int, args.k_range.split(','))
    k_range = (k_min, k_max)

    # Parse n_clusters
    if args.n_clusters == 'auto':
        n_clusters = 'auto'
    else:
        n_clusters = int(args.n_clusters)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("TCR-pMHC STATIC ANGLE CLUSTERING ANALYSIS (STEP 1)")
    print("="*70)
    print(f"Input file: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"Algorithm: {args.algorithm}")
    print(f"Number of clusters: {n_clusters}")
    print("="*70 + "\n")

    # Load data
    print("Loading angle data...")
    df = pd.read_csv(args.input)

    # Validate required columns
    required_cols = ['PDB ID', 'Docking angle', 'Incident angle']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Remove rows with NaN values
    df_clean = df[required_cols].dropna()
    print(f"Loaded {len(df_clean)} structures (removed {len(df) - len(df_clean)} with missing data)")

    # Extract data
    pdb_ids = df_clean['PDB ID'].values
    docking_angles = df_clean['Docking angle'].values
    incident_angles = df_clean['Incident angle'].values

    # Initialize clusterer
    clusterer = AngleClusterer(random_state=args.random_seed)
    clusterer.load_data(pdb_ids, docking_angles, incident_angles)

    # Determine optimal k if needed
    if n_clusters == 'auto' or args.algorithm == 'kmeans':
        print(f"\nDetermining optimal number of clusters (k range: {k_range})...")
        metrics_file = output_dir / 'elbow_silhouette.png'
        metrics = clusterer.determine_optimal_clusters(
            k_range=k_range,
            plot_metrics=True,
            output_file=str(metrics_file)
        )

        if n_clusters == 'auto':
            n_clusters = metrics['optimal_k']

    # Perform clustering
    if args.algorithm == 'kmeans':
        clusterer.fit_kmeans(n_clusters=n_clusters, k_range=k_range)
    elif args.algorithm == 'hierarchical':
        clusterer.fit_hierarchical(n_clusters=n_clusters)

        # Plot dendrogram
        dendrogram_file = output_dir / 'dendrogram.png'
        clusterer.plot_dendrogram(output_file=str(dendrogram_file))

    # Generate outputs
    print("\nGenerating outputs...")

    # 1. Scatter plot
    scatter_file = output_dir / 'clustering_scatter.png'
    clusterer.plot_clusters(output_file=str(scatter_file))

    # 2. Cluster assignments
    assignments_file = output_dir / 'cluster_assignments.csv'
    assignments_df = clusterer.export_assignments(output_file=str(assignments_file))

    # 3. Cluster centers
    centers_file = output_dir / 'cluster_centers.csv'
    centers_df = clusterer.export_cluster_centers(output_file=str(centers_file))

    # 4. Text report
    report_file = output_dir / 'clustering_report.txt'
    clusterer.generate_report(output_file=str(report_file))

    # Print summary
    print("\n" + "="*70)
    print("CLUSTERING COMPLETE")
    print("="*70)
    print(f"Number of clusters identified: {n_clusters}")
    print(f"Silhouette score: {clusterer.silhouette:.4f}")
    print(f"\nOutput files:")
    print(f"  - Cluster assignments: {assignments_file}")
    print(f"  - Cluster centers: {centers_file}")
    print(f"  - Scatter plot: {scatter_file}")
    print(f"  - Text report: {report_file}")

    if args.algorithm == 'kmeans':
        print(f"  - Metrics plot: {output_dir / 'elbow_silhouette.png'}")
    if args.algorithm == 'hierarchical':
        print(f"  - Dendrogram: {output_dir / 'dendrogram.png'}")

    print("\n" + "="*70)
    print("NEXT STEP:")
    print("Use cluster_assignments.csv and cluster_centers.csv for Step 2")
    print("(MD trajectory stability analysis)")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
