#!/usr/bin/env python3
"""
Angle-based clustering and stability analysis for TCR-pMHC complexes

This module provides:
1. AngleClusterer: Static angle clustering for identifying binding modes
2. TrajectoryStabilityAnalyzer: MD trajectory stability analysis
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from sklearn.cluster import KMeans, DBSCAN
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')


class AngleClusterer:
    """
    Clustering analysis for TCR-pMHC binding angles

    Identifies distinct binding modes based on Docking angle and Incident angle.
    """

    def __init__(self, random_state=42):
        """
        Initialize AngleClusterer

        Parameters:
        -----------
        random_state : int
            Random seed for reproducibility
        """
        self.random_state = random_state
        self.scaler = MinMaxScaler()

        # Cluster results
        self.cluster_labels = None
        self.cluster_centers = None
        self.n_clusters = None

        # Input data
        self.pdb_ids = None
        self.docking_angles = None
        self.incident_angles = None
        self.features_scaled = None

        # Metrics
        self.silhouette = None
        self.davies_bouldin = None
        self.calinski_harabasz = None

    def load_data(self, pdb_ids: List[str], docking_angles: np.ndarray,
                  incident_angles: np.ndarray):
        """
        Load angle data

        Parameters:
        -----------
        pdb_ids : list
            PDB identifiers
        docking_angles : array
            Docking angles in degrees
        incident_angles : array
            Incident angles in degrees
        """
        self.pdb_ids = np.array(pdb_ids)
        self.docking_angles = np.array(docking_angles)
        self.incident_angles = np.array(incident_angles)

        # Create feature matrix and normalize
        features = np.column_stack([self.docking_angles, self.incident_angles])
        self.features_scaled = self.scaler.fit_transform(features)

        print(f"Loaded {len(pdb_ids)} structures")
        print(f"Docking angle range: {self.docking_angles.min():.1f}° - {self.docking_angles.max():.1f}°")
        print(f"Incident angle range: {self.incident_angles.min():.1f}° - {self.incident_angles.max():.1f}°")

    def determine_optimal_clusters(self, k_range=(2, 11),
                                   plot_metrics=True,
                                   output_file='elbow_silhouette.png'):
        """
        Determine optimal number of clusters using Elbow method and Silhouette score

        Parameters:
        -----------
        k_range : tuple
            Range of k values to test (min, max)
        plot_metrics : bool
            Whether to plot metrics
        output_file : str
            Output file for metrics plot

        Returns:
        --------
        dict : Metrics for each k value
        """
        if self.features_scaled is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        k_values = range(k_range[0], k_range[1])
        inertias = []
        silhouette_scores = []
        davies_bouldin_scores = []
        calinski_harabasz_scores = []

        print(f"\nEvaluating k from {k_range[0]} to {k_range[1]-1}...")

        for k in k_values:
            kmeans = KMeans(n_clusters=k, random_state=self.random_state, n_init=20)
            labels = kmeans.fit_predict(self.features_scaled)

            inertias.append(kmeans.inertia_)
            silhouette_scores.append(silhouette_score(self.features_scaled, labels))
            davies_bouldin_scores.append(davies_bouldin_score(self.features_scaled, labels))
            calinski_harabasz_scores.append(calinski_harabasz_score(self.features_scaled, labels))

            print(f"  k={k}: Silhouette={silhouette_scores[-1]:.3f}, "
                  f"Davies-Bouldin={davies_bouldin_scores[-1]:.3f}, "
                  f"Calinski-Harabasz={calinski_harabasz_scores[-1]:.1f}")

        # Find optimal k based on silhouette score
        optimal_k = k_values[np.argmax(silhouette_scores)]
        print(f"\nOptimal k based on Silhouette score: {optimal_k}")

        if plot_metrics:
            self._plot_cluster_metrics(k_values, inertias, silhouette_scores,
                                      davies_bouldin_scores, optimal_k, output_file)

        return {
            'k_values': list(k_values),
            'inertias': inertias,
            'silhouette_scores': silhouette_scores,
            'davies_bouldin_scores': davies_bouldin_scores,
            'calinski_harabasz_scores': calinski_harabasz_scores,
            'optimal_k': optimal_k
        }

    def _plot_cluster_metrics(self, k_values, inertias, silhouette_scores,
                             davies_bouldin_scores, optimal_k, output_file):
        """Plot elbow curve and silhouette scores"""
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        # Elbow curve
        axes[0].plot(k_values, inertias, 'bo-', linewidth=2, markersize=8)
        axes[0].set_xlabel('Number of clusters (k)', fontsize=12)
        axes[0].set_ylabel('Inertia', fontsize=12)
        axes[0].set_title('Elbow Method', fontsize=13, fontweight='bold')
        axes[0].grid(True, alpha=0.3)

        # Silhouette scores
        axes[1].plot(k_values, silhouette_scores, 'go-', linewidth=2, markersize=8)
        axes[1].axvline(optimal_k, color='red', linestyle='--', linewidth=2,
                        label=f'Optimal k={optimal_k}')
        axes[1].set_xlabel('Number of clusters (k)', fontsize=12)
        axes[1].set_ylabel('Silhouette Score', fontsize=12)
        axes[1].set_title('Silhouette Analysis', fontsize=13, fontweight='bold')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        # Davies-Bouldin index (lower is better)
        axes[2].plot(k_values, davies_bouldin_scores, 'ro-', linewidth=2, markersize=8)
        axes[2].set_xlabel('Number of clusters (k)', fontsize=12)
        axes[2].set_ylabel('Davies-Bouldin Index', fontsize=12)
        axes[2].set_title('Davies-Bouldin Index (lower better)', fontsize=13, fontweight='bold')
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Metrics plot saved: {output_file}")
        plt.close()

    def fit_kmeans(self, n_clusters='auto', k_range=(2, 11)):
        """
        Perform K-means clustering

        Parameters:
        -----------
        n_clusters : int or 'auto'
            Number of clusters or 'auto' to determine automatically
        k_range : tuple
            Range for auto selection

        Returns:
        --------
        self
        """
        if self.features_scaled is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        if n_clusters == 'auto':
            metrics = self.determine_optimal_clusters(k_range=k_range, plot_metrics=False)
            n_clusters = metrics['optimal_k']

        self.n_clusters = n_clusters

        print(f"\nPerforming K-means clustering with k={n_clusters}...")

        kmeans = KMeans(n_clusters=n_clusters, random_state=self.random_state, n_init=20)
        self.cluster_labels = kmeans.fit_predict(self.features_scaled)

        # Get cluster centers in original scale
        centers_scaled = kmeans.cluster_centers_
        centers_original = self.scaler.inverse_transform(centers_scaled)
        self.cluster_centers = centers_original

        # Calculate metrics
        self.silhouette = silhouette_score(self.features_scaled, self.cluster_labels)
        self.davies_bouldin = davies_bouldin_score(self.features_scaled, self.cluster_labels)
        self.calinski_harabasz = calinski_harabasz_score(self.features_scaled, self.cluster_labels)

        print(f"Clustering complete!")
        print(f"  Silhouette score: {self.silhouette:.3f}")
        print(f"  Davies-Bouldin index: {self.davies_bouldin:.3f}")
        print(f"  Calinski-Harabasz score: {self.calinski_harabasz:.1f}")

        # Print cluster summary
        self._print_cluster_summary()

        return self

    def fit_hierarchical(self, n_clusters=4, method='ward', metric='euclidean'):
        """
        Perform hierarchical clustering

        Parameters:
        -----------
        n_clusters : int
            Number of clusters to form
        method : str
            Linkage method ('ward', 'complete', 'average', 'single')
        metric : str
            Distance metric

        Returns:
        --------
        self
        """
        if self.features_scaled is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        self.n_clusters = n_clusters

        print(f"\nPerforming hierarchical clustering ({method} linkage)...")

        # Compute linkage
        self.linkage_matrix = linkage(self.features_scaled, method=method, metric=metric)

        # Cut dendrogram
        self.cluster_labels = fcluster(self.linkage_matrix, n_clusters, criterion='maxclust') - 1

        # Compute cluster centers
        centers = []
        for i in range(n_clusters):
            mask = self.cluster_labels == i
            center = self.features_scaled[mask].mean(axis=0)
            centers.append(center)

        self.cluster_centers = self.scaler.inverse_transform(np.array(centers))

        # Calculate metrics
        self.silhouette = silhouette_score(self.features_scaled, self.cluster_labels)
        self.davies_bouldin = davies_bouldin_score(self.features_scaled, self.cluster_labels)
        self.calinski_harabasz = calinski_harabasz_score(self.features_scaled, self.cluster_labels)

        print(f"Clustering complete!")
        print(f"  Silhouette score: {self.silhouette:.3f}")

        self._print_cluster_summary()

        return self

    def _print_cluster_summary(self):
        """Print summary statistics for each cluster"""
        print(f"\n{'='*70}")
        print("CLUSTER SUMMARY")
        print(f"{'='*70}")
        print(f"{'Cluster':<10} {'Size':<8} {'Docking (mean±std)':<25} {'Incident (mean±std)':<25}")
        print(f"{'-'*70}")

        for i in range(self.n_clusters):
            mask = self.cluster_labels == i
            size = mask.sum()

            dock_mean = self.docking_angles[mask].mean()
            dock_std = self.docking_angles[mask].std()
            inc_mean = self.incident_angles[mask].mean()
            inc_std = self.incident_angles[mask].std()

            print(f"{i:<10} {size:<8} {dock_mean:6.1f}° ± {dock_std:5.1f}°{'':<9} "
                  f"{inc_mean:6.1f}° ± {inc_std:5.1f}°")

        print(f"{'='*70}\n")

    def plot_clusters(self, output_file='clustering_scatter.png',
                     figsize=(10, 8), show_labels=False):
        """
        Plot clustering results

        Parameters:
        -----------
        output_file : str
            Output file path
        figsize : tuple
            Figure size
        show_labels : bool
            Whether to show PDB labels
        """
        if self.cluster_labels is None:
            raise ValueError("Clustering not performed. Call fit_kmeans() or fit_hierarchical() first.")

        fig, ax = plt.subplots(figsize=figsize)

        # Plot each cluster with different color
        colors = sns.color_palette('tab10', self.n_clusters)

        for i in range(self.n_clusters):
            mask = self.cluster_labels == i
            ax.scatter(self.docking_angles[mask], self.incident_angles[mask],
                      c=[colors[i]], label=f'Cluster {i} (n={mask.sum()})',
                      s=60, alpha=0.7, edgecolors='black', linewidth=0.5)

            # Optionally show labels for some points
            if show_labels and mask.sum() <= 10:
                for idx in np.where(mask)[0]:
                    ax.annotate(self.pdb_ids[idx],
                               (self.docking_angles[idx], self.incident_angles[idx]),
                               fontsize=7, alpha=0.7)

        # Plot cluster centers
        ax.scatter(self.cluster_centers[:, 0], self.cluster_centers[:, 1],
                  marker='X', s=300, c='red', edgecolors='black', linewidth=2,
                  label='Cluster Centers', zorder=10)

        # Add center labels
        for i, center in enumerate(self.cluster_centers):
            ax.annotate(f'C{i}', (center[0], center[1]),
                       fontsize=12, fontweight='bold', ha='center', va='center',
                       bbox=dict(boxstyle='circle', facecolor='white', alpha=0.8))

        ax.set_xlabel('Docking Angle (degrees)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Incident Angle (degrees)', fontsize=14, fontweight='bold')
        ax.set_title(f'TCR-pMHC Angle Clustering ({self.n_clusters} clusters)',
                    fontsize=15, fontweight='bold', pad=15)
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')

        # Add metrics text box
        textstr = f'Silhouette: {self.silhouette:.3f}\nDavies-Bouldin: {self.davies_bouldin:.3f}'
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Clustering scatter plot saved: {output_file}")
        plt.close()

    def plot_dendrogram(self, output_file='dendrogram.png', figsize=(12, 6)):
        """
        Plot dendrogram for hierarchical clustering

        Parameters:
        -----------
        output_file : str
            Output file path
        figsize : tuple
            Figure size
        """
        if not hasattr(self, 'linkage_matrix'):
            raise ValueError("Hierarchical clustering not performed.")

        fig, ax = plt.subplots(figsize=figsize)

        dendrogram(self.linkage_matrix, ax=ax,
                   color_threshold=self.linkage_matrix[-self.n_clusters, 2],
                   above_threshold_color='gray')

        ax.set_xlabel('Sample Index', fontsize=12)
        ax.set_ylabel('Distance', fontsize=12)
        ax.set_title('Hierarchical Clustering Dendrogram', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Dendrogram saved: {output_file}")
        plt.close()

    def export_assignments(self, output_file='cluster_assignments.csv'):
        """
        Export cluster assignments to CSV

        Parameters:
        -----------
        output_file : str
            Output CSV file path

        Returns:
        --------
        pd.DataFrame : Cluster assignments
        """
        if self.cluster_labels is None:
            raise ValueError("Clustering not performed.")

        df = pd.DataFrame({
            'PDB_ID': self.pdb_ids,
            'Docking_Angle': self.docking_angles,
            'Incident_Angle': self.incident_angles,
            'Cluster': self.cluster_labels
        })

        df = df.sort_values('Cluster')
        df.to_csv(output_file, index=False)

        print(f"Cluster assignments saved: {output_file}")
        return df

    def export_cluster_centers(self, output_file='cluster_centers.csv'):
        """
        Export cluster centers to CSV

        Parameters:
        -----------
        output_file : str
            Output CSV file path

        Returns:
        --------
        pd.DataFrame : Cluster centers
        """
        if self.cluster_centers is None:
            raise ValueError("Clustering not performed.")

        df = pd.DataFrame({
            'Cluster': range(self.n_clusters),
            'Center_Docking_Angle': self.cluster_centers[:, 0],
            'Center_Incident_Angle': self.cluster_centers[:, 1]
        })

        # Add cluster sizes
        sizes = [np.sum(self.cluster_labels == i) for i in range(self.n_clusters)]
        df['Cluster_Size'] = sizes

        df.to_csv(output_file, index=False)

        print(f"Cluster centers saved: {output_file}")
        return df

    def generate_report(self, output_file='clustering_report.txt'):
        """
        Generate detailed clustering report

        Parameters:
        -----------
        output_file : str
            Output text file path
        """
        if self.cluster_labels is None:
            raise ValueError("Clustering not performed.")

        with open(output_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("TCR-pMHC ANGLE CLUSTERING REPORT\n")
            f.write("="*70 + "\n\n")

            f.write(f"Total structures: {len(self.pdb_ids)}\n")
            f.write(f"Number of clusters: {self.n_clusters}\n")
            f.write(f"Clustering algorithm: K-means\n\n")

            f.write("QUALITY METRICS:\n")
            f.write(f"  Silhouette score: {self.silhouette:.4f}\n")
            f.write(f"  Davies-Bouldin index: {self.davies_bouldin:.4f}\n")
            f.write(f"  Calinski-Harabasz score: {self.calinski_harabasz:.2f}\n\n")

            f.write("="*70 + "\n")
            f.write("CLUSTER DETAILS\n")
            f.write("="*70 + "\n\n")

            for i in range(self.n_clusters):
                mask = self.cluster_labels == i
                size = mask.sum()

                f.write(f"Cluster {i}:\n")
                f.write(f"  Size: {size} structures ({size/len(self.pdb_ids)*100:.1f}%)\n")
                f.write(f"  Center: Docking {self.cluster_centers[i, 0]:.2f}°, "
                       f"Incident {self.cluster_centers[i, 1]:.2f}°\n")

                dock_mean = self.docking_angles[mask].mean()
                dock_std = self.docking_angles[mask].std()
                inc_mean = self.incident_angles[mask].mean()
                inc_std = self.incident_angles[mask].std()

                f.write(f"  Docking angle: {dock_mean:.2f}° ± {dock_std:.2f}°\n")
                f.write(f"  Incident angle: {inc_mean:.2f}° ± {inc_std:.2f}°\n")

                # List some example PDBs
                examples = self.pdb_ids[mask][:5]
                f.write(f"  Example PDBs: {', '.join(examples)}\n")

                # Check for obtuse docking angles
                obtuse_count = np.sum(self.docking_angles[mask] > 90)
                if obtuse_count > 0:
                    f.write(f"  WARNING: {obtuse_count} structures with obtuse docking angles (>90°)\n")

                f.write("\n")

        print(f"Clustering report saved: {output_file}")


# Placeholder for TrajectoryStabilityAnalyzer (Step 2)
class TrajectoryStabilityAnalyzer:
    """
    Analyze angle stability in MD trajectories

    This class will be implemented for Step 2 when MD trajectory data is available.
    """

    def __init__(self, cluster_assignments, cluster_centers):
        """
        Initialize analyzer

        Parameters:
        -----------
        cluster_assignments : pd.DataFrame
            PDB to cluster mapping
        cluster_centers : pd.DataFrame
            Cluster center angles
        """
        self.cluster_assignments = cluster_assignments
        self.cluster_centers = cluster_centers

        print("TrajectoryStabilityAnalyzer initialized (Step 2 - requires MD data)")

    def calculate_trajectory_angles(self, topology, trajectory, output_file):
        """Calculate angles for each frame in trajectory (to be implemented)"""
        raise NotImplementedError("Step 2 implementation requires MD trajectory data")

    def compute_stability_metrics(self, angle_trajectory, cluster_center):
        """Compute stability metrics (to be implemented)"""
        raise NotImplementedError("Step 2 implementation requires MD trajectory data")
