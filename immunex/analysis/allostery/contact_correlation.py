import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Set
from collections import defaultdict, Counter
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from networkx.algorithms import community
import logging

logger = logging.getLogger(__name__)

try:
    import sys
    mosaic_path = Path(__file__).parent.parent.parent.parent / "allosteric_workspace" / "MoSAIC" / "src"
    if mosaic_path.exists():
        sys.path.insert(0, str(mosaic_path))
    import mosaic
    MOSAIC_AVAILABLE = True
except ImportError:
    MOSAIC_AVAILABLE = False
    logger.warning("MoSAIC not available. Install with: pip install mosaic-clustering")


class ContactCorrelationAnalyzer:
    """Dynamic Cross-Correlation Analysis of Residue Contacts

    Analyzes cooperative motion within a protein by computing Pearson
    correlation between contact time series. This method identifies
    residue pairs whose contacts form/break in a coordinated manner,
    revealing allosteric communication pathways.

    Attributes:
        universe (MDAnalysis.Universe): Trajectory data container
        topology (str): Path to topology file
        trajectory (str): Path to trajectory file

    Example:
        >>> analyzer = ContactCorrelationAnalyzer("md.tpr", "md.xtc")
        >>> results = analyzer.run_full_analysis(
        ...     selection="protein",
        ...     output_dir="./allostery_analysis"
        ... )
    """

    def __init__(self, topology: str, trajectory: str):
        """
        Initialize ContactCorrelationAnalyzer.

        Args:
            topology: Topology file (.tpr, .gro, .pdb)
            trajectory: Trajectory file (.xtc, .trr)

        Raises:
            ValueError: If topology or trajectory cannot be loaded
        """
        self.topology = topology
        self.trajectory = trajectory

        try:
            self.universe = mda.Universe(topology, trajectory)
            logger.info(
                f"Contact Correlation Analyzer initialized: "
                f"{len(self.universe.trajectory)} frames, "
                f"{self.universe.atoms.n_atoms} atoms"
            )
        except Exception as e:
            logger.error(f"Failed to load trajectory: {e}")
            raise ValueError(f"Cannot initialize analyzer: {e}")

    def calculate_contact_matrix(
        self,
        selection: str = "protein",
        cutoff: float = 4.5,
        seq_dist_cutoff: int = 3,
        min_frequency: float = 0.1
    ) -> Tuple[np.ndarray, List[Tuple[int, str, int, str]], np.ndarray]:
        """Calculate binary contact time series matrix (OPTIMIZED)

        Identifies persistent residue-residue contacts across trajectory
        and constructs a binary matrix indicating contact states.

        Uses optimized single-pass algorithm with vectorized distance calculations:
        - Computes all atom-atom distances once per frame using MDAnalysis.lib.distances
        - Extracts residue-residue minimum distances from atom distance matrix
        - Filters for persistent contacts after trajectory scan

        Args:
            selection: MDAnalysis selection string (default: "protein")
            cutoff: Contact distance threshold in Angstroms (default: 4.5)
            seq_dist_cutoff: Exclude pairs with |i-j| < cutoff (default: 3)
            min_frequency: Minimum fraction of frames for persistent contact (default: 0.1)

        Returns:
            Tuple containing:
                - contact_matrix: Binary array of shape (N_contacts, T_frames), dtype=uint8
                - contact_labels: List of tuples (res_i_id, res_i_name, res_j_id, res_j_name)
                - times: Array of time stamps (ps) for each frame

        Raises:
            ValueError: If selection returns no atoms or no persistent contacts found

        Note:
            Only heavy atoms (non-hydrogen) are considered for distance calculations.
            This optimized version is ~10-50x faster than the naive implementation.
        """
        logger.info(f"Calculating contact matrix with selection: {selection}")
        logger.info(f"Parameters: cutoff={cutoff} A, seq_dist_cutoff={seq_dist_cutoff}, min_frequency={min_frequency}")

        # Select heavy atoms only
        try:
            atoms = self.universe.select_atoms(f"({selection}) and not name H*")
            residues = atoms.residues
        except Exception as e:
            logger.error(f"Selection failed: {e}")
            raise ValueError(f"Invalid selection: {e}")

        if len(residues) == 0:
            raise ValueError(f"Selection '{selection}' returned no residues")

        logger.info(f"Analyzing {len(residues)} residues, {len(atoms)} heavy atoms")

        # Build residue-to-atom index mapping for fast lookup
        logger.info("Building residue-atom mapping...")
        res_atom_indices = []
        atom_to_local_idx = {global_idx: local_idx for local_idx, global_idx in enumerate(atoms.indices)}

        for res in residues:
            # Get indices of heavy atoms in this residue
            res_heavy_atoms = res.atoms.select_atoms("not name H*")
            # Map global atom indices to local indices in the 'atoms' selection
            local_indices = [atom_to_local_idx[atom.index] for atom in res_heavy_atoms]
            res_atom_indices.append(local_indices)

        # Pre-compute residue pairs to check (respecting seq_dist_cutoff)
        logger.info("Generating residue pair list...")
        residue_pairs = []
        for i in range(len(residues)):
            for j in range(i + seq_dist_cutoff, len(residues)):
                residue_pairs.append((i, j))

        logger.info(f"Will check {len(residue_pairs)} residue pairs per frame")

        # Single pass: Compute contact counts for all residue pairs
        logger.info("Computing contact time series...")
        contact_counts = defaultdict(int)
        total_frames = len(self.universe.trajectory)
        times = np.zeros(total_frames)

        # Store contact states for all pairs across all frames
        # This allows us to build the matrix directly without second pass
        all_contacts_per_frame = []

        for frame_idx, ts in enumerate(self.universe.trajectory):
            if (frame_idx + 1) % 50 == 0:
                logger.info(f"  Frame {frame_idx + 1}/{total_frames}")

            times[frame_idx] = ts.time

            # Get all atom positions for this frame
            all_positions = atoms.positions

            # Compute full distance matrix once per frame (OPTIMIZED)
            dist_matrix = distance_array(all_positions, all_positions)

            # Extract minimum distances for each residue pair
            frame_contacts = []
            for res_i_idx, res_j_idx in residue_pairs:
                res_i = residues[res_i_idx]
                res_j = residues[res_j_idx]

                # Get atom indices for this residue pair
                atoms_i = res_atom_indices[res_i_idx]
                atoms_j = res_atom_indices[res_j_idx]

                # Extract sub-matrix for this residue pair
                # dist_matrix[atoms_i][:, atoms_j] gives all atom-atom distances
                sub_dist = dist_matrix[np.ix_(atoms_i, atoms_j)]
                min_dist = np.min(sub_dist)

                # Check contact
                if min_dist < cutoff:
                    contact_pair = (res_i.resid, res_i.resname, res_j.resid, res_j.resname)
                    contact_counts[contact_pair] += 1
                    frame_contacts.append(contact_pair)

            all_contacts_per_frame.append(set(frame_contacts))

        # Filter by minimum frequency
        logger.info("Filtering persistent contacts...")
        persistent_contacts = [
            pair for pair, count in contact_counts.items()
            if count / total_frames >= min_frequency
        ]

        if len(persistent_contacts) == 0:
            raise ValueError(
                f"No persistent contacts found with cutoff={cutoff} A and min_frequency={min_frequency}. "
                f"Try increasing cutoff or decreasing min_frequency."
            )

        logger.info(
            f"Identified {len(persistent_contacts)} persistent contacts "
            f"(frequency >= {min_frequency})"
        )

        # Build binary matrix for persistent contacts
        logger.info("Building contact matrix from stored data...")
        n_contacts = len(persistent_contacts)
        contact_matrix = np.zeros((n_contacts, total_frames), dtype=np.uint8)

        # Create lookup dictionary
        contact_to_idx = {pair: idx for idx, pair in enumerate(persistent_contacts)}

        # Fill matrix from stored contact data
        for frame_idx, frame_contacts in enumerate(all_contacts_per_frame):
            for contact_pair in frame_contacts:
                if contact_pair in contact_to_idx:
                    contact_idx = contact_to_idx[contact_pair]
                    contact_matrix[contact_idx, frame_idx] = 1

        logger.info(f"Contact matrix shape: {contact_matrix.shape}")
        logger.info(f"Matrix sparsity: {np.sum(contact_matrix) / contact_matrix.size:.2%}")

        return contact_matrix, persistent_contacts, times

    def calculate_correlation_matrix(
        self,
        contact_matrix: np.ndarray
    ) -> np.ndarray:
        """Calculate Pearson correlation matrix between contact time series

        Args:
            contact_matrix: Binary contact matrix of shape (N_contacts, T_frames)

        Returns:
            Correlation matrix of shape (N_contacts, N_contacts)
            Values range from -1 (perfect anti-correlation) to +1 (perfect correlation)

        Note:
            Uses NumPy's corrcoef for efficient correlation computation.
            Diagonal elements are always 1.0 (self-correlation).
        """
        logger.info(f"Calculating correlation matrix for {contact_matrix.shape[0]} contacts...")

        # Compute Pearson correlation coefficient matrix
        correlation_matrix = np.corrcoef(contact_matrix)

        # Log statistics
        logger.info(f"Correlation matrix shape: {correlation_matrix.shape}")

        # Get off-diagonal correlations (exclude self-correlation)
        n = correlation_matrix.shape[0]
        mask = ~np.eye(n, dtype=bool)
        off_diag = correlation_matrix[mask]

        logger.info(f"Correlation statistics:")
        logger.info(f"  Mean: {np.mean(off_diag):.3f}")
        logger.info(f"  Std: {np.std(off_diag):.3f}")
        logger.info(f"  Min: {np.min(off_diag):.3f}")
        logger.info(f"  Max: {np.max(off_diag):.3f}")

        # Count high correlations
        high_pos = np.sum(off_diag > 0.7)
        high_neg = np.sum(off_diag < -0.7)
        logger.info(f"  High positive correlation (>0.7): {high_pos} pairs")
        logger.info(f"  High negative correlation (<-0.7): {high_neg} pairs")

        return correlation_matrix

    def plot_correlation_heatmap(
        self,
        correlation_matrix: np.ndarray,
        contact_labels: List[Tuple[int, str, int, str]],
        output_file: str,
        figsize: Tuple[int, int] = (14, 12)
    ) -> None:
        """Generate correlation heatmap visualization (reference figure style)

        Args:
            correlation_matrix: Correlation matrix (N_contacts × N_contacts)
            contact_labels: List of contact pair identifiers
            output_file: Path to save the heatmap image
            figsize: Figure size in inches (default: (14, 12))

        Note:
            Uses Yellow-Green-Blue colormap similar to reference figure.
            Shows correlation values from 0 (weak) to 1 (strong).
        """
        logger.info(f"Generating correlation heatmap: {output_file}")

        # Use seaborn-v0_8 style (consistent with PlotManager)
        plt.style.use('seaborn-v0_8')

        fig, ax = plt.subplots(figsize=figsize)

        # Convert to absolute correlation for visualization (like reference)
        # Reference figure shows 0-1 range (all positive)
        abs_correlation = np.abs(correlation_matrix)

        # Handle NaN values: replace with 0 (no correlation)
        # NaN occurs when a contact has zero variance (always present or absent)
        abs_correlation = np.nan_to_num(abs_correlation, nan=0.0)

        # Yellow-Green-Blue colormap (like reference figure)
        sns.heatmap(
            abs_correlation,
            cmap='YlGnBu',  # Yellow-Green-Blue: weak to strong correlation
            vmin=0, vmax=1,
            cbar_kws={'label': 'Correlation Strength', 'shrink': 0.8},
            square=True,
            linewidths=0.5,    # Add gridlines like reference
            linecolor='gray',  # Gray gridlines
            ax=ax
        )

        ax.set_title(
            'Contact Correlation Matrix',
            fontsize=16, fontweight='bold', pad=20
        )
        ax.set_xlabel('Contact Pair Index', fontsize=12, fontweight='bold')
        ax.set_ylabel('Contact Pair Index', fontsize=12, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Correlation heatmap saved: {output_file}")

    def save_results(
        self,
        contact_matrix: np.ndarray,
        correlation_matrix: np.ndarray,
        contact_labels: List[Tuple[int, str, int, str]],
        times: np.ndarray,
        output_dir: str
    ) -> Dict[str, str]:
        """Save analysis results to files

        Args:
            contact_matrix: Binary contact matrix (N_contacts × T_frames)
            correlation_matrix: Correlation matrix (N_contacts × N_contacts)
            contact_labels: List of contact pair identifiers
            times: Time stamps for each frame (ps)
            output_dir: Directory to save output files

        Returns:
            Dictionary with paths to saved files
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        files_saved = {}

        # 1. Save contact labels
        logger.info("Saving contact labels...")
        labels_df = pd.DataFrame(contact_labels, columns=['ResA_ID', 'ResA_Name', 'ResB_ID', 'ResB_Name'])
        labels_df['ContactID'] = np.arange(len(contact_labels))
        labels_df['Frequency'] = np.sum(contact_matrix, axis=1) / contact_matrix.shape[1]
        labels_df = labels_df[['ContactID', 'ResA_ID', 'ResA_Name', 'ResB_ID', 'ResB_Name', 'Frequency']]

        labels_file = output_path / "contact_labels.csv"
        labels_df.to_csv(labels_file, index=False)
        files_saved['contact_labels'] = str(labels_file)
        logger.info(f"  Saved: {labels_file}")

        # 2. Save correlation matrix (CSV)
        logger.info("Saving correlation matrix (CSV)...")
        corr_df = pd.DataFrame(
            correlation_matrix,
            columns=[f"Contact_{i}" for i in range(len(contact_labels))],
            index=[f"Contact_{i}" for i in range(len(contact_labels))]
        )
        corr_csv_file = output_path / "correlation_matrix.csv"
        corr_df.to_csv(corr_csv_file)
        files_saved['correlation_matrix_csv'] = str(corr_csv_file)
        logger.info(f"  Saved: {corr_csv_file}")

        # 3. Save correlation matrix (NPZ - compressed binary)
        logger.info("Saving correlation matrix (NPZ)...")
        corr_npz_file = output_path / "correlation_matrix.npz"
        np.savez_compressed(
            corr_npz_file,
            correlation_matrix=correlation_matrix,
            contact_labels=np.array(contact_labels, dtype=object)
        )
        files_saved['correlation_matrix_npz'] = str(corr_npz_file)
        logger.info(f"  Saved: {corr_npz_file}")

        # 4. Save contact time series (NPZ - compressed binary)
        logger.info("Saving contact time series matrix (NPZ)...")
        timeseries_file = output_path / "contact_time_series.npz"
        np.savez_compressed(
            timeseries_file,
            contact_matrix=contact_matrix,
            times=times,
            contact_labels=np.array(contact_labels, dtype=object)
        )
        files_saved['contact_time_series'] = str(timeseries_file)
        logger.info(f"  Saved: {timeseries_file}")

        return files_saved

    def run_full_analysis(
        self,
        selection: str = "protein",
        output_dir: str = "./allostery_analysis",
        cutoff: float = 4.5,
        seq_dist_cutoff: int = 3,
        min_frequency: float = 0.1,
        plot_heatmap: bool = True,
        heatmap_figsize: Tuple[int, int] = (14, 12)
    ) -> Dict:
        """Run complete contact correlation analysis workflow

        This convenience method performs all analysis steps:
        1. Calculate contact matrix
        2. Calculate correlation matrix
        3. Save results to files
        4. Generate heatmap visualization (optional)

        Args:
            selection: MDAnalysis selection string (default: "protein")
            output_dir: Directory to save all output files (default: "./allostery_analysis")
            cutoff: Contact distance threshold in Angstroms (default: 4.5)
            seq_dist_cutoff: Exclude pairs with |i-j| < cutoff (default: 3)
            min_frequency: Minimum fraction of frames for persistent contact (default: 0.1)
            plot_heatmap: Whether to generate heatmap visualization (default: True)
            heatmap_figsize: Figure size for heatmap (default: (14, 12))

        Returns:
            Dictionary with analysis results and file paths:
                - n_contacts: Number of persistent contacts
                - n_frames: Number of trajectory frames
                - mean_correlation: Mean off-diagonal correlation
                - n_high_pos_correlation: Number of highly positively correlated pairs (>0.7)
                - n_high_neg_correlation: Number of highly negatively correlated pairs (<-0.7)
                - files: Dictionary of saved file paths

        Example:
            >>> analyzer = ContactCorrelationAnalyzer("md.tpr", "md.xtc")
            >>> results = analyzer.run_full_analysis(
            ...     selection="protein",
            ...     output_dir="./results",
            ...     cutoff=4.5
            ... )
            >>> print(f"Found {results['n_contacts']} persistent contacts")
        """
        logger.info("=" * 60)
        logger.info("Starting full contact correlation analysis")
        logger.info("=" * 60)

        # Step 1: Calculate contact matrix
        contact_matrix, contact_labels, times = self.calculate_contact_matrix(
            selection=selection,
            cutoff=cutoff,
            seq_dist_cutoff=seq_dist_cutoff,
            min_frequency=min_frequency
        )

        # Step 2: Calculate correlation matrix
        correlation_matrix = self.calculate_correlation_matrix(contact_matrix)

        # Step 3: Save results
        saved_files = self.save_results(
            contact_matrix=contact_matrix,
            correlation_matrix=correlation_matrix,
            contact_labels=contact_labels,
            times=times,
            output_dir=output_dir
        )

        # Step 4: Generate heatmap (optional)
        if plot_heatmap:
            heatmap_file = str(Path(output_dir) / "correlation_heatmap.png")
            self.plot_correlation_heatmap(
                correlation_matrix=correlation_matrix,
                contact_labels=contact_labels,
                output_file=heatmap_file,
                figsize=heatmap_figsize
            )
            saved_files['correlation_heatmap'] = heatmap_file

        # Compute summary statistics
        n = correlation_matrix.shape[0]
        mask = ~np.eye(n, dtype=bool)
        off_diag = correlation_matrix[mask]

        results = {
            'n_contacts': len(contact_labels),
            'n_frames': len(times),
            'mean_correlation': float(np.mean(off_diag)),
            'std_correlation': float(np.std(off_diag)),
            'min_correlation': float(np.min(off_diag)),
            'max_correlation': float(np.max(off_diag)),
            'n_high_pos_correlation': int(np.sum(off_diag > 0.7)),
            'n_high_neg_correlation': int(np.sum(off_diag < -0.7)),
            'files': saved_files
        }

        logger.info("=" * 60)
        logger.info("Analysis completed successfully!")
        logger.info("=" * 60)
        logger.info(f"Summary:")
        logger.info(f"  Persistent contacts: {results['n_contacts']}")
        logger.info(f"  Trajectory frames: {results['n_frames']}")
        logger.info(f"  Mean correlation: {results['mean_correlation']:.3f}")
        logger.info(f"  High positive correlations (>0.7): {results['n_high_pos_correlation']}")
        logger.info(f"  High negative correlations (<-0.7): {results['n_high_neg_correlation']}")
        logger.info(f"  Output directory: {output_dir}")

        return results

    # ============================================================================
    # Clustering Methods
    # ============================================================================

    def cluster_hierarchical(
        self,
        correlation_matrix: np.ndarray = None,
        method: str = 'average',
        threshold: float = 0.7,
        correlation_threshold: float = 0.0,
        min_cluster_size: int = 3
    ) -> Dict:
        """Perform hierarchical clustering on correlation matrix

        Uses distance metric: d = 1 - |correlation|

        Args:
            correlation_matrix: Correlation matrix (uses self if None)
            method: Linkage method (single, complete, average, ward)
            threshold: Distance threshold for cutting tree
            correlation_threshold: Minimum correlation to consider (default 0.0)
                Correlations with |r| < threshold are set to 0 before clustering
            min_cluster_size: Minimum contacts per cluster

        Returns:
            Dict with cluster_labels, linkage_matrix, n_clusters, distance_matrix
        """
        logger.info("Performing hierarchical clustering...")
        logger.info(f"  Method: {method}")
        logger.info(f"  Distance threshold: {threshold}")
        logger.info(f"  Correlation threshold: {correlation_threshold}")
        logger.info(f"  Min cluster size: {min_cluster_size}")

        # Use provided matrix or self.correlation_matrix
        if correlation_matrix is None:
            if not hasattr(self, 'correlation_matrix'):
                raise ValueError("No correlation matrix available")
            correlation_matrix = self.correlation_matrix

        # Apply correlation threshold filtering if specified
        if correlation_threshold > 0:
            logger.info(f"Applying correlation threshold filter...")
            n_edges_before = np.sum(np.abs(correlation_matrix) > 0)
            correlation_matrix_filtered = correlation_matrix.copy()
            correlation_matrix_filtered[np.abs(correlation_matrix_filtered) < correlation_threshold] = 0.0
            n_edges_after = np.sum(np.abs(correlation_matrix_filtered) > 0)
            logger.info(f"  Edges before filtering: {n_edges_before}")
            logger.info(f"  Edges after filtering: {n_edges_after} ({100*n_edges_after/n_edges_before:.1f}%)")
            correlation_matrix = correlation_matrix_filtered

        # Convert correlation to distance: 1 - |corr|
        distance_matrix = 1 - np.abs(correlation_matrix)
        np.fill_diagonal(distance_matrix, 0)

        # Handle NaN values
        nan_mask = np.isnan(distance_matrix)
        if np.any(nan_mask):
            n_nan = np.sum(nan_mask)
            logger.warning(f"Found {n_nan} NaN values, replacing with 1.0")
            distance_matrix[nan_mask] = 1.0

        # Condensed distance matrix for scipy
        condensed_dist = squareform(distance_matrix, checks=False)

        logger.info(f"Distance matrix stats: min={np.min(condensed_dist):.3f}, "
                   f"max={np.max(condensed_dist):.3f}, mean={np.mean(condensed_dist):.3f}")

        # Hierarchical clustering
        Z = linkage(condensed_dist, method=method)

        # Cut tree
        clusters = fcluster(Z, t=threshold, criterion='distance')

        # Filter small clusters
        cluster_sizes = Counter(clusters)
        valid_clusters = {c for c, size in cluster_sizes.items()
                         if size >= min_cluster_size}

        logger.info(f"Initial clustering: {len(set(clusters))} clusters")
        logger.info(f"After size filtering: {len(valid_clusters)} valid clusters")

        # Remap to 0, 1, 2, ... and -1 for noise
        cluster_mapping = {old_c: new_c for new_c, old_c in enumerate(sorted(valid_clusters))}
        cluster_labels = np.array([cluster_mapping.get(c, -1) for c in clusters])

        n_noise = np.sum(cluster_labels == -1)
        logger.info(f"Noise contacts (not assigned): {n_noise}")

        return {
            'cluster_labels': cluster_labels,
            'linkage_matrix': Z,
            'n_clusters': len(valid_clusters),
            'distance_matrix': distance_matrix
        }

    def cluster_network(
        self,
        correlation_matrix: np.ndarray = None,
        threshold: float = 0.7,
        min_cluster_size: int = 3,
        method: str = 'louvain'
    ) -> Dict:
        """Perform network-based clustering on correlation matrix

        Builds a graph where nodes are contacts and edges represent
        high correlations, then applies community detection.

        Args:
            correlation_matrix: Correlation matrix (uses self if None)
            threshold: Correlation threshold for edge creation
            min_cluster_size: Minimum contacts per cluster
            method: Community detection method (louvain, greedy, label_prop)

        Returns:
            Dict with cluster_labels, graph, n_clusters, modularity
        """
        logger.info("Performing network clustering...")
        logger.info(f"  Method: {method}")
        logger.info(f"  Correlation threshold: {threshold}")
        logger.info(f"  Min cluster size: {min_cluster_size}")

        # Use provided matrix or self
        if correlation_matrix is None:
            if not hasattr(self, 'correlation_matrix'):
                raise ValueError("No correlation matrix available")
            correlation_matrix = self.correlation_matrix

        # Build graph
        G = nx.Graph()
        n = correlation_matrix.shape[0]

        for i in range(n):
            for j in range(i+1, n):
                corr = abs(correlation_matrix[i, j])
                # Skip NaN values
                if np.isnan(corr):
                    continue
                if corr > threshold:
                    G.add_edge(i, j, weight=corr)

        logger.info(f"Graph built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

        if G.number_of_edges() == 0:
            logger.warning(f"No edges found with threshold {threshold}")
            return {
                'cluster_labels': np.full(n, -1),
                'graph': G,
                'n_clusters': 0,
                'modularity': 0.0
            }

        # Community detection
        if method == 'louvain':
            communities = community.louvain_communities(G, seed=42)
        elif method == 'greedy':
            communities = community.greedy_modularity_communities(G)
        elif method == 'label_prop':
            communities = community.label_propagation_communities(G)
        else:
            raise ValueError(f"Unknown method: {method}")

        logger.info(f"Initial communities: {len(communities)}")

        # Convert to cluster labels
        cluster_labels = np.full(n, -1)
        cluster_id = 0
        for comm in communities:
            if len(comm) >= min_cluster_size:
                for node in comm:
                    cluster_labels[node] = cluster_id
                cluster_id += 1

        logger.info(f"After size filtering: {cluster_id} valid clusters")

        # Calculate modularity
        communities_for_mod = []
        for c in range(cluster_id):
            communities_for_mod.append(set(np.where(cluster_labels == c)[0]))
        noise_nodes = set(np.where(cluster_labels == -1)[0]) & set(G.nodes())
        if noise_nodes:
            communities_for_mod.append(noise_nodes)

        mod = community.modularity(G, communities_for_mod) if len(communities_for_mod) > 0 else 0.0

        logger.info(f"Modularity: {mod:.3f}")

        return {
            'cluster_labels': cluster_labels,
            'graph': G,
            'n_clusters': cluster_id,
            'modularity': mod
        }

    def plot_cluster_dendrogram(
        self,
        linkage_matrix: np.ndarray,
        output_file: str,
        threshold: float = 0.7,
        max_items: int = 500
    ) -> None:
        """Plot hierarchical clustering dendrogram

        Args:
            linkage_matrix: Linkage matrix from hierarchical clustering
            output_file: Path to save the dendrogram image
            threshold: Threshold line to show on dendrogram
            max_items: Maximum items to display (uses truncation for larger)
        """
        logger.info(f"Generating dendrogram: {output_file}")

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        n_samples = linkage_matrix.shape[0] + 1

        plt.style.use('seaborn-v0_8')
        fig, ax = plt.subplots(figsize=(16, 8))

        # Use truncation for large dendrograms
        if n_samples > max_items:
            logger.warning(f"Large dendrogram ({n_samples} samples), using truncated display")
            dendrogram_kwargs = {
                'truncate_mode': 'lastp',
                'p': min(50, n_samples // 20),
                'color_threshold': threshold,
                'above_threshold_color': 'gray'
            }
        else:
            dendrogram_kwargs = {
                'color_threshold': threshold,
                'above_threshold_color': 'gray',
                'no_labels': n_samples > 100
            }

        dendrogram(linkage_matrix, ax=ax, **dendrogram_kwargs)

        ax.axhline(y=threshold, color='r', linestyle='--', linewidth=1.5,
                   label=f'Threshold = {threshold}')
        ax.set_xlabel('Contact Pair Index (or Merged Cluster)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Distance (1 - |correlation|)', fontsize=12, fontweight='bold')

        if n_samples > max_items:
            ax.set_title(f'Hierarchical Clustering Dendrogram (Truncated View)',
                        fontsize=14, fontweight='bold')
        else:
            ax.set_title('Hierarchical Clustering Dendrogram',
                        fontsize=14, fontweight='bold')

        ax.legend(loc='upper right')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Dendrogram saved: {output_file}")

    def plot_cluster_network(
        self,
        graph: nx.Graph,
        cluster_labels: np.ndarray,
        output_file: str,
        layout: str = 'spring'
    ) -> None:
        """Plot network clustering graph

        Args:
            graph: NetworkX graph (nodes=contacts, edges=correlations)
            cluster_labels: Array of cluster assignments
            output_file: Path to save the network image
            layout: Layout algorithm (spring, kamada)
        """
        logger.info(f"Generating network graph: {output_file}")

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        plt.style.use('seaborn-v0_8')
        fig, ax = plt.subplots(figsize=(14, 14))

        # Layout
        logger.info(f"Computing {layout} layout...")
        if layout == 'spring':
            pos = nx.spring_layout(graph, k=0.5, iterations=50, seed=42)
        elif layout == 'kamada':
            pos = nx.kamada_kawai_layout(graph)
        else:
            pos = nx.spring_layout(graph, seed=42)

        # Node colors by cluster
        unique_clusters = sorted([c for c in set(cluster_labels) if c >= 0])
        n_clusters = len(unique_clusters)

        if n_clusters > 0:
            colors = plt.cm.Set3(np.linspace(0, 1, max(n_clusters, 3)))
            node_colors = [
                colors[cluster_labels[node]] if cluster_labels[node] >= 0 else 'lightgray'
                for node in graph.nodes()
            ]
        else:
            node_colors = 'lightgray'

        # Draw network
        nx.draw_networkx_nodes(
            graph, pos,
            node_color=node_colors,
            node_size=50,
            alpha=0.8,
            ax=ax
        )

        nx.draw_networkx_edges(
            graph, pos,
            alpha=0.2,
            width=0.5,
            ax=ax
        )

        ax.set_title(
            f'Contact Correlation Network ({n_clusters} clusters)',
            fontsize=14, fontweight='bold'
        )
        ax.axis('off')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Network graph saved: {output_file}")

    def save_cluster_results(
        self,
        cluster_labels: np.ndarray,
        contact_labels: List[Tuple],
        correlation_matrix: np.ndarray,
        output_dir: str
    ) -> Dict[str, str]:
        """Save clustering results to files

        Args:
            cluster_labels: Array of cluster assignments (-1 for noise)
            contact_labels: List of contact pair identifiers
            correlation_matrix: Correlation matrix (for intra-cluster stats)
            output_dir: Directory to save output files

        Returns:
            Dictionary with paths to saved files
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        logger.info("Saving clustering results...")

        # Cluster assignments
        labels_df = pd.DataFrame(contact_labels, columns=['ResA_ID', 'ResA_Name', 'ResB_ID', 'ResB_Name'])
        labels_df['Cluster'] = cluster_labels

        assignments_file = output_path / "cluster_assignments.csv"
        labels_df.to_csv(assignments_file, index=False)
        logger.info(f"  Saved: {assignments_file}")

        # Cluster summary
        summary = []
        for cluster_id in sorted(set(cluster_labels)):
            if cluster_id < 0:
                continue  # Skip noise

            mask = cluster_labels == cluster_id
            indices = np.where(mask)[0]
            intra_corr = []
            for i in indices:
                for j in indices:
                    if i < j:
                        corr_val = abs(correlation_matrix[i, j])
                        # Skip NaN values
                        if not np.isnan(corr_val):
                            intra_corr.append(corr_val)

            summary.append({
                'Cluster': int(cluster_id),
                'Size': int(np.sum(mask)),
                'Mean_Intra_Correlation': np.mean(intra_corr) if intra_corr else 0.0,
                'Std_Intra_Correlation': np.std(intra_corr) if intra_corr else 0.0,
                'Min_Intra_Correlation': np.min(intra_corr) if intra_corr else 0.0,
                'Max_Intra_Correlation': np.max(intra_corr) if intra_corr else 0.0
            })

        summary_df = pd.DataFrame(summary)
        summary_file = output_path / "cluster_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"  Saved: {summary_file}")

        return {
            'cluster_assignments': str(assignments_file),
            'cluster_summary': str(summary_file)
        }

    def get_cluster_residues(
        self,
        cluster_labels: np.ndarray,
        contact_labels: List[Tuple]
    ) -> Dict[int, Set[int]]:
        """Extract residues involved in each cluster (for PyMOL visualization)

        Args:
            cluster_labels: Array of cluster assignments
            contact_labels: List of contact pair identifiers

        Returns:
            Dictionary mapping cluster_id -> set of residue IDs
        """
        cluster_residues = defaultdict(set)

        for idx, contact in enumerate(contact_labels):
            cluster_id = int(cluster_labels[idx])
            if cluster_id < 0:
                continue  # Skip noise

            res_a = int(contact[0])  # ResA_ID
            res_b = int(contact[2])  # ResB_ID

            cluster_residues[cluster_id].add(res_a)
            cluster_residues[cluster_id].add(res_b)

        return dict(cluster_residues)

    # ============================================================================
    # MoSAIC Clustering Methods
    # ============================================================================

    def cluster_mosaic(
        self,
        correlation_matrix: np.ndarray = None,
        mode: str = 'CPM',
        resolution_parameter: Optional[float] = None,
        correlation_threshold: float = 0.3,
        weighted: bool = True,
        n_neighbors: Optional[int] = None,
        seed: Optional[int] = 42
    ) -> Dict:
        """Perform MoSAIC clustering on correlation matrix

        MoSAIC (Molecular Systems Automated Identification of Cooperativity)
        uses Leiden algorithm to identify collective motion modules while
        automatically detecting uncorrelated coordinates as noise.

        Args:
            correlation_matrix: Correlation matrix (uses self if None)
            mode: Clustering mode ('CPM', 'modularity', or 'linkage')
                - 'CPM': Constant Potts Model (full weighted graph)
                - 'modularity': Modularity optimization (knn-graph)
                - 'linkage': Complete-linkage hierarchical clustering
            resolution_parameter: Resolution parameter for CPM/linkage modes
                If None, uses 3rd quartile (0.75) of matrix values
            correlation_threshold: Minimum correlation to keep edge (default 0.3)
                Correlations with |r| < threshold are set to 0 (removed)
            weighted: Use weighted edges in graph construction
            n_neighbors: Number of neighbors for knn-graph (None = full graph)
            seed: Random seed for reproducibility

        Returns:
            Dict with:
                - cluster_labels: Array of cluster assignments
                - mosaic_clusters: List of arrays with contact indices per cluster
                - n_clusters: Number of valid clusters
                - matrix_reordered: Block-diagonal reordered matrix
                - permutation: Permutation indices
                - ticks: Cumulative cluster sizes (for plotting)
                - silhouette_score: Clustering quality metric

        Raises:
            ImportError: If MoSAIC is not installed
            ValueError: If correlation matrix is invalid

        Example:
            >>> analyzer = ContactCorrelationAnalyzer("md.tpr", "md.xtc")
            >>> contact_matrix, contact_labels, times = analyzer.calculate_contact_matrix()
            >>> correlation_matrix = analyzer.calculate_correlation_matrix(contact_matrix)
            >>> mosaic_result = analyzer.cluster_mosaic(
            ...     correlation_matrix,
            ...     mode='CPM',
            ...     resolution_parameter=0.3,
            ...     correlation_threshold=0.3
            ... )
            >>> print(f"Identified {mosaic_result['n_clusters']} cooperative modules")
        """
        if not MOSAIC_AVAILABLE:
            raise ImportError(
                "MoSAIC is not installed. Install with: pip install mosaic-clustering "
                "or from source: https://github.com/moldyn/MoSAIC"
            )

        logger.info("Performing MoSAIC clustering...")
        logger.info(f"  Mode: {mode}")
        logger.info(f"  Resolution parameter: {resolution_parameter}")
        logger.info(f"  Correlation threshold: {correlation_threshold}")
        logger.info(f"  Weighted: {weighted}")
        logger.info(f"  KNN neighbors: {n_neighbors}")
        logger.info(f"  Seed: {seed}")

        # Use provided matrix or self
        if correlation_matrix is None:
            if not hasattr(self, 'correlation_matrix'):
                raise ValueError("No correlation matrix available")
            correlation_matrix = self.correlation_matrix

        # Handle NaN values: replace with 0 (no correlation)
        correlation_matrix_clean = np.nan_to_num(correlation_matrix, nan=0.0)

        # Threshold filtering: set weak correlations to 0
        logger.info(f"Applying correlation threshold filter...")
        n_edges_before = np.sum(np.abs(correlation_matrix_clean) > 0)
        correlation_matrix_filtered = correlation_matrix_clean.copy()
        correlation_matrix_filtered[np.abs(correlation_matrix_filtered) < correlation_threshold] = 0.0
        n_edges_after = np.sum(np.abs(correlation_matrix_filtered) > 0)
        logger.info(f"  Edges before filtering: {n_edges_before}")
        logger.info(f"  Edges after filtering: {n_edges_after} ({100*n_edges_after/n_edges_before:.1f}%)")

        # Convert to absolute correlation (MoSAIC expects similarity 0-1)
        similarity_matrix = np.abs(correlation_matrix_filtered)

        # Ensure diagonal is 1
        np.fill_diagonal(similarity_matrix, 1.0)

        # Initialize MoSAIC clustering
        try:
            clust = mosaic.Clustering(
                mode=mode,
                weighted=weighted,
                n_neighbors=n_neighbors,
                resolution_parameter=resolution_parameter,
                seed=seed
            )
        except Exception as e:
            logger.error(f"Failed to initialize MoSAIC Clustering: {e}")
            raise

        # Fit clustering
        logger.info("Fitting MoSAIC clustering...")
        try:
            clust.fit(similarity_matrix)
        except Exception as e:
            logger.error(f"MoSAIC clustering failed: {e}")
            raise

        # Extract results
        cluster_labels = clust.labels_
        mosaic_clusters = clust.clusters_
        n_clusters = len(mosaic_clusters)
        matrix_reordered = clust.matrix_
        permutation = clust.permutation_
        ticks = clust.ticks_

        logger.info(f"MoSAIC clustering completed: {n_clusters} clusters")
        logger.info(f"  Cluster sizes: {[len(c) for c in mosaic_clusters]}")

        # Calculate silhouette score
        try:
            silhouette = clust.score(similarity_matrix)
            logger.info(f"  Silhouette score: {silhouette:.3f}")
        except Exception as e:
            logger.warning(f"Could not calculate silhouette score: {e}")
            silhouette = None

        # Log resolution parameter used
        if hasattr(clust, 'resolution_param_'):
            logger.info(f"  Resolution parameter used: {clust.resolution_param_:.3f}")

        return {
            'cluster_labels': cluster_labels,
            'mosaic_clusters': mosaic_clusters,
            'n_clusters': n_clusters,
            'matrix_reordered': matrix_reordered,
            'permutation': permutation,
            'ticks': ticks,
            'silhouette_score': silhouette,
            'resolution_parameter': clust.resolution_param_ if hasattr(clust, 'resolution_param_') else resolution_parameter
        }

    def plot_mosaic_matrix(
        self,
        matrix_reordered: np.ndarray,
        ticks: np.ndarray,
        output_file: str,
        mosaic_result: Dict = None,
        min_cluster_size: int = 5,
        figsize: Tuple[int, int] = (12, 10),
        title: str = 'MoSAIC Clustered Correlation Matrix'
    ) -> None:
        """Plot block-diagonal matrix from MoSAIC clustering

        Shows only valid clusters (size >= min_cluster_size) with cluster IDs as axis labels.
        Filters out singletons and small clusters for clearer visualization.

        Args:
            matrix_reordered: Reordered correlation matrix from MoSAIC
            ticks: Cumulative cluster sizes (for plotting cluster boundaries)
            output_file: Path to save the plot
            mosaic_result: Complete MoSAIC result dict (required for filtering)
            min_cluster_size: Minimum cluster size to include (default: 5)
            figsize: Figure size in inches
            title: Plot title
        """
        logger.info(f"Generating MoSAIC matrix plot: {output_file}")

        # Step 1: Identify valid clusters (size >= min_cluster_size)
        valid_cluster_info = []

        if mosaic_result is not None and 'mosaic_clusters' in mosaic_result:
            mosaic_clusters = mosaic_result['mosaic_clusters']

            for cluster_id, members in enumerate(mosaic_clusters):
                cluster_size = len(members)
                if cluster_size >= min_cluster_size:
                    valid_cluster_info.append((cluster_id, cluster_size, members))

            logger.info(f"Found {len(valid_cluster_info)} valid clusters (size >= {min_cluster_size})")

            if len(valid_cluster_info) == 0:
                logger.warning(f"No clusters with size >= {min_cluster_size}, cannot plot filtered matrix")
                raise ValueError(f"No valid clusters found (min_size={min_cluster_size})")

            # Extract cluster IDs and sizes
            valid_cluster_ids = [info[0] for info in valid_cluster_info]
            valid_cluster_sizes = [info[1] for info in valid_cluster_info]

            # Extract all contact indices from valid clusters
            valid_contact_indices = []
            for _, _, members in valid_cluster_info:
                valid_contact_indices.extend(members)
            valid_contact_indices = sorted(valid_contact_indices)

            logger.info(f"Total contacts in valid clusters: {len(valid_contact_indices)}")

            # Step 2: Slice matrix to extract only valid cluster submatrix
            contact_mask = np.isin(np.arange(len(matrix_reordered)), valid_contact_indices)
            matrix_filtered = matrix_reordered[contact_mask, :][:, contact_mask]
            logger.info(f"Filtered matrix shape: {matrix_filtered.shape} (from {matrix_reordered.shape})")

            # Recalculate ticks for filtered clusters
            ticks_filtered = np.cumsum([0] + valid_cluster_sizes)
            logger.info(f"Filtered ticks: {ticks_filtered}")

            # Step 3: Calculate cluster center positions for axis labels
            cluster_centers = []
            for i in range(len(valid_cluster_sizes)):
                start = ticks_filtered[i]
                end = ticks_filtered[i + 1]
                center = (start + end) / 2
                cluster_centers.append(center)

            # Use original cluster IDs as labels
            cluster_labels_display = [str(cid) for cid in valid_cluster_ids]
            logger.info(f"Cluster centers: {cluster_centers}")
            logger.info(f"Cluster labels: {cluster_labels_display}")

            # Use filtered data
            abs_matrix = np.abs(matrix_filtered)
            ticks_to_plot = ticks_filtered

        else:
            # Fallback: no filtering, use original data
            logger.warning("No MoSAIC result provided, plotting full matrix without filtering")
            abs_matrix = np.abs(matrix_reordered)
            ticks_to_plot = ticks
            cluster_centers = None
            cluster_labels_display = None

        # Step 4: Plot the (filtered) matrix
        abs_matrix = np.nan_to_num(abs_matrix, nan=0.0)

        plt.style.use('seaborn-v0_8')
        fig, ax = plt.subplots(figsize=figsize)

        # Plot heatmap
        im = ax.imshow(
            abs_matrix,
            cmap='YlGnBu',
            aspect='auto',
            interpolation='nearest',
            vmin=0,
            vmax=1
        )

        # Add cluster boundary lines
        for tick in ticks_to_plot[1:-1]:
            ax.axhline(y=tick - 0.5, color='red', linewidth=1.5, linestyle='--', alpha=0.7)
            ax.axvline(x=tick - 0.5, color='red', linewidth=1.5, linestyle='--', alpha=0.7)

        # Colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Correlation Strength', fontsize=12, fontweight='bold')

        # Labels and ticks
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

        if cluster_centers is not None and cluster_labels_display is not None:
            # Use cluster IDs as axis labels (filtered mode)
            ax.set_xlabel('Cluster', fontsize=12, fontweight='bold')
            ax.set_ylabel('Cluster', fontsize=12, fontweight='bold')

            # Adaptive font size based on number of clusters
            if len(valid_cluster_ids) > 20:
                label_fontsize = 7
            elif len(valid_cluster_ids) > 10:
                label_fontsize = 9
            else:
                label_fontsize = 10

            ax.set_xticks(cluster_centers)
            ax.set_xticklabels(cluster_labels_display, fontsize=label_fontsize, rotation=45)
            ax.set_yticks(cluster_centers)
            ax.set_yticklabels(cluster_labels_display, fontsize=label_fontsize)
        else:
            # Fallback: use original contact index labels
            ax.set_xlabel('Contact Pair Index (Reordered)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Contact Pair Index (Reordered)', fontsize=12, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"MoSAIC matrix plot saved: {output_file}")

    def save_mosaic_results(
        self,
        mosaic_result: Dict,
        contact_labels: List[Tuple],
        correlation_matrix: np.ndarray,
        output_dir: str
    ) -> Dict[str, str]:
        """Save MoSAIC clustering results to files

        Args:
            mosaic_result: Dictionary returned by cluster_mosaic()
            contact_labels: List of contact pair identifiers
            correlation_matrix: Original correlation matrix
            output_dir: Directory to save output files

        Returns:
            Dictionary with paths to saved files
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        logger.info("Saving MoSAIC clustering results...")

        files_saved = {}

        # Extract results
        cluster_labels = mosaic_result['cluster_labels']
        mosaic_clusters = mosaic_result['mosaic_clusters']

        # 1. Cluster assignments
        labels_df = pd.DataFrame(contact_labels, columns=['ResA_ID', 'ResA_Name', 'ResB_ID', 'ResB_Name'])
        labels_df['Cluster'] = cluster_labels

        assignments_file = output_path / "mosaic_cluster_assignments.csv"
        labels_df.to_csv(assignments_file, index=False)
        files_saved['cluster_assignments'] = str(assignments_file)
        logger.info(f"  Saved: {assignments_file}")

        # 2. Cluster summary
        summary = []
        for cluster_id, cluster_indices in enumerate(mosaic_clusters):
            # Convert to list if needed
            if isinstance(cluster_indices, np.ndarray):
                cluster_indices = cluster_indices.tolist()
            elif not isinstance(cluster_indices, list):
                cluster_indices = list(cluster_indices)

            # Calculate intra-cluster correlation
            if len(cluster_indices) > 1:
                intra_corr = []
                for i in cluster_indices:
                    for j in cluster_indices:
                        if i < j:
                            corr_val = abs(correlation_matrix[i, j])
                            if not np.isnan(corr_val):
                                intra_corr.append(corr_val)

                summary.append({
                    'Cluster': int(cluster_id),
                    'Size': len(cluster_indices),
                    'Mean_Intra_Correlation': np.mean(intra_corr) if intra_corr else 0.0,
                    'Std_Intra_Correlation': np.std(intra_corr) if intra_corr else 0.0,
                    'Min_Intra_Correlation': np.min(intra_corr) if intra_corr else 0.0,
                    'Max_Intra_Correlation': np.max(intra_corr) if intra_corr else 0.0
                })
            else:
                summary.append({
                    'Cluster': int(cluster_id),
                    'Size': 1,
                    'Mean_Intra_Correlation': 0.0,
                    'Std_Intra_Correlation': 0.0,
                    'Min_Intra_Correlation': 0.0,
                    'Max_Intra_Correlation': 0.0
                })

        summary_df = pd.DataFrame(summary)
        summary_file = output_path / "mosaic_cluster_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        files_saved['cluster_summary'] = str(summary_file)
        logger.info(f"  Saved: {summary_file}")

        # 3. Reordered matrix (NPZ)
        matrix_file = output_path / "mosaic_matrix_reordered.npz"
        np.savez_compressed(
            matrix_file,
            matrix=mosaic_result['matrix_reordered'],
            permutation=mosaic_result['permutation'],
            ticks=mosaic_result['ticks']
        )
        files_saved['reordered_matrix'] = str(matrix_file)
        logger.info(f"  Saved: {matrix_file}")

        # 4. Metadata
        metadata = {
            'n_clusters': int(mosaic_result['n_clusters']),
            'n_contacts': len(contact_labels),
            'silhouette_score': float(mosaic_result['silhouette_score']) if mosaic_result['silhouette_score'] is not None else None,
            'resolution_parameter': float(mosaic_result['resolution_parameter']) if mosaic_result['resolution_parameter'] is not None else None
        }

        metadata_file = output_path / "mosaic_metadata.json"
        import json
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        files_saved['metadata'] = str(metadata_file)
        logger.info(f"  Saved: {metadata_file}")

        return files_saved
