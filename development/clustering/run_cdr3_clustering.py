#!/usr/bin/env python3
"""
Perform CDR3β cluster analysis using GROMACS gmx cluster.
Align to pHLA, cluster CDR3β region.
Exclude tasks with abnormal TCR overall RMSD.
"""

import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import logging
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_full_structures(topology, trajectory, middle_times, output_pdb, work_dir):
    """
    Extract full complex structures for cluster middle structures.

    Args:
        topology: Path to topology file (md.tpr)
        trajectory: Path to processed trajectory
        middle_times: List of time points (ps) for cluster middle structures
        output_pdb: Output PDB file path
        work_dir: Working directory for temporary files

    Returns:
        bool: True if successful
    """
    import tempfile
    import os

    try:
        temp_pdbs = []

        # Extract each middle structure
        for i, time_ps in enumerate(middle_times, 1):
            temp_pdb = work_dir / f"temp_cluster_{i}.pdb"
            temp_pdbs.append(temp_pdb)

            cmd = [
                'gmx', 'trjconv',
                '-s', str(topology),
                '-f', str(trajectory),
                '-o', str(temp_pdb),
                '-dump', str(time_ps),
                '-pbc', 'mol'
            ]

            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            stdout, stderr = process.communicate(input="System\n", timeout=60)

            if process.returncode != 0 or not temp_pdb.exists():
                logger.warning(f"Failed to extract structure at {time_ps} ps")
                continue

        # Combine all temp PDB files into multi-model PDB
        if temp_pdbs:
            with open(output_pdb, 'w') as outfile:
                for model_num, temp_pdb in enumerate(temp_pdbs, 1):
                    if temp_pdb.exists():
                        with open(temp_pdb, 'r') as infile:
                            outfile.write(f"MODEL     {model_num:4d}\n")
                            for line in infile:
                                if line.startswith(('ATOM', 'HETATM', 'TER')):
                                    outfile.write(line)
                            outfile.write("ENDMDL\n")

                        # Clean up temp file
                        temp_pdb.unlink()

            return True

        return False

    except Exception as e:
        logger.error(f"Error extracting full structures: {e}")
        return False


def run_cluster_analysis(task_info):
    """Run gmx cluster for a single task."""
    task_name, task_dir, output_dir = task_info

    result = {
        'task_name': task_name,
        'status': 'failed',
        'error': None
    }

    try:
        topology = task_dir / "md.tpr"
        trajectory = list(task_dir.glob("*_processed.xtc"))[0]
        index_file = task_dir / "combined_cdr3_index.ndx"

        task_output_dir = output_dir / task_name
        task_output_dir.mkdir(exist_ok=True)

        # Output files
        cluster_log = task_output_dir / "cluster.log"
        cluster_pdb = task_output_dir / "clusters.pdb"
        rmsd_dist = task_output_dir / "rmsd-dist.xvg"
        rmsd_matrix = task_output_dir / "rmsd-matrix.xpm"

        # Check files exist
        if not all([topology.exists(), trajectory.exists(), index_file.exists()]):
            result['error'] = 'Required files not found'
            return result

        # Run gmx cluster
        # Fit: pHLA (fixed reference)
        # Cluster: CDR3_TCR_beta_CA (measure CDR3β conformational changes)
        # Method: gromos (RMSD-based, most commonly used)
        # Cutoff: 0.15 nm (typical for CDR3 regions)
        cmd = [
            'gmx', 'cluster',
            '-s', str(topology),
            '-f', str(trajectory),
            '-n', str(index_file),
            '-method', 'gromos',
            '-cutoff', '0.15',
            '-cl', str(cluster_pdb),
            '-o', str(rmsd_matrix),
            '-dist', str(rmsd_dist),
            '-g', str(cluster_log)
        ]

        # Provide selections: pHLA for fitting, CDR3_TCR_beta_CA for clustering
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=str(task_output_dir)
        )

        stdout, stderr = process.communicate(
            input="pHLA\nCDR3_TCR_beta_CA\n",
            timeout=300
        )

        if process.returncode != 0:
            result['error'] = f'gmx cluster failed: {stderr}'
            return result

        # Parse cluster.log for cluster information
        if cluster_log.exists():
            with open(cluster_log, 'r') as f:
                log_content = f.read()

            # Extract number of clusters and largest cluster size
            num_clusters = None
            largest_cluster_size = None
            middle_times = []

            for line in log_content.split('\n'):
                if 'Found' in line and 'cluster' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'Found' and i+1 < len(parts):
                            try:
                                num_clusters = int(parts[i+1])
                                break
                            except:
                                pass

                # Look for cluster info: "cl. | #st  rmsd | middle rmsd | cluster members"
                if 'cl.' in line and '|' in line:
                    parts = line.split('|')
                    if len(parts) >= 3 and parts[0].strip().isdigit():
                        # First data line: extract largest cluster size
                        if largest_cluster_size is None:
                            try:
                                largest_cluster_size = int(parts[1].strip().split()[0])
                            except:
                                pass
                        # Extract middle structure time (ps)
                        try:
                            middle_time_str = parts[2].strip().split()[0]
                            middle_times.append(float(middle_time_str))
                        except:
                            pass

            # Extract full complex structures for cluster middle structures
            cluster_full_pdb = task_output_dir / "clusters_full.pdb"
            if middle_times:
                extract_success = extract_full_structures(
                    topology, trajectory, middle_times, cluster_full_pdb, task_output_dir
                )
                if not extract_success:
                    logger.warning(f"  {task_name}: Failed to extract full structures")

            result['status'] = 'success'
            result['num_clusters'] = num_clusters if num_clusters else 'unknown'
            result['largest_cluster_size'] = largest_cluster_size if largest_cluster_size else 'unknown'
            result['cluster_pdb'] = str(cluster_pdb)
            result['cluster_full_pdb'] = str(cluster_full_pdb) if middle_times else None
            result['cluster_log'] = str(cluster_log)

            logger.info(f"✓ {task_name}: {num_clusters} clusters (largest: {largest_cluster_size} frames)")
        else:
            result['error'] = 'cluster.log not generated'
            return result

        return result

    except subprocess.TimeoutExpired:
        result['error'] = 'gmx cluster timed out'
        logger.error(f"✗ {task_name}: Timeout")
        return result
    except Exception as e:
        result['error'] = str(e)
        logger.error(f"✗ {task_name}: {e}")
        return result


def main():
    """Main function to run CDR3β cluster analysis."""

    logger.info("=" * 80)
    logger.info("CDR3β Cluster Analysis (Align pHLA, Cluster CDR3β)")
    logger.info("=" * 80)

    # Paths
    trajectory_base = Path("/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step")
    analysis_output = Path("/home/xumy/work/development/AfterMD/output/cdr3_analysis")
    cluster_output = analysis_output / "cluster_analysis"
    cluster_output.mkdir(exist_ok=True)

    # Load CDR3β summary
    cdr3_summary = analysis_output / "cdr3_beta_rmsd_complete_summary.csv"
    df_cdr3 = pd.read_csv(cdr3_summary)

    # Load extreme outlier list
    outlier_file = analysis_output / "tcr_rmsd_extreme_outliers.txt"
    excluded_tasks = []
    if outlier_file.exists():
        with open(outlier_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    excluded_tasks.append(line)

    logger.info(f"\nExcluding {len(excluded_tasks)} tasks with abnormal TCR RMSD:")
    for task in excluded_tasks:
        logger.info(f"  - {task}")

    # Filter tasks
    df_filtered = df_cdr3[~df_cdr3['task_name'].isin(excluded_tasks)]
    logger.info(f"\nTotal tasks for clustering: {len(df_filtered)}")

    # Prepare task list
    tasks_to_process = []
    for _, row in df_filtered.iterrows():
        task_name = row['task_name']
        task_dir = trajectory_base / task_name
        if task_dir.exists():
            tasks_to_process.append((task_name, task_dir, cluster_output))

    logger.info(f"Processing {len(tasks_to_process)} tasks")
    logger.info(f"\nClustering parameters:")
    logger.info(f"  Fitting group:     pHLA (fixed reference)")
    logger.info(f"  Clustering group:  CDR3_TCR_beta_CA")
    logger.info(f"  Method:            gromos")
    logger.info(f"  RMSD cutoff:       0.15 nm\n")

    # Process in parallel (use fewer workers for cluster analysis)
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(run_cluster_analysis, tasks_to_process))

    # Summary
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    logger.info("\n" + "=" * 80)
    logger.info("Cluster Analysis Summary")
    logger.info("=" * 80)
    logger.info(f"Total tasks: {len(results)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"Success rate: {len(successful)/len(results)*100:.1f}%")

    if successful:
        # Save summary
        df_results = pd.DataFrame(successful)
        summary_csv = cluster_output / "cluster_summary.csv"
        df_results[['task_name', 'num_clusters', 'largest_cluster_size', 'cluster_pdb', 'cluster_log']].to_csv(
            summary_csv, index=False
        )

        logger.info(f"\n✓ Cluster summary saved to: {summary_csv}")

        # Statistics
        cluster_counts = [r['num_clusters'] for r in successful if isinstance(r['num_clusters'], int)]
        if cluster_counts:
            avg_clusters = sum(cluster_counts) / len(cluster_counts)
            logger.info(f"\nCluster Statistics:")
            logger.info(f"  Tasks analyzed: {len(cluster_counts)}")
            logger.info(f"  Average clusters per task: {avg_clusters:.1f}")
            logger.info(f"  Range: {min(cluster_counts)} - {max(cluster_counts)} clusters")

            # Largest cluster sizes
            largest_sizes = [r['largest_cluster_size'] for r in successful if isinstance(r['largest_cluster_size'], int)]
            if largest_sizes:
                avg_largest = sum(largest_sizes) / len(largest_sizes)
                logger.info(f"\nLargest Cluster Size:")
                logger.info(f"  Average: {avg_largest:.1f} frames")
                logger.info(f"  Range: {min(largest_sizes)} - {max(largest_sizes)} frames")

    if failed:
        logger.info(f"\n✗ Failed tasks: {len(failed)}")
        for r in failed[:10]:
            logger.info(f"  {r['task_name']}: {r['error']}")

    logger.info("\n" + "=" * 80)
    logger.info("Output Directory Structure:")
    logger.info("=" * 80)
    logger.info(f"  {cluster_output}/")
    logger.info(f"    ├── cluster_summary.csv          (summary table)")
    logger.info(f"    └── <task_name>/")
    logger.info(f"        ├── clusters.pdb             (CDR3β representative structures)")
    logger.info(f"        ├── clusters_full.pdb        (full complex representative structures)")
    logger.info(f"        ├── cluster.log              (detailed cluster info)")
    logger.info(f"        ├── rmsd-matrix.xpm          (RMSD matrix)")
    logger.info(f"        └── rmsd-dist.xvg            (RMSD distribution)")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
