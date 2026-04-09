#!/usr/bin/env python3
"""
Extract full complex structures from existing cluster analysis results.
Based on cluster middle structures, extract complete complex conformations.

Optimization:
- Only extract Top 80 largest clusters (by member count)
- Only include protein atoms (exclude water and ions)

This reduces storage from ~230 GB to ~25-35 GB.
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


def extract_full_structures_for_task(task_info):
    """Extract full complex structures for one task."""
    task_name, cluster_dir, trajectory_base = task_info

    result = {
        'task_name': task_name,
        'status': 'failed',
        'error': None
    }

    try:
        # Input files
        cluster_log = cluster_dir / "cluster.log"
        topology = trajectory_base / task_name / "md.tpr"
        trajectory = list((trajectory_base / task_name).glob("*_processed.xtc"))[0]

        # Output file
        cluster_full_pdb = cluster_dir / "clusters_full.pdb"

        # Check if already exists
        if cluster_full_pdb.exists():
            logger.info(f"  {task_name}: clusters_full.pdb already exists, skipping")
            result['status'] = 'skipped'
            return result

        # Check files exist
        if not all([cluster_log.exists(), topology.exists(), trajectory.exists()]):
            result['error'] = 'Required files not found'
            return result

        # Parse cluster.log to get middle structure time points
        # Only extract Top 80 largest clusters
        with open(cluster_log, 'r') as f:
            log_content = f.read()

        cluster_info = []  # [(cluster_id, size, middle_time), ...]
        in_cluster_section = False

        for line in log_content.split('\n'):
            # Detect cluster section header
            if 'cl.' in line and 'middle' in line and '|' in line:
                in_cluster_section = True
                continue

            # Parse cluster data lines
            if in_cluster_section and '|' in line:
                parts = line.split('|')
                if len(parts) >= 3:
                    # Check if first column is a cluster ID number
                    first_col = parts[0].strip()
                    if first_col.isdigit():
                        try:
                            cluster_id = int(first_col)
                            # Extract cluster size (# of structures)
                            size_str = parts[1].strip().split()[0]
                            cluster_size = int(size_str)
                            # Extract middle time
                            middle_time_str = parts[2].strip().split()[0]
                            middle_time = float(middle_time_str)
                            cluster_info.append((cluster_id, cluster_size, middle_time))
                        except:
                            pass
                    elif not first_col:  # Continuation line
                        continue
                    else:  # End of cluster section
                        break

        if not cluster_info:
            result['error'] = 'No cluster information found'
            return result

        # Sort by cluster size (descending) and select Top 80
        cluster_info.sort(key=lambda x: x[1], reverse=True)
        top_clusters = cluster_info[:80]

        # Extract middle times for Top 80 clusters
        middle_times = [info[2] for info in top_clusters]

        logger.debug(f"  {task_name}: Total {len(cluster_info)} clusters, extracting Top {len(top_clusters)} (sizes: {top_clusters[0][1]}-{top_clusters[-1][1]} frames)")

        # Extract each middle structure
        temp_pdbs = []
        for i, time_ps in enumerate(middle_times, 1):
            temp_pdb = cluster_dir / f"temp_cluster_{i}.pdb"
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

            # Select "Protein" group (excludes water and ions)
            stdout, stderr = process.communicate(input="Protein\n", timeout=60)

            if process.returncode != 0 or not temp_pdb.exists():
                logger.warning(f"  {task_name}: Failed to extract frame {i} at {time_ps} ps")
                continue

        # Combine all temp PDB files into multi-model PDB
        extracted_count = 0
        with open(cluster_full_pdb, 'w') as outfile:
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
                    extracted_count += 1

        if extracted_count > 0:
            result['status'] = 'success'
            result['num_structures'] = extracted_count
            logger.info(f"✓ {task_name}: Extracted {extracted_count} full complex structures")
        else:
            result['error'] = 'No structures extracted'

        return result

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"✗ {task_name}: {e}")
        return result


def main():
    """Main function to extract full structures."""

    logger.info("=" * 80)
    logger.info("Extract Full Complex Structures from Cluster Results")
    logger.info("Strategy: Top 80 largest clusters + Protein only (no water)")
    logger.info("=" * 80)

    # Paths
    trajectory_base = Path("/home/xumy/work/development/Immunex/input/pbc_1000frames_2step")
    cluster_output = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/cluster_analysis")

    # Get all task directories
    task_dirs = [d for d in cluster_output.iterdir() if d.is_dir()]

    logger.info(f"\nFound {len(task_dirs)} task directories with cluster results")
    logger.info(f"Extracting Top 80 largest clusters per task (protein atoms only)\n")

    # Prepare task list
    tasks_to_process = []
    for task_dir in task_dirs:
        task_name = task_dir.name
        if (task_dir / "cluster.log").exists():
            tasks_to_process.append((task_name, task_dir, trajectory_base))

    logger.info(f"Processing {len(tasks_to_process)} tasks with valid cluster.log files\n")

    # Process in parallel (use fewer workers due to I/O intensive operations)
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(extract_full_structures_for_task, tasks_to_process))

    # Summary
    successful = [r for r in results if r['status'] == 'success']
    skipped = [r for r in results if r['status'] == 'skipped']
    failed = [r for r in results if r['status'] == 'failed']

    logger.info("\n" + "=" * 80)
    logger.info("Extraction Summary")
    logger.info("=" * 80)
    logger.info(f"Total tasks: {len(results)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Skipped (already exists): {len(skipped)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"Success rate: {len(successful)/(len(results)-len(skipped))*100:.1f}% (excluding skipped)")

    if successful:
        total_structures = sum(r['num_structures'] for r in successful)
        avg_structures = total_structures / len(successful)
        logger.info(f"\n✓ Total structures extracted: {total_structures}")
        logger.info(f"  Average per task: {avg_structures:.1f} structures")

    if failed:
        logger.info(f"\n✗ Failed tasks: {len(failed)}")
        for r in failed[:10]:
            logger.info(f"  {r['task_name']}: {r['error']}")

    logger.info("\n" + "=" * 80)
    logger.info("Output Files")
    logger.info("=" * 80)
    logger.info(f"  Location: {cluster_output}/<task_name>/clusters_full.pdb")
    logger.info(f"  Format: Multi-model PDB with full complex structures")
    logger.info(f"  Content: Top 80 largest clusters, protein atoms only (no water/ions)")
    logger.info(f"  Space saving: ~85-90% compared to full extraction")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
