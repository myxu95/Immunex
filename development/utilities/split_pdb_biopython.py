#!/usr/bin/env python3
"""
Split all multi-model PDB files using BioPython.
Save all structures to a unified training_data directory.
Naming format: {pdb_id}_{model_num:03d}.pdb
"""

from pathlib import Path
from Bio import PDB
import logging
from concurrent.futures import ProcessPoolExecutor

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def split_multimodel_pdb(task_info):
    """
    Split a multi-model PDB using BioPython.

    Args:
        task_info: Tuple of (input_pdb, output_dir, pdb_id)

    Returns:
        dict: Result with status and count
    """
    input_pdb, output_dir, pdb_id = task_info

    result = {
        'pdb_id': pdb_id,
        'status': 'failed',
        'count': 0,
        'error': None
    }

    try:
        # Parse the structure
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(input_pdb))

        # Get all models
        models = list(structure.get_models())
        n_models = len(models)

        # Save each model
        io = PDB.PDBIO()
        for i, model in enumerate(models, 1):
            # Create a new structure with single model
            output_file = output_dir / f"{pdb_id}_{i:03d}.pdb"

            # Save the model
            io.set_structure(model)
            io.save(str(output_file))

        result['status'] = 'success'
        result['count'] = n_models
        logger.info(f"✓ {pdb_id}: Extracted {n_models} structures")

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"✗ {pdb_id}: {e}")

    return result


def main():
    """Main function to split all clusters_full.pdb files."""

    logger.info("=" * 80)
    logger.info("Split Multi-model PDB Files to Training Data")
    logger.info("=" * 80)

    # Paths
    cluster_dir = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/cluster_analysis")
    output_dir = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/training_data")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all clusters_full.pdb files
    pdb_files = sorted(cluster_dir.glob("*/clusters_full.pdb"))

    logger.info(f"\nFound {len(pdb_files)} multi-model PDB files")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Naming format: {{pdb_id}}_{{num:03d}}.pdb")
    logger.info(f"Using BioPython for parsing\n")

    # Prepare tasks
    tasks = []
    for pdb_file in pdb_files:
        task_name = pdb_file.parent.name
        tasks.append((pdb_file, output_dir, task_name))

    logger.info(f"Processing {len(tasks)} tasks in parallel (4 workers)...\n")

    # Process in parallel
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(split_multimodel_pdb, tasks))

    # Summary
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']
    total_structures = sum(r['count'] for r in successful)

    logger.info("\n" + "=" * 80)
    logger.info("Split Summary")
    logger.info("=" * 80)
    logger.info(f"Total tasks: {len(results)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"\nTotal structures extracted: {total_structures:,}")
    logger.info(f"Average per task: {total_structures/len(successful):.1f}")

    logger.info(f"\nOutput directory: {output_dir}")
    logger.info(f"Naming format: {{pdb_id}}_{{num:03d}}.pdb")

    if failed:
        logger.info(f"\nFailed tasks ({len(failed)}):")
        for r in failed[:10]:
            logger.info(f"  {r['pdb_id']}: {r['error']}")

    # Check disk space
    logger.info("\n" + "=" * 80)
    logger.info("Storage Info")
    logger.info("=" * 80)
    logger.info(f"Total PDB files created: {total_structures:,}")
    logger.info(f"Estimated size: ~{total_structures * 0.98:.1f} MB (~{total_structures * 0.98 / 1024:.1f} GB)")
    logger.info("=" * 80)

    return 0 if len(failed) == 0 else 1


if __name__ == "__main__":
    exit(main())
