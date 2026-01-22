#!/usr/bin/env python3
"""
Extract protein chains from standardized PDB files.
Remove water chains to create clean PDB files for index generation.
"""

import os
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import logging
import MDAnalysis as mda

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_protein_chains(pdb_info):
    """Extract protein chains (A-E) from a standardized PDB file."""
    pdb_file, output_file = pdb_info

    result = {
        'pdb_id': pdb_file.stem.replace('_standardized', ''),
        'status': 'failed',
        'error': None
    }

    try:
        if not pdb_file.exists():
            result['error'] = f'PDB file not found: {pdb_file}'
            return result

        # Load structure
        u = mda.Universe(str(pdb_file))

        # Select protein chains A, B, C, D, E
        # Using protein selection to ensure we only get protein atoms
        protein = u.select_atoms("protein and (segid A or segid B or segid C or segid D or segid E)")

        if len(protein) == 0:
            # Try chainID instead of segid
            protein = u.select_atoms("protein and (chainID A or chainID B or chainID C or chainID D or chainID E)")

        if len(protein) == 0:
            result['error'] = 'No protein atoms found in chains A-E'
            return result

        # Write protein-only PDB
        protein.write(str(output_file))

        # Verify output
        if not output_file.exists():
            result['error'] = 'Output file not created'
            return result

        # Get chain info
        chains = set(protein.chainIDs)

        result['status'] = 'success'
        result['output_file'] = str(output_file)
        result['n_atoms'] = len(protein)
        result['chains'] = sorted(chains)

        logger.info(f"✓ {result['pdb_id']}: {result['n_atoms']} atoms, chains {chains}")

        return result

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"✗ {result['pdb_id']}: {e}")
        return result


def main():
    """Main function to extract protein chains."""

    logger.info("=" * 60)
    logger.info("Extracting Protein Chains from Standardized PDB Files")
    logger.info("=" * 60)

    # Paths
    standardized_pdb_dir = Path("/home/xumy/work/development/AfterMD/input/standardized_pdbs")
    output_dir = Path("/home/xumy/work/development/AfterMD/input/protein_only_pdbs")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get failed tasks
    failed_log = Path("/home/xumy/work/development/AfterMD/output/cdr3_analysis/cdr3_beta_rmsd_failed.log")

    failed_tasks = []
    with open(failed_log) as f:
        for line in f:
            if line.startswith("Task: "):
                task_name = line.strip().split(": ")[1]
                failed_tasks.append(task_name)

    logger.info(f"Found {len(failed_tasks)} failed tasks")

    # Prepare PDB files to process
    pdbs_to_process = []
    for task_name in failed_tasks:
        pdb_id_base = task_name.replace('_1', '').replace('_run1', '')

        # Try both uppercase and lowercase
        pdb_file = None
        for pdb_id in [pdb_id_base.upper(), pdb_id_base.lower()]:
            candidate = standardized_pdb_dir / f"{pdb_id}_standardized.pdb"
            if candidate.exists():
                pdb_file = candidate
                break

        if pdb_file:
            output_file = output_dir / f"{pdb_id_base.lower()}_protein_only.pdb"
            pdbs_to_process.append((pdb_file, output_file))
        else:
            logger.warning(f"Skipping {task_name}: PDB file not found")

    logger.info(f"Processing {len(pdbs_to_process)} PDB files")

    # Process in parallel
    logger.info("\nExtracting protein chains...")

    with ProcessPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(extract_protein_chains, pdbs_to_process))

    # Summarize results
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    logger.info("\n" + "=" * 60)
    logger.info("Protein Chain Extraction Summary")
    logger.info("=" * 60)
    logger.info(f"Total PDB files: {len(pdbs_to_process)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"Success rate: {len(successful)/len(pdbs_to_process)*100:.1f}%")

    if successful:
        logger.info("\nChain statistics:")
        total_atoms = sum(r['n_atoms'] for r in successful)
        logger.info(f"  Total protein atoms: {total_atoms:,}")
        logger.info(f"  Average atoms per structure: {total_atoms/len(successful):.0f}")

    if failed:
        logger.info("\nFailed PDB files:")
        for r in failed:
            logger.info(f"  {r['pdb_id']}: {r['error']}")

    logger.info("=" * 60)


if __name__ == "__main__":
    main()
