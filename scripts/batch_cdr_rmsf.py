#!/usr/bin/env python3
"""
Batch RMSF calculation for CDR loops (CDR1-3, Alpha and Beta)
Uses exact sequence matching via CDRSelector
"""

import argparse
import logging
import subprocess
from pathlib import Path
from typing import Dict, List
import json
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

from aftermd.utils import CDRSelector

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CDRRMSFAnalyzer:
    """Batch RMSF analysis for CDR loops"""

    def __init__(self, input_dirs: List[str], output_dir: str,
                 cdr_csv: str, max_workers: int = 8):
        self.input_dirs = [Path(d) for d in input_dirs]
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.cdr_selector = CDRSelector(cdr_csv)
        self.max_workers = max_workers
        self.gmx = 'gmx'

    def discover_tasks(self) -> List[Dict]:
        """Discover valid MD tasks"""
        tasks = []

        for input_dir in self.input_dirs:
            if not input_dir.exists():
                logger.warning(f"Input directory not found: {input_dir}")
                continue

            for task_dir in sorted(input_dir.iterdir()):
                if not task_dir.is_dir():
                    continue

                task_name = task_dir.name
                pdb_id = task_name[:4].upper()

                # Required files
                topology = task_dir / "md.tpr"
                trajectory = task_dir / f"{task_name}_2step_processed.xtc"
                pdb_file = task_dir / "md_converted.pdb"

                if not all([topology.exists(), trajectory.exists(), pdb_file.exists()]):
                    continue

                tasks.append({
                    'name': task_name,
                    'pdb_id': pdb_id,
                    'topology': str(topology),
                    'trajectory': str(trajectory),
                    'pdb_file': str(pdb_file)
                })

        logger.info(f"Discovered {len(tasks)} valid MD tasks")
        return tasks

    def generate_cdr_index(self, pdb_file: str, pdb_id: str, output_file: str) -> Dict:
        """Generate GROMACS index file for all CDRs"""
        import MDAnalysis as mda

        # Get CDR regions
        cdr_info = self.cdr_selector.get_cdr_regions(pdb_file, pdb_id)

        # Load structure
        u = mda.Universe(pdb_file)

        with open(output_file, 'w') as f:
            cdr_count = 0
            cdr_map = {}

            for chain_type in ['alpha', 'beta']:
                for cdr_num in [1, 2, 3]:
                    cdr_key = f'cdr{cdr_num}'

                    if cdr_key not in cdr_info[chain_type]:
                        continue

                    cdr_data = cdr_info[chain_type][cdr_key]
                    chain_id = cdr_data['chain_id']
                    resid_start, resid_end = cdr_data['residue_range']

                    # Select CA atoms for this CDR
                    selection = (f"chainID {chain_id} and "
                               f"resid {resid_start}:{resid_end} and "
                               f"name CA")

                    try:
                        atoms = u.select_atoms(selection)
                        if len(atoms) == 0:
                            logger.warning(f"No atoms found for {chain_type} CDR{cdr_num}")
                            continue

                        # Write index group
                        group_name = f"CDR{cdr_num}_{chain_type}_CA"
                        f.write(f"[ {group_name} ]\n")

                        indices = atoms.indices + 1  # 1-based indexing
                        for i in range(0, len(indices), 10):
                            f.write(" ".join(map(str, indices[i:i+10])) + "\n")
                        f.write("\n")

                        cdr_map[f"{chain_type}_cdr{cdr_num}"] = {
                            'group_name': group_name,
                            'n_atoms': len(atoms),
                            'chain': chain_id,
                            'residues': f"{resid_start}-{resid_end}"
                        }
                        cdr_count += 1
                        logger.info(f"  {group_name}: {len(atoms)} CA atoms")

                    except Exception as e:
                        logger.error(f"Failed to select {chain_type} CDR{cdr_num}: {e}")

        logger.info(f"Index file generated: {output_file} ({cdr_count} CDRs)")
        return cdr_map

    def calculate_rmsf(self, topology: str, trajectory: str, index_file: str,
                      group_name: str, output_file: str) -> Dict:
        """Calculate RMSF using gmx rmsf"""
        cmd = [
            self.gmx, 'rmsf',
            '-s', topology,
            '-f', trajectory,
            '-n', index_file,
            '-o', output_file,
            '-res'  # Output per-residue RMSF
        ]

        try:
            # Provide group selection via stdin
            input_text = f"{group_name}\n"

            result = subprocess.run(
                cmd,
                input=input_text,
                capture_output=True,
                text=True,
                timeout=300,
                check=True
            )

            # Parse XVG file
            residue_ids = []
            rmsf_values = []

            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith(('#', '@')):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        residue_ids.append(int(parts[0]))
                        rmsf_values.append(float(parts[1]))

            if len(rmsf_values) == 0:
                return {'status': 'failed', 'error': 'No RMSF values parsed'}

            rmsf_values = np.array(rmsf_values)

            return {
                'status': 'success',
                'mean_rmsf_nm': float(np.mean(rmsf_values)),
                'std_rmsf_nm': float(np.std(rmsf_values)),
                'min_rmsf_nm': float(np.min(rmsf_values)),
                'max_rmsf_nm': float(np.max(rmsf_values)),
                'n_residues': len(rmsf_values)
            }

        except subprocess.CalledProcessError as e:
            logger.error(f"gmx rmsf failed: {e.stderr}")
            return {'status': 'failed', 'error': e.stderr}
        except subprocess.TimeoutExpired:
            return {'status': 'failed', 'error': 'Timeout'}
        except Exception as e:
            return {'status': 'failed', 'error': str(e)}

    def process_single_task(self, task: Dict) -> Dict:
        """Process a single task"""
        task_name = task['name']
        pdb_id = task['pdb_id']
        task_output_dir = self.output_dir / task_name
        task_output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"\n{'='*80}")
        logger.info(f"Processing task: {task_name} (PDB: {pdb_id})")
        logger.info(f"{'='*80}")

        try:
            # Step 1: Generate CDR index file
            logger.info(f"[{task_name}] Step 1: Generating CDR index file...")
            index_file = task_output_dir / "cdr_index.ndx"
            cdr_map = self.generate_cdr_index(
                task['pdb_file'],
                pdb_id,
                str(index_file)
            )

            if len(cdr_map) == 0:
                logger.warning(f"[{task_name}] No CDRs found")
                return {
                    'task': task_name,
                    'status': 'failed',
                    'error': 'No CDRs found',
                    'results': []
                }

            # Save CDR mapping
            mapping_file = task_output_dir / "cdr_mapping.json"
            with open(mapping_file, 'w') as f:
                json.dump(cdr_map, f, indent=2)

            # Step 2: Calculate RMSF for each CDR
            logger.info(f"[{task_name}] Step 2: Calculating RMSF...")
            results = []

            for cdr_key, cdr_info in cdr_map.items():
                chain_type, cdr_region = cdr_key.split('_')
                cdr_num = cdr_region[-1]
                group_name = cdr_info['group_name']

                output_xvg = task_output_dir / f"rmsf_{cdr_key}.xvg"

                rmsf_result = self.calculate_rmsf(
                    task['topology'],
                    task['trajectory'],
                    str(index_file),
                    group_name,
                    str(output_xvg)
                )

                if rmsf_result['status'] == 'success':
                    logger.info(f"[{task_name}] {chain_type.upper()} CDR{cdr_num}: "
                              f"RMSF = {rmsf_result['mean_rmsf_nm']:.4f} ± "
                              f"{rmsf_result['std_rmsf_nm']:.4f} nm")

                    results.append({
                        'task': task_name,
                        'pdb_id': pdb_id,
                        'chain': f"TCR_{chain_type}",
                        'cdr_region': f"CDR{cdr_num}",
                        'mean_rmsf_nm': rmsf_result['mean_rmsf_nm'],
                        'std_rmsf_nm': rmsf_result['std_rmsf_nm'],
                        'min_rmsf_nm': rmsf_result['min_rmsf_nm'],
                        'max_rmsf_nm': rmsf_result['max_rmsf_nm'],
                        'n_residues': rmsf_result['n_residues'],
                        'output_file': str(output_xvg),
                        'status': 'success'
                    })
                else:
                    logger.warning(f"[{task_name}] {chain_type.upper()} CDR{cdr_num}: "
                                 f"Failed - {rmsf_result.get('error', 'Unknown error')}")
                    results.append({
                        'task': task_name,
                        'pdb_id': pdb_id,
                        'chain': f"TCR_{chain_type}",
                        'cdr_region': f"CDR{cdr_num}",
                        'status': 'failed',
                        'error': rmsf_result.get('error', 'Unknown error')
                    })

            logger.info(f"[{task_name}] Completed: {len([r for r in results if r.get('status') == 'success'])}/{len(results)} CDRs")

            return {
                'task': task_name,
                'status': 'success',
                'results': results
            }

        except Exception as e:
            logger.error(f"[{task_name}] Task failed: {e}")
            import traceback
            traceback.print_exc()
            return {
                'task': task_name,
                'status': 'failed',
                'error': str(e),
                'results': []
            }

    def run_batch(self) -> pd.DataFrame:
        """Run batch analysis"""
        tasks = self.discover_tasks()

        logger.info("=" * 80)
        logger.info(f"Starting CDR RMSF analysis for {len(tasks)} tasks")
        logger.info(f"Max workers: {self.max_workers}")
        logger.info("=" * 80)

        all_results = []
        completed = 0
        failed = 0

        if self.max_workers == 1:
            # Serial processing
            for i, task in enumerate(tasks, 1):
                logger.info(f"\nProcessing task {i}/{len(tasks)}: {task['name']}")
                task_result = self.process_single_task(task)

                if task_result['status'] == 'success':
                    completed += 1
                    all_results.extend(task_result['results'])
                else:
                    failed += 1
                    logger.error(f"[{i}/{len(tasks)}] {task['name']}: Failed")

        else:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_task = {executor.submit(self.process_single_task, task): task
                                for task in tasks}

                for i, future in enumerate(as_completed(future_to_task), 1):
                    task = future_to_task[future]
                    try:
                        task_result = future.result()
                        if task_result['status'] == 'success':
                            completed += 1
                            all_results.extend(task_result['results'])
                            logger.info(f"[{i}/{len(tasks)}] {task['name']}: Completed "
                                      f"({len([r for r in task_result['results'] if r.get('status') == 'success'])} CDRs)")
                        else:
                            failed += 1
                            logger.error(f"[{i}/{len(tasks)}] {task['name']}: Failed")

                    except Exception as e:
                        failed += 1
                        logger.error(f"[{i}/{len(tasks)}] {task['name']}: Exception - {e}")

        # Convert to DataFrame
        df = pd.DataFrame(all_results)

        # Save summary
        summary_file = self.output_dir / "batch_cdr_rmsf_summary.csv"
        df.to_csv(summary_file, index=False)
        logger.info(f"Summary saved to: {summary_file}")

        # Generate statistics
        self.generate_statistics(df)

        logger.info("\n" + "=" * 80)
        logger.info(f"Batch RMSF analysis completed")
        logger.info(f"Total tasks: {len(tasks)}")
        logger.info(f"Completed: {completed}")
        logger.info(f"Failed: {failed}")
        logger.info(f"Success rate: {completed/len(tasks)*100:.1f}%")
        logger.info(f"Total CDR RMSF records: {len(df[df['status'] == 'success'])}")
        logger.info("=" * 80)

        return df

    def generate_statistics(self, df: pd.DataFrame):
        """Generate statistics report"""
        stats_file = self.output_dir / "batch_cdr_rmsf_statistics.txt"

        with open(stats_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CDR RMSF Analysis Statistics\n")
            f.write("=" * 80 + "\n\n")

            # Check if DataFrame is empty
            if len(df) == 0:
                f.write("No RMSF data available.\n")
                logger.warning("No RMSF data to generate statistics")
                return

            # Filter successful results
            df_success = df[df['status'] == 'success']

            f.write(f"Total RMSF records: {len(df_success)}\n")
            f.write(f"Total tasks analyzed: {df_success['task'].nunique()}\n\n")

            # Statistics by CDR region
            for chain in ['TCR_alpha', 'TCR_beta']:
                for cdr_num in ['1', '2', '3']:
                    region_df = df_success[
                        (df_success['chain'] == chain) &
                        (df_success['cdr_region'] == f'CDR{cdr_num}')
                    ]

                    if len(region_df) > 0:
                        f.write(f"{chain.replace('TCR_', '').upper()} CDR{cdr_num}:\n")
                        f.write(f"  N: {len(region_df)}\n")
                        f.write(f"  Mean RMSF: {region_df['mean_rmsf_nm'].mean():.4f} ± {region_df['mean_rmsf_nm'].std():.4f} nm\n")
                        f.write(f"  Median RMSF: {region_df['mean_rmsf_nm'].median():.4f} nm\n")
                        f.write(f"  Range: {region_df['mean_rmsf_nm'].min():.4f} - {region_df['mean_rmsf_nm'].max():.4f} nm\n\n")

        logger.info(f"Statistics saved to: {stats_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Batch RMSF calculation for CDR loops'
    )
    parser.add_argument(
        '--input-dirs',
        nargs='+',
        default=[
            '/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step',
            '/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step_patch'
        ],
        help='Input directories containing MD trajectories'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='output/rmsf_cdr',
        help='Output directory for RMSF results'
    )
    parser.add_argument(
        '--cdr-csv',
        type=str,
        default='/home/xumy/work/development/AfterMD/input/standardizedpdbs/tcr_cdrs_output.csv',
        help='CDR reference CSV file'
    )
    parser.add_argument(
        '--max-workers',
        type=int,
        default=8,
        help='Maximum number of parallel workers'
    )

    args = parser.parse_args()

    analyzer = CDRRMSFAnalyzer(
        input_dirs=args.input_dirs,
        output_dir=args.output_dir,
        cdr_csv=args.cdr_csv,
        max_workers=args.max_workers
    )

    try:
        analyzer.run_batch()
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
