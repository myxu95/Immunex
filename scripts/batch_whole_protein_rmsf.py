#!/usr/bin/env python3
"""
Calculate whole protein RMSF and extract CDR loop regions
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


class WholeProteinRMSFAnalyzer:
    """Calculate whole protein RMSF and extract CDR regions"""

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

    def calculate_whole_protein_rmsf(self, topology: str, trajectory: str,
                                     output_file: str) -> Dict:
        """Calculate RMSF for whole protein (all CA atoms)"""
        cmd = [
            self.gmx, 'rmsf',
            '-s', topology,
            '-f', trajectory,
            '-o', output_file,
            '-res'  # Per-residue RMSF
        ]

        try:
            # Select C-alpha atoms (group 3 in standard GROMACS groups)
            input_text = "C-alpha\n"

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

            return {
                'status': 'success',
                'residue_ids': residue_ids,
                'rmsf_values': rmsf_values,
                'n_residues': len(rmsf_values)
            }

        except subprocess.CalledProcessError as e:
            logger.error(f"gmx rmsf failed: {e.stderr}")
            return {'status': 'failed', 'error': e.stderr}
        except subprocess.TimeoutExpired:
            return {'status': 'failed', 'error': 'Timeout'}
        except Exception as e:
            return {'status': 'failed', 'error': str(e)}

    def extract_cdr_rmsf(self, residue_ids: List[int], rmsf_values: List[float],
                        cdr_info: Dict) -> Dict:
        """Extract RMSF values for CDR regions from whole protein RMSF"""
        cdr_rmsf = {}

        # Convert to numpy arrays for easier indexing
        residue_ids = np.array(residue_ids)
        rmsf_values = np.array(rmsf_values)

        for chain_type in ['alpha', 'beta']:
            for cdr_num in [1, 2, 3]:
                cdr_key = f'cdr{cdr_num}'

                if cdr_key not in cdr_info[chain_type]:
                    continue

                cdr_data = cdr_info[chain_type][cdr_key]
                resid_start, resid_end = cdr_data['residue_range']

                # Find RMSF values for this residue range
                mask = (residue_ids >= resid_start) & (residue_ids <= resid_end)
                cdr_residues = residue_ids[mask]
                cdr_rmsf_vals = rmsf_values[mask]

                if len(cdr_rmsf_vals) > 0:
                    cdr_rmsf[f"{chain_type}_cdr{cdr_num}"] = {
                        'residue_ids': cdr_residues.tolist(),
                        'rmsf_values': cdr_rmsf_vals.tolist(),
                        'mean_rmsf': float(np.mean(cdr_rmsf_vals)),
                        'std_rmsf': float(np.std(cdr_rmsf_vals)),
                        'min_rmsf': float(np.min(cdr_rmsf_vals)),
                        'max_rmsf': float(np.max(cdr_rmsf_vals)),
                        'n_residues': len(cdr_rmsf_vals),
                        'residue_range': [int(resid_start), int(resid_end)]
                    }

        return cdr_rmsf

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
            # Step 1: Get CDR regions
            logger.info(f"[{task_name}] Step 1: Identifying CDR regions...")
            cdr_info = self.cdr_selector.get_cdr_regions(
                task['pdb_file'],
                pdb_id
            )

            # Count CDRs
            n_cdrs = sum(1 for chain in ['alpha', 'beta']
                        for cdr in ['cdr1', 'cdr2', 'cdr3']
                        if cdr in cdr_info[chain])

            if n_cdrs == 0:
                logger.warning(f"[{task_name}] No CDRs found")
                return {
                    'task': task_name,
                    'status': 'no_cdrs',
                    'results': []
                }

            logger.info(f"[{task_name}] Found {n_cdrs} CDR regions")

            # Step 2: Calculate whole protein RMSF
            logger.info(f"[{task_name}] Step 2: Calculating whole protein RMSF...")
            rmsf_file = task_output_dir / "rmsf_whole_protein.xvg"

            rmsf_result = self.calculate_whole_protein_rmsf(
                task['topology'],
                task['trajectory'],
                str(rmsf_file)
            )

            if rmsf_result['status'] != 'success':
                logger.error(f"[{task_name}] RMSF calculation failed: {rmsf_result.get('error')}")
                return {
                    'task': task_name,
                    'status': 'failed',
                    'error': rmsf_result.get('error'),
                    'results': []
                }

            logger.info(f"[{task_name}] Calculated RMSF for {rmsf_result['n_residues']} residues")

            # Step 3: Extract CDR RMSF from whole protein data
            logger.info(f"[{task_name}] Step 3: Extracting CDR RMSF...")
            cdr_rmsf = self.extract_cdr_rmsf(
                rmsf_result['residue_ids'],
                rmsf_result['rmsf_values'],
                cdr_info
            )

            # Save CDR RMSF data
            cdr_rmsf_file = task_output_dir / "cdr_rmsf_extracted.json"
            with open(cdr_rmsf_file, 'w') as f:
                json.dump(cdr_rmsf, f, indent=2)

            # Step 4: Generate summary results
            results = []
            for cdr_key, cdr_data in cdr_rmsf.items():
                chain_type, cdr_region = cdr_key.split('_')
                cdr_num = cdr_region[-1]

                logger.info(f"[{task_name}] {chain_type.upper()} CDR{cdr_num}: "
                          f"RMSF = {cdr_data['mean_rmsf']:.4f} ± {cdr_data['std_rmsf']:.4f} nm")

                results.append({
                    'task': task_name,
                    'pdb_id': pdb_id,
                    'chain': f"TCR_{chain_type}",
                    'cdr_region': f"CDR{cdr_num}",
                    'mean_rmsf_nm': cdr_data['mean_rmsf'],
                    'std_rmsf_nm': cdr_data['std_rmsf'],
                    'min_rmsf_nm': cdr_data['min_rmsf'],
                    'max_rmsf_nm': cdr_data['max_rmsf'],
                    'n_residues': cdr_data['n_residues'],
                    'residue_range': f"{cdr_data['residue_range'][0]}-{cdr_data['residue_range'][1]}",
                    'whole_protein_rmsf_file': str(rmsf_file),
                    'status': 'success'
                })

            logger.info(f"[{task_name}] Completed: {len(results)}/{n_cdrs} CDRs extracted")

            return {
                'task': task_name,
                'status': 'success',
                'n_total_residues': rmsf_result['n_residues'],
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
        logger.info(f"Starting whole protein RMSF analysis for {len(tasks)} tasks")
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
                                      f"({len(task_result['results'])} CDRs)")
                        else:
                            failed += 1
                            logger.error(f"[{i}/{len(tasks)}] {task['name']}: Failed")

                    except Exception as e:
                        failed += 1
                        logger.error(f"[{i}/{len(tasks)}] {task['name']}: Exception - {e}")

        # Convert to DataFrame
        df = pd.DataFrame(all_results)

        # Save summary
        if len(df) > 0:
            summary_file = self.output_dir / "batch_whole_protein_rmsf_summary.csv"
            df.to_csv(summary_file, index=False)
            logger.info(f"Summary saved to: {summary_file}")

            # Generate statistics
            self.generate_statistics(df)

        logger.info("\n" + "=" * 80)
        logger.info(f"Batch analysis completed")
        logger.info(f"Total tasks: {len(tasks)}")
        logger.info(f"Completed: {completed}")
        logger.info(f"Failed: {failed}")
        logger.info(f"Success rate: {completed/len(tasks)*100:.1f}%")
        logger.info(f"Total CDR RMSF records: {len(df)}")
        logger.info("=" * 80)

        return df

    def generate_statistics(self, df: pd.DataFrame):
        """Generate statistics report"""
        stats_file = self.output_dir / "batch_whole_protein_rmsf_statistics.txt"

        with open(stats_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CDR RMSF Statistics (Extracted from Whole Protein RMSF)\n")
            f.write("=" * 80 + "\n\n")

            if len(df) == 0:
                f.write("No RMSF data available.\n")
                return

            # Filter successful results
            df_success = df[df['status'] == 'success']

            f.write(f"Total CDR RMSF records: {len(df_success)}\n")
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
        description='Calculate whole protein RMSF and extract CDR regions'
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

    analyzer = WholeProteinRMSFAnalyzer(
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
