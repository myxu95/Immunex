#!/usr/bin/env python3
"""
Batch TCR RMSD Analysis - Aligned to pHLA
Calculate RMSD for whole TCR (alpha + beta chains) with pHLA alignment
"""
import sys
from pathlib import Path
import json
import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import MDAnalysis as mda

sys.path.insert(0, str(Path(__file__).parent.parent))

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class TCRRMSDAnalyzer:
    """
    Analyze whole TCR RMSD with pHLA alignment
    """

    def __init__(self, input_dirs: List[str], output_dir: str, max_workers: int = 8):
        self.input_dirs = [Path(d) for d in input_dirs]
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.gmx = 'gmx'

    def discover_tasks(self) -> List[Dict]:
        """Discover all valid MD tasks"""
        tasks = []

        for input_dir in self.input_dirs:
            if not input_dir.exists():
                logger.warning(f"Directory not found: {input_dir}")
                continue

            for task_dir in sorted(input_dir.iterdir()):
                if not task_dir.is_dir():
                    continue

                task_name = task_dir.name

                # Find required files
                trajectory_files = list(task_dir.glob("*_processed.xtc"))
                pdb_file = task_dir / "md_converted.pdb"
                topology = task_dir / "md.tpr"

                if not trajectory_files or not pdb_file.exists() or not topology.exists():
                    continue

                trajectory = trajectory_files[0]

                tasks.append({
                    'name': task_name,
                    'pdb_id': task_name[:4].upper(),
                    'topology': str(topology),
                    'trajectory': str(trajectory),
                    'pdb_file': str(pdb_file),
                    'task_dir': str(task_dir)
                })

        logger.info(f"Discovered {len(tasks)} valid MD tasks")
        return tasks

    def identify_chains(self, pdb_file: str) -> Dict[str, str]:
        """Identify chain IDs based on residue counts"""
        try:
            u = mda.Universe(pdb_file)
            chains = [c for c in set(u.atoms.chainIDs) if c]

            chain_info = {}
            for chain_id in chains:
                chain_atoms = u.select_atoms(f'chainID {chain_id}')
                ca_atoms = chain_atoms.select_atoms('name CA')
                chain_info[chain_id] = ca_atoms.n_atoms

            # Sort chains by residue count
            sorted_chains = sorted(chain_info.items(), key=lambda x: x[1])

            result = {}
            # Assign chains based on residue count
            for chain_id, res_count in sorted_chains:
                if res_count < 15 and 'peptide' not in result:
                    result['peptide'] = chain_id
                elif 90 < res_count < 110 and 'b2m' not in result:
                    result['b2m'] = chain_id
                elif 100 < res_count < 150 and 'tcr_alpha' not in result:
                    result['tcr_alpha'] = chain_id
                elif 150 <= res_count < 250 and 'tcr_beta' not in result:
                    result['tcr_beta'] = chain_id
                elif res_count >= 250 and 'mhc_alpha' not in result:
                    result['mhc_alpha'] = chain_id

            return result

        except Exception as e:
            logger.error(f"Chain identification failed: {e}")
            return {
                'mhc_alpha': 'A',
                'b2m': 'B',
                'peptide': 'C',
                'tcr_alpha': 'D',
                'tcr_beta': 'E'
            }

    def generate_index_file(self, pdb_file: str, chain_mapping: Dict, output_file: str):
        """Generate GROMACS index file for pHLA and TCR"""
        u = mda.Universe(pdb_file)

        with open(output_file, 'w') as f:
            # Group 1: pHLA (alignment group)
            try:
                mhc_chain = chain_mapping.get('mhc_alpha', 'A')
                b2m_chain = chain_mapping.get('b2m', 'B')
                peptide_chain = chain_mapping.get('peptide', 'C')

                phla_selection = f'chainID {mhc_chain} {b2m_chain} {peptide_chain}'
                phla_atoms = u.select_atoms(phla_selection)

                f.write("[ pHLA ]\n")
                indices = [a.index + 1 for a in phla_atoms]
                for i in range(0, len(indices), 10):
                    f.write(" ".join(map(str, indices[i:i+10])) + "\n")
                f.write("\n")

                logger.info(f"  pHLA: {len(indices)} atoms")

            except Exception as e:
                logger.error(f"Failed to select pHLA: {e}")
                raise

            # Group 2: TCR (calculation group)
            try:
                tcr_alpha_chain = chain_mapping.get('tcr_alpha', 'D')
                tcr_beta_chain = chain_mapping.get('tcr_beta', 'E')

                tcr_selection = f'chainID {tcr_alpha_chain} {tcr_beta_chain}'
                tcr_atoms = u.select_atoms(tcr_selection)

                f.write("[ TCR ]\n")
                indices = [a.index + 1 for a in tcr_atoms]
                for i in range(0, len(indices), 10):
                    f.write(" ".join(map(str, indices[i:i+10])) + "\n")
                f.write("\n")

                logger.info(f"  TCR: {len(indices)} atoms")

            except Exception as e:
                logger.error(f"Failed to select TCR: {e}")
                raise

            # Group 3: TCR CA atoms only
            try:
                tcr_ca_selection = f'chainID {tcr_alpha_chain} {tcr_beta_chain} and name CA'
                tcr_ca_atoms = u.select_atoms(tcr_ca_selection)

                f.write("[ TCR_CA ]\n")
                indices = [a.index + 1 for a in tcr_ca_atoms]
                for i in range(0, len(indices), 10):
                    f.write(" ".join(map(str, indices[i:i+10])) + "\n")
                f.write("\n")

                logger.info(f"  TCR_CA: {len(indices)} atoms")

            except Exception as e:
                logger.error(f"Failed to select TCR CA: {e}")
                raise

        logger.info(f"Index file generated: {output_file}")

    def calculate_rmsd(self, topology: str, trajectory: str, index_file: str,
                      align_group: str, calc_group: str, output_file: str) -> Dict:
        """Calculate RMSD using gmx rms"""
        cmd = [
            self.gmx, 'rms',
            '-s', topology,
            '-f', trajectory,
            '-n', index_file,
            '-o', output_file,
            '-tu', 'ns'
        ]

        try:
            # Provide group selections via stdin
            input_text = f"{align_group}\n{calc_group}\n"

            result = subprocess.run(
                cmd,
                input=input_text,
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0:
                logger.error(f"gmx rms failed: {result.stderr}")
                return {'status': 'failed', 'error': result.stderr}

            # Parse XVG file
            times = []
            rmsd_values = []

            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('@'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        times.append(float(parts[0]))
                        rmsd_values.append(float(parts[1]))

            if len(rmsd_values) == 0:
                return {'status': 'failed', 'error': 'No RMSD values parsed'}

            # Convert to numpy arrays
            times = np.array(times)
            rmsd_values = np.array(rmsd_values)

            # Exclude first 20 ns equilibration period
            equilibration_time_ns = 20.0
            mask = times >= equilibration_time_ns

            if np.sum(mask) > 0:
                equilibrated_rmsd = rmsd_values[mask]
                mean_rmsd = np.mean(equilibrated_rmsd)
                std_rmsd = np.std(equilibrated_rmsd)
                n_frames_used = len(equilibrated_rmsd)
            else:
                # If trajectory is shorter than 20 ns, use all data
                mean_rmsd = np.mean(rmsd_values)
                std_rmsd = np.std(rmsd_values)
                n_frames_used = len(rmsd_values)

            return {
                'status': 'success',
                'mean_rmsd_nm': float(mean_rmsd),
                'std_rmsd_nm': float(std_rmsd),
                'n_frames_total': len(rmsd_values),
                'n_frames_used': n_frames_used,
                'equilibration_time_ns': equilibration_time_ns
            }

        except subprocess.TimeoutExpired:
            return {'status': 'failed', 'error': 'Timeout'}
        except Exception as e:
            return {'status': 'failed', 'error': str(e)}

    def process_single_task(self, task: Dict) -> Dict:
        """Process a single task"""
        task_name = task['name']
        task_output_dir = self.output_dir / task_name
        task_output_dir.mkdir(exist_ok=True)

        logger.info(f"\n{'='*80}")
        logger.info(f"Processing task: {task_name} (PDB: {task['pdb_id']})")
        logger.info(f"{'='*80}")

        try:
            # Step 1: Identify chains
            logger.info(f"[{task_name}] Step 1: Identifying chains...")
            chain_mapping = self.identify_chains(task['pdb_file'])

            chain_file = task_output_dir / "chain_mapping.json"
            with open(chain_file, 'w') as f:
                json.dump(chain_mapping, f, indent=2)

            # Step 2: Generate index file
            logger.info(f"[{task_name}] Step 2: Generating index file...")
            index_file = task_output_dir / "tcr_phla_index.ndx"
            self.generate_index_file(
                task['pdb_file'],
                chain_mapping,
                str(index_file)
            )

            # Step 3: Calculate RMSD for whole TCR
            logger.info(f"[{task_name}] Step 3: Calculating TCR RMSD...")

            results = []

            # All atoms RMSD
            output_file = task_output_dir / "rmsd_tcr_all.xvg"
            rmsd_result = self.calculate_rmsd(
                task['topology'],
                task['trajectory'],
                str(index_file),
                align_group='pHLA',
                calc_group='TCR',
                output_file=str(output_file)
            )

            if rmsd_result['status'] == 'success':
                logger.info(f"[{task_name}] TCR (all atoms): RMSD = {rmsd_result['mean_rmsd_nm']:.4f} ± {rmsd_result['std_rmsd_nm']:.4f} nm (used {rmsd_result['n_frames_used']}/{rmsd_result['n_frames_total']} frames after {rmsd_result['equilibration_time_ns']} ns equilibration)")
                results.append({
                    'task': task_name,
                    'pdb_id': task['pdb_id'],
                    'selection': 'TCR_all_atoms',
                    'mean_rmsd_nm': rmsd_result['mean_rmsd_nm'],
                    'std_rmsd_nm': rmsd_result['std_rmsd_nm'],
                    'n_frames_total': rmsd_result['n_frames_total'],
                    'n_frames_used': rmsd_result['n_frames_used'],
                    'equilibration_time_ns': rmsd_result['equilibration_time_ns'],
                    'status': 'success'
                })
            else:
                logger.warning(f"[{task_name}] TCR (all atoms): Failed - {rmsd_result.get('error', 'Unknown error')}")

            # CA atoms only RMSD
            output_file = task_output_dir / "rmsd_tcr_ca.xvg"
            rmsd_result = self.calculate_rmsd(
                task['topology'],
                task['trajectory'],
                str(index_file),
                align_group='pHLA',
                calc_group='TCR_CA',
                output_file=str(output_file)
            )

            if rmsd_result['status'] == 'success':
                logger.info(f"[{task_name}] TCR (CA atoms): RMSD = {rmsd_result['mean_rmsd_nm']:.4f} ± {rmsd_result['std_rmsd_nm']:.4f} nm (used {rmsd_result['n_frames_used']}/{rmsd_result['n_frames_total']} frames after {rmsd_result['equilibration_time_ns']} ns equilibration)")
                results.append({
                    'task': task_name,
                    'pdb_id': task['pdb_id'],
                    'selection': 'TCR_CA_atoms',
                    'mean_rmsd_nm': rmsd_result['mean_rmsd_nm'],
                    'std_rmsd_nm': rmsd_result['std_rmsd_nm'],
                    'n_frames_total': rmsd_result['n_frames_total'],
                    'n_frames_used': rmsd_result['n_frames_used'],
                    'equilibration_time_ns': rmsd_result['equilibration_time_ns'],
                    'status': 'success'
                })
            else:
                logger.warning(f"[{task_name}] TCR (CA atoms): Failed - {rmsd_result.get('error', 'Unknown error')}")

            logger.info(f"[{task_name}] Completed: {len(results)}/2 calculations")

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
        logger.info(f"Starting TCR RMSD analysis for {len(tasks)} tasks")
        logger.info(f"Max workers: {self.max_workers}")
        logger.info("=" * 80)

        all_results = []
        completed = 0
        failed = 0

        if self.max_workers == 1:
            # Serial processing
            for i, task in enumerate(tasks, 1):
                logger.info(f"\nProcessing task {i}/{len(tasks)}: {task['name']}")
                try:
                    task_result = self.process_single_task(task)

                    if task_result['status'] == 'success':
                        all_results.extend(task_result['results'])
                        completed += 1
                    else:
                        failed += 1

                except Exception as e:
                    logger.error(f"Task {task['name']} failed: {e}")
                    failed += 1
        else:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_task = {executor.submit(self.process_single_task, task): task for task in tasks}

                for i, future in enumerate(as_completed(future_to_task), 1):
                    task = future_to_task[future]
                    task_name = task['name']

                    try:
                        task_result = future.result()
                        if task_result['status'] == 'success':
                            all_results.extend(task_result['results'])
                            completed += 1
                            logger.info(f"[{i}/{len(tasks)}] {task_name}: Completed ({len(task_result['results'])} calculations)")
                        else:
                            failed += 1
                            logger.warning(f"[{i}/{len(tasks)}] {task_name}: Failed")
                    except Exception as e:
                        logger.error(f"[{i}/{len(tasks)}] {task_name}: Exception - {e}")
                        failed += 1

        logger.info("=" * 80)
        logger.info(f"Batch processing complete!")
        logger.info(f"  Completed: {completed}/{len(tasks)} tasks")
        logger.info(f"  Failed: {failed}/{len(tasks)} tasks")
        logger.info("=" * 80)

        # Save results
        df = pd.DataFrame(all_results)

        if len(df) > 0:
            summary_file = self.output_dir / "batch_tcr_rmsd_summary.csv"
            df.to_csv(summary_file, index=False)
            logger.info(f"Summary saved to: {summary_file}")

            # Generate statistics report
            self._generate_statistics_report(df)
        else:
            logger.warning("No results to save!")

        return df

    def _generate_statistics_report(self, df: pd.DataFrame):
        """Generate statistics report"""
        report_file = self.output_dir / "batch_tcr_rmsd_statistics.txt"
        success_df = df[df['status'] == 'success']

        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("Whole TCR RMSD Analysis - Statistics Report\n")
            f.write("Alignment: pHLA | Analysis: TCR (alpha + beta chains)\n")
            f.write("=" * 80 + "\n\n")

            total = len(df)
            successful = len(success_df)
            failed = len(df[df['status'] == 'failed'])

            f.write(f"Total calculations: {total}\n")
            if total > 0:
                f.write(f"Successful: {successful} ({successful/total*100:.1f}%)\n")
                f.write(f"Failed: {failed} ({failed/total*100:.1f}%)\n\n")

            if successful > 0:
                f.write("-" * 80 + "\n")
                f.write("Statistics by Selection\n")
                f.write("-" * 80 + "\n\n")

                for selection in ['TCR_all_atoms', 'TCR_CA_atoms']:
                    sel_df = success_df[success_df['selection'] == selection]

                    if len(sel_df) > 0:
                        f.write(f"{selection}:\n")
                        f.write(f"  N = {len(sel_df)} tasks\n")
                        f.write(f"  Mean RMSD: {sel_df['mean_rmsd_nm'].mean():.4f} ± {sel_df['mean_rmsd_nm'].std():.4f} nm\n")
                        f.write(f"  Median RMSD: {sel_df['mean_rmsd_nm'].median():.4f} nm\n")
                        f.write(f"  Range: {sel_df['mean_rmsd_nm'].min():.4f} - {sel_df['mean_rmsd_nm'].max():.4f} nm\n\n")

        logger.info(f"Statistics report saved to: {report_file}")


def main():
    """Main function"""
    import argparse

    parser = argparse.ArgumentParser(description='Batch TCR RMSD Analysis (aligned to pHLA)')
    parser.add_argument('--max-workers', type=int, default=8, help='Number of parallel workers (default: 8)')
    args = parser.parse_args()

    input_dirs = [
        "/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step",
        "/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step_patch"
    ]

    output_dir = "/home/xumy/work/development/AfterMD/output/tcr_rmsd_phla_align"

    analyzer = TCRRMSDAnalyzer(
        input_dirs=input_dirs,
        output_dir=output_dir,
        max_workers=args.max_workers
    )

    df_results = analyzer.run_batch()

    logger.info(f"\nBatch processing complete!")
    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"Total successful calculations: {len(df_results[df_results['status'] == 'success'])}")


if __name__ == "__main__":
    main()
