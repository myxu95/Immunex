#!/usr/bin/env python3
"""
CDR RMSD Batch Analysis - Using Exact CDR Sequences from CSV
Based on tcr_cdrs_output.csv for precise CDR region identification
"""
import sys
from pathlib import Path
import json
import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import MDAnalysis as mda

sys.path.insert(0, str(Path(__file__).parent.parent))

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ExactCDRAnalyzer:
    """
    CDR RMSD Analyzer using exact CDR sequences from CSV
    """

    def __init__(self,
                 cdr_csv_file: str,
                 input_dirs: List[str],
                 output_dir: str,
                 max_workers: int = 1,
                 gmx_executable: str = "gmx"):

        self.cdr_csv = Path(cdr_csv_file)
        self.input_dirs = [Path(d) for d in input_dirs]
        self.output_dir = Path(output_dir)
        self.max_workers = max_workers
        self.gmx = gmx_executable

        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load CDR reference data
        self.cdr_data = self._load_cdr_reference()
        logger.info(f"Loaded CDR data for {len(self.cdr_data)} PDB structures")

    def _load_cdr_reference(self) -> Dict:
        """Load CDR sequences from CSV"""
        df = pd.read_csv(self.cdr_csv)

        cdr_dict = {}
        for _, row in df.iterrows():
            pdb_id = row['PDB_ID'].upper()

            cdr_dict[pdb_id] = {
                'alpha': {
                    'sequence': row['alpha_variable_domain'],
                    'cdr1': row['alpha_cdr1'] if pd.notna(row['alpha_cdr1']) else None,
                    'cdr2': row['alpha_cdr2'] if pd.notna(row['alpha_cdr2']) else None,
                    'cdr3': row['alpha_cdr3'] if pd.notna(row['alpha_cdr3']) else None,
                },
                'beta': {
                    'sequence': row['beta_variable_domain'],
                    'cdr1': row['beta_cdr1'] if pd.notna(row['beta_cdr1']) else None,
                    'cdr2': row['beta_cdr2'] if pd.notna(row['beta_cdr2']) else None,
                    'cdr3': row['beta_cdr3'] if pd.notna(row['beta_cdr3']) else None,
                }
            }

        return cdr_dict

    def discover_tasks(self) -> List[Dict]:
        """Discover all valid MD tasks"""
        tasks = []

        for input_dir in self.input_dirs:
            if not input_dir.exists():
                logger.warning(f"Input directory not found: {input_dir}")
                continue

            logger.info(f"Scanning directory: {input_dir}")

            for task_dir in sorted(input_dir.iterdir()):
                if not task_dir.is_dir():
                    continue

                task_name = task_dir.name

                # Extract PDB ID (first 4 characters)
                pdb_id = task_name[:4].upper()

                # Check if we have CDR data for this PDB
                if pdb_id not in self.cdr_data:
                    logger.warning(f"No CDR data for {task_name} (PDB: {pdb_id})")
                    continue

                # Find required files
                trajectory_files = list(task_dir.glob("*_processed.xtc"))
                pdb_file = task_dir / "md_converted.pdb"
                topology = task_dir / "md.tpr"

                if not trajectory_files or not pdb_file.exists() or not topology.exists():
                    continue

                trajectory = trajectory_files[0]

                tasks.append({
                    'name': task_name,
                    'pdb_id': pdb_id,
                    'topology': str(topology),
                    'trajectory': str(trajectory),
                    'pdb_file': str(pdb_file),
                    'task_dir': str(task_dir),
                    'cdr_reference': self.cdr_data[pdb_id]
                })

        logger.info(f"Discovered {len(tasks)} valid MD tasks with CDR data")
        return tasks

    def find_cdr_in_sequence(self, full_sequence: str, cdr_sequence: str) -> Optional[Tuple[int, int]]:
        """
        Find CDR sequence position in full sequence
        Returns: (start_residue, end_residue) in 1-based numbering
        """
        if not cdr_sequence or pd.isna(cdr_sequence):
            return None

        # Find CDR sequence in full sequence
        pos = full_sequence.find(cdr_sequence)

        if pos == -1:
            return None

        # Convert to 1-based residue numbering
        start = pos + 1
        end = pos + len(cdr_sequence)

        return (start, end)

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

    def map_cdr_regions(self, task: Dict, chain_mapping: Dict) -> Dict:
        """
        Map CDR regions using exact sequences from CSV
        Search all chains to find CDR sequences
        """
        pdb_file = task['pdb_file']
        cdr_ref = task['cdr_reference']

        u = mda.Universe(pdb_file)

        # Build AA code dictionary
        aa_code = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }

        # Get all chains with their sequences
        chain_sequences = {}
        for chain_id in set(u.atoms.chainIDs):
            if not chain_id:
                continue
            chain_atoms = u.select_atoms(f'chainID {chain_id} and name CA')
            if len(chain_atoms) > 0:
                seq = ''.join([aa_code.get(r.resname, 'X') for r in chain_atoms.residues])
                chain_sequences[chain_id] = {
                    'sequence': seq,
                    'residues': chain_atoms.residues
                }

        cdr_regions = {}

        for chain_type in ['alpha', 'beta']:
            chain_data = {}

            # For each CDR, search in all chains
            for cdr_num in [1, 2, 3]:
                cdr_key = f'cdr{cdr_num}'
                cdr_seq = cdr_ref[chain_type][cdr_key]

                if not cdr_seq or pd.isna(cdr_seq):
                    continue

                # Search for CDR sequence in all chains
                found = False
                for chain_id, chain_info in chain_sequences.items():
                    pdb_sequence = chain_info['sequence']
                    cdr_range = self.find_cdr_in_sequence(pdb_sequence, cdr_seq)

                    if cdr_range:
                        start, end = cdr_range

                        # Get atom indices for this range
                        cdr_atoms = u.select_atoms(
                            f'chainID {chain_id} and name CA and resid {start}:{end}'
                        )

                        chain_data[cdr_key] = {
                            'chain_id': chain_id,
                            'sequence': cdr_seq,
                            'residue_range': [start, end],
                            'residue_count': end - start + 1,
                            'ca_count': len(cdr_atoms),
                            'ca_indices': [a.index + 1 for a in cdr_atoms]
                        }

                        logger.info(f"  {chain_type.upper()} CDR{cdr_num}: {cdr_seq}")
                        logger.info(f"    Chain {chain_id}, Residues: {start}-{end} ({len(cdr_atoms)} CA atoms)")
                        found = True
                        break

                if not found:
                    logger.warning(f"  {chain_type.upper()} CDR{cdr_num}: Sequence '{cdr_seq}' not found in any chain")

            cdr_regions[chain_type] = chain_data

        return cdr_regions

    def generate_index_file(self, pdb_file: str, chain_mapping: Dict,
                          cdr_regions: Dict, output_file: str):
        """Generate GROMACS index file"""
        u = mda.Universe(pdb_file)

        with open(output_file, 'w') as f:
            # Group 1: pMHC
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
            except Exception as e:
                logger.error(f"Failed to select pMHC: {e}")
                raise

            # CDR groups
            for chain_type in ['alpha', 'beta']:
                chain_data = cdr_regions.get(chain_type, {})

                for cdr_num in [1, 2, 3]:
                    cdr_key = f'cdr{cdr_num}'
                    if cdr_key not in chain_data:
                        continue

                    cdr_info = chain_data[cdr_key]
                    group_name = f"CDR{cdr_num}_{chain_type}_CA"

                    f.write(f"[ {group_name} ]\n")
                    ca_indices = cdr_info['ca_indices']
                    for i in range(0, len(ca_indices), 10):
                        f.write(" ".join(map(str, ca_indices[i:i+10])) + "\n")
                    f.write("\n")

                    logger.info(f"  {group_name}: {len(ca_indices)} CA atoms")

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

        stdin_input = f"{align_group}\n{calc_group}\n"

        subprocess.run(
            cmd,
            input=stdin_input,
            text=True,
            capture_output=True,
            timeout=300,
            check=True
        )

        return self._parse_rmsd_xvg(output_file)

    def _parse_rmsd_xvg(self, xvg_file: str) -> Dict:
        """Parse RMSD XVG file and exclude first 20 ns equilibration"""
        data = []
        with open(xvg_file, 'r') as f:
            for line in f:
                if line.startswith(('#', '@')):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    data.append([float(parts[0]), float(parts[1])])

        data = np.array(data)
        times = data[:, 0]
        rmsd_values = data[:, 1]

        # Exclude first 20 ns equilibration period
        equilibration_time_ns = 20.0
        mask = times >= equilibration_time_ns

        if np.sum(mask) > 0:
            equilibrated_rmsd = rmsd_values[mask]
            n_frames_used = len(equilibrated_rmsd)
        else:
            # If trajectory is shorter than 20 ns, use all data
            equilibrated_rmsd = rmsd_values
            n_frames_used = len(rmsd_values)

        return {
            'mean': float(np.mean(equilibrated_rmsd)),
            'std': float(np.std(equilibrated_rmsd)),
            'min': float(np.min(equilibrated_rmsd)),
            'max': float(np.max(equilibrated_rmsd)),
            'n_frames_total': len(rmsd_values),
            'n_frames_used': n_frames_used,
            'equilibration_time_ns': equilibration_time_ns
        }

    def process_single_task(self, task: Dict) -> Dict:
        """Process single task"""
        task_name = task['name']
        pdb_id = task['pdb_id']
        task_output_dir = self.output_dir / task_name
        task_output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"\n{'='*80}")
        logger.info(f"Processing task: {task_name} (PDB: {pdb_id})")
        logger.info(f"{'='*80}")

        try:
            # Step 1: Identify chains
            logger.info(f"[{task_name}] Step 1: Identifying chains...")
            chain_mapping = self.identify_chains(task['pdb_file'])

            chain_file = task_output_dir / "chain_mapping.json"
            with open(chain_file, 'w') as f:
                json.dump(chain_mapping, f, indent=2)

            # Step 2: Map CDR regions using exact sequences
            logger.info(f"[{task_name}] Step 2: Mapping CDR regions from CSV...")
            cdr_regions = self.map_cdr_regions(task, chain_mapping)

            ranges_file = task_output_dir / "cdr_ranges.json"
            # Convert numpy types to native Python types for JSON serialization
            def convert_numpy_types(obj):
                if isinstance(obj, dict):
                    return {k: convert_numpy_types(v) for k, v in obj.items()}
                elif isinstance(obj, list):
                    return [convert_numpy_types(item) for item in obj]
                elif isinstance(obj, (np.integer, np.int64)):
                    return int(obj)
                elif isinstance(obj, (np.floating, np.float64)):
                    return float(obj)
                return obj

            with open(ranges_file, 'w') as f:
                json.dump(convert_numpy_types(cdr_regions), f, indent=2)

            # Step 3: Generate index file
            logger.info(f"[{task_name}] Step 3: Generating index file...")
            index_file = task_output_dir / "cdr_phla_index.ndx"
            self.generate_index_file(
                task['pdb_file'],
                chain_mapping,
                cdr_regions,
                str(index_file)
            )

            # Step 4: Calculate RMSD
            logger.info(f"[{task_name}] Step 4: Calculating RMSD...")
            results = []

            for chain_type in ['alpha', 'beta']:
                chain_name = f'TCR_{chain_type}'
                chain_data = cdr_regions.get(chain_type, {})

                for cdr_num in [1, 2, 3]:
                    cdr_key = f'cdr{cdr_num}'
                    if cdr_key not in chain_data:
                        logger.info(f"[{task_name}] {chain_name} CDR{cdr_num}: Not available")
                        continue

                    cdr_info = chain_data[cdr_key]
                    group_name = f"CDR{cdr_num}_{chain_type}_CA"
                    output_xvg = task_output_dir / f"rmsd_pHLA_to_CDR{cdr_num}_{chain_type}.xvg"

                    try:
                        rmsd_stats = self.calculate_rmsd(
                            topology=task['topology'],
                            trajectory=task['trajectory'],
                            index_file=str(index_file),
                            align_group="pHLA",
                            calc_group=group_name,
                            output_file=str(output_xvg)
                        )

                        results.append({
                            'task': task_name,
                            'pdb_id': pdb_id,
                            'chain': chain_name,
                            'cdr_region': f'CDR{cdr_num}',
                            'cdr_sequence': cdr_info['sequence'],
                            'residue_range': f"{cdr_info['residue_range'][0]}-{cdr_info['residue_range'][1]}",
                            'mean_rmsd_nm': rmsd_stats['mean'],
                            'std_rmsd_nm': rmsd_stats['std'],
                            'min_rmsd_nm': rmsd_stats['min'],
                            'max_rmsd_nm': rmsd_stats['max'],
                            'n_frames_total': rmsd_stats['n_frames_total'],
                            'n_frames_used': rmsd_stats['n_frames_used'],
                            'equilibration_time_ns': rmsd_stats['equilibration_time_ns'],
                            'output_file': str(output_xvg),
                            'status': 'success'
                        })

                        logger.info(f"[{task_name}] {chain_name} CDR{cdr_num}: "
                                  f"RMSD = {rmsd_stats['mean']:.4f} ± {rmsd_stats['std']:.4f} nm "
                                  f"(used {rmsd_stats['n_frames_used']}/{rmsd_stats['n_frames_total']} frames after {rmsd_stats['equilibration_time_ns']} ns)")

                    except Exception as e:
                        logger.error(f"[{task_name}] {chain_name} CDR{cdr_num} failed: {e}")
                        results.append({
                            'task': task_name,
                            'pdb_id': pdb_id,
                            'chain': chain_name,
                            'cdr_region': f'CDR{cdr_num}',
                            'status': 'failed',
                            'error': str(e)
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
        logger.info(f"Starting exact CDR RMSD analysis for {len(tasks)} tasks")
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
                            logger.info(f"[{i}/{len(tasks)}] {task_name}: Completed ({len(task_result['results'])} CDRs)")
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
            summary_file = self.output_dir / "batch_cdr_rmsd_summary.csv"
            df.to_csv(summary_file, index=False)
            logger.info(f"Summary saved to: {summary_file}")

            # Generate statistics report
            self._generate_statistics_report(df)
        else:
            logger.warning("No results to save!")

        return df

    def _generate_statistics_report(self, df: pd.DataFrame):
        """Generate statistics report"""
        report_file = self.output_dir / "batch_cdr_rmsd_statistics.txt"
        success_df = df[df['status'] == 'success']

        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("Exact CDR RMSD Analysis - Statistics Report\n")
            f.write("CDR regions identified from tcr_cdrs_output.csv\n")
            f.write("Alignment: pMHC | Analysis: CDR Loops (CA atoms)\n")
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
                f.write("Statistics by CDR Region\n")
                f.write("-" * 80 + "\n\n")

                for chain in ['TCR_alpha', 'TCR_beta']:
                    for cdr_num in [1, 2, 3]:
                        cdr_region = f'CDR{cdr_num}'
                        region_df = success_df[(success_df['chain'] == chain) &
                                              (success_df['cdr_region'] == cdr_region)]

                        if len(region_df) > 0:
                            f.write(f"{chain} - {cdr_region}:\n")
                            f.write(f"  N = {len(region_df)} tasks\n")
                            f.write(f"  Mean RMSD: {region_df['mean_rmsd_nm'].mean():.4f} ± {region_df['mean_rmsd_nm'].std():.4f} nm\n")
                            f.write(f"  Median RMSD: {region_df['mean_rmsd_nm'].median():.4f} nm\n")
                            f.write(f"  Range: {region_df['mean_rmsd_nm'].min():.4f} - {region_df['mean_rmsd_nm'].max():.4f} nm\n\n")

        logger.info(f"Statistics report saved to: {report_file}")


def main():
    """Main function"""
    import argparse

    parser = argparse.ArgumentParser(description='CDR RMSD Batch Analysis using exact sequences')
    parser.add_argument('--max-workers', type=int, default=8, help='Number of parallel workers (default: 8)')
    parser.add_argument('--test-single', type=str, help='Test on a single task (provide task name)')
    args = parser.parse_args()

    cdr_csv = "/home/xumy/work/development/AfterMD/input/standardizedpdbs/tcr_cdrs_output.csv"

    input_dirs = [
        "/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step",
        "/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step_patch"
    ]

    output_dir = "/home/xumy/work/development/AfterMD/output/cdr_rmsd_exact_analysis"

    analyzer = ExactCDRAnalyzer(
        cdr_csv_file=cdr_csv,
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
