#!/usr/bin/env python3
"""
Full batch CDR RMSD analysis for all 241 tasks
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.batch_cdr_rmsd_phla_align import BatchCDRRMSDAnalyzer
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    # Two input directories
    input_dirs = [
        "/home/xumy/work/development/Immunex/input/pbc_1000frames_2step",
        "/home/xumy/work/development/Immunex/input/pbc_1000frames_2step_patch"
    ]

    output_dir = "/home/xumy/work/development/Immunex/output/cdr_rmsd_phla_align_full"
    cdr3_ranges_file = "/home/xumy/work/development/Immunex/cdr3_analysis_toolkit/cdr3_ranges.json"

    print("="*80)
    print("Full Batch CDR RMSD Analysis")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Input directories:")
    for d in input_dirs:
        print(f"  - {d}")
    print(f"Output directory: {output_dir}")
    print(f"CDR3 ranges file: {cdr3_ranges_file}")
    print(f"Analysis: All CDR loops (CDR1/2/3) with ANARCI integration")
    print("="*80)
    print()

    # Initialize analyzer
    analyzer = BatchCDRRMSDAnalyzer(
        input_dirs=input_dirs,
        output_dir=output_dir,
        cdr3_ranges_file=cdr3_ranges_file,
        max_workers=1,  # Sequential processing for stability
        analyze_cdr123=True  # Enable all CDR loops
    )

    # Run batch analysis
    print("Starting batch analysis...")
    df_results = analyzer.run_batch()

    print()
    print("="*80)
    print("Batch Analysis Complete!")
    print("="*80)
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total calculations: {len(df_results)}")
    print(f"Successful: {len(df_results[df_results['status'] == 'success'])}")
    print(f"Failed: {len(df_results[df_results['status'] != 'success'])}")
    print()

    if len(df_results) > 0:
        print("Statistics by CDR region:")
        for cdr_region in ['CDR1', 'CDR2', 'CDR3']:
            cdr_data = df_results[df_results['cdr_region'] == cdr_region]
            if len(cdr_data) > 0:
                print(f"\n{cdr_region}:")
                print(f"  N = {len(cdr_data)} calculations")
                print(f"  Mean RMSD = {cdr_data['mean_rmsd_nm'].mean():.4f} +/- {cdr_data['mean_rmsd_nm'].std():.4f} nm")
                print(f"  Range: {cdr_data['mean_rmsd_nm'].min():.4f} - {cdr_data['mean_rmsd_nm'].max():.4f} nm")

    print()
    print("="*80)
    print(f"Results saved to: {output_dir}")
    print("="*80)

if __name__ == "__main__":
    main()
