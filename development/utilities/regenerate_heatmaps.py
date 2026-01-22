#!/usr/bin/env python3
"""Regenerate correlation heatmaps with NaN fix"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.allostery import ContactCorrelationAnalyzer

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def regenerate_heatmap(task_dir):
    """Regenerate heatmap for a single task"""
    task_name = task_dir.name

    try:
        npz_file = task_dir / "correlation_matrix.npz"
        labels_file = task_dir / "contact_labels.csv"

        if not npz_file.exists() or not labels_file.exists():
            return f"{task_name}: SKIP (missing files)"

        # Load data
        data = np.load(npz_file)
        correlation_matrix = data['correlation_matrix']

        labels_df = pd.read_csv(labels_file)
        contact_labels = list(labels_df[['ResA_ID', 'ResA_Name', 'ResB_ID', 'ResB_Name']].itertuples(index=False, name=None))

        # Create temporary analyzer
        analyzer = ContactCorrelationAnalyzer.__new__(ContactCorrelationAnalyzer)

        # Regenerate heatmap
        output_file = task_dir / "correlation_heatmap.png"
        analyzer.plot_correlation_heatmap(
            correlation_matrix=correlation_matrix,
            contact_labels=contact_labels,
            output_file=str(output_file)
        )

        return f"{task_name}: SUCCESS"

    except Exception as e:
        return f"{task_name}: FAILED ({str(e)})"


def main():
    output_dir = Path("output/allostery_analysis")
    task_dirs = [d for d in sorted(output_dir.iterdir()) if d.is_dir()]

    logger.info(f"Found {len(task_dirs)} tasks")
    logger.info("Regenerating heatmaps with NaN fix...")

    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(regenerate_heatmap, task_dir) for task_dir in task_dirs]

        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            if (i-1) % 20 == 0:
                logger.info(f"Progress: {i}/{len(task_dirs)}")

    logger.info("Completed!")


if __name__ == "__main__":
    main()
