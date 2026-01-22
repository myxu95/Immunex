#!/usr/bin/env python3
"""
Rescue script to supplement missing md.tpr files in processed trajectories.

This script scans processed trajectory directories and copies md.tpr files
from the original MD production directories.
"""

import sys
from pathlib import Path
import shutil
import logging
import argparse
from typing import List, Tuple, Dict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class TPRRescuer:
    """Rescue missing md.tpr files in processed trajectory directories."""

    def __init__(self, original_md_root: str, processed_root: str, dry_run: bool = False):
        """
        Initialize TPR rescuer.

        Args:
            original_md_root: Root directory of original MD production results
            processed_root: Root directory of processed trajectories
            dry_run: If True, only show what would be copied without actually copying
        """
        self.original_md_root = Path(original_md_root)
        self.processed_root = Path(processed_root)
        self.dry_run = dry_run

        if not self.original_md_root.exists():
            raise ValueError(f"Original MD root not found: {self.original_md_root}")
        if not self.processed_root.exists():
            raise ValueError(f"Processed root not found: {self.processed_root}")

    def find_processed_dirs(self) -> List[Path]:
        """
        Find all processed trajectory directories.

        Returns:
            List of paths to processed directories
        """
        processed_dirs = []

        for item in self.processed_root.rglob('*'):
            if item.is_dir():
                # Check if directory contains processed trajectory
                has_processed_xtc = any(item.glob('*_processed.xtc'))
                has_gro = any(item.glob('md.gro'))

                if has_processed_xtc or has_gro:
                    processed_dirs.append(item)

        logger.info(f"Found {len(processed_dirs)} processed trajectory directories")
        return processed_dirs

    def find_original_tpr(self, processed_dir: Path) -> Path:
        """
        Find the original md.tpr file corresponding to a processed directory.

        Args:
            processed_dir: Path to processed directory

        Returns:
            Path to original md.tpr file

        Raises:
            FileNotFoundError: If original tpr file cannot be found
        """
        # Extract job identifier from directory name
        # e.g., 5w1w_fixed_1_processed -> 5w1w_fixed_1
        # e.g., 5bs0_run1_processed -> 5bs0_run1
        # e.g., 6eqb_1_processed -> 6eqb_1
        dir_name = processed_dir.name

        # Remove '_processed' suffix
        if dir_name.endswith('_processed'):
            job_id = dir_name.replace('_processed', '')
        else:
            job_id = dir_name

        # Search for original md.tpr file in the job directory
        # Priority 1: original_md_root/job_id/prod/md.tpr
        # Priority 2: original_md_root/job_id/md.tpr
        # Priority 3: Recursive search in job_id directory
        tpr_candidates = [
            self.original_md_root / job_id / "prod" / "md.tpr",
            self.original_md_root / job_id / "production" / "md.tpr",
            self.original_md_root / job_id / "md.tpr",
        ]

        for tpr_file in tpr_candidates:
            if tpr_file.exists():
                return tpr_file

        # Fallback: Recursive search in the job directory
        job_dir = self.original_md_root / job_id
        if job_dir.exists():
            for tpr_file in job_dir.rglob("md.tpr"):
                return tpr_file

        raise FileNotFoundError(f"Cannot find md.tpr for {job_id} in {self.original_md_root}")

    def copy_tpr_file(self, source: Path, dest_dir: Path) -> bool:
        """
        Copy md.tpr file to destination directory.

        Args:
            source: Source tpr file path
            dest_dir: Destination directory

        Returns:
            True if copied successfully, False otherwise
        """
        dest_file = dest_dir / "md.tpr"

        if dest_file.exists():
            logger.info(f"  TPR already exists: {dest_file}")
            return False

        if self.dry_run:
            logger.info(f"  [DRY RUN] Would copy: {source} -> {dest_file}")
            return True

        try:
            shutil.copy2(source, dest_file)
            logger.info(f"  Copied: {source} -> {dest_file}")
            return True
        except Exception as e:
            logger.error(f"  Failed to copy {source}: {e}")
            return False

    def rescue_all(self) -> Dict[str, int]:
        """
        Rescue all missing md.tpr files.

        Returns:
            Statistics dictionary
        """
        processed_dirs = self.find_processed_dirs()

        stats = {
            'total': len(processed_dirs),
            'copied': 0,
            'already_exists': 0,
            'failed': 0
        }

        logger.info(f"\nStarting TPR rescue operation...")
        logger.info(f"Mode: {'DRY RUN' if self.dry_run else 'LIVE'}\n")

        for i, proc_dir in enumerate(processed_dirs, 1):
            logger.info(f"[{i}/{len(processed_dirs)}] Processing: {proc_dir.name}")

            try:
                original_tpr = self.find_original_tpr(proc_dir)
                logger.info(f"  Found original: {original_tpr}")

                if self.copy_tpr_file(original_tpr, proc_dir):
                    if (proc_dir / "md.tpr").exists():
                        stats['already_exists'] += 1
                    else:
                        stats['copied'] += 1

            except FileNotFoundError as e:
                logger.warning(f"  {e}")
                stats['failed'] += 1
            except Exception as e:
                logger.error(f"  Unexpected error: {e}")
                stats['failed'] += 1

        return stats

    def print_summary(self, stats: Dict[str, int]):
        """Print summary statistics."""
        logger.info("\n" + "="*60)
        logger.info("RESCUE OPERATION SUMMARY")
        logger.info("="*60)
        logger.info(f"Total directories processed: {stats['total']}")
        logger.info(f"TPR files copied: {stats['copied']}")
        logger.info(f"Already exists: {stats['already_exists']}")
        logger.info(f"Failed: {stats['failed']}")
        logger.info("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Rescue missing md.tpr files in processed trajectories'
    )
    parser.add_argument(
        '--original-md-root',
        required=True,
        help='Root directory of original MD production results'
    )
    parser.add_argument(
        '--processed-root',
        required=True,
        help='Root directory of processed trajectories'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be copied without actually copying'
    )

    args = parser.parse_args()

    try:
        rescuer = TPRRescuer(
            original_md_root=args.original_md_root,
            processed_root=args.processed_root,
            dry_run=args.dry_run
        )

        stats = rescuer.rescue_all()
        rescuer.print_summary(stats)

        if args.dry_run:
            logger.info("\nThis was a DRY RUN. Use without --dry-run to actually copy files.")

    except Exception as e:
        logger.error(f"Fatal error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
