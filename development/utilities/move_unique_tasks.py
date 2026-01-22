#!/usr/bin/env python3
"""
Move Unique Tasks Script

Move unique incomplete tasks (those without completed versions) to a separate directory.

Usage:
    python scripts/move_unique_tasks.py /path/to/incomplete /path/to/unique --task-list unique_tasks.txt
    python scripts/move_unique_tasks.py /path/to/incomplete /path/to/unique --auto-detect
"""

import argparse
import shutil
import re
from pathlib import Path
from typing import Set, List


def extract_pdb_id(task_name: str) -> str:
    """Extract PDB ID from task name (first 4 characters, case-insensitive)."""
    match = re.match(r'^([a-zA-Z0-9]{4})', task_name)
    if match:
        return match.group(1).lower()
    else:
        return task_name[:4].lower()


def load_task_list(task_list_file: Path) -> Set[str]:
    """Load task names from a text file."""
    task_names = set()

    try:
        with open(task_list_file, 'r') as f:
            for line in f:
                task_name = line.strip()
                if task_name:
                    task_names.add(task_name)

        print(f"📋 Loaded {len(task_names)} task names from: {task_list_file}")
        return task_names

    except Exception as e:
        print(f"❌ Failed to load task list: {e}")
        return set()


def move_tasks(source_dir: Path,
               dest_dir: Path,
               task_names: Set[str],
               dry_run: bool = False) -> dict:
    """Move specified tasks from source to destination directory."""

    if not dry_run:
        dest_dir.mkdir(parents=True, exist_ok=True)

    stats = {"moved": 0, "not_found": 0, "failed": 0, "skipped": 0}

    for task_name in sorted(task_names):
        source_path = source_dir / task_name
        dest_path = dest_dir / task_name

        # Check if source exists
        if not source_path.exists():
            print(f"⚠️  Task not found: {task_name}")
            stats["not_found"] += 1
            continue

        # Check if destination already exists
        if dest_path.exists():
            print(f"⚠️  Destination exists, skipping: {task_name}")
            stats["skipped"] += 1
            continue

        # Move the task
        try:
            if dry_run:
                print(f"[DRY RUN] Would move: {source_path} -> {dest_path}")
            else:
                shutil.move(str(source_path), str(dest_path))
                print(f"✅ Moved: {task_name}")

            stats["moved"] += 1

        except Exception as e:
            print(f"❌ Failed to move {task_name}: {e}")
            stats["failed"] += 1

    return stats


def auto_detect_unique_tasks(incomplete_dir: Path, completed_dir: Path) -> Set[str]:
    """Auto-detect unique tasks by comparing PDB IDs."""

    # Get PDB IDs from completed directory
    completed_pdbs = set()
    if completed_dir.exists():
        for item in completed_dir.iterdir():
            if item.is_dir():
                pdb_id = extract_pdb_id(item.name)
                completed_pdbs.add(pdb_id)

    # Find unique tasks in incomplete directory
    unique_tasks = set()
    for item in incomplete_dir.iterdir():
        if item.is_dir():
            task_name = item.name
            pdb_id = extract_pdb_id(task_name)

            if pdb_id not in completed_pdbs:
                unique_tasks.add(task_name)

    return unique_tasks


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Move unique incomplete tasks to a separate directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Move tasks from a task list file
    python scripts/move_unique_tasks.py /data/incomplete /data/unique --task-list unique_tasks.txt

    # Auto-detect unique tasks by comparing with completed directory
    python scripts/move_unique_tasks.py /data/incomplete /data/unique --completed-dir /data/completed

    # Dry run to preview operations
    python scripts/move_unique_tasks.py /data/incomplete /data/unique --task-list unique_tasks.txt --dry-run

    # Real example:
    python scripts/move_unique_tasks.py \\
        /public/home/xmy/dataset/Standard_Stimulation/total_traj/Human_ClassI_imcomplete \\
        /public/home/xmy/dataset/Standard_Stimulation/total_traj/Human_ClassI_unique \\
        --completed-dir /public/home/xmy/dataset/Standard_Stimulation/total_traj/Human_ClassI
        """
    )

    # Required arguments
    parser.add_argument(
        "source_dir",
        help="Source directory containing incomplete tasks"
    )

    parser.add_argument(
        "dest_dir",
        help="Destination directory for unique tasks"
    )

    # Task selection methods (mutually exclusive)
    task_group = parser.add_mutually_exclusive_group(required=True)

    task_group.add_argument(
        "--task-list",
        type=Path,
        help="Text file containing task names to move (one per line)"
    )

    task_group.add_argument(
        "--completed-dir",
        type=Path,
        help="Completed tasks directory for auto-detection of unique tasks"
    )

    # Options
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview operations without making changes"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show verbose output"
    )

    args = parser.parse_args()

    # Validate paths
    source_dir = Path(args.source_dir)
    dest_dir = Path(args.dest_dir)

    if not source_dir.exists():
        print(f"❌ Source directory does not exist: {source_dir}")
        return 1

    if not source_dir.is_dir():
        print(f"❌ Source path is not a directory: {source_dir}")
        return 1

    print("📦 Move Unique Tasks")
    print("=" * 50)
    print(f"📂 Source: {source_dir}")
    print(f"📁 Destination: {dest_dir}")

    if args.dry_run:
        print("🧪 DRY RUN MODE - No changes will be made")

    # Determine which tasks to move
    if args.task_list:
        print(f"📋 Using task list: {args.task_list}")
        task_names = load_task_list(args.task_list)

        if not task_names:
            print("❌ No valid task names found in task list")
            return 1

    elif args.completed_dir:
        print(f"🔍 Auto-detecting unique tasks (comparing with: {args.completed_dir})")

        if not args.completed_dir.exists():
            print(f"❌ Completed directory does not exist: {args.completed_dir}")
            return 1

        task_names = auto_detect_unique_tasks(source_dir, args.completed_dir)

        if not task_names:
            print("✨ No unique tasks found - all incomplete tasks have completed versions")
            return 0

        print(f"🆕 Found {len(task_names)} unique tasks")

    # Show tasks to be moved (if verbose or dry run)
    if args.verbose or args.dry_run:
        print(f"\n📋 Tasks to be moved ({len(task_names)}):")
        for i, task_name in enumerate(sorted(task_names), 1):
            print(f"   {i:3d}. {task_name}")

    if not args.dry_run:
        # Confirm operation
        response = input(f"\n❓ Move {len(task_names)} tasks to {dest_dir}? (y/N): ")
        if response.lower() not in ['y', 'yes']:
            print("❌ Operation cancelled")
            return 0

    # Move tasks
    print(f"\n📦 Moving {len(task_names)} unique tasks...")
    stats = move_tasks(source_dir, dest_dir, task_names, args.dry_run)

    # Show results
    print(f"\n📊 Operation Results:")
    print(f"   ✅ Moved: {stats['moved']}")
    print(f"   ❌ Failed: {stats['failed']}")
    print(f"   ⚠️  Not found: {stats['not_found']}")
    print(f"   ⏭️  Skipped: {stats['skipped']}")

    if stats["moved"] > 0 and not args.dry_run:
        print(f"\n✨ Successfully moved {stats['moved']} unique tasks to: {dest_dir}")
    elif args.dry_run:
        print(f"\n💡 Run without --dry-run to execute the move operation")

    return 0 if stats["failed"] == 0 else 1


if __name__ == "__main__":
    exit(main())