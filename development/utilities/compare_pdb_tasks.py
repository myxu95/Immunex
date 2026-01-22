#!/usr/bin/env python3
"""
PDB Task Comparison Script

Compare PDB IDs between completed and incomplete task directories
to find unique incomplete tasks that need to be processed.

Usage:
    python scripts/compare_pdb_tasks.py /path/to/completed /path/to/incomplete
    python scripts/compare_pdb_tasks.py /path/to/completed /path/to/incomplete --output unique_tasks.txt
"""

import argparse
import re
from pathlib import Path
from typing import Set, List, Dict


def extract_pdb_id(task_name: str) -> str:
    """
    Extract PDB ID from task name (first 4 characters, case-insensitive).

    Args:
        task_name: Task directory name

    Returns:
        PDB ID in lowercase
    """
    match = re.match(r'^([a-zA-Z0-9]{4})', task_name)
    if match:
        return match.group(1).lower()
    else:
        return task_name[:4].lower()


def get_pdb_ids_from_directory(directory: Path) -> Dict[str, List[str]]:
    """
    Get all PDB IDs from task directories.

    Args:
        directory: Directory containing task folders

    Returns:
        Dictionary mapping PDB IDs to list of task names
    """
    pdb_mapping = {}

    if not directory.exists():
        print(f"❌ Directory does not exist: {directory}")
        return pdb_mapping

    for item in directory.iterdir():
        if item.is_dir():
            task_name = item.name
            pdb_id = extract_pdb_id(task_name)

            if pdb_id not in pdb_mapping:
                pdb_mapping[pdb_id] = []
            pdb_mapping[pdb_id].append(task_name)

    return pdb_mapping


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Compare PDB IDs between completed and incomplete task directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic comparison
    python scripts/compare_pdb_tasks.py /data/completed_tasks /data/incomplete_tasks

    # Save unique incomplete tasks to file
    python scripts/compare_pdb_tasks.py /data/completed /data/incomplete --output unique.txt

    # Show detailed mapping
    python scripts/compare_pdb_tasks.py /data/completed /data/incomplete --verbose
        """
    )

    parser.add_argument(
        "completed_dir",
        help="Directory containing completed tasks"
    )

    parser.add_argument(
        "incomplete_dir",
        help="Directory containing incomplete tasks"
    )

    parser.add_argument(
        "--output", "-o",
        help="Output file to save unique incomplete task names"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed information"
    )

    args = parser.parse_args()

    # Validate paths
    completed_dir = Path(args.completed_dir)
    incomplete_dir = Path(args.incomplete_dir)

    if not completed_dir.exists():
        print(f"❌ Completed directory does not exist: {completed_dir}")
        return 1

    if not incomplete_dir.exists():
        print(f"❌ Incomplete directory does not exist: {incomplete_dir}")
        return 1

    print("🔍 PDB Task Comparison")
    print("=" * 50)
    print(f"✅ Completed tasks: {completed_dir}")
    print(f"❌ Incomplete tasks: {incomplete_dir}")
    print()

    # Get PDB mappings
    print("📊 Analyzing directories...")
    completed_pdbs = get_pdb_ids_from_directory(completed_dir)
    incomplete_pdbs = get_pdb_ids_from_directory(incomplete_dir)

    print(f"✅ Completed PDB IDs: {len(completed_pdbs)}")
    print(f"❌ Incomplete PDB IDs: {len(incomplete_pdbs)}")

    # Find unique incomplete tasks
    completed_pdb_set = set(completed_pdbs.keys())
    incomplete_pdb_set = set(incomplete_pdbs.keys())

    unique_incomplete_pdbs = incomplete_pdb_set - completed_pdb_set
    overlapping_pdbs = incomplete_pdb_set & completed_pdb_set

    print()
    print("📊 Comparison Results:")
    print(f"🆕 Unique incomplete PDBs: {len(unique_incomplete_pdbs)}")
    print(f"🔄 Overlapping PDBs: {len(overlapping_pdbs)}")

    # Show unique incomplete tasks
    if unique_incomplete_pdbs:
        print()
        print("🆕 Unique Incomplete Tasks (need processing):")
        print("-" * 60)

        unique_task_names = []
        for pdb_id in sorted(unique_incomplete_pdbs):
            task_names = incomplete_pdbs[pdb_id]
            for task_name in task_names:
                unique_task_names.append(task_name)
                print(f"   {pdb_id.upper()}: {task_name}")

        print(f"\n📋 Total unique incomplete tasks: {len(unique_task_names)}")

        # Save to file if requested
        if args.output:
            output_file = Path(args.output)
            try:
                with open(output_file, 'w') as f:
                    for task_name in sorted(unique_task_names):
                        f.write(f"{task_name}\n")
                print(f"💾 Saved unique task names to: {output_file}")
            except Exception as e:
                print(f"❌ Failed to save output file: {e}")

    else:
        print("\n✨ No unique incomplete tasks found!")
        print("All incomplete tasks have corresponding completed versions.")

    # Show overlapping tasks if verbose
    if args.verbose and overlapping_pdbs:
        print()
        print("🔄 Overlapping PDBs (already completed):")
        print("-" * 60)
        for pdb_id in sorted(overlapping_pdbs):
            completed_tasks = completed_pdbs[pdb_id]
            incomplete_tasks = incomplete_pdbs[pdb_id]

            print(f"   {pdb_id.upper()}:")
            print(f"      ✅ Completed: {', '.join(completed_tasks)}")
            print(f"      ❌ Incomplete: {', '.join(incomplete_tasks)}")

    # Show summary
    print()
    print("=" * 50)
    print("📋 Summary:")
    print(f"   📁 Total directories scanned: {len(completed_pdbs) + len(incomplete_pdbs)}")
    print(f"   ✅ Completed PDBs: {len(completed_pdbs)}")
    print(f"   ❌ Incomplete PDBs: {len(incomplete_pdbs)}")
    print(f"   🆕 Unique incomplete: {len(unique_incomplete_pdbs)}")
    print(f"   🔄 Already completed: {len(overlapping_pdbs)}")

    if unique_incomplete_pdbs:
        print()
        print("💡 Next steps:")
        print("   1. Process the unique incomplete tasks listed above")
        print("   2. Consider removing overlapping incomplete tasks (duplicates)")
        if args.output:
            print(f"   3. Use the task list saved in {args.output}")

    return 0


if __name__ == "__main__":
    exit(main())