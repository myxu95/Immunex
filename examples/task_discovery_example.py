#!/usr/bin/env python3
"""
Example: Task Discovery and Manifest Generation

This script demonstrates how to use the Immunex task discovery module
to scan a batch root directory and generate manifest files.
"""

from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import discover_tasks, write_manifest, discovery_report_to_contexts


def main():
    # Example usage
    print("=" * 80)
    print("Immunex Task Discovery Example")
    print("=" * 80)

    # For this example, create a mock batch structure
    batch_root = Path("./test_batch")

    if not batch_root.exists():
        print(f"\nCreating example batch structure at: {batch_root}")
        create_example_batch(batch_root)

    # Step 1: Discover tasks
    print(f"\nDiscovering tasks in: {batch_root}")
    report = discover_tasks(batch_root)

    # Step 2: Print summary
    report.print_summary()

    # Step 3: Show detailed results
    print("\nDetailed Results:")
    print("-" * 80)

    if report.valid_tasks:
        print(f"\nValid Tasks ({len(report.valid_tasks)}):")
        for task in report.valid_tasks:
            print(f"  - {task.task_id}")
            print(f"    Structure: {Path(task.input_files.structure_path).name}")
            print(f"    Topology:  {Path(task.input_files.topology_path).name}")
            print(f"    Trajectory: {Path(task.input_files.trajectory_path).name}")

    if report.invalid_tasks:
        print(f"\nInvalid Tasks ({len(report.invalid_tasks)}):")
        for task in report.invalid_tasks:
            print(f"  - {task.task_id}")
            for msg in task.validation_messages:
                print(f"    ERROR: {msg}")

    if report.ambiguous_tasks:
        print(f"\nAmbiguous Tasks ({len(report.ambiguous_tasks)}):")
        for task in report.ambiguous_tasks:
            print(f"  - {task.task_id}")
            for msg in task.validation_messages:
                print(f"    WARNING: {msg}")

    # Step 4: Convert valid tasks for legacy PipelineContext-based pipelines
    contexts = discovery_report_to_contexts(report)
    print(f"\nConverted {len(contexts)} valid tasks to PipelineContext objects")
    if contexts:
        print(f"  - First context system_id: {contexts[0].system_id}")
        print(f"  - First context output_dir: {contexts[0].output_dir}")

    # Step 5: Write manifests
    print("\nWriting manifests...")
    manifest_dir = Path("./manifests")
    manifest_dir.mkdir(exist_ok=True)

    write_manifest(report, manifest_dir / "tasks.jsonl", format="jsonl")
    write_manifest(report, manifest_dir / "tasks.csv", format="csv")

    print(f"  - JSONL: {manifest_dir / 'tasks.jsonl'}")
    print(f"  - CSV:   {manifest_dir / 'tasks.csv'}")

    print("\nDone!")
    print("=" * 80)


def create_example_batch(batch_root: Path):
    """Create example batch structure for demonstration."""
    batch_root.mkdir(exist_ok=True)

    # Task 1: Valid task with input/ directory
    task1 = batch_root / "system_001"
    task1_input = task1 / "input"
    task1_input.mkdir(parents=True, exist_ok=True)
    (task1_input / "structure.pdb").touch()
    (task1_input / "topology.tpr").touch()
    (task1_input / "trajectory.xtc").touch()

    # Task 2: Valid task with files in root
    task2 = batch_root / "system_002"
    task2.mkdir(exist_ok=True)
    (task2 / "protein.pdb").touch()
    (task2 / "md.tpr").touch()
    (task2 / "md.xtc").touch()

    # Task 3: Invalid task (missing trajectory)
    task3 = batch_root / "system_003"
    task3.mkdir(exist_ok=True)
    (task3 / "structure.pdb").touch()
    (task3 / "topology.tpr").touch()

    # Task 4: Ambiguous task (multiple trajectories)
    task4 = batch_root / "system_004"
    task4.mkdir(exist_ok=True)
    (task4 / "complex.pdb").touch()
    (task4 / "system.tpr").touch()
    (task4 / "run1.xtc").touch()
    (task4 / "run2.xtc").touch()

    # Task 5: Task with task.yaml
    task5 = batch_root / "system_005"
    task5_data = task5 / "data"
    task5_data.mkdir(parents=True, exist_ok=True)
    (task5_data / "my_structure.pdb").touch()
    (task5_data / "my_topology.tpr").touch()
    (task5_data / "my_trajectory.xtc").touch()

    # Write task.yaml
    import yaml
    task_config = {
        'task_id': 'custom_system',
        'input_files': {
            'structure': 'data/my_structure.pdb',
            'topology': 'data/my_topology.tpr',
            'trajectory': 'data/my_trajectory.xtc'
        },
        'metadata': {
            'system_type': 'protein-ligand',
            'temperature': 300,
            'pressure': 1.0
        },
        'tags': ['example', 'custom']
    }
    with open(task5 / "task.yaml", 'w') as f:
        yaml.dump(task_config, f)

    print(f"Created example batch with 5 tasks:")
    print(f"  - system_001: Valid (input/ directory)")
    print(f"  - system_002: Valid (root directory)")
    print(f"  - system_003: Invalid (missing trajectory)")
    print(f"  - system_004: Ambiguous (multiple trajectories)")
    print(f"  - system_005: Valid (with task.yaml)")


if __name__ == "__main__":
    main()
