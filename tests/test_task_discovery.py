"""
Tests for Immunex Task Discovery and Manifest Generation

Tests cover:
1. Valid task discovery (input/ directory and root directory)
2. Invalid task discovery (missing files)
3. Ambiguous task discovery (multiple candidates)
4. Task.yaml parsing and file resolution
5. Manifest generation (JSONL and CSV)
6. Error handling
"""

import pytest
from pathlib import Path
import tempfile
import shutil
import yaml
import json
import csv

from immunex.core import (
    TaskInputFiles,
    DiscoveredTask,
    DiscoveryReport,
    TaskDiscoverer,
    ManifestWriter,
    discover_tasks,
    discover_tasks_from_list,
    write_manifest,
    discovered_task_to_context,
    discovery_report_to_contexts,
    TaskDiscoveryError,
    TaskFileParseError
)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def temp_batch_root(tmp_path):
    """Create a temporary batch root directory."""
    batch_root = tmp_path / "batch_root"
    batch_root.mkdir()
    return batch_root


@pytest.fixture
def valid_task_standard(temp_batch_root):
    """Create a valid task with input/ directory."""
    task_dir = temp_batch_root / "task_001"
    input_dir = task_dir / "input"
    input_dir.mkdir(parents=True)

    # Create required files
    (input_dir / "structure.pdb").touch()
    (input_dir / "topology.tpr").touch()
    (input_dir / "trajectory.xtc").touch()

    return task_dir


@pytest.fixture
def valid_task_root(temp_batch_root):
    """Create a valid task with files in root directory."""
    task_dir = temp_batch_root / "task_002"
    task_dir.mkdir()

    # Create required files in root
    (task_dir / "protein.pdb").touch()
    (task_dir / "md.tpr").touch()
    (task_dir / "md.xtc").touch()

    return task_dir


@pytest.fixture
def invalid_task_missing_trajectory(temp_batch_root):
    """Create an invalid task missing trajectory."""
    task_dir = temp_batch_root / "task_003"
    input_dir = task_dir / "input"
    input_dir.mkdir(parents=True)

    # Only create structure and topology
    (input_dir / "structure.pdb").touch()
    (input_dir / "topology.tpr").touch()

    return task_dir


@pytest.fixture
def ambiguous_task_multiple_trajectories(temp_batch_root):
    """Create an ambiguous task with multiple trajectory candidates."""
    task_dir = temp_batch_root / "task_004"
    input_dir = task_dir / "input"
    input_dir.mkdir(parents=True)

    # Create required files
    (input_dir / "structure.pdb").touch()
    (input_dir / "topology.tpr").touch()

    # Create multiple trajectory candidates
    (input_dir / "md.xtc").touch()
    (input_dir / "prod.xtc").touch()

    return task_dir


@pytest.fixture
def task_with_yaml(temp_batch_root):
    """Create a task with task.yaml specifying explicit paths."""
    task_dir = temp_batch_root / "task_005"
    task_dir.mkdir()

    # Create task.yaml
    task_config = {
        'task_id': 'custom_task_id',
        'input_files': {
            'structure': 'data/my_structure.pdb',
            'topology': 'data/my_topology.tpr',
            'trajectory': 'data/my_trajectory.xtc'
        },
        'metadata': {
            'system_type': 'protein-ligand',
            'temperature': 300
        },
        'tags': ['test', 'example']
    }

    # Create data directory and files
    data_dir = task_dir / "data"
    data_dir.mkdir()
    (data_dir / "my_structure.pdb").touch()
    (data_dir / "my_topology.tpr").touch()
    (data_dir / "my_trajectory.xtc").touch()

    # Write task.yaml
    with open(task_dir / "task.yaml", 'w') as f:
        yaml.dump(task_config, f)

    return task_dir


# ============================================================================
# Test Task Discovery
# ============================================================================

def test_discover_valid_task_input_dir(valid_task_standard, temp_batch_root):
    """Test discovery of valid task with input/ directory."""
    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 1
    assert report.num_valid == 1
    assert report.num_invalid == 0
    assert report.num_ambiguous == 0

    task = report.valid_tasks[0]
    assert task.task_id == "task_001"
    assert task.validation_status == "valid"
    assert task.input_files.structure_path is not None
    assert task.input_files.topology_path is not None
    assert task.input_files.trajectory_path is not None
    assert len(task.validation_messages) == 0


def test_discover_valid_task_root_dir(valid_task_root, temp_batch_root):
    """Test discovery of valid task with files in root directory."""
    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 1
    assert report.num_valid == 1

    task = report.valid_tasks[0]
    assert task.task_id == "task_002"
    assert task.validation_status == "valid"
    assert "protein.pdb" in task.input_files.structure_path
    assert "md.tpr" in task.input_files.topology_path
    assert "md.xtc" in task.input_files.trajectory_path


def test_discover_invalid_task_missing_trajectory(
    invalid_task_missing_trajectory,
    temp_batch_root
):
    """Test discovery of invalid task missing trajectory."""
    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 1
    assert report.num_valid == 0
    assert report.num_invalid == 1
    assert report.num_ambiguous == 0

    task = report.invalid_tasks[0]
    assert task.task_id == "task_003"
    assert task.validation_status == "invalid"
    assert "Missing trajectory file" in task.validation_messages[0]


def test_discover_ambiguous_task_multiple_trajectories(
    ambiguous_task_multiple_trajectories,
    temp_batch_root
):
    """Test discovery of ambiguous task with multiple trajectory candidates."""
    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 1
    assert report.num_valid == 0
    assert report.num_invalid == 0
    assert report.num_ambiguous == 1

    task = report.ambiguous_tasks[0]
    assert task.task_id == "task_004"
    assert task.validation_status == "ambiguous"
    assert any("Ambiguous trajectory" in msg for msg in task.validation_messages)


def test_discover_task_with_yaml(task_with_yaml, temp_batch_root):
    """Test discovery of task with task.yaml."""
    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 1
    assert report.num_valid == 1

    task = report.valid_tasks[0]
    assert task.task_id == "custom_task_id"  # Overridden by task.yaml
    assert task.validation_status == "valid"
    assert task.task_file_path is not None
    assert "my_structure.pdb" in task.input_files.structure_path
    assert "my_topology.tpr" in task.input_files.topology_path
    assert "my_trajectory.xtc" in task.input_files.trajectory_path

    # Check metadata and tags
    assert task.metadata['system_type'] == 'protein-ligand'
    assert task.metadata['temperature'] == 300
    assert 'test' in task.tags
    assert 'example' in task.tags


def test_discover_multiple_tasks(temp_batch_root):
    """Test discovery of multiple tasks with mixed statuses."""
    # Create valid task
    valid_dir = temp_batch_root / "valid_task"
    valid_input = valid_dir / "input"
    valid_input.mkdir(parents=True)
    (valid_input / "structure.pdb").touch()
    (valid_input / "topology.tpr").touch()
    (valid_input / "trajectory.xtc").touch()

    # Create invalid task
    invalid_dir = temp_batch_root / "invalid_task"
    invalid_dir.mkdir()
    (invalid_dir / "structure.pdb").touch()
    # Missing topology and trajectory

    # Create ambiguous task
    ambiguous_dir = temp_batch_root / "ambiguous_task"
    ambiguous_dir.mkdir()
    (ambiguous_dir / "structure.pdb").touch()
    (ambiguous_dir / "topology.tpr").touch()
    (ambiguous_dir / "traj1.xtc").touch()
    (ambiguous_dir / "traj2.xtc").touch()

    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 3
    assert report.num_valid == 1
    assert report.num_invalid == 1
    assert report.num_ambiguous == 1
    assert report.success_rate == pytest.approx(1.0 / 3.0)


# ============================================================================
# Test Error Handling
# ============================================================================

def test_discover_nonexistent_root():
    """Test that discovering nonexistent root raises TaskDiscoveryError."""
    nonexistent = Path("/nonexistent/path/to/batch")

    with pytest.raises(TaskDiscoveryError) as exc_info:
        discover_tasks(nonexistent)

    assert "does not exist" in str(exc_info.value)


def test_discover_file_instead_of_directory(tmp_path):
    """Test that discovering a file instead of directory raises TaskDiscoveryError."""
    file_path = tmp_path / "file.txt"
    file_path.touch()

    with pytest.raises(TaskDiscoveryError) as exc_info:
        discover_tasks(file_path)

    assert "not a directory" in str(exc_info.value)


def test_parse_invalid_yaml(temp_batch_root):
    """Test that invalid YAML in task.yaml is handled gracefully."""
    task_dir = temp_batch_root / "bad_yaml_task"
    task_dir.mkdir()

    # Write invalid YAML
    with open(task_dir / "task.yaml", 'w') as f:
        f.write("invalid: yaml: content: [unclosed")

    # Create valid input files
    (task_dir / "structure.pdb").touch()
    (task_dir / "topology.tpr").touch()
    (task_dir / "trajectory.xtc").touch()

    # Should not crash - should fall back to auto-discovery
    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 1
    # Task should still be valid due to auto-discovery fallback
    assert report.num_valid == 1


# ============================================================================
# Test Manifest Generation
# ============================================================================

def test_write_manifest_jsonl(valid_task_standard, temp_batch_root, tmp_path):
    """Test writing manifest in JSONL format."""
    report = discover_tasks(temp_batch_root)

    manifest_path = tmp_path / "manifest.jsonl"
    write_manifest(report, manifest_path, format="jsonl")

    assert manifest_path.exists()

    # Read and verify
    with open(manifest_path, 'r') as f:
        lines = f.readlines()

    assert len(lines) == 1

    task_data = json.loads(lines[0])
    assert task_data['task_id'] == "task_001"
    assert task_data['validation_status'] == "valid"
    assert 'input_files' in task_data
    assert 'structure_path' in task_data['input_files']


def test_write_manifest_csv(valid_task_standard, temp_batch_root, tmp_path):
    """Test writing manifest in CSV format."""
    report = discover_tasks(temp_batch_root)

    manifest_path = tmp_path / "manifest.csv"
    write_manifest(report, manifest_path, format="csv")

    assert manifest_path.exists()

    # Read and verify
    with open(manifest_path, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]['task_id'] == "task_001"
    assert rows[0]['validation_status'] == "valid"
    assert rows[0]['structure_path'] != ''


def test_write_manifest_unsupported_format(valid_task_standard, temp_batch_root, tmp_path):
    """Test that unsupported format raises ValueError."""
    report = discover_tasks(temp_batch_root)

    manifest_path = tmp_path / "manifest.xml"

    with pytest.raises(ValueError) as exc_info:
        write_manifest(report, manifest_path, format="xml")

    assert "Unsupported format" in str(exc_info.value)


def test_write_manifest_creates_parent_dirs(valid_task_standard, temp_batch_root, tmp_path):
    """Test that write_manifest creates parent directories."""
    report = discover_tasks(temp_batch_root)

    nested_path = tmp_path / "deep" / "nested" / "path" / "manifest.jsonl"

    # Should not raise error
    write_manifest(report, nested_path, format="jsonl")

    assert nested_path.exists()


# ============================================================================
# Test Data Models
# ============================================================================

def test_task_input_files_completeness():
    """Test TaskInputFiles completeness checks."""
    # Complete
    complete_files = TaskInputFiles(
        structure_path="/path/to/structure.pdb",
        topology_path="/path/to/topology.tpr",
        trajectory_path="/path/to/trajectory.xtc"
    )
    assert complete_files.is_complete()
    assert len(complete_files.missing_files()) == 0

    # Incomplete
    incomplete_files = TaskInputFiles(
        structure_path="/path/to/structure.pdb"
    )
    assert not incomplete_files.is_complete()
    missing = incomplete_files.missing_files()
    assert "topology" in missing
    assert "trajectory" in missing


def test_discovered_task_status_checks():
    """Test DiscoveredTask status check methods."""
    valid_task = DiscoveredTask(
        task_id="test_valid",
        task_root="/path/to/task",
        source_root="/path/to/batch",
        input_files=TaskInputFiles(),
        validation_status="valid"
    )
    assert valid_task.is_valid()
    assert not valid_task.is_invalid()
    assert not valid_task.is_ambiguous()

    invalid_task = DiscoveredTask(
        task_id="test_invalid",
        task_root="/path/to/task",
        source_root="/path/to/batch",
        input_files=TaskInputFiles(),
        validation_status="invalid"
    )
    assert invalid_task.is_invalid()
    assert not invalid_task.is_valid()

    ambiguous_task = DiscoveredTask(
        task_id="test_ambiguous",
        task_root="/path/to/task",
        source_root="/path/to/batch",
        input_files=TaskInputFiles(),
        validation_status="ambiguous"
    )
    assert ambiguous_task.is_ambiguous()
    assert not ambiguous_task.is_valid()


def test_discovery_report_properties():
    """Test DiscoveryReport computed properties."""
    valid_task = DiscoveredTask(
        task_id="valid",
        task_root="/path/to/valid",
        source_root="/batch",
        input_files=TaskInputFiles(),
        validation_status="valid"
    )

    invalid_task = DiscoveredTask(
        task_id="invalid",
        task_root="/path/to/invalid",
        source_root="/batch",
        input_files=TaskInputFiles(),
        validation_status="invalid"
    )

    report = DiscoveryReport(
        source_root="/batch",
        discovered_at="2024-01-01T00:00:00",
        total_tasks=2,
        valid_tasks=[valid_task],
        invalid_tasks=[invalid_task],
        ambiguous_tasks=[],
        all_tasks=[valid_task, invalid_task]
    )

    assert report.num_valid == 1
    assert report.num_invalid == 1
    assert report.num_ambiguous == 0
    assert report.success_rate == 0.5


def test_discovered_task_serialization():
    """Test DiscoveredTask to_dict and from_dict."""
    task = DiscoveredTask(
        task_id="test_task",
        task_root="/path/to/task",
        source_root="/path/to/batch",
        input_files=TaskInputFiles(
            structure_path="/path/to/structure.pdb",
            topology_path="/path/to/topology.tpr",
            trajectory_path="/path/to/trajectory.xtc"
        ),
        validation_status="valid",
        tags=["test", "example"],
        metadata={"key": "value"}
    )

    # Serialize
    task_dict = task.to_dict()
    assert task_dict['task_id'] == "test_task"
    assert task_dict['input_files']['structure_path'] == "/path/to/structure.pdb"

    # Deserialize
    reconstructed = DiscoveredTask.from_dict(task_dict)
    assert reconstructed.task_id == task.task_id
    assert reconstructed.input_files.structure_path == task.input_files.structure_path
    assert reconstructed.tags == task.tags
    assert reconstructed.metadata == task.metadata


# ============================================================================
# Integration Tests
# ============================================================================

def test_full_workflow(temp_batch_root, tmp_path):
    """Test complete workflow: discover -> validate -> write manifest."""
    # Setup: Create 3 tasks with different statuses
    # Valid task
    valid_dir = temp_batch_root / "project_001"
    valid_input = valid_dir / "input"
    valid_input.mkdir(parents=True)
    (valid_input / "complex.pdb").touch()
    (valid_input / "system.tpr").touch()
    (valid_input / "traj.xtc").touch()

    # Invalid task
    invalid_dir = temp_batch_root / "project_002"
    invalid_dir.mkdir()
    (invalid_dir / "structure.pdb").touch()

    # Ambiguous task
    ambiguous_dir = temp_batch_root / "project_003"
    ambiguous_dir.mkdir()
    (ambiguous_dir / "struct.pdb").touch()
    (ambiguous_dir / "topol.tpr").touch()
    (ambiguous_dir / "run1.xtc").touch()
    (ambiguous_dir / "run2.xtc").touch()

    # Step 1: Discover
    report = discover_tasks(temp_batch_root)

    assert report.total_tasks == 3
    assert report.num_valid == 1
    assert report.num_invalid == 1
    assert report.num_ambiguous == 1

    # Step 2: Print summary
    report.print_summary()  # Should not crash

    # Step 3: Write manifests
    jsonl_path = tmp_path / "tasks.jsonl"
    csv_path = tmp_path / "tasks.csv"

    write_manifest(report, jsonl_path, format="jsonl")
    write_manifest(report, csv_path, format="csv")

    assert jsonl_path.exists()
    assert csv_path.exists()

    # Verify JSONL content
    with open(jsonl_path, 'r') as f:
        lines = f.readlines()
    assert len(lines) == 3

    # Verify CSV content
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert len(rows) == 3

    # Check that validation statuses are correct
    statuses = {row['task_id']: row['validation_status'] for row in rows}
    assert statuses['project_001'] == 'valid'
    assert statuses['project_002'] == 'invalid'
    assert statuses['project_003'] == 'ambiguous'


def test_discovered_task_to_context(valid_task_standard, temp_batch_root):
    """Valid discovered task should convert to PipelineContext."""
    report = discover_tasks(temp_batch_root)
    task = report.valid_tasks[0]

    context = discovered_task_to_context(task, output_base_dir='./custom_results')

    assert context.system_id == task.task_id
    assert context.topology == task.input_files.topology_path
    assert context.trajectory_raw == task.input_files.trajectory_path
    assert context.structure_pdb == task.input_files.structure_path
    assert context.output_dir.endswith(f"custom_results/{task.task_id}")
    assert context.metadata['discovery_contract'] == 'manifest'


def test_discovery_report_to_contexts_filters_valid_tasks(temp_batch_root, valid_task_standard, valid_task_root, invalid_task_missing_trajectory):
    """Only valid tasks should be adapted by default."""
    report = discover_tasks(temp_batch_root)

    contexts = discovery_report_to_contexts(report)

    assert len(contexts) == report.num_valid
    assert all(context.system_id in {task.task_id for task in report.valid_tasks} for context in contexts)


def test_discover_tasks_with_task_depth(tmp_path):
    root = tmp_path / 'nested_root'
    task_dir = root / '1ao7' / 'standard'
    task_dir.mkdir(parents=True)
    (task_dir / 'structure.pdb').touch()
    (task_dir / 'md.gro').touch()
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks(root, task_depth=2)

    assert report.total_tasks == 1
    assert report.num_valid == 1
    assert report.valid_tasks[0].task_id == 'standard'


def test_discover_tasks_from_list():
    report = discover_tasks_from_list([
        {
            'task_id': 'manual_task',
            'topology': 'missing_topology.tpr',
            'trajectory_raw': 'missing_trajectory.xtc',
        },
        {
            'task_id': 'incomplete_task',
            'topology': 'only_topology.tpr',
        },
    ])

    assert report.total_tasks == 2
    assert report.num_invalid == 2
    assert report.num_valid == 0



def test_discover_tasks_from_list_resolves_relative_paths(tmp_path):
    source_root = tmp_path / 'manual_batch'
    task_dir = source_root / 'manual_task'
    task_dir.mkdir(parents=True)
    (task_dir / 'md.gro').touch()
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks_from_list([
        {
            'task_id': 'manual_task',
            'structure': 'manual_task/md.gro',
            'topology': 'manual_task/md.tpr',
            'trajectory_raw': 'manual_task/md.xtc',
        },
    ], source_root=str(source_root))

    assert report.total_tasks == 1
    assert report.num_valid == 1
    task = report.valid_tasks[0]
    assert task.input_files.structure_path == str((task_dir / 'md.gro').absolute())
    assert task.input_files.topology_path == str((task_dir / 'md.tpr').absolute())
    assert task.input_files.trajectory_path == str((task_dir / 'md.xtc').absolute())
    assert task.task_root == str(task_dir.absolute())


def test_discover_tasks_from_list_bridges_to_contexts(tmp_path):
    source_root = tmp_path / 'manual_contexts'
    task_dir = source_root / 'task_a'
    task_dir.mkdir(parents=True)
    (task_dir / 'md.gro').touch()
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks_from_list([
        {
            'task_id': 'task_a',
            'structure': 'task_a/md.gro',
            'topology': 'task_a/md.tpr',
            'trajectory_raw': 'task_a/md.xtc',
            'metadata': {'batch': 'manual'},
        },
    ], source_root=str(source_root))
    contexts = discovery_report_to_contexts(report, output_base_dir=str(tmp_path / 'results'))

    assert len(contexts) == 1
    assert contexts[0].system_id == 'task_a'
    assert contexts[0].metadata['discovery_contract'] == 'manifest'
    assert contexts[0].metadata['batch'] == 'manual'



def test_discover_tasks_with_stage_specific_required_files(tmp_path):
    root = tmp_path / 'stage_specific'
    task_dir = root / 'task_001'
    task_dir.mkdir(parents=True)
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks(root, required_files=['topology', 'trajectory'])

    assert report.total_tasks == 1
    assert report.num_valid == 1
    assert report.num_invalid == 0



def test_discover_tasks_from_list_with_stage_specific_required_files(tmp_path):
    source_root = tmp_path / 'manual_stage_specific'
    task_dir = source_root / 'task_001'
    task_dir.mkdir(parents=True)
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks_from_list([
        {
            'task_id': 'task_001',
            'topology': 'task_001/md.tpr',
            'trajectory_raw': 'task_001/md.xtc',
        },
    ], source_root=str(source_root), required_files=['topology', 'trajectory'])

    assert report.total_tasks == 1
    assert report.num_valid == 1
    assert report.valid_tasks[0].task_id == 'task_001'



def test_discover_tasks_supports_prod_subdirectories(tmp_path):
    root = tmp_path / 'prod_layout'
    task_dir = root / 'task_001' / 'prod'
    task_dir.mkdir(parents=True)
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks(root, required_files=['topology', 'trajectory'])

    assert report.total_tasks == 1
    assert report.num_valid == 1
    assert report.valid_tasks[0].task_id == 'task_001'
    assert Path(report.valid_tasks[0].input_files.trajectory_path).parent.name == 'prod'
