"""Tests for BatchExecutor discovery-native execution paths."""

from pathlib import Path

from immunex.core import PipelineContext, discover_tasks
from immunex.pipeline import BatchExecutor, Pipeline


class DummyPipeline(Pipeline):
    """Simple pipeline used to verify BatchExecutor wiring."""

    def execute(self, context: PipelineContext) -> PipelineContext:
        context.results['dummy'] = {'seen': context.system_id}
        return context



def test_prepare_tasks_accepts_discovery_report(tmp_path):
    batch_root = tmp_path / 'batch'
    task_dir = batch_root / 'task_001'
    task_dir.mkdir(parents=True)
    (task_dir / 'md.gro').touch()
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks(batch_root)
    executor = BatchExecutor(max_workers=1)

    tasks = executor.prepare_tasks(report, output_base_dir=str(tmp_path / 'results'))

    assert len(tasks) == 1
    assert tasks[0].system_id == 'task_001'
    assert tasks[0].output_dir.endswith('results/task_001')
    assert tasks[0].metadata['discovery_contract'] == 'manifest'



def test_execute_pipeline_accepts_discovery_report_serially(tmp_path):
    batch_root = tmp_path / 'batch'
    task_dir = batch_root / 'task_001'
    task_dir.mkdir(parents=True)
    (task_dir / 'md.gro').touch()
    (task_dir / 'md.tpr').touch()
    (task_dir / 'md.xtc').touch()

    report = discover_tasks(batch_root)
    executor = BatchExecutor(max_workers=1)
    results = executor.execute_pipeline(
        report,
        DummyPipeline(),
        show_progress=False,
        output_base_dir=str(tmp_path / 'results'),
    )

    assert len(results) == 1
    assert results[0].system_id == 'task_001'
    assert results[0].results['dummy']['seen'] == 'task_001'
    assert not results[0].has_errors()
