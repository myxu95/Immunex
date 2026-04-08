"""Tests for the batch workflow compatibility wrapper."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import immunex.pipeline.batch_workflow as batch_workflow
from immunex.core import PipelineContext


class DummyPipeline:
    def __init__(self, *args, **kwargs):
        pass

    def execute(self, context: PipelineContext):
        context.trajectory_processed = str(Path(context.output_dir) / 'md_processed.xtc')
        return context



def test_process_md_tasks_uses_discovery_wrapper(tmp_path, monkeypatch):
    base = tmp_path / 'raw_md' / 'task_001' / 'prod'
    base.mkdir(parents=True)
    (base / 'md.tpr').touch()
    (base / 'md.xtc').touch()

    import immunex.pipeline as pipeline_module
    monkeypatch.setattr(pipeline_module, 'PreprocessOnlyPipeline', DummyPipeline)

    results = batch_workflow.process_md_tasks(str(base.parent.parent), max_workers=1)

    assert results['total_tasks'] == 1
    assert results['successful'] == 1
    assert results['failed'] == 0
    assert results['results'][0]['task_name'] == 'task_001'



def test_discover_md_tasks_uses_stage_specific_discovery(tmp_path):
    base = tmp_path / 'raw_md' / 'task_001' / 'prod'
    base.mkdir(parents=True)
    (base / 'md.tpr').touch()
    (base / 'md.xtc').touch()

    discovered = batch_workflow.discover_md_tasks(str(base.parent.parent))

    assert set(discovered) == {'task_001'}
    trajectory, topology = discovered['task_001']
    assert trajectory.endswith('md.xtc')
    assert topology.endswith('md.tpr')
