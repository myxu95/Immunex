"""
Quick Start Example - New Architecture

This example demonstrates the new four-layer architecture:
1. Core Modules (Layer 1): RMSDCalculator
2. Pipeline Nodes (Layer 2): PreprocessNode, RMSDNode
3. Pipeline Orchestration (Layer 3): StandardTrajectoryPipeline
4. Configuration/UI (Layer 4): Python API

Example usage of the new Immunex pipeline architecture.
"""

import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PipelineContext, discover_tasks
from immunex.pipeline import StandardTrajectoryPipeline, BatchExecutor


def example_single_task():
    """Example: Process a single task."""
    print('=' * 60)
    print('Example 1: Single Task Processing')
    print('=' * 60)

    context = PipelineContext(
        system_id='1ao7_example',
        topology='data/1ao7/md.tpr',
        trajectory_raw='data/1ao7/md.xtc',
        output_dir='./results/example_single'
    )

    pipeline = StandardTrajectoryPipeline()

    print(f'\nProcessing task: {context.system_id}')
    result = pipeline.execute(context)

    if result.has_errors():
        print(f'\nTask failed with {len(result.errors)} errors:')
        for error in result.errors:
            print(f'  - {error}')
    else:
        print('\nTask completed successfully')
        print(f'  Processed trajectory: {result.trajectory_processed}')
        if 'rmsd' in result.results:
            rmsd_data = result.results['rmsd']
            print(f"  Mean RMSD: {rmsd_data['mean_rmsd']:.3f} nm")
            print(f"  Std RMSD: {rmsd_data['std_rmsd']:.3f} nm")


def example_batch_processing():
    """Example: Batch processing with public task discovery."""
    print('\n\n' + '=' * 60)
    print('Example 2: Batch Processing with DiscoveryReport')
    print('=' * 60)

    report = discover_tasks(Path('./test_batch'))
    print(f'\nDiscovered {report.num_valid} valid tasks from {report.total_tasks} total:')
    for task in report.valid_tasks:
        print(f'  - {task.task_id}')

    pipeline = StandardTrajectoryPipeline()
    executor = BatchExecutor(max_workers=2)
    print(f'\nExecuting pipeline with {executor.max_workers} workers...')

    results = executor.execute_pipeline(
        report,
        pipeline,
        show_progress=True,
    )

    summary = executor.summarize_results(results)

    print('\nBatch Summary:')
    print(f"  Total tasks: {summary['total_tasks']}")
    print(f"  Successful: {summary['successful']}")
    print(f"  Failed: {summary['failed']}")
    print(f"  Success rate: {summary['success_rate']:.1f}%")


def example_custom_pipeline():
    """Example: Create a custom pipeline."""
    print('\n\n' + '=' * 60)
    print('Example 3: Custom Pipeline')
    print('=' * 60)

    from immunex.pipeline import Pipeline, PreprocessNode, RMSDNode

    class MyCustomPipeline(Pipeline):
        def __init__(self):
            nodes = [
                PreprocessNode(method='3step', dt=10.0),
                RMSDNode(selection='protein and name CA', reference_frame=0),
                RMSDNode(selection='backbone', reference_frame=0),
            ]
            super().__init__(nodes=nodes)

    context = PipelineContext(
        system_id='custom_example',
        topology='data/example/md.tpr',
        trajectory_raw='data/example/md.xtc'
    )

    pipeline = MyCustomPipeline()
    print(f'\nCustom pipeline with {len(pipeline.nodes)} nodes:')
    for i, node in enumerate(pipeline.nodes, start=1):
        print(f'  {i}. {node.name}')

    result = pipeline.execute(context)

    if not result.has_errors():
        print('\nCustom pipeline completed successfully')
    else:
        print(f'\nCustom pipeline failed: {result.errors}')


def example_context_management():
    """Example: Context data passing."""
    print('\n\n' + '=' * 60)
    print('Example 4: Context Data Passing')
    print('=' * 60)

    context = PipelineContext(
        system_id='context_example',
        topology='data/example/md.tpr',
        trajectory_raw='data/example/md.xtc'
    )

    context.set_selection('protein', 'protein')
    context.set_selection('backbone', 'name CA C N O')
    context.set_selection('tcr', 'chainID D E')

    context.metadata['description'] = 'TCR-pMHC complex'
    context.metadata['temperature'] = 300.0

    output_file = './results/example_context.json'
    context.save(output_file)
    print(f'\nContext saved to: {output_file}')


if __name__ == '__main__':
    example_single_task()
    example_batch_processing()
    example_custom_pipeline()
    example_context_management()
