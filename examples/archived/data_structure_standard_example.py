#!/usr/bin/env python3
"""
Data Structure Standard Example

Demonstrates the unified data structure standard for Immunex:
1. Standardized input directory structures (nested and flat)
2. Standardized output directory organization
3. Enhanced PipelineContext path management
4. Manifest-based task discovery with multi-level and list-based support
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import (
    PipelineContext,
    discover_tasks,
    discover_tasks_from_list,
    discovery_report_to_contexts,
)
from immunex.pipeline import StandardTrajectoryPipeline


def example_context_path_management():
    print('=' * 60)
    print('Example 1: PipelineContext Path Management')
    print('=' * 60)

    context = PipelineContext(
        system_id='1ao7_standard',
        topology='data/1ao7/standard/md.tpr',
        trajectory_raw='data/1ao7/standard/md.xtc',
        output_dir='./results/1ao7_standard'
    )

    print('\nStandardized output paths:')
    print(f"  Preprocessing: {context.get_preprocessing_path('processed.xtc')}")
    print(f"  RMSD analysis: {context.get_analysis_path('rmsd', 'rmsd_protein.xvg')}")
    print(f"  Angle analysis: {context.get_analysis_path('angles', 'docking_angles.csv')}")
    print(f"  Interface analysis: {context.get_analysis_path('interface', 'contact_heatmap.png')}")
    print(f"  Plot: {context.get_plot_path('rmsd_overview.png')}")
    print(f"  Quality: {context.get_quality_path('energy_quality.json')}")
    print(f"  Index: {context.get_index_path('protein.ndx')}")
    print(f"  Temp: {context.get_temp_path('intermediate.xtc')}")


def example_nested_structure_discovery():
    print('\n\n' + '=' * 60)
    print('Example 2: Nested Structure Discovery')
    print('=' * 60)

    print('\nExpected nested structure:')
    print('  input/')
    print('    {pdb_id}/')
    print('      {variant}/')
    print('        md.tpr, md.xtc, md.gro')

    print('\nManifest-based discovery:')
    print('  report = discover_tasks(base_dir, task_depth=2)')
    print('\nWould discover tasks like:')
    print('  - 1ao7_standard')
    print('  - 1ao7_rest2')
    print('  - 1bd2_standard')


def example_flat_structure_discovery():
    print('\n\n' + '=' * 60)
    print('Example 3: Flat Structure Discovery')
    print('=' * 60)

    print('\nExpected flat structure:')
    print('  input/')
    print('    {task_name}/')
    print('      md.tpr, md.xtc, md.gro')
    print('    OR')
    print('    {task_name}/')
    print('      prod/')
    print('        md.tpr, md.xtc, md.gro')

    print('\nManifest-based discovery:')
    print('  report = discover_tasks(base_dir, task_depth=1)')


def example_custom_output_directory():
    print('\n\n' + '=' * 60)
    print('Example 4: Custom Output Directory')
    print('=' * 60)

    report = discover_tasks_from_list([
        {
            'task_id': '1ao7_standard',
            'structure': 'data/1ao7/standard/md.gro',
            'topology': 'data/1ao7/standard/md.tpr',
            'trajectory_raw': 'data/1ao7/standard/md.xtc'
        }
    ], source_root='.')
    contexts = discovery_report_to_contexts(report, output_base_dir='/custom/output/batch1')

    print('\nCustom output directory:')
    print('  Base: /custom/output/batch1/')
    print(f'  Task output: {contexts[0].output_dir}')


def example_output_directory_structure():
    print('\n\n' + '=' * 60)
    print('Example 5: Standard Output Directory Structure')
    print('=' * 60)

    context = PipelineContext(
        system_id='1ao7_standard',
        topology='data/1ao7/standard/md.tpr',
        trajectory_raw='data/1ao7/standard/md.xtc'
    )

    print("\nStandard output structure for task '1ao7_standard':")
    print(context.get_preprocessing_path('processed.xtc'))
    print(context.get_analysis_path('rmsd', 'rmsd_protein.xvg'))
    print(context.get_quality_path('quality_report.txt'))


def example_report_to_pipeline_bridge():
    print('\n\n' + '=' * 60)
    print('Example 6: DiscoveryReport To Pipeline Bridge')
    print('=' * 60)

    report = discover_tasks(Path('./test_batch'))
    pipeline = StandardTrajectoryPipeline()
    contexts = discovery_report_to_contexts(report)

    print(f'\nDiscovered {report.num_valid} valid tasks')
    print(f'Adapted {len(contexts)} contexts for pipeline execution')
    print(f'Pipeline type: {type(pipeline).__name__}')


if __name__ == '__main__':
    example_context_path_management()
    example_nested_structure_discovery()
    example_flat_structure_discovery()
    example_custom_output_directory()
    example_output_directory_structure()
    example_report_to_pipeline_bridge()
