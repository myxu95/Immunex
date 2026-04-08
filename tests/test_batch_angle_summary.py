import tempfile
import unittest
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.cli.commands.batch import _load_angle_timeseries, _write_angle_summary
from immunex.core.context import PipelineContext
from immunex.pipeline.nodes.docking_angle_node import DockingAngleNode


class TestBatchAngleSummary(unittest.TestCase):
    def test_load_angle_timeseries_computes_tail90_variation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_file = Path(tmpdir) / 'docking_angles.csv'
            pd.DataFrame(
                {
                    'Time(ps)': [0, 1, 2, 3, 4],
                    'Crossing(deg)': [10.0, 12.0, 11.0, 14.0, 13.0],
                    'Incident(deg)': [30.0, 31.0, 29.0, 33.0, 32.0],
                }
            ).to_csv(csv_file, index=False)

            metrics = _load_angle_timeseries(str(csv_file))

            self.assertEqual(metrics['n_frames'], 5)
            self.assertEqual(metrics['trimmed_start_frame'], 0)
            self.assertAlmostEqual(metrics['crossing_full_variation_deg'], 4.0)
            self.assertAlmostEqual(metrics['incident_tail90_variation_deg'], 4.0)

    def test_write_angle_summary_outputs_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            task_dir = tmp_path / 'task1'
            angles_dir = task_dir / 'analysis' / 'angles'
            angles_dir.mkdir(parents=True)
            csv_file = angles_dir / 'docking_angles.csv'
            pd.DataFrame(
                {
                    'Time(ps)': [0, 1, 2],
                    'Crossing(deg)': [40.0, 45.0, 43.0],
                    'Incident(deg)': [12.0, 16.0, 13.0],
                }
            ).to_csv(csv_file, index=False)

            results = {
                'results': [
                    {
                        'task_name': 'task1',
                        'status': 'success',
                        'output_dir': str(task_dir),
                        'errors': [],
                        'docking_angles': {
                            'statistics': {
                                'crossing_mean': 42.67,
                                'crossing_std': 2.05,
                                'incident_mean': 13.67,
                                'incident_std': 1.70,
                            },
                            'output_files': [str(csv_file)],
                        },
                    }
                ]
            }

            summary_files = _write_angle_summary(results, str(tmp_path / 'summary'))

            self.assertTrue(Path(summary_files['csv']).exists())
            self.assertTrue(Path(summary_files['json']).exists())
            self.assertTrue(Path(summary_files['markdown']).exists())
            self.assertTrue(Path(summary_files['crossing_plot']).exists())
            self.assertTrue(Path(summary_files['incident_plot']).exists())

    def test_docking_angle_node_uses_matching_topology_for_trajectory(self):
        context = PipelineContext(
            system_id='demo',
            topology='/tmp/demo.tpr',
            trajectory_raw='/tmp/demo.xtc',
            structure_pdb='/tmp/demo.pdb',
            trajectory_processed='/tmp/demo_processed.xtc',
        )
        node = DockingAngleNode(stride=5)

        input_params = node._build_input_params(context, context.trajectory_processed)

        self.assertEqual(input_params.topology, '/tmp/demo.tpr')
        self.assertEqual(input_params.trajectory, '/tmp/demo_processed.xtc')
        self.assertEqual(input_params.stride, 5)

    def test_docking_angle_node_uses_chain_mapping_as_manual_selections(self):
        context = PipelineContext(
            system_id='demo',
            topology='/tmp/demo.tpr',
            trajectory_raw='/tmp/demo.xtc',
            structure_pdb='/tmp/demo.pdb',
            trajectory_processed='/tmp/demo_processed.xtc',
            metadata={
                'chain_mapping': {
                    'mhc_alpha': 'A',
                    'tcr_alpha': 'D',
                    'tcr_beta': 'E',
                }
            }
        )
        node = DockingAngleNode(stride=5)

        input_params = node._build_input_params(context, context.trajectory_processed)

        self.assertFalse(input_params.auto_identify_chains)
        self.assertEqual(input_params.mhc_selection, 'chainID A')
        self.assertEqual(input_params.tcr_alpha_selection, 'chainID D')
        self.assertEqual(input_params.tcr_beta_selection, 'chainID E')


if __name__ == '__main__':
    unittest.main()
