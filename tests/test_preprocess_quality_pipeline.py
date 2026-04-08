import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.context import PipelineContext
from immunex.pipeline import PreprocessQualityPipeline
from immunex.pipeline.nodes.rmsd_quality_node import RMSDQualityNode


class TestPreprocessQualityPipeline(unittest.TestCase):
    def test_pipeline_contains_preprocess_rmsd_and_quality_nodes(self):
        pipeline = PreprocessQualityPipeline(method='2step')
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ['PreprocessNode', 'RMSDNode', 'RMSDPlotNode', 'RMSDQualityNode'])

    def test_rmsd_quality_node_generates_reports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            rmsd_dir = tmp_path / 'analysis' / 'rmsd'
            rmsd_dir.mkdir(parents=True)
            rmsd_file = rmsd_dir / 'rmsd.xvg'
            rmsd_file.write_text(
                "\n".join(
                    [
                        '# RMSD calculation',
                        '@    title "RMSD"',
                        '0.0 0.10',
                        '1.0 0.11',
                        '2.0 0.09',
                        '3.0 0.10',
                        '4.0 0.11',
                    ]
                ),
                encoding='utf-8',
            )

            context = PipelineContext(
                system_id='demo',
                topology=str(tmp_path / 'md.tpr'),
                trajectory_raw=str(tmp_path / 'md.xtc'),
                output_dir=str(tmp_path),
            )
            context.add_result(
                'rmsd',
                {
                    'success': True,
                    'output_file': str(rmsd_file),
                    'mean_rmsd': 0.102,
                    'std_rmsd': 0.008,
                }
            )

            result = RMSDQualityNode().execute(context)
            quality = result.get_result('preprocess_quality')

            self.assertIsNotNone(quality)
            self.assertIn('tail90_variation_nm', quality['metrics'])
            self.assertTrue(Path(quality['reports']['markdown']).exists())
            self.assertTrue(Path(quality['reports']['json']).exists())


if __name__ == '__main__':
    unittest.main()
