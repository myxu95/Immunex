import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.interface import BuriedSurfaceAreaResult
from immunex.core.context import PipelineContext
from immunex.pipeline import BSAPipeline
from immunex.pipeline.nodes.bsa_node import BSAAnalysisNode


class TestBSAPipeline(unittest.TestCase):
    def test_pipeline_with_auto_chain_identification_contains_two_nodes(self):
        pipeline = BSAPipeline()
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ["ChainIdentificationNode", "BSAAnalysisNode"])

    def test_pipeline_with_manual_selections_skips_chain_identification(self):
        pipeline = BSAPipeline(
            selection_a="chainID A or chainID B or chainID C",
            selection_b="chainID D or chainID E",
            auto_identify_chains=False,
        )
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ["BSAAnalysisNode"])

    @patch("immunex.pipeline.nodes.bsa_node.BuriedSurfaceAreaAnalyzer")
    def test_bsa_node_writes_timeseries_summary_and_plot(self, analyzer_cls):
        result = BuriedSurfaceAreaResult(
            times=np.array([0.0, 10.0, 20.0], dtype=float),
            sasa_a=np.array([100.0, 101.0, 102.0], dtype=float),
            sasa_b=np.array([80.0, 81.0, 82.0], dtype=float),
            sasa_ab=np.array([150.0, 151.0, 152.0], dtype=float),
            buried_surface_area=np.array([15.0, 15.5, 16.0], dtype=float),
            interface_ratio=np.array([0.1667, 0.1703, 0.1739], dtype=float),
            selection_a="chainID A or chainID B or chainID C",
            selection_b="chainID D or chainID E",
            probe_radius=1.4,
            stride=5,
            time_unit="ps",
        )
        analyzer_instance = analyzer_cls.return_value
        analyzer_instance.calculate.return_value = result

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            context = PipelineContext(
                system_id="demo_bsa",
                topology=str(tmp_path / "md.tpr"),
                trajectory_raw=str(tmp_path / "md.xtc"),
                output_dir=str(tmp_path),
            )
            node = BSAAnalysisNode(
                selection_a="chainID A or chainID B or chainID C",
                selection_b="chainID D or chainID E",
                stride=5,
            )

            updated = node.execute(context)
            bsa_result = updated.get_result("bsa")

            self.assertIsNotNone(bsa_result)
            self.assertTrue(Path(bsa_result["timeseries_file"]).exists())
            self.assertTrue(Path(bsa_result["summary_file"]).exists())
            self.assertTrue(Path(bsa_result["plot_file"]).exists())
            self.assertEqual(bsa_result["selection_a"], "chainID A or chainID B or chainID C")
            self.assertEqual(bsa_result["selection_b"], "chainID D or chainID E")

            analyzer_cls.assert_called_once_with(str(tmp_path / "md.tpr"), str(tmp_path / "md.xtc"))
            analyzer_instance.calculate.assert_called_once()


if __name__ == "__main__":
    unittest.main()
