import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.context import PipelineContext
from immunex.pipeline import AnnotatedRMSFPipeline
from immunex.pipeline.nodes.residue_rmsf_node import ResidueRMSFNode
from immunex.analysis.trajectory import ResidueRMSFResult


class TestAnnotatedRMSFPipeline(unittest.TestCase):
    def test_pipeline_contains_chain_cdr_and_rmsf_nodes(self):
        pipeline = AnnotatedRMSFPipeline()
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ["ChainIdentificationNode", "CDRDetectionNode", "ResidueRMSFNode"])

    def test_pipeline_with_manual_metadata_skips_detection_nodes(self):
        pipeline = AnnotatedRMSFPipeline(auto_identify_chains=False, auto_detect_cdr=False)
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ["ResidueRMSFNode"])

    @patch("immunex.pipeline.nodes.residue_rmsf_node.ResidueRMSFAnalyzer")
    def test_residue_rmsf_node_writes_expected_outputs(self, analyzer_cls):
        residue_frame = pd.DataFrame(
            [
                {
                    "chain_id": "D",
                    "resid": 93,
                    "resname": "ASP",
                    "component": "TCR_alpha",
                    "complex_side": "tcr",
                    "phla_region": None,
                    "mhc_region": None,
                    "mhc_subregion": None,
                    "tcr_chain": "alpha",
                    "tcr_region": "CDR3",
                    "tcr_region_detailed": "CDR3_alpha",
                    "region_group": "CDR3_alpha",
                    "rmsf_angstrom": 1.25,
                },
                {
                    "chain_id": "C",
                    "resid": 5,
                    "resname": "GLY",
                    "component": "peptide",
                    "complex_side": "phla",
                    "phla_region": "peptide",
                    "mhc_region": "peptide",
                    "mhc_subregion": "peptide",
                    "tcr_chain": None,
                    "tcr_region": None,
                    "tcr_region_detailed": None,
                    "region_group": "peptide",
                    "rmsf_angstrom": 0.82,
                },
            ]
        )
        region_summary = pd.DataFrame(
            [
                {"region_group": "CDR3_alpha", "n_residues": 1, "mean_rmsf_angstrom": 1.25, "median_rmsf_angstrom": 1.25, "max_rmsf_angstrom": 1.25, "min_rmsf_angstrom": 1.25},
                {"region_group": "peptide", "n_residues": 1, "mean_rmsf_angstrom": 0.82, "median_rmsf_angstrom": 0.82, "max_rmsf_angstrom": 0.82, "min_rmsf_angstrom": 0.82},
            ]
        )
        result = ResidueRMSFResult(
            residue_frame=residue_frame,
            region_summary=region_summary,
            summary={
                "selection": "name CA and (chainID C or chainID D)",
                "stride": 10,
                "n_residues": 2,
                "n_frames": 50,
                "time_unit": "ps",
                "time_start": 0.0,
                "time_end": 5000.0,
                "mean_rmsf_angstrom": 1.035,
                "max_rmsf_angstrom": 1.25,
                "tcr_mean_rmsf_angstrom": 1.25,
                "phla_mean_rmsf_angstrom": 0.82,
            },
            stride=10,
            time_unit="ps",
        )
        analyzer_instance = analyzer_cls.return_value
        analyzer_instance.calculate.return_value = result

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            context = PipelineContext(
                system_id="demo_rmsf",
                topology=str(tmp_path / "md.tpr"),
                trajectory_raw=str(tmp_path / "md.xtc"),
                structure_pdb=str(tmp_path / "md_processed_converted.pdb"),
                output_dir=str(tmp_path),
                metadata={
                    "chain_mapping": {"mhc_alpha": "A", "b2m": "B", "peptide": "C", "tcr_alpha": "D", "tcr_beta": "E"},
                    "cdr_detection": {"chains": {}},
                },
            )
            node = ResidueRMSFNode(stride=10)
            updated = node.execute(context)
            rmsf_result = updated.get_result("residue_rmsf")

            self.assertTrue(Path(rmsf_result["residue_csv"]).exists())
            self.assertTrue(Path(rmsf_result["region_summary_csv"]).exists())
            self.assertTrue(Path(rmsf_result["summary_json"]).exists())
            self.assertTrue(Path(rmsf_result["tcr_plot"]).exists())
            self.assertTrue(Path(rmsf_result["phla_plot"]).exists())
            self.assertTrue(Path(rmsf_result["region_plot"]).exists())

            summary = json.loads(Path(rmsf_result["summary_json"]).read_text(encoding="utf-8"))
            self.assertEqual(summary["n_residues"], 2)
            self.assertEqual(summary["stride"], 10)


if __name__ == "__main__":
    unittest.main()
