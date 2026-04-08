import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.context import PipelineContext
from immunex.pipeline import NormalModePipeline
from immunex.pipeline.nodes.normal_mode_node import NormalModeNode
from immunex.analysis.allostery import NormalModeResult


class TestNormalModePipeline(unittest.TestCase):
    def test_pipeline_contains_expected_nodes(self):
        pipeline = NormalModePipeline()
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ["ChainIdentificationNode", "CDRDetectionNode", "NormalModeNode"])

    def test_pipeline_without_auto_nodes(self):
        pipeline = NormalModePipeline(auto_identify_chains=False, auto_detect_cdr=False)
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ["NormalModeNode"])

    @patch("immunex.pipeline.nodes.normal_mode_node.NormalModeAnalyzer")
    def test_normal_mode_node_writes_outputs(self, analyzer_cls):
        residue_frame = pd.DataFrame(
            [
                {
                    "chain_id": "D",
                    "resid": 93,
                    "resname": "ASP",
                    "residue_label": "ASP93",
                    "component": "TCR_alpha",
                    "complex_side": "tcr",
                    "phla_region": None,
                    "mhc_region": None,
                    "mhc_subregion": None,
                    "tcr_chain": "alpha",
                    "tcr_region": "CDR3",
                    "tcr_region_detailed": "CDR3_alpha",
                    "region_group": "CDR3_alpha",
                    "gnm_mobility": 1.2,
                    "gnm_mobility_norm": 0.8,
                    "anm_mobility": 1.0,
                    "anm_mobility_norm": 0.7,
                    "hinge_score": 0.9,
                    "prs_effectiveness": 1.3,
                    "prs_effectiveness_norm": 0.85,
                    "prs_sensitivity": 1.1,
                    "prs_sensitivity_norm": 0.75,
                    "combined_key_residue_score": 0.88,
                },
                {
                    "chain_id": "A",
                    "resid": 154,
                    "resname": "GLU",
                    "residue_label": "GLU154",
                    "component": "HLA_alpha",
                    "complex_side": "phla",
                    "phla_region": "HLA_alpha",
                    "mhc_region": "groove",
                    "mhc_subregion": "alpha2_helix",
                    "tcr_chain": None,
                    "tcr_region": None,
                    "tcr_region_detailed": None,
                    "region_group": "alpha2_helix",
                    "gnm_mobility": 0.6,
                    "gnm_mobility_norm": 0.3,
                    "anm_mobility": 0.4,
                    "anm_mobility_norm": 0.2,
                    "hinge_score": 0.4,
                    "prs_effectiveness": 0.5,
                    "prs_effectiveness_norm": 0.25,
                    "prs_sensitivity": 0.6,
                    "prs_sensitivity_norm": 0.35,
                    "combined_key_residue_score": 0.31,
                },
            ]
        )
        region_summary = pd.DataFrame(
            [
                {
                    "region_group": "CDR3_alpha",
                    "n_residues": 1,
                    "mean_gnm_mobility": 1.2,
                    "mean_anm_mobility": 1.0,
                    "mean_hinge_score": 0.9,
                    "mean_prs_effectiveness": 1.3,
                    "max_combined_score": 0.88,
                },
                {
                    "region_group": "alpha2_helix",
                    "n_residues": 1,
                    "mean_gnm_mobility": 0.6,
                    "mean_anm_mobility": 0.4,
                    "mean_hinge_score": 0.4,
                    "mean_prs_effectiveness": 0.5,
                    "max_combined_score": 0.31,
                },
            ]
        )
        result = NormalModeResult(
            residue_frame=residue_frame,
            region_summary=region_summary,
            summary={
                "n_residues": 2,
                "network_edges": 1,
                "cutoff_angstrom": 10.0,
                "n_low_modes": 10,
                "prs_force_directions": 8,
                "top_key_residues": [{"residue_label": "ASP93"}],
            },
            mode_count=10,
            cutoff_angstrom=10.0,
        )
        analyzer_cls.return_value.calculate.return_value = result

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            context = PipelineContext(
                system_id="demo_nma",
                topology=str(tmp_path / "md.tpr"),
                trajectory_raw=str(tmp_path / "md.xtc"),
                structure_pdb=str(tmp_path / "md_processed_converted.pdb"),
                output_dir=str(tmp_path),
                metadata={
                    "chain_mapping": {"mhc_alpha": "A", "b2m": "B", "peptide": "C", "tcr_alpha": "D", "tcr_beta": "E"},
                    "cdr_detection": {"chains": {}},
                },
            )
            node = NormalModeNode(cutoff_angstrom=10.0, n_low_modes=10, prs_force_directions=8)
            updated = node.execute(context)
            nma_result = updated.results["normal_mode"]

            self.assertTrue(Path(nma_result["residue_csv"]).exists())
            self.assertTrue(Path(nma_result["region_summary_csv"]).exists())
            self.assertTrue(Path(nma_result["summary_json"]).exists())
            self.assertTrue(Path(nma_result["mobility_plot"]).exists())
            self.assertTrue(Path(nma_result["hinge_plot"]).exists())
            self.assertTrue(Path(nma_result["prs_plot"]).exists())

            summary = json.loads(Path(nma_result["summary_json"]).read_text(encoding="utf-8"))
            self.assertEqual(summary["n_residues"], 2)
            self.assertEqual(summary["network_edges"], 1)


if __name__ == "__main__":
    unittest.main()
