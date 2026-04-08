import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.conformation import InterfaceClusteringResult
from immunex.analysis.conformation.interface_clustering import _derive_contact_stride
from immunex.core.context import PipelineContext
from immunex.pipeline import InterfaceClusteringPipeline
from immunex.pipeline.nodes.interface_clustering_node import InterfaceClusteringNode


class TestInterfaceClusteringPipeline(unittest.TestCase):
    def test_contact_stride_rule(self):
        self.assertEqual(_derive_contact_stride(1), 1)
        self.assertEqual(_derive_contact_stride(2), 1)
        self.assertEqual(_derive_contact_stride(3), 2)
        self.assertEqual(_derive_contact_stride(5), 3)
        self.assertEqual(_derive_contact_stride(10), 5)

    def test_pipeline_contains_expected_nodes(self):
        pipeline = InterfaceClusteringPipeline()
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(
            node_names,
            ["ChainIdentificationNode", "CDRDetectionNode", "InterfaceClusteringNode"],
        )

    def test_pipeline_without_auto_nodes(self):
        pipeline = InterfaceClusteringPipeline(auto_identify_chains=False, auto_detect_cdr=False)
        node_names = [node.__class__.__name__ for node in pipeline.nodes]
        self.assertEqual(node_names, ["InterfaceClusteringNode"])

    @patch("immunex.pipeline.nodes.interface_clustering_node.InterfaceClusteringAnalyzer")
    def test_interface_clustering_node_writes_outputs(self, analyzer_cls):
        frame_assignments = pd.DataFrame(
            [
                {"frame": 0, "time_ps": 0.0, "cluster_id": 0},
                {"frame": 10, "time_ps": 20.0, "cluster_id": 1},
            ]
        )
        cluster_summary = pd.DataFrame(
            [
                {"cluster_id": 0, "population_percent": 50.0, "frames": 1, "representative_frame": 0, "medoid_frame": 0, "first_appearance_time_ps": 0.0, "mean_dwell_time_ps": 20.0, "transitions_in": 0, "transitions_out": 1, "main_structural_descriptor": "E:95__C:5[peptide]", "mean_distance_to_medoid": 0.0},
                {"cluster_id": 1, "population_percent": 50.0, "frames": 1, "representative_frame": 10, "medoid_frame": 10, "first_appearance_time_ps": 20.0, "mean_dwell_time_ps": 20.0, "transitions_in": 1, "transitions_out": 0, "main_structural_descriptor": "E:100__A:154[alpha2_helix]", "mean_distance_to_medoid": 0.0},
            ]
        )
        cluster_difference_matrix = pd.DataFrame(
            [[0.0, 0.6], [0.6, 0.0]],
            index=["cluster_0", "cluster_1"],
            columns=["cluster_0", "cluster_1"],
        )
        cluster_feature_summary = pd.DataFrame(
            [
                {"cluster_id": 0, "pair_label": "E:95__C:5[peptide]", "contact_occupancy": 1.0, "mean_pair_distance_angstrom": 3.2, "rank": 1},
                {"cluster_id": 1, "pair_label": "E:100__A:154[alpha2_helix]", "contact_occupancy": 0.7, "mean_pair_distance_angstrom": 4.1, "rank": 1},
            ]
        )
        cluster_signatures = pd.DataFrame(
            [
                {"cluster_id": 0, "pair_label": "E:95__C:5[peptide]", "contact_occupancy": 1.0, "mean_pair_distance_angstrom": 3.2, "occupancy_delta_vs_others": 0.6, "signature_score": 0.6},
                {"cluster_id": 1, "pair_label": "E:100__A:154[alpha2_helix]", "contact_occupancy": 0.7, "mean_pair_distance_angstrom": 4.1, "occupancy_delta_vs_others": 0.3, "signature_score": 0.21},
            ]
        )
        cluster_feature_digest = pd.DataFrame(
            [
                {"cluster_id": 0, "digest_type": "key_feature", "pair_label": "E:95__C:5[peptide]", "contact_occupancy": 1.0, "mean_pair_distance_angstrom": 3.2, "score": 1.0, "rank": 1},
                {"cluster_id": 1, "digest_type": "state_signature", "pair_label": "E:100__A:154[alpha2_helix]", "contact_occupancy": 0.7, "mean_pair_distance_angstrom": 4.1, "score": 0.21, "rank": 1},
            ]
        )
        result = InterfaceClusteringResult(
            frame_assignments=frame_assignments,
            cluster_summary=cluster_summary,
            cluster_difference_matrix=cluster_difference_matrix,
            cluster_feature_summary=cluster_feature_summary,
            cluster_signatures=cluster_signatures,
            cluster_feature_digest=cluster_feature_digest,
            summary={"n_frames": 2, "n_clusters": 2, "largest_cluster": {"cluster_id": 0}},
        )
        analyzer_cls.return_value.calculate.return_value = result
        analyzer_cls.return_value.export_cluster_structures.return_value = [
            {"cluster_id": 0, "representative_pdb": "rep0.pdb", "medoid_pdb": "med0.pdb", "average_pdb": "avg0.pdb"},
            {"cluster_id": 1, "representative_pdb": "rep1.pdb", "medoid_pdb": "med1.pdb", "average_pdb": "avg1.pdb"},
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            context = PipelineContext(
                system_id="demo_interface_cluster",
                topology=str(tmp_path / "md.tpr"),
                trajectory_raw=str(tmp_path / "md.xtc"),
                trajectory_processed=str(tmp_path / "md.xtc"),
                structure_pdb=str(tmp_path / "md_processed_converted.pdb"),
                output_dir=str(tmp_path),
                metadata={
                    "chain_mapping": {"mhc_alpha": "A", "b2m": "B", "peptide": "C", "tcr_alpha": "D", "tcr_beta": "E"},
                    "cdr_detection": {
                        "chains": {
                            "TCR_alpha": {"chain_id": "D", "cdrs": {3: {"residue_range": [90, 100]}}},
                            "TCR_beta": {"chain_id": "E", "cdrs": {3: {"residue_range": [95, 105]}}},
                        }
                    },
                },
            )
            node = InterfaceClusteringNode()
            updated = node.execute(context)
            cluster_result = updated.results["interface_clustering"]

            self.assertTrue(Path(cluster_result["frame_assignments_csv"]).exists())
            self.assertTrue(Path(cluster_result["summary_table_csv"]).exists())
            self.assertTrue(Path(cluster_result["cluster_difference_matrix_csv"]).exists())
            self.assertTrue(Path(cluster_result["cluster_feature_summary_csv"]).exists())
            self.assertTrue(Path(cluster_result["cluster_signatures_csv"]).exists())
            self.assertTrue(Path(cluster_result["cluster_feature_digest_csv"]).exists())
            self.assertTrue(Path(cluster_result["summary_json"]).exists())
            self.assertTrue(Path(cluster_result["cluster_id_vs_time_plot"]).exists())
            self.assertTrue(Path(cluster_result["state_population_plot"]).exists())
            self.assertEqual(len(cluster_result["structures"]), 2)

            summary = json.loads(Path(cluster_result["summary_json"]).read_text(encoding="utf-8"))
            self.assertEqual(summary["n_frames"], 2)
            self.assertEqual(summary["n_clusters"], 2)


if __name__ == "__main__":
    unittest.main()
