#!/usr/bin/env python3
"""Hydrogen-bond interaction pipeline smoke tests."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.context import PipelineContext
from immunex.pipeline.nodes.hydrogen_bond_pair_node import HydrogenBondPairNode
from immunex.pipeline.nodes.hydrogen_bond_annotation_node import HydrogenBondAnnotationNode


class TestHydrogenBondPairNode(unittest.TestCase):
    def test_hbond_pair_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            topology = temp_path / "test.pdb"
            trajectory = temp_path / "traj.xtc"
            topology.touch()
            trajectory.touch()

            context = PipelineContext(
                system_id="demo",
                topology=str(topology),
                trajectory_raw=str(trajectory),
                structure_pdb=str(topology),
                trajectory_processed=str(trajectory),
                output_dir=str(temp_path / "results"),
            )
            context.metadata["chain_mapping"] = {
                "mhc_alpha": "A",
                "b2m": "B",
                "peptide": "C",
                "tcr_alpha": "D",
                "tcr_beta": "E",
            }

            mock_pairs = pd.DataFrame(
                [
                    {
                        "chain_id_1": "A",
                        "resid_1": 65,
                        "resname_1": "ASP",
                        "residue_label_1": "ASP65",
                        "chain_id_2": "E",
                        "resid_2": 55,
                        "resname_2": "SER",
                        "residue_label_2": "SER55",
                        "contact_frames": 10,
                        "total_frames": 20,
                        "contact_frequency": 0.5,
                        "min_distance_observed": 2.9,
                        "mean_angle_observed": 162.0,
                        "interaction_family": "hbond",
                        "distance_cutoff_angstrom": 3.5,
                        "angle_cutoff_deg": 150.0,
                    }
                ]
            )

            with patch(
                "immunex.pipeline.nodes.hydrogen_bond_pair_node.HydrogenBondPairAnalyzer"
            ) as analyzer_cls:
                analyzer = analyzer_cls.return_value
                analyzer.calculate_interface_pairs.return_value = mock_pairs

                result = HydrogenBondPairNode().execute(context)

            self.assertIn("hbond_pairs", result.results)
            raw_file = Path(result.results["hbond_pairs"]["raw_pairs_file"])
            self.assertTrue(raw_file.exists())
            written = pd.read_csv(raw_file)
            self.assertEqual(len(written), 1)

    def test_hbond_pair_node_prefers_topology_over_structure_pdb(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            topology = temp_path / "test.tpr"
            structure = temp_path / "test.pdb"
            trajectory = temp_path / "traj.xtc"
            topology.touch()
            structure.touch()
            trajectory.touch()

            context = PipelineContext(
                system_id="demo",
                topology=str(topology),
                trajectory_raw=str(trajectory),
                structure_pdb=str(structure),
                trajectory_processed=str(trajectory),
                output_dir=str(temp_path / "results"),
            )
            context.metadata["chain_mapping"] = {
                "mhc_alpha": "A",
                "b2m": "B",
                "peptide": "C",
                "tcr_alpha": "D",
                "tcr_beta": "E",
            }

            mock_pairs = pd.DataFrame(columns=[
                "chain_id_1", "resid_1", "resname_1", "residue_label_1",
                "chain_id_2", "resid_2", "resname_2", "residue_label_2",
                "contact_frames", "total_frames", "contact_frequency",
                "min_distance_observed", "mean_angle_observed",
                "interaction_family", "distance_cutoff_angstrom", "angle_cutoff_deg",
            ])

            with patch(
                "immunex.pipeline.nodes.hydrogen_bond_pair_node.HydrogenBondPairAnalyzer"
            ) as analyzer_cls:
                analyzer = analyzer_cls.return_value
                analyzer.calculate_interface_pairs.return_value = mock_pairs

                HydrogenBondPairNode().execute(context)

            analyzer_cls.assert_called_once_with(
                str(topology),
                str(trajectory),
                reference_structure=str(structure),
            )


class TestHydrogenBondAnnotationNode(unittest.TestCase):
    def test_hbond_annotation_node_writes_region_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            raw_file = temp_path / "raw_hbonds.csv"
            pd.DataFrame(
                [
                    {
                        "chain_id_1": "A",
                        "resid_1": 65,
                        "resname_1": "ASP",
                        "residue_label_1": "ASP65",
                        "chain_id_2": "E",
                        "resid_2": 55,
                        "resname_2": "SER",
                        "residue_label_2": "SER55",
                        "contact_frames": 10,
                        "total_frames": 20,
                        "contact_frequency": 0.5,
                        "min_distance_observed": 2.9,
                        "mean_angle_observed": 162.0,
                        "interaction_family": "hbond",
                        "distance_cutoff_angstrom": 3.5,
                        "angle_cutoff_deg": 150.0,
                    }
                ]
            ).to_csv(raw_file, index=False)

            context = PipelineContext(
                system_id="demo",
                topology=str(temp_path / "test.pdb"),
                trajectory_raw=str(temp_path / "traj.xtc"),
                output_dir=str(temp_path / "results"),
            )
            context.results["hbond_pairs"] = {"raw_pairs_file": str(raw_file)}
            context.metadata["chain_mapping"] = {
                "mhc_alpha": "A",
                "b2m": "B",
                "peptide": "C",
                "tcr_alpha": "D",
                "tcr_beta": "E",
            }
            context.metadata["cdr_detection"] = {
                "chains": {
                    "TCR_beta": {"chain_id": "E", "cdrs": {2: {"residue_range": [50, 60]}}},
                }
            }

            result = HydrogenBondAnnotationNode().execute(context)
            outputs = result.results["hbond_annotation"]
            self.assertTrue(Path(outputs["annotated_pairs_file"]).exists())
            self.assertTrue(Path(outputs["hbond_report_file"]).exists())
            self.assertTrue(Path(outputs["groove_tcr_hbonds_file"]).exists())
            self.assertTrue(Path(outputs["cdr2"]["groove_tcr_contacts_file"]).exists())


if __name__ == "__main__":
    unittest.main()
