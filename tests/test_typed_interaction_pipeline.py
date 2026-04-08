#!/usr/bin/env python3
"""Salt bridge / hydrophobic interaction pipeline smoke tests."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.context import PipelineContext
from immunex.pipeline.nodes.salt_bridge_pair_node import SaltBridgePairNode
from immunex.pipeline.nodes.salt_bridge_annotation_node import SaltBridgeAnnotationNode
from immunex.pipeline.nodes.hydrophobic_contact_pair_node import HydrophobicContactPairNode
from immunex.pipeline.nodes.hydrophobic_annotation_node import HydrophobicAnnotationNode
from immunex.pipeline.nodes.pi_stacking_pair_node import PiStackingPairNode
from immunex.pipeline.nodes.pi_stacking_annotation_node import PiStackingAnnotationNode
from immunex.pipeline.nodes.cation_pi_pair_node import CationPiPairNode
from immunex.pipeline.nodes.cation_pi_annotation_node import CationPiAnnotationNode


def _chain_mapping():
    return {
        "mhc_alpha": "A",
        "b2m": "B",
        "peptide": "C",
        "tcr_alpha": "D",
        "tcr_beta": "E",
    }


def _pair_frame(interaction_family: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "chain_id_1": "A",
                "resid_1": 65,
                "resname_1": "ASP",
                "residue_label_1": "ASP65",
                "chain_id_2": "E",
                "resid_2": 55,
                "resname_2": "LYS",
                "residue_label_2": "LYS55",
                "contact_frames": 10,
                "total_frames": 20,
                "contact_frequency": 0.5,
                "min_distance_observed": 3.2,
                "interaction_family": interaction_family,
                "distance_cutoff_angstrom": 4.0,
            }
        ]
    )


def _pi_pair_frame(interaction_family: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "chain_id_1": "A",
                "resid_1": 65,
                "resname_1": "TYR",
                "residue_label_1": "TYR65",
                "chain_id_2": "E",
                "resid_2": 95,
                "resname_2": "TRP",
                "residue_label_2": "TRP95",
                "contact_frames": 8,
                "total_frames": 20,
                "contact_frequency": 0.4,
                "min_distance_observed": 4.8,
                "min_geometry_angle_observed": 18.0,
                "interaction_family": interaction_family,
                "distance_cutoff_angstrom": 6.5,
            }
        ]
    )


class TestSaltBridgeNodes(unittest.TestCase):
    def test_pair_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            topology = temp / "test.tpr"
            structure = temp / "test.pdb"
            trajectory = temp / "traj.xtc"
            topology.touch()
            structure.touch()
            trajectory.touch()
            context = PipelineContext(
                system_id="demo",
                topology=str(topology),
                trajectory_raw=str(trajectory),
                structure_pdb=str(structure),
                trajectory_processed=str(trajectory),
                output_dir=str(temp / "results"),
            )
            context.metadata["chain_mapping"] = _chain_mapping()

            with patch("immunex.pipeline.nodes.salt_bridge_pair_node.SaltBridgePairAnalyzer") as analyzer_cls:
                analyzer = analyzer_cls.return_value
                analyzer.calculate_interface_pairs.return_value = _pair_frame("salt_bridge")
                result = SaltBridgePairNode().execute(context)

            self.assertIn("salt_bridge_pairs", result.results)
            self.assertTrue(Path(result.results["salt_bridge_pairs"]["raw_pairs_file"]).exists())

    def test_annotation_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            raw_file = temp / "raw_salt.csv"
            _pair_frame("salt_bridge").to_csv(raw_file, index=False)
            context = PipelineContext(
                system_id="demo",
                topology=str(temp / "test.pdb"),
                trajectory_raw=str(temp / "traj.xtc"),
                output_dir=str(temp / "results"),
            )
            context.results["salt_bridge_pairs"] = {"raw_pairs_file": str(raw_file)}
            context.metadata["chain_mapping"] = _chain_mapping()
            context.metadata["cdr_detection"] = {
                "chains": {"TCR_beta": {"chain_id": "E", "cdrs": {2: {"residue_range": [50, 60]}}}}
            }
            result = SaltBridgeAnnotationNode().execute(context)
            outputs = result.results["salt_bridge_annotation"]
            self.assertTrue(Path(outputs["salt_bridge_report_file"]).exists())
            self.assertTrue(Path(outputs["groove_tcr_salt_bridges_file"]).exists())


class TestHydrophobicNodes(unittest.TestCase):
    def test_pair_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            topology = temp / "test.tpr"
            structure = temp / "test.pdb"
            trajectory = temp / "traj.xtc"
            topology.touch()
            structure.touch()
            trajectory.touch()
            context = PipelineContext(
                system_id="demo",
                topology=str(topology),
                trajectory_raw=str(trajectory),
                structure_pdb=str(structure),
                trajectory_processed=str(trajectory),
                output_dir=str(temp / "results"),
            )
            context.metadata["chain_mapping"] = _chain_mapping()

            with patch("immunex.pipeline.nodes.hydrophobic_contact_pair_node.HydrophobicContactPairAnalyzer") as analyzer_cls:
                analyzer = analyzer_cls.return_value
                analyzer.calculate_interface_pairs.return_value = _pair_frame("hydrophobic_contact")
                result = HydrophobicContactPairNode().execute(context)

            self.assertIn("hydrophobic_pairs", result.results)
            self.assertTrue(Path(result.results["hydrophobic_pairs"]["raw_pairs_file"]).exists())

    def test_annotation_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            raw_file = temp / "raw_hydrophobic.csv"
            _pair_frame("hydrophobic_contact").to_csv(raw_file, index=False)
            context = PipelineContext(
                system_id="demo",
                topology=str(temp / "test.pdb"),
                trajectory_raw=str(temp / "traj.xtc"),
                output_dir=str(temp / "results"),
            )
            context.results["hydrophobic_pairs"] = {"raw_pairs_file": str(raw_file)}
            context.metadata["chain_mapping"] = _chain_mapping()
            context.metadata["cdr_detection"] = {
                "chains": {"TCR_beta": {"chain_id": "E", "cdrs": {2: {"residue_range": [50, 60]}}}}
            }
            result = HydrophobicAnnotationNode().execute(context)
            outputs = result.results["hydrophobic_annotation"]
            self.assertTrue(Path(outputs["hydrophobic_report_file"]).exists())
            self.assertTrue(Path(outputs["groove_tcr_hydrophobic_contacts_file"]).exists())


class TestPiNodes(unittest.TestCase):
    def test_pair_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            topology = temp / "test.tpr"
            structure = temp / "test.pdb"
            trajectory = temp / "traj.xtc"
            topology.touch()
            structure.touch()
            trajectory.touch()
            context = PipelineContext(
                system_id="demo",
                topology=str(topology),
                trajectory_raw=str(trajectory),
                structure_pdb=str(structure),
                trajectory_processed=str(trajectory),
                output_dir=str(temp / "results"),
            )
            context.metadata["chain_mapping"] = _chain_mapping()

            with patch("immunex.pipeline.nodes.pi_stacking_pair_node.PiStackingPairAnalyzer") as analyzer_cls:
                analyzer = analyzer_cls.return_value
                analyzer.calculate_interface_pairs.return_value = _pi_pair_frame("pi_pi")
                result = PiStackingPairNode().execute(context)

            self.assertIn("pi_pi_pairs", result.results)
            self.assertTrue(Path(result.results["pi_pi_pairs"]["raw_pairs_file"]).exists())

    def test_annotation_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            raw_file = temp / "raw_pi.csv"
            _pi_pair_frame("pi_pi").to_csv(raw_file, index=False)
            context = PipelineContext(
                system_id="demo",
                topology=str(temp / "test.pdb"),
                trajectory_raw=str(temp / "traj.xtc"),
                output_dir=str(temp / "results"),
            )
            context.results["pi_pi_pairs"] = {"raw_pairs_file": str(raw_file)}
            context.metadata["chain_mapping"] = _chain_mapping()
            context.metadata["cdr_detection"] = {
                "chains": {"TCR_beta": {"chain_id": "E", "cdrs": {3: {"residue_range": [90, 100]}}}}
            }
            result = PiStackingAnnotationNode().execute(context)
            outputs = result.results["pi_pi_annotation"]
            self.assertTrue(Path(outputs["pi_pi_report_file"]).exists())
            self.assertTrue(Path(outputs["groove_tcr_pi_pi_file"]).exists())


class TestCationPiNodes(unittest.TestCase):
    def test_pair_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            topology = temp / "test.tpr"
            structure = temp / "test.pdb"
            trajectory = temp / "traj.xtc"
            topology.touch()
            structure.touch()
            trajectory.touch()
            context = PipelineContext(
                system_id="demo",
                topology=str(topology),
                trajectory_raw=str(trajectory),
                structure_pdb=str(structure),
                trajectory_processed=str(trajectory),
                output_dir=str(temp / "results"),
            )
            context.metadata["chain_mapping"] = _chain_mapping()

            with patch("immunex.pipeline.nodes.cation_pi_pair_node.CationPiPairAnalyzer") as analyzer_cls:
                analyzer = analyzer_cls.return_value
                analyzer.calculate_interface_pairs.return_value = _pi_pair_frame("cation_pi")
                result = CationPiPairNode().execute(context)

            self.assertIn("cation_pi_pairs", result.results)
            self.assertTrue(Path(result.results["cation_pi_pairs"]["raw_pairs_file"]).exists())

    def test_annotation_node_writes_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            raw_file = temp / "raw_cation_pi.csv"
            _pi_pair_frame("cation_pi").to_csv(raw_file, index=False)
            context = PipelineContext(
                system_id="demo",
                topology=str(temp / "test.pdb"),
                trajectory_raw=str(temp / "traj.xtc"),
                output_dir=str(temp / "results"),
            )
            context.results["cation_pi_pairs"] = {"raw_pairs_file": str(raw_file)}
            context.metadata["chain_mapping"] = _chain_mapping()
            context.metadata["cdr_detection"] = {
                "chains": {"TCR_beta": {"chain_id": "E", "cdrs": {3: {"residue_range": [90, 100]}}}}
            }
            result = CationPiAnnotationNode().execute(context)
            outputs = result.results["cation_pi_annotation"]
            self.assertTrue(Path(outputs["cation_pi_report_file"]).exists())
            self.assertTrue(Path(outputs["groove_tcr_cation_pi_file"]).exists())


if __name__ == "__main__":
    unittest.main()
