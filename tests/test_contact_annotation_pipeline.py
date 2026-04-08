#!/usr/bin/env python3
"""
Smoke tests for residue contact annotation pipeline pieces.
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.topology import (
    ContactAnnotationAnnotator,
    ResiduePairAnnotationAnnotator,
    TCRResidueSemanticAnnotator,
    ComplexResidueSemanticAnnotator,
)
from immunex.core.context import PipelineContext
from immunex.pipeline.nodes.contact_annotation_node import ContactAnnotationNode
from immunex.pipeline.nodes.contact_frequency_node import ContactFrequencyNode


class TestContactAnnotationAnnotator(unittest.TestCase):
    def test_annotate_contacts_assigns_semantics(self):
        contacts = pd.DataFrame(
            [
                {
                    "chain_id_1": "C",
                    "resid_1": 7,
                    "resname_1": "GLY",
                    "residue_label_1": "GLY7",
                    "chain_id_2": "D",
                    "resid_2": 110,
                    "resname_2": "TYR",
                    "residue_label_2": "TYR110",
                    "contact_frames": 50,
                    "total_frames": 100,
                    "contact_frequency": 0.5,
                    "min_distance_observed": 3.2,
                }
            ]
        )
        chain_mapping = {
            "mhc_alpha": "A",
            "b2m": "B",
            "peptide": "C",
            "tcr_alpha": "D",
            "tcr_beta": "E",
        }
        cdr_detection = {
            "chains": {
                "TCR_alpha": {
                    "chain_id": "D",
                    "cdrs": {
                        3: {"residue_range": [105, 117]}
                    },
                }
            }
        }

        annotator = ContactAnnotationAnnotator()
        annotated = annotator.annotate_contacts(contacts, chain_mapping, cdr_detection)

        row = annotated.iloc[0]
        self.assertEqual(row["interaction_class"], "peptide_tcr")
        self.assertEqual(row["tcr_chain"], "alpha")
        self.assertEqual(row["tcr_region"], "CDR3")
        self.assertEqual(row["phla_region"], "peptide")
        self.assertEqual(row["mhc_region"], "peptide")
        self.assertEqual(row["partner_component"], "peptide")

    def test_pair_annotator_marks_hla_groove(self):
        pairs = pd.DataFrame(
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
                    "contact_frames": 25,
                    "total_frames": 100,
                    "contact_frequency": 0.25,
                    "min_distance_observed": 3.8,
                }
            ]
        )
        annotator = ResiduePairAnnotationAnnotator()
        annotated = annotator.annotate_pairs(
            pairs,
            chain_mapping={
                "mhc_alpha": "A",
                "b2m": "B",
                "peptide": "C",
                "tcr_alpha": "D",
                "tcr_beta": "E",
            },
            cdr_detection={
                "chains": {
                    "TCR_beta": {"chain_id": "E", "cdrs": {2: {"residue_range": [50, 60]}}},
                }
            },
        )
        row = annotated.iloc[0]
        self.assertEqual(row["interaction_class"], "hla_tcr")
        self.assertEqual(row["mhc_region_1"], "groove")
        self.assertEqual(row["mhc_region"], "groove")


class TestTCRResidueSemanticAnnotator(unittest.TestCase):
    def test_annotate_residue_assigns_component_and_cdr(self):
        annotator = TCRResidueSemanticAnnotator(
            chain_mapping={
                "mhc_alpha": "A",
                "b2m": "B",
                "peptide": "C",
                "tcr_alpha": "D",
                "tcr_beta": "E",
            },
            cdr_detection={
                "chains": {
                    "TCR_alpha": {"chain_id": "D", "cdrs": {1: {"residue_range": [27, 38]}}},
                    "TCR_beta": {"chain_id": "E", "cdrs": {3: {"residue_range": [105, 117]}}},
                }
            },
        )

        alpha_cdr1 = annotator.annotate_residue("D", 31)
        beta_non_cdr = annotator.annotate_residue("E", 90)
        peptide = annotator.annotate_residue("C", 5)

        self.assertEqual(alpha_cdr1.component, "TCR_alpha")
        self.assertEqual(alpha_cdr1.tcr_chain, "alpha")
        self.assertEqual(alpha_cdr1.tcr_region, "CDR1")
        self.assertEqual(alpha_cdr1.tcr_region_detailed, "CDR1_alpha")

        self.assertEqual(beta_non_cdr.component, "TCR_beta")
        self.assertEqual(beta_non_cdr.tcr_chain, "beta")
        self.assertEqual(beta_non_cdr.tcr_region, "non_cdr")

        self.assertEqual(peptide.component, "peptide")
        self.assertIsNone(peptide.tcr_chain)
        self.assertIsNone(peptide.tcr_region)


class TestComplexResidueSemanticAnnotator(unittest.TestCase):
    def test_annotate_residue_assigns_complex_side_and_phla_region(self):
        annotator = ComplexResidueSemanticAnnotator(
            chain_mapping={
                "mhc_alpha": "A",
                "b2m": "B",
                "peptide": "C",
                "tcr_alpha": "D",
                "tcr_beta": "E",
            },
            cdr_detection={
                "chains": {
                    "TCR_beta": {"chain_id": "E", "cdrs": {2: {"residue_range": [56, 65]}}},
                }
            },
        )

        hla = annotator.annotate_residue("A", 70)
        peptide = annotator.annotate_residue("C", 5)
        tcr_beta = annotator.annotate_residue("E", 60)

        self.assertEqual(hla.component, "HLA_alpha")
        self.assertEqual(hla.complex_side, "phla")
        self.assertEqual(hla.phla_region, "HLA_alpha")
        self.assertEqual(hla.mhc_region, "groove")

        self.assertEqual(peptide.component, "peptide")
        self.assertEqual(peptide.complex_side, "phla")
        self.assertEqual(peptide.phla_region, "peptide")
        self.assertEqual(peptide.mhc_region, "peptide")

        self.assertEqual(tcr_beta.component, "TCR_beta")
        self.assertEqual(tcr_beta.complex_side, "tcr")
        self.assertIsNone(tcr_beta.phla_region)
        self.assertIsNone(tcr_beta.mhc_region)
        self.assertEqual(tcr_beta.tcr_region, "CDR2")


class TestContactFrequencyNode(unittest.TestCase):
    def test_contact_frequency_node_writes_outputs(self):
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

            mock_contacts = pd.DataFrame(
                [
                    {
                        "chain_id_1": "C",
                        "resid_1": 7,
                        "resname_1": "GLY",
                        "residue_label_1": "GLY7",
                        "chain_id_2": "D",
                        "resid_2": 110,
                        "resname_2": "TYR",
                        "residue_label_2": "TYR110",
                        "contact_frames": 50,
                        "total_frames": 100,
                        "contact_frequency": 0.5,
                        "min_distance_observed": 3.2,
                        "cutoff_angstrom": 4.5,
                        "selection_1": "chainID A or chainID B or chainID C",
                        "selection_2": "chainID D or chainID E",
                    }
                ]
            )

            with patch(
                "immunex.pipeline.nodes.contact_frequency_node.ResidueContactFrequencyAnalyzer"
            ) as analyzer_cls:
                analyzer = analyzer_cls.return_value
                analyzer.calculate_residue_contact_frequencies.return_value = mock_contacts

                node = ContactFrequencyNode(cutoff=4.5)
                result = node.execute(context)

            self.assertIn("contact_frequency", result.results)
            raw_file = Path(result.results["contact_frequency"]["raw_contacts_file"])
            self.assertTrue(raw_file.exists())
            written = pd.read_csv(raw_file)
            self.assertEqual(len(written), 1)


class TestContactAnnotationNode(unittest.TestCase):
    def test_contact_annotation_node_writes_filtered_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            raw_file = temp_path / "raw_contacts.csv"
            pd.DataFrame(
                [
                    {
                        "chain_id_1": "C",
                        "resid_1": 7,
                        "resname_1": "GLY",
                        "residue_label_1": "GLY7",
                        "chain_id_2": "D",
                        "resid_2": 110,
                        "resname_2": "TYR",
                        "residue_label_2": "TYR110",
                        "contact_frames": 50,
                        "total_frames": 100,
                        "contact_frequency": 0.5,
                        "min_distance_observed": 3.2,
                    },
                    {
                        "chain_id_1": "A",
                        "resid_1": 65,
                        "resname_1": "ASP",
                        "residue_label_1": "ASP65",
                        "chain_id_2": "E",
                        "resid_2": 55,
                        "resname_2": "SER",
                        "residue_label_2": "SER55",
                        "contact_frames": 25,
                        "total_frames": 100,
                        "contact_frequency": 0.25,
                        "min_distance_observed": 3.8,
                    },
                ]
            ).to_csv(raw_file, index=False)

            context = PipelineContext(
                system_id="demo",
                topology=str(temp_path / "test.pdb"),
                trajectory_raw=str(temp_path / "traj.xtc"),
                output_dir=str(temp_path / "results"),
            )
            context.results["contact_frequency"] = {"raw_contacts_file": str(raw_file)}
            context.metadata["chain_mapping"] = {
                "mhc_alpha": "A",
                "b2m": "B",
                "peptide": "C",
                "tcr_alpha": "D",
                "tcr_beta": "E",
            }
            context.metadata["cdr_detection"] = {
                "chains": {
                    "TCR_alpha": {"chain_id": "D", "cdrs": {3: {"residue_range": [105, 117]}}},
                    "TCR_beta": {"chain_id": "E", "cdrs": {2: {"residue_range": [50, 60]}}},
                }
            }

            node = ContactAnnotationNode()
            result = node.execute(context)

            contact_outputs = result.results["contact_annotation"]
            annotated_file = Path(contact_outputs["annotated_contacts_file"])
            summary_file = Path(contact_outputs["summary_file"])
            groove_file = Path(contact_outputs["groove_tcr_contacts_file"])
            cdr3_dir = contact_outputs["cdr3"]
            cdr2_dir = contact_outputs["cdr2"]

            self.assertTrue(annotated_file.exists())
            self.assertTrue(summary_file.exists())
            self.assertTrue(groove_file.exists())
            self.assertTrue(Path(cdr3_dir["contacts_file"]).exists())
            self.assertTrue(Path(cdr3_dir["peptide_tcr_contacts_file"]).exists())
            self.assertTrue(Path(cdr2_dir["hla_tcr_contacts_file"]).exists())
            self.assertTrue(Path(cdr2_dir["groove_tcr_contacts_file"]).exists())

            cdr3_contacts = pd.read_csv(cdr3_dir["contacts_file"])
            cdr3_peptide_contacts = pd.read_csv(cdr3_dir["peptide_tcr_contacts_file"])
            cdr2_hla_contacts = pd.read_csv(cdr2_dir["hla_tcr_contacts_file"])
            cdr2_groove_contacts = pd.read_csv(cdr2_dir["groove_tcr_contacts_file"])
            global_groove_contacts = pd.read_csv(groove_file)
            self.assertEqual(len(cdr3_contacts), 1)
            self.assertEqual(len(cdr3_peptide_contacts), 1)
            self.assertEqual(len(cdr2_hla_contacts), 1)
            self.assertEqual(len(cdr2_groove_contacts), 1)
            self.assertEqual(len(global_groove_contacts), 1)


if __name__ == "__main__":
    unittest.main()
