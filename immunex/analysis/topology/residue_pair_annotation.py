"""
通用 residue-pair 语义注释器。
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import pandas as pd

from .complex_residue_semantics import ComplexResidueSemanticAnnotator


class ResiduePairAnnotationAnnotator:
    """
    为任意 residue-pair 结果表添加 pHLA/TCR 语义标签。

    当前主用于 contact pairs，后续 typed interactions 也应复用同一套字段。
    """

    EMPTY_COLUMNS = [
        "component_1",
        "complex_side_1",
        "phla_region_1",
        "mhc_region_1",
        "mhc_subregion_1",
        "tcr_chain_1",
        "tcr_region_1",
        "tcr_region_detailed_1",
        "component_2",
        "complex_side_2",
        "phla_region_2",
        "mhc_region_2",
        "mhc_subregion_2",
        "tcr_chain_2",
        "tcr_region_2",
        "tcr_region_detailed_2",
        "interaction_class",
        "tcr_chain",
        "tcr_region",
        "tcr_region_detailed",
        "phla_region",
        "mhc_region",
        "mhc_subregion",
        "tcr_residue_label",
        "partner_residue_label",
        "partner_component",
    ]

    def annotate_pairs(
        self,
        pairs: pd.DataFrame,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict] = None,
    ) -> pd.DataFrame:
        annotated = pairs.copy()
        if "interaction_family" not in annotated.columns:
            annotated["interaction_family"] = "contact"
        if annotated.empty:
            for column in self.EMPTY_COLUMNS:
                annotated[column] = pd.Series(dtype="object")
            return annotated

        semantic_annotator = ComplexResidueSemanticAnnotator(
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

        side1 = annotated.apply(
            lambda row: self._annotate_residue(
                chain_id=row["chain_id_1"],
                resid=int(row["resid_1"]),
                semantic_annotator=semantic_annotator,
            ),
            axis=1,
            result_type="expand",
        )
        side1.columns = [
            "component_1",
            "complex_side_1",
            "phla_region_1",
            "mhc_region_1",
            "mhc_subregion_1",
            "tcr_chain_1",
            "tcr_region_1",
            "tcr_region_detailed_1",
        ]

        side2 = annotated.apply(
            lambda row: self._annotate_residue(
                chain_id=row["chain_id_2"],
                resid=int(row["resid_2"]),
                semantic_annotator=semantic_annotator,
            ),
            axis=1,
            result_type="expand",
        )
        side2.columns = [
            "component_2",
            "complex_side_2",
            "phla_region_2",
            "mhc_region_2",
            "mhc_subregion_2",
            "tcr_chain_2",
            "tcr_region_2",
            "tcr_region_detailed_2",
        ]

        annotated = pd.concat([annotated, side1, side2], axis=1)

        interaction_info = annotated.apply(
            self._derive_interaction_semantics,
            axis=1,
            result_type="expand",
        )
        interaction_info.columns = [
            "interaction_class",
            "tcr_chain",
            "tcr_region",
            "tcr_region_detailed",
            "phla_region",
            "mhc_region",
            "mhc_subregion",
            "tcr_residue_label",
            "partner_residue_label",
            "partner_component",
        ]

        annotated = pd.concat([annotated, interaction_info], axis=1)
        return annotated

    @staticmethod
    def filter_pairs(
        pairs: pd.DataFrame,
        interaction_class: Optional[str] = None,
        tcr_region: Optional[str] = None,
        partner_component: Optional[str] = None,
        mhc_region: Optional[str] = None,
        mhc_subregion: Optional[str] = None,
    ) -> pd.DataFrame:
        filtered = pairs.copy()
        if interaction_class is not None:
            filtered = filtered[filtered["interaction_class"] == interaction_class]
        if tcr_region is not None:
            filtered = filtered[filtered["tcr_region"] == tcr_region]
        if partner_component is not None:
            filtered = filtered[filtered["partner_component"] == partner_component]
        if mhc_region is not None:
            filtered = filtered[filtered["mhc_region"] == mhc_region]
        if mhc_subregion is not None:
            filtered = filtered[filtered["mhc_subregion"] == mhc_subregion]
        return filtered.reset_index(drop=True)

    @staticmethod
    def summarize_pair_partners(
        pairs: pd.DataFrame,
        label_column: str,
    ) -> pd.DataFrame:
        if pairs.empty:
            return pd.DataFrame(
                columns=[label_column, "contact_frequency_sum", "contact_frequency_max", "n_pairs"]
            )

        return (
            pairs.groupby(label_column)["contact_frequency"]
            .agg(["sum", "max", "count"])
            .reset_index()
            .rename(
                columns={
                    "sum": "contact_frequency_sum",
                    "max": "contact_frequency_max",
                    "count": "n_pairs",
                }
            )
            .sort_values(
                by=["contact_frequency_sum", "contact_frequency_max"],
                ascending=[False, False],
            )
            .reset_index(drop=True)
        )

    @staticmethod
    def _annotate_residue(
        chain_id: str,
        resid: int,
        semantic_annotator: ComplexResidueSemanticAnnotator,
    ) -> Tuple[str, str, Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
        semantics = semantic_annotator.annotate_residue(chain_id=chain_id, resid=resid)
        return (
            semantics.component,
            semantics.complex_side,
            semantics.phla_region,
            semantics.mhc_region,
            semantics.mhc_subregion,
            semantics.tcr_chain,
            semantics.tcr_region,
            semantics.tcr_region_detailed,
        )

    @staticmethod
    def _derive_interaction_semantics(row) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
        side1_is_tcr = row["complex_side_1"] == "tcr"
        side2_is_tcr = row["complex_side_2"] == "tcr"

        if side1_is_tcr and not side2_is_tcr:
            tcr_side = 1
            partner_side = 2
        elif side2_is_tcr and not side1_is_tcr:
            tcr_side = 2
            partner_side = 1
        else:
            return (None, None, None, None, None, None, None, None, None, None)

        partner_component = row[f"component_{partner_side}"]
        partner_phla_region = row[f"phla_region_{partner_side}"]
        partner_mhc_region = row[f"mhc_region_{partner_side}"]
        partner_mhc_subregion = row[f"mhc_subregion_{partner_side}"]

        if partner_component == "peptide":
            interaction_class = "peptide_tcr"
            phla_region = "peptide"
        elif partner_component == "HLA_alpha":
            interaction_class = "hla_tcr"
            phla_region = "HLA_alpha"
        elif partner_component == "beta2m":
            interaction_class = "beta2m_tcr"
            phla_region = "beta2m"
        else:
            interaction_class = "other_tcr"
            phla_region = partner_component

        return (
            interaction_class,
            row[f"tcr_chain_{tcr_side}"],
            row[f"tcr_region_{tcr_side}"],
            row[f"tcr_region_detailed_{tcr_side}"],
            partner_phla_region or phla_region,
            partner_mhc_region,
            partner_mhc_subregion,
            row[f"residue_label_{tcr_side}"],
            row[f"residue_label_{partner_side}"],
            partner_component,
        )


class ContactAnnotationAnnotator(ResiduePairAnnotationAnnotator):
    """兼容旧名称；后续 typed interactions 与 contact 统一走 pair annotator。"""

    def annotate_contacts(
        self,
        contacts: pd.DataFrame,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict] = None,
    ) -> pd.DataFrame:
        return self.annotate_pairs(contacts, chain_mapping=chain_mapping, cdr_detection=cdr_detection)

    def filter_contacts(
        self,
        contacts: pd.DataFrame,
        interaction_class: Optional[str] = None,
        tcr_region: Optional[str] = None,
        partner_component: Optional[str] = None,
    ) -> pd.DataFrame:
        return self.filter_pairs(
            contacts,
            interaction_class=interaction_class,
            tcr_region=tcr_region,
            partner_component=partner_component,
        )

    def summarize_contact_partners(self, contacts: pd.DataFrame, label_column: str) -> pd.DataFrame:
        return self.summarize_pair_partners(contacts, label_column)
