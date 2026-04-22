"""兼容导出：contact 专用注释器已收敛到通用 pair annotator。"""

from .residue_pair_annotation import ContactAnnotationAnnotator, ResiduePairAnnotationAnnotator

__all__ = ["ResiduePairAnnotationAnnotator", "ContactAnnotationAnnotator"]
