# 1OGA Query Examples

These examples show how `reporter-query` should answer lightweight lookup questions using the available `1OGA` outputs.

## Example 1

### User question
Which cluster is dominant?

### Recommended answer shape
The dominant cluster is cluster 1, based on the interface clustering summary.

- Supported by `interface_clustering_summary.json` and `summary_table.csv`.
- Cluster 1 has the largest population share among all identified interface states.
- The system is best described as having a dominant state with smaller secondary states.

## Example 2

### User question
Which RRCS hotspot is strongest?

### Recommended answer shape
The strongest RRCS hotspot pair is the top-ranked annotated pair in the RRCS summary output.

- Supported by `rrcs_summary.json` and `annotated_rrcs_pair_summary.csv`.
- The answer should report the first-ranked hotspot pair by the RRCS ranking metric.
- The pair label should be returned in annotated form when available, not just raw chain and residue identifiers.

## Example 3

### User question
Which region is most active at the interface?

### Recommended answer shape
The dominant interface region is the top-ranked region in the region-level interaction or RRCS digest.

- Supported by `interaction_overview.json`, region-level contact digest, and `rrcs_region_summary.csv`.
- The answer should name the region directly, for example peptide-facing, groove-facing, or a TCR-region to pHLA-region pairing.
- If contact and RRCS point to the same region, say that the signal is consistently supported across modules.

## Example 4

### User question
Which residues are the most important?

### Recommended answer shape
The most important residues are the top hotspot-supported residues in the RRCS pair summary, cross-checked against region-level interface activity when available.

- Supported by `annotated_rrcs_pair_summary.csv` and `rrcs_summary.json`.
- State clearly that the current ranking is derived from RRCS hotspot support instead of a dedicated residue-ranking module.
- Return at most the top 5 residues by default.

## Example 5

### User question
Is this trajectory stable?

### Recommended answer shape
The trajectory appears stable, or should be treated with caution, based on the quality summary.

- Supported by `quality_summary`-style outputs such as RMSD digest or tail-90 variation.
- The answer should use explicit wording: stable, caution, or unstable.
- The strongest supporting metric should be named directly.
