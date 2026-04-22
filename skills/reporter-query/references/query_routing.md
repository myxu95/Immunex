# Query Routing

Use the smallest source that can answer the question.

## Residue Importance Questions

Examples:
- Which residues are most important?
- Which residues look like hotspots?

Preferred sources:
1. `annotated_rrcs_pair_summary.csv`
2. `rrcs_summary.json`
3. `rrcs_region_summary.csv`
4. contact landscape digest or persistence digest

Answer focus:
- top RRCS-supported residues
- region labels
- supporting module names
- make it explicit when the answer is hotspot-support driven rather than a dedicated residue-ranking module

## Region Importance Questions

Examples:
- Which region is most active?
- Which interface region is dominant?

Preferred sources:
1. region summary outputs
2. `rrcs_region_summary.csv`
3. `cluster_feature_digest.csv`
4. contact landscape digest

Answer focus:
- one dominant region or top 3 regions
- metric used for ranking
- whether the signal is peptide-facing or groove-facing

## RRCS Questions

Examples:
- Which RRCS hotspot is strongest?
- How many RRCS pairs were identified?

Preferred sources:
1. `rrcs_summary.json`
2. `annotated_rrcs_pair_summary.csv`
3. `rrcs_region_summary.csv`

Answer focus:
- identified pair count
- top hotspot pair
- dominant active region

## Cluster Questions

Examples:
- Which cluster is dominant?
- Is there a single dominant state?

Preferred sources:
1. `interface_clustering_summary.json`
2. `summary_table.csv`
3. `cluster_feature_digest.csv`

Answer focus:
- dominant cluster id
- population share
- whether the distribution is concentrated or distributed

## Persistence / Occupancy Questions

Examples:
- Which pair is most stable?
- Which interface is most persistent?

Preferred sources:
1. occupancy or persistence summary JSON
2. top persistent pair digest
3. occupancy region summary

Answer focus:
- top stable pair
- peptide-facing vs groove-facing balance
- persistence strength

## Quality / Stability Questions

Examples:
- Is the trajectory stable?
- Which quality metric is strongest?

Preferred sources:
1. quality summary JSON
2. RMSD digest
3. tail-90 summary outputs

Answer focus:
- overall stable / caution / unstable label
- strongest supporting metric
