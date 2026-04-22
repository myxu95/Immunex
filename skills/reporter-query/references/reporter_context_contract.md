# Reporter Context Contract

`reporter-query` should prefer a prebuilt `reporter_context.json` when available.

## Preferred Contract

The following top-level keys are preferred:
- `system_id`
- `identity_summary`
- `quality_summary`
- `interface_summary`
- `contact_digest`
- `persistence_digest`
- `rrcs_digest`
- `cluster_digest`
- `nma_digest`

## Minimal Query-Safe Fields

### Residue importance
- `top_key_residues`
- `top_hinges`
- `top_hotspot_pairs`
- `top_rrcs_pairs`

### Region importance
- `top_regions`
- `region_summary`
- `dominant_interface_region`

### Cluster lookup
- `dominant_cluster`
- `cluster_population`
- `is_single_dominant_state`

### RRCS lookup
- `identified_pairs`
- `top_rrcs_pair`
- `dominant_rrcs_region`

### Persistence lookup
- `top_persistent_pair`
- `peptide_facing_fraction`
- `groove_facing_fraction`

### Quality lookup
- `quality_label`
- `tail90_variation`
- `mean_rmsd`

## Fallback Rule

If `reporter_context.json` is missing, query the smallest summary or digest file that can answer the question.
