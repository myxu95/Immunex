# Answer Templates

## Template: ranked answer

Direct answer:
`The strongest <entity> is <item>, based on <primary metric>.`

Evidence bullets:
- `Supported by <summary or digest source>.`
- `<item> ranks first by <metric>.`
- `<item> is annotated as <region or component>.`

## Template: region answer

Direct answer:
`The dominant region is <region>, which carries the strongest <signal>.`

Evidence bullets:
- `Supported by <region summary source>.`
- `<region> ranks first by <metric>.`
- `This signal is primarily <peptide-facing / groove-facing / cluster-specific>.`

## Template: cluster answer

Direct answer:
`Cluster <id> is dominant and accounts for <share>% of sampled frames.`

Evidence bullets:
- `Supported by interface clustering summary and population table.`
- `The system is best described as <single-dominant / distributed-state>.`
- `The dominant state is associated with <signature or hotspot region>, when available.`

## Template: uncertainty answer

Direct answer:
`The current answer is provisional because the preferred summary source is missing.`

Evidence bullets:
- `Missing source: <source name>.`
- `Fallback source used: <fallback source>.`
- `Interpretation should be treated as lower confidence.`
