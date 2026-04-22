# Evidence Rules

## Answer Shape

Always use this order:
1. One direct answer sentence.
2. 3 to 5 evidence bullets.
3. Optional top-5 ranking only when ranking is central to the question.

## Evidence Priority

Prefer evidence in this order:
1. explicit summary fields
2. digest tables
3. region-level summaries
4. top-pair tables
5. full raw tables

## Evidence Wording

Use explicit support wording, for example:
- Supported by RRCS summary and annotated hotspot pairs.
- Supported by cluster summary and population table.
- Supported by NMA residue ranking and region summary.

## Uncertainty Handling

If the ideal source is missing:
- say which source is missing
- answer from the next-best source
- label the answer as provisional when needed

## Ranking Rules

When ranking residues or pairs:
- return at most top 5 by default
- prefer meaningful labels over raw identifiers when available
- include the ranking metric in the wording

## Forbidden Behavior

- Do not invent missing evidence.
- Do not turn a single weak signal into a strong conclusion.
- Do not expand into design recommendations.
- Do not rewrite raw tables into long prose.
