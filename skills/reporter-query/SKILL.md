---
name: reporter-query
description: Use this skill for lightweight lookup questions over existing Immunex analysis outputs, such as identifying the most important residues, the most active regions, the dominant clusters, the strongest RRCS hotspots, or the most persistent interface pairs.
---

# Reporter Query

Use this skill when the user asks short, direct, single-system questions that can be answered from existing Immunex results without generating a full diagnosis or recommendation.

## Use This Skill For

Examples:
- Which residues are the most important?
- Which hotspot pairs are strongest?
- Which region has the strongest contact signal?
- Which cluster is dominant?
- Which RRCS hotspot is strongest?
- Which persistent pair is most stable?

Do not use this skill for:
- full single-system diagnosis
- two-system comparison
- mutation or experiment recommendation
- open-ended mechanistic interpretation

## Workflow

1. Identify the query type.
2. Read the smallest relevant summary or digest source.
3. Prefer ranked or region-level outputs over full raw tables.
4. Return one direct answer sentence first.
5. Support the answer with 3 to 5 evidence bullets.
6. If ranking is involved, provide at most top 5 items.

## Source Selection Order

1. If `reporter_context.json` exists, use it first.
2. Otherwise prefer summary and digest files over full outputs.
3. Only fall back to full tables if the answer cannot be recovered from summary or digest files.

## Read These References

- For question-to-source routing: `references/query_routing.md`
- For evidence ordering and answer style: `references/evidence_rules.md`

## Output Constraints

- Keep the answer short and direct.
- Do not generate long narrative summaries.
- Do not expand into design suggestions.
- Do not introduce claims unsupported by analysis results.
- State uncertainty explicitly when the available outputs are incomplete.
