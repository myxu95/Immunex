---
name: reporter-compare
description: Use this skill for two-condition or two-complex comparison based on Immunex results, including WT vs mutant, binder vs weak binder, standard MD vs REST2, or different binding modes.
---

# Reporter Compare

## Scope

Use this skill when the task is to compare two systems or conditions, for example:
- WT vs mutant
- binder vs weak binder
- standard MD vs REST2
- different binding modes
- replicate A vs replicate B

Do not use this skill for:
- single-system diagnosis
- simple query tasks
- design recommendation

## Comparison Order

Organize the comparison in this order:
1. `Comparability`
2. `Identity / context`
3. `Quality & interface`
4. `Flexibility shift`
5. `Landscape / persistence shift`
6. `RRCS shift`
7. `Cluster state shift`
8. `Mechanical hotspot shift`

## Working Principles

1. Use neutral comparative language throughout:
   - `Condition A`
   - `Condition B`
2. Assess comparability before interpreting differences.
3. Do not turn unsupported differences into causal claims.
4. Prioritize:
   - strongest differences
   - region-level shifts
   - hotspot gain/loss

## Recommended Output Structure

- `Comparison takeaway`
- `Major shifts`
- `Hotspot differences`
- `Interpretation and caution`

## Output Constraints

- Do not rewrite full comparison tables as long prose.
- Every difference statement should resolve to a concrete region or hotspot.
- If the inputs are not strictly comparable, state the limitation explicitly.
