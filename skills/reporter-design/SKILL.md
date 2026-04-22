---
name: reporter-design
description: Use this skill for advanced design-oriented recommendations based on Immunex outputs, including candidate mutation sites, validation targets, and regions that deserve focused visualization or follow-up analysis.
---

# Reporter Design

## Scope

Use this skill when the task enters a design-oriented stage, for example:
- recommend candidate mutation sites
- recommend experimental validation sites
- recommend priority visualization regions
- recommend regions for deeper enhanced-sampling follow-up

This skill depends on prior evidence, such as:
- diagnostic conclusions
- hotspot candidates
- comparison results when applicable

Do not use this skill when evidence is too weak to support actionable recommendations.

## Design Order

1. Confirm that candidate residues truly participate in the interface.
2. Check RRCS, persistence, and top-pair support.
3. Check cluster, NMA, and PRS support.
4. Evaluate penalties:
   - overly core-like
   - non-interface
   - excessively rigid
5. Then output recommendations.

## Recommendation Types

### Candidate Mutation Sites
Prefer residues that are:
- clearly interface-engaged
- strongly supported by hotspot evidence
- not obvious structural core residues
- interpretable in terms of a plausible engineering direction

### Candidate Validation Sites
Prefer residues that are:
- supported by multiple modules
- suitable for alanine scanning, focused mutation, or binding assays

### Priority Visualization Regions
Prefer regions such as:
- CDR3-peptide hotspots
- groove-facing hotspots
- cluster-discriminating signatures
- NMA/PRS-supported mechanical hotspots

## Recommended Output Structure

- `Candidate mutation sites`
- `Candidate validation sites`
- `Priority visualization regions`
- `Rationale and risk`

## Output Constraints

- Every recommendation must cite evidence sources.
- Recommendations must be labeled as `high priority`, `medium priority`, or `exploratory`.
- Do not present pure speculation as a confirmed conclusion.
