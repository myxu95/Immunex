---
name: reporter-diagnostic
description: Use this skill to diagnose a single Immunex system based on existing analysis outputs, focusing on stability, interface organization, flexibility, mechanical hotspots, and state distribution.
---

# Reporter Diagnostic

## Scope

Use this skill for single-system diagnostic interpretation, for example:
- Is the system stable?
- Which regions dominate the interface?
- Are there clear hotspot residues or pairs?
- Does clustering indicate one dominant state or multiple states?
- Do NMA or PRS highlight important control residues?

Do not use this skill for:
- two-system comparison
- mutation or experiment design recommendation
- simple lookup questions

## Diagnostic Order

Organize conclusions in the following order:
1. `Quality`
2. `Interface / BSA`
3. `Landscape / Persistence`
4. `RRCS`
5. `Flexibility / RMSF`
6. `Cluster`
7. `NMA / Mechanical hotspots`

## Working Principles

1. State diagnostic conclusions first, then give evidence.
2. Consume structured summary and digest outputs instead of full raw tables.
3. Prioritize identifying:
   - stable dominant interfaces
   - peptide-facing versus groove-facing bias
   - strong hotspot pairs
   - dominant and minor states
   - mechanically important residues
4. Clearly separate:
   - evidence-backed findings
   - interpretation or inference

## Recommended Output Structure

- `Overall diagnosis`
- `Evidence highlights`
- `Potential hotspots`
- `Open questions / caution`

## Output Constraints

- Each diagnostic conclusion must point back to at least one analysis module.
- Do not make two-system difference claims.
- Do not jump directly to intervention design; only mark positions worth further attention.
