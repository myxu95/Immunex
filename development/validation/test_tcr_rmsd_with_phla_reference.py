#!/usr/bin/env python3
"""
Calculate TCR RMSD with pHLA as reference frame

This script demonstrates:
1. Standard PBC processing (no modification)
2. Generate index file with pHLA and TCR groups
3. Calculate TCR RMSD using pHLA as alignment reference

Key concept: PBC processing uses default settings, but RMSD calculation
uses pHLA as the reference frame to measure TCR movement.

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.index_generation import (
    IndexGenerator, IndexGenerationInput,
    IndexGenerationMethod, ComponentDefinition
)
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def main():
    """
    Workflow:
    1. Generate index file with pHLA and TCR groups
    2. (PBC processing - use standard workflow, not shown here)
    3. Calculate TCR RMSD with pHLA alignment reference
    """

    # Test PDB file
    structure_pdb = "input/standardizedpdbs/raw_data/1AO7_SOL.pdb"

    if not Path(structure_pdb).exists():
        logger.error(f"Test PDB file not found: {structure_pdb}")
        return 1

    # Output directory
    output_dir = Path("development/test_output/tcr_rmsd_phla_ref")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("TCR RMSD Calculation with pHLA as Reference Frame")
    print("=" * 70)
    print()

    # Step 1: Generate index file
    print("Step 1: Generating index file with pHLA and TCR groups")
    print("-" * 70)

    index_file = str(output_dir / "index.ndx")

    components = [
        ComponentDefinition(
            name="pHLA",
            selection="chainID A B C"  # MHC alpha + beta2m + peptide
        ),
        ComponentDefinition(
            name="pHLA_backbone",
            selection="chainID A B C and (name CA C N O)"
        ),
        ComponentDefinition(
            name="TCR",
            selection="chainID D E"  # TCR alpha + beta
        ),
        ComponentDefinition(
            name="TCR_backbone",
            selection="chainID D E and (name CA C N O)"
        ),
    ]

    index_gen = IndexGenerator()
    index_input = IndexGenerationInput(
        topology=structure_pdb,
        method=IndexGenerationMethod.CUSTOM,
        standardized_pdb=structure_pdb,
        components=components,
        output_file=index_file
    )

    result = index_gen.generate(index_input)

    if not result.success:
        logger.error(f"Index generation failed: {result.error_message}")
        return 1

    logger.info(f"Generated index file: {index_file}")
    for comp in result.components:
        logger.info(f"  - {comp.name}: {comp.atom_count} atoms")

    print()

    # Step 2: PBC processing (standard workflow - no changes needed)
    print("Step 2: PBC processing (standard workflow)")
    print("-" * 70)
    print("Use normal PBC processing:")
    print("""
from immunex.core.pbc_processor import PBCProcessor

pbc_processor = PBCProcessor()
pbc_processor.remove_pbc_2step(
    trajectory="md.xtc",
    topology="md.tpr",
    output="md_pbc.xtc"
)
# This uses default Backbone alignment - no changes needed!
    """)
    print()

    # Step 3: Calculate TCR RMSD with pHLA reference
    print("Step 3: Calculate TCR RMSD (using pHLA as reference)")
    print("-" * 70)
    print()
    print("Key command:")
    print(f"gmx rms -f md_pbc.xtc -s md.tpr -n {index_file} -o tcr_rmsd.xvg")
    print()
    print("When prompted:")
    print("  Select group 1 (pHLA_backbone) for least squares fit")
    print("  Select group 2 (TCR_backbone) for RMSD calculation")
    print()
    print("Or use stdin input:")
    print('subprocess.run([...], input="pHLA_backbone\\nTCR_backbone\\n", text=True)')
    print()

    # Step 4: Interpretation
    print("=" * 70)
    print("Biological Interpretation")
    print("=" * 70)
    print()
    print("What does this RMSD measure?")
    print("  - Fit on pHLA: The system is aligned using pHLA as reference")
    print("  - Calculate TCR: RMSD measures how much TCR moves relative to pHLA")
    print()
    print("Biological meaning:")
    print("  - pHLA is the 'fixed' reference (stable MHC-peptide complex)")
    print("  - TCR RMSD shows receptor dynamics during recognition")
    print("  - Lower RMSD = stable binding, higher RMSD = dynamic recognition")
    print()
    print("Comparison:")
    print("  - Standard RMSD: Fit on Backbone, calc Backbone → overall stability")
    print("  - This approach: Fit on pHLA, calc TCR → TCR-pMHC dynamics")
    print()

    # Summary
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print()
    print("Generated files:")
    print(f"  - Index file: {index_file}")
    print()
    print("Complete workflow:")
    print("  1. Generate index with pHLA and TCR groups ✓")
    print("  2. Run standard PBC processing (no changes)")
    print("  3. Run: gmx rms -n index.ndx")
    print("     - Select pHLA_backbone for fit")
    print("     - Select TCR_backbone for RMSD")
    print()
    print("Key point:")
    print("  PBC processing unchanged - customization is in RMSD calculation step!")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
