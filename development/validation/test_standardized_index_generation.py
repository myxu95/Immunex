#!/usr/bin/env python3
"""
Test Standardized Index Generation with Real Data

This script demonstrates the new standardized IndexGenerator API
by generating actual index files.

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.index_generation import (
    IndexGenerator,
    IndexGenerationInput,
    IndexGenerationMethod,
    ComponentDefinition
)


def test_pdb_based_generation():
    """Test PDB-based index generation with real PDB file."""
    print("=" * 70)
    print("Test 1: PDB-Based Index Generation")
    print("=" * 70)
    
    # Use available test PDB
    pdb_file = "test_batch/system_005/data/my_structure.pdb"
    output_dir = Path("development/test_output")
    output_dir.mkdir(exist_ok=True)
    
    # Check if file exists
    if not Path(pdb_file).exists():
        print(f"PDB file not found: {pdb_file}")
        print("Skipping PDB-based test")
        return False
    
    print(f"\nInput: {pdb_file}")
    print(f"Output: {output_dir / 'test_pdb_based.ndx'}")
    
    # Initialize generator
    generator = IndexGenerator()
    
    # Set progress callback
    def show_progress(progress, message):
        bar_width = 40
        filled = int(bar_width * progress)
        bar = '#' * filled + '-' * (bar_width - filled)
        print(f"\r[{bar}] {progress*100:3.0f}% - {message}", end='', flush=True)
    
    generator.set_progress_callback(show_progress)
    
    # Define input parameters
    input_params = IndexGenerationInput(
        topology=pdb_file,
        method=IndexGenerationMethod.PDB_BASED,
        components=[
            ComponentDefinition(name="pHLA", chains=['A', 'B', 'C']),
            ComponentDefinition(name="TCR", chains=['D', 'E']),
            ComponentDefinition(name="peptide", chains=['C']),
            ComponentDefinition(name="HLA_alpha", chains=['A']),
            ComponentDefinition(name="beta2m", chains=['B'])
        ],
        output_file=str(output_dir / "test_pdb_based.ndx"),
        auto_standardize=False  # Use PDB as-is
    )
    
    # Generate index
    print("\nGenerating index...")
    result = generator.generate(input_params)
    print()  # New line after progress bar
    
    # Display results
    if result.success:
        print("\n✓ SUCCESS!")
        print(f"\nOutput file: {result.output_file}")
        print(f"\nComponents generated: {len(result.components)}")
        for comp in result.components:
            print(f"  • {comp.name:12s} - Group {comp.group_id:2d} ({comp.atom_count:5d} atoms)")
        
        print(f"\nProcessing stats:")
        print(f"  • Time: {result.processing_stats.get('processing_time_sec', 0):.3f}s")
        print(f"  • Total atoms: {result.processing_stats.get('total_atoms', 0)}")
        
        print(f"\nMetadata:")
        print(f"  • Method: {result.metadata.get('method')}")
        print(f"  • Timestamp: {result.metadata.get('timestamp')}")
        
        # Show first few lines of generated file
        print(f"\nGenerated index file preview:")
        with open(result.output_file, 'r') as f:
            lines = f.readlines()[:20]
            for line in lines:
                print(f"  {line.rstrip()}")
        
        return True
    else:
        print("\n✗ FAILED!")
        print(f"Error: {result.error_message}")
        return False


def test_topology_based_generation():
    """Test topology-based index generation with TPR file."""
    print("\n\n" + "=" * 70)
    print("Test 2: Topology-Based Index Generation")
    print("=" * 70)
    
    # Use available test TPR
    tpr_file = "test_batch/system_005/data/my_topology.tpr"
    output_dir = Path("development/test_output")
    
    if not Path(tpr_file).exists():
        print(f"TPR file not found: {tpr_file}")
        print("Skipping topology-based test")
        return False
    
    print(f"\nInput: {tpr_file}")
    
    # First, identify chains from topology
    print("\nStep 1: Identifying chains from topology...")
    
    from immunex.utils.topology_chain_identifier import TopologyChainIdentifier
    
    identifier = TopologyChainIdentifier()
    chain_result = identifier.identify_chains_from_topology(
        topology_file=tpr_file
    )
    
    if not chain_result.success:
        print(f"✗ Chain identification failed: {chain_result.error_message}")
        return False
    
    print("✓ Chain identification successful!")
    print(f"\nDetected components:")
    for comp_name, group_id in sorted(chain_result.component_map.items()):
        chain = chain_result.chains[group_id]
        print(f"  • {comp_name:12s} -> Group {group_id:2d} "
              f"({chain.residue_count:3d} residues, {chain.atom_count:5d} atoms)")
    
    # Now generate index using identified group IDs
    print("\nStep 2: Generating index from topology groups...")
    
    generator = IndexGenerator()
    
    # Build components from identified chains
    component_map = chain_result.component_map
    
    input_params = IndexGenerationInput(
        topology=tpr_file,
        method=IndexGenerationMethod.TOPOLOGY_BASED,
        components=[
            ComponentDefinition(
                name="pHLA",
                group_ids=[
                    component_map.get('HLA_alpha'),
                    component_map.get('beta2m'),
                    component_map.get('peptide')
                ]
            ),
            ComponentDefinition(
                name="TCR",
                group_ids=[
                    component_map.get('TCR_alpha'),
                    component_map.get('TCR_beta')
                ]
            ),
            ComponentDefinition(
                name="peptide",
                group_ids=[component_map.get('peptide')]
            )
        ],
        output_file=str(output_dir / "test_topology_based.ndx")
    )
    
    result = generator.generate(input_params)
    
    if result.success:
        print("\n✓ SUCCESS!")
        print(f"\nOutput file: {result.output_file}")
        print(f"Components: {len(result.components)}")
        
        # Show file preview
        print(f"\nGenerated index file preview:")
        with open(result.output_file, 'r') as f:
            lines = f.readlines()[:15]
            for line in lines:
                print(f"  {line.rstrip()}")
        
        return True
    else:
        print(f"\n✗ FAILED: {result.error_message}")
        return False


def test_custom_selections():
    """Test custom MDAnalysis-based selections."""
    print("\n\n" + "=" * 70)
    print("Test 3: Custom Selection-Based Index Generation")
    print("=" * 70)
    
    pdb_file = "test_batch/system_005/data/my_structure.pdb"
    output_dir = Path("development/test_output")
    
    if not Path(pdb_file).exists():
        print(f"PDB file not found: {pdb_file}")
        print("Skipping custom selection test")
        return False
    
    print(f"\nInput: {pdb_file}")
    print(f"\nDefining custom selections...")
    
    generator = IndexGenerator()
    
    input_params = IndexGenerationInput(
        topology=pdb_file,
        method=IndexGenerationMethod.CUSTOM,
        components=[
            ComponentDefinition(
                name="backbone",
                selection="protein and name CA C N O"
            ),
            ComponentDefinition(
                name="CA_atoms",
                selection="protein and name CA"
            ),
            ComponentDefinition(
                name="hydrophobic",
                selection="protein and (resname ALA VAL LEU ILE PHE TRP MET)"
            )
        ],
        output_file=str(output_dir / "test_custom.ndx")
    )
    
    result = generator.generate(input_params)
    
    if result.success:
        print("\n✓ SUCCESS!")
        print(f"\nOutput file: {result.output_file}")
        print(f"\nCustom selections generated:")
        for comp in result.components:
            print(f"  • {comp.name:15s}: {comp.atom_count:5d} atoms, {comp.residue_count:4d} residues")
        
        return True
    else:
        print(f"\n✗ FAILED: {result.error_message}")
        return False


def main():
    """Run all tests."""
    print("\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + " " * 10 + "Testing Standardized Index Generation" + " " * 20 + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "═" * 68 + "╝")
    print()
    
    results = {}
    
    # Test 1: PDB-based
    try:
        results['pdb_based'] = test_pdb_based_generation()
    except Exception as e:
        print(f"\n✗ PDB-based test crashed: {e}")
        import traceback
        traceback.print_exc()
        results['pdb_based'] = False
    
    # Test 2: Topology-based
    try:
        results['topology_based'] = test_topology_based_generation()
    except Exception as e:
        print(f"\n✗ Topology-based test crashed: {e}")
        import traceback
        traceback.print_exc()
        results['topology_based'] = False
    
    # Test 3: Custom selections
    try:
        results['custom'] = test_custom_selections()
    except Exception as e:
        print(f"\n✗ Custom selection test crashed: {e}")
        import traceback
        traceback.print_exc()
        results['custom'] = False
    
    # Summary
    print("\n\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    total = len(results)
    passed = sum(1 for v in results.values() if v)
    failed = total - passed
    
    print(f"\nTotal tests: {total}")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    
    for test_name, success in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"  {status} - {test_name}")
    
    print("\n" + "=" * 70)
    
    if failed == 0:
        print("All tests passed!")
        return 0
    else:
        print(f"{failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
