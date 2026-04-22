"""
PDB Downloader Usage Examples

This script demonstrates various use cases of the PDBDownloader class
for downloading PDB structures from RCSB PDB database.

Examples include:
1. Basic single PDB download
2. Batch download from a list
3. Query assembly information
4. Complete processing pipeline integration
5. Error handling patterns
"""

import logging
from pathlib import Path

from immunex.utils import (
    PDBDownloader,
    PDBStructureFixer,
    PDBChainStandardizer,
    PDBDistanceTrimmer
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def example1_basic_download():
    """
    Example 1: Basic single PDB download with auto assembly selection.
    """
    print("\n" + "="*60)
    print("Example 1: Basic Single PDB Download")
    print("="*60)

    # Initialize downloader with default settings
    downloader = PDBDownloader(
        download_dir="input/standardizedpdbs/pdbs_raw",
        assembly_preference='auto'  # Recommended default
    )

    # Download a single PDB
    pdb_id = '1ao7'
    result = downloader.download_pdb(pdb_id)

    if result['success']:
        print(f"\nSuccess! Downloaded {pdb_id}")
        print(f"  File: {result['file_path']}")
        print(f"  Assembly: {result['assembly_number']}")
    else:
        print(f"\nFailed to download {pdb_id}")
        print(f"  Reason: {result['message']}")


def example2_batch_download():
    """
    Example 2: Batch download multiple PDB files in parallel.
    """
    print("\n" + "="*60)
    print("Example 2: Batch Download")
    print("="*60)

    # Initialize downloader
    downloader = PDBDownloader(
        download_dir="input/standardizedpdbs/pdbs_raw",
        assembly_preference='auto'
    )

    # List of PDB IDs to download
    pdb_ids = ['1ao7', '1bd2', '1fo0', '1g4m', '1hhh']

    # Batch download with 4 parallel workers
    results = downloader.batch_download(
        pdb_ids,
        max_workers=4,
        force=False  # Skip existing files
    )

    # Print summary
    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]

    print(f"\nBatch download completed:")
    print(f"  Total: {len(results)}")
    print(f"  Successful: {len(successful)}")
    print(f"  Failed: {len(failed)}")

    if failed:
        print("\nFailed downloads:")
        for r in failed:
            print(f"  {r['pdb_id']}: {r['message']}")


def example3_query_assembly_info():
    """
    Example 3: Query assembly information without downloading.
    """
    print("\n" + "="*60)
    print("Example 3: Query Assembly Information")
    print("="*60)

    downloader = PDBDownloader()

    pdb_id = '1ao7'
    print(f"\nQuerying assembly information for {pdb_id}...")

    # Get entry metadata
    metadata = downloader.get_entry_metadata(pdb_id)
    if metadata:
        entry_info = metadata.get('rcsb_entry_info', {})
        print(f"  Title: {metadata.get('struct', {}).get('title', 'N/A')}")
        print(f"  Resolution: {entry_info.get('resolution_combined', ['N/A'])[0]} Å")

    # Get assembly information
    assemblies = downloader.get_assembly_info(pdb_id)
    if assemblies:
        print(f"\nAvailable biological assemblies: {len(assemblies)}")
        for asm in assemblies:
            preferred = " (PREFERRED)" if asm['preferred'] else ""
            print(
                f"  Assembly {asm['assembly_id']}: "
                f"{asm['oligomeric_count']}{preferred}"
            )
    else:
        print(f"No assembly information found for {pdb_id}")


def example4_different_assembly_strategies():
    """
    Example 4: Different assembly selection strategies.
    """
    print("\n" + "="*60)
    print("Example 4: Different Assembly Selection Strategies")
    print("="*60)

    strategies = {
        'auto': 'Auto-select (prefer assembly 1)',
        'first': 'Fast mode (always assembly 1)',
        'asymmetric': 'Asymmetric unit only',
        1: 'Specific assembly number (1)',
        2: 'Specific assembly number (2)'
    }

    pdb_id = '1ao7'

    for strategy, description in strategies.items():
        print(f"\nStrategy: {strategy} - {description}")

        downloader = PDBDownloader(
            download_dir="input/standardizedpdbs/pdbs_raw",
            assembly_preference=strategy
        )

        # Just show what would be selected (don't actually download)
        try:
            url, assembly_num = downloader.select_assembly(pdb_id)
            print(f"  Would download: {url}")
            print(f"  Assembly number: {assembly_num}")
        except Exception as e:
            print(f"  Error: {e}")


def example5_complete_pipeline():
    """
    Example 5: Complete processing pipeline integration.

    Download → Fix → Standardize → Trim
    """
    print("\n" + "="*60)
    print("Example 5: Complete Processing Pipeline")
    print("="*60)

    pdb_id = '1ao7'

    # Setup directories
    base_dir = Path("input/standardizedpdbs")
    raw_dir = base_dir / "pdbs_raw"
    fix_dir = base_dir / "fix_pdbs"
    protein_dir = base_dir / "protein_only"
    trimmed_dir = base_dir / "trimmed_sd_pdbs"

    # Create directories
    for d in [raw_dir, fix_dir, protein_dir, trimmed_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Step 1: Download biological assembly
    print("\n[Step 1] Downloading biological assembly...")
    downloader = PDBDownloader(download_dir=str(raw_dir))
    result = downloader.download_pdb(pdb_id, force=False)

    if not result['success']:
        print(f"Download failed: {result['message']}")
        return

    raw_file = Path(result['file_path'])
    print(f"  Downloaded: {raw_file}")

    # Step 2: Fix structure (add missing atoms, etc.)
    print("\n[Step 2] Fixing structure...")
    fixed_file = fix_dir / f"{pdb_id}_fixed.pdb"

    fixer = PDBStructureFixer(
        skip_terminal_residues=True,
        remove_heterogens=False
    )
    fixer.fix_structure(str(raw_file), str(fixed_file))
    print(f"  Fixed: {fixed_file}")

    # Step 3: Standardize chain IDs
    print("\n[Step 3] Standardizing chain IDs...")
    std_file = protein_dir / f"{pdb_id}_std.pdb"

    standardizer = PDBChainStandardizer(
        use_intelligent_identification=True
    )
    std_result = standardizer.process_single(str(fixed_file), str(std_file))
    print(f"  Standardized: {std_file}")
    print(f"  Chain mapping: {std_result.chain_mapping}")

    # Step 4: Trim by distance from peptide
    print("\n[Step 4] Trimming by distance...")
    trimmed_file = trimmed_dir / f"{pdb_id}_sd.pdb"

    trimmer = PDBDistanceTrimmer(distance_cutoff_nm=0.8)
    stats = trimmer.trim_structure(str(std_file), str(trimmed_file))
    print(f"  Trimmed: {trimmed_file}")
    print(f"  Kept {stats['atoms_kept']}/{stats['total_atoms']} atoms")

    print("\n[Pipeline Complete]")
    print(f"  Final file: {trimmed_file}")


def example6_error_handling():
    """
    Example 6: Error handling patterns.
    """
    print("\n" + "="*60)
    print("Example 6: Error Handling")
    print("="*60)

    downloader = PDBDownloader()

    # Test cases
    test_cases = [
        ('1ao7', 'Valid PDB ID'),
        ('XXXX', 'Invalid/non-existent PDB'),
        ('abc', 'Too short PDB ID'),
        ('12345', 'Too long PDB ID')
    ]

    for pdb_id, description in test_cases:
        print(f"\nTest: {description} ({pdb_id})")
        result = downloader.download_pdb(pdb_id)

        if result['success']:
            print(f"  Status: Success")
            print(f"  File: {result['file_path']}")
        else:
            print(f"  Status: Failed")
            print(f"  Message: {result['message']}")


def example7_batch_from_file():
    """
    Example 7: Batch download from a text file.
    """
    print("\n" + "="*60)
    print("Example 7: Batch Download from File")
    print("="*60)

    # Create a sample PDB ID list file
    list_file = Path("example_pdb_list.txt")
    pdb_ids = ['1ao7', '1bd2', '1fo0', '1g4m', '1hhh']

    print(f"\nCreating list file: {list_file}")
    with open(list_file, 'w') as f:
        for pdb_id in pdb_ids:
            f.write(f"{pdb_id}\n")
    print(f"  Added {len(pdb_ids)} PDB IDs")

    # Read from file and download
    print(f"\nReading PDB IDs from {list_file}...")
    with open(list_file, 'r') as f:
        pdb_list = [line.strip() for line in f if line.strip()]

    print(f"  Found {len(pdb_list)} PDB IDs")

    downloader = PDBDownloader()
    results = downloader.batch_download(pdb_list, max_workers=4)

    successful = sum(1 for r in results if r['success'])
    print(f"\nDownload completed: {successful}/{len(results)} successful")

    # Clean up
    list_file.unlink()
    print(f"\nCleaned up temporary file: {list_file}")


def example8_verify_assembly():
    """
    Example 8: Verify that biological assembly was downloaded.
    """
    print("\n" + "="*60)
    print("Example 8: Verify Biological Assembly")
    print("="*60)

    downloader = PDBDownloader(assembly_preference='auto')

    pdb_id = '1ao7'
    result = downloader.download_pdb(pdb_id)

    if not result['success']:
        print(f"Download failed: {result['message']}")
        return

    pdb_file = Path(result['file_path'])
    print(f"\nVerifying assembly information in {pdb_file}...")

    # Check for REMARK 350 (biological assembly information)
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    remark_350 = [line for line in lines if line.startswith('REMARK 350')]

    if remark_350:
        print(f"  Found {len(remark_350)} REMARK 350 lines")
        print(f"  This confirms biological assembly was downloaded")
        print(f"\nFirst few lines:")
        for line in remark_350[:5]:
            print(f"  {line.rstrip()}")
    else:
        print(f"  No REMARK 350 found")
        print(f"  This may be an asymmetric unit")


def main():
    """Run all examples."""
    print("="*60)
    print("PDB Downloader Usage Examples")
    print("="*60)

    examples = [
        ("Basic Download", example1_basic_download),
        ("Batch Download", example2_batch_download),
        ("Query Assembly Info", example3_query_assembly_info),
        ("Assembly Strategies", example4_different_assembly_strategies),
        ("Complete Pipeline", example5_complete_pipeline),
        ("Error Handling", example6_error_handling),
        ("Batch from File", example7_batch_from_file),
        ("Verify Assembly", example8_verify_assembly)
    ]

    print("\nAvailable examples:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")

    print("\nRunning all examples...\n")

    for name, func in examples:
        try:
            func()
        except Exception as e:
            logger.error(f"Error in {name}: {e}", exc_info=True)

    print("\n" + "="*60)
    print("All examples completed!")
    print("="*60)


if __name__ == '__main__':
    # You can run specific examples or all examples
    import sys

    if len(sys.argv) > 1:
        example_num = int(sys.argv[1])
        examples = [
            example1_basic_download,
            example2_batch_download,
            example3_query_assembly_info,
            example4_different_assembly_strategies,
            example5_complete_pipeline,
            example6_error_handling,
            example7_batch_from_file,
            example8_verify_assembly
        ]

        if 1 <= example_num <= len(examples):
            examples[example_num - 1]()
        else:
            print(f"Invalid example number. Choose 1-{len(examples)}")
    else:
        # Run all examples
        main()
