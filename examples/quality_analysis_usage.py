#!/usr/bin/env python3
"""
Immunex Quality Analysis Usage Example

Demonstrates how to use the MD Production quality analysis module
for comprehensive MD simulation quality control.
"""

import sys
from pathlib import Path
import logging

# Add the parent directory to Python path to import immunex
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.quality import (
    MDCompletenessChecker,
    StructureValidator,
    BatchTracker,
    QualityReporter
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def example_md_completeness_check():
    """Example: Check MD simulation completeness."""
    print("\n" + "="*60)
    print("MD Completeness Check Example")
    print("="*60)
    
    # Initialize checker with custom parameters
    checker = MDCompletenessChecker(
        min_trajectory_size_mb=2.0,     # Minimum 2MB trajectory
        min_simulation_time_ps=5000.0,  # Minimum 5ns simulation
        required_files=['md.gro', 'md.xtc', 'md.log']
    )
    
    # Example directory paths (replace with actual paths)
    example_dirs = [
        "/path/to/md/simulation1",
        "/path/to/md/simulation2", 
        "/path/to/md/simulation3"
    ]
    
    print(f"Checking {len(example_dirs)} MD simulation directories...")
    
    # Check individual directory
    if Path(example_dirs[0]).exists():
        result = checker.check_directory(example_dirs[0])
        print(f"\nSample result for {example_dirs[0]}:")
        print(f"  Status: {result['status']}")
        print(f"  Completeness Score: {result.get('completeness_score', 0):.1f}/100")
        print(f"  Issues: {result.get('issues', [])}")
    
    # Batch check all directories
    results = []
    for md_dir in example_dirs:
        if Path(md_dir).exists():
            result = checker.check_directory(md_dir)
            results.append(result)
    
    if results:
        # Generate summary statistics
        summary = checker.get_summary_statistics(results)
        print(f"\nBatch Summary:")
        print(f"  Total directories: {summary['total_directories']}")
        print(f"  Complete: {summary['complete']}")
        print(f"  Incomplete: {summary['incomplete']}")
        print(f"  Success rate: {summary['success_rate']:.1f}%")
        print(f"  Average score: {summary['average_completeness_score']:.1f}")
    
    return results


def example_structure_validation():
    """Example: Validate PDB structures for chain analysis."""
    print("\n" + "="*60)
    print("Structure Validation Example (PDB Files Only)")
    print("="*60)
    
    # Initialize validator
    validator = StructureValidator(
        expected_chain_count=5,           # Expecting 5 protein chains
        max_coordinate_range=500.0,       # Max 500 Angstrom coordinate range
        max_missing_residue_ratio=0.1     # Max 10% missing residues
    )
    
    # Example structure files - NOTE: Only PDB files are supported for chain analysis
    example_structures = [
        "/path/to/structure1.pdb",        # PDB format - contains chain information
        "/path/to/structure2.pdb",        # PDB format - contains chain information  
        "/path/to/structure3.ent",        # ENT format - same as PDB
        "/path/to/structure4.gro"         # GRO format - will be skipped (no chain info)
    ]
    
    print(f"Processing {len(example_structures)} structure files...")
    print("Note: Only PDB/ENT files will be validated (GRO files lack chain information)")
    
    results = []
    for structure_file in example_structures:
        if Path(structure_file).exists():
            result = validator.validate_structure(structure_file)
            results.append(result)
            
            print(f"\nValidation result for {Path(structure_file).name}:")
            print(f"  Status: {result['status']}")
            
            if result['status'] == 'skipped':
                print(f"  Reason: {result.get('message', 'Unknown')}")
                print(f"  File format: {result.get('file_format', 'Unknown')}")
            else:
                print(f"  Validation Score: {result.get('validation_score', 0):.1f}/100")
                
                # Chain analysis details (only for PDB files)
                chain_analysis = result.get('chain_analysis', {})
                if 'protein_chains' in chain_analysis:
                    print(f"  Protein chains: {chain_analysis['protein_chains']}")
                    print(f"  Chain count normal: {chain_analysis.get('chain_count_normal', False)}")
                    if 'chain_ids' in chain_analysis:
                        print(f"  Chain IDs: {chain_analysis['chain_ids']}")
    
    if results:
        # Get validation summary
        summary = validator.get_validation_summary(results)
        print(f"\nValidation Summary:")
        print(f"  Total files processed: {summary['total_structures']}")
        print(f"  Valid structures: {summary['valid']}")
        print(f"  Acceptable structures: {summary.get('acceptable', 0)}")
        print(f"  Problematic structures: {summary['problematic']}")
        print(f"  Validation rate: {summary['validation_rate']:.1f}%")
        
        # Chain count distribution
        chain_dist = summary.get('chain_count_distribution', {})
        if 'mean' in chain_dist:
            print(f"  Average chain count: {chain_dist['mean']:.1f}")
            print(f"  Expected chain count: {chain_dist['expected']}")
    
    print(f"\nImportant Note:")
    print(f"- StructureValidator is designed specifically for PDB files")
    print(f"- GRO files do not contain chain information and will be skipped")
    print(f"- For GRO file validation, use coordinate and atom count checks only")
    
    return results


def example_batch_tracking():
    """Example: Track batch processing and detect duplicates."""
    print("\n" + "="*60)
    print("Batch Tracking Example")
    print("="*60)
    
    # Initialize batch tracker
    tracker = BatchTracker(
        pdb_id_pattern=r'^([A-Za-z0-9]{4})',  # Extract first 4 characters as PDB ID
        tracking_db_path="./batch_tracking.json"
    )
    
    # Example batch directories
    example_batches = [
        "/path/to/batch1",
        "/path/to/batch2",
        "/path/to/batch3"
    ]
    
    batch_results = []
    
    for batch_dir in example_batches:
        if Path(batch_dir).exists():
            print(f"\nAnalyzing batch: {batch_dir}")
            
            # Analyze batch
            batch_data = tracker.analyze_batch(batch_dir)
            batch_results.append(batch_data)
            
            stats = batch_data['statistics']
            print(f"  Total directories: {stats['total_directories']}")
            print(f"  Identified PDbs: {stats['identified_pdbs']}")
            print(f"  Unique PDbs: {stats['unique_pdbs']}")
            print(f"  Duplicates in batch: {stats['duplicate_pdbs']}")
    
    # Cross-batch analysis
    print(f"\nCross-batch Analysis:")
    cross_duplicates = tracker.find_duplicates_across_batches()
    print(f"  PDB IDs in multiple batches: {len(cross_duplicates)}")
    
    # Global statistics
    global_stats = tracker.get_global_statistics()
    print(f"  Total batches tracked: {global_stats['total_batches']}")
    print(f"  Total unique PDB IDs: {global_stats['unique_pdbs']}")
    print(f"  Identification rate: {global_stats['identification_rate']:.1f}%")
    
    # Example: Check for missing PDB IDs
    expected_pdbs = {"1ABC", "2DEF", "3GHI", "4JKL", "5MNO"}  # Example expected set
    missing_analysis = tracker.find_missing_pdbs(expected_pdbs)
    print(f"\nMissing PDB Analysis:")
    print(f"  Expected: {missing_analysis['expected_total']}")
    print(f"  Processed: {missing_analysis['processed_total']}")
    print(f"  Missing: {missing_analysis['missing_count']}")
    print(f"  Completion rate: {missing_analysis['completion_rate']:.1f}%")
    
    return batch_results, global_stats


def example_comprehensive_reporting():
    """Example: Generate comprehensive quality report."""
    print("\n" + "="*60)
    print("Comprehensive Quality Reporting Example")
    print("="*60)
    
    # Initialize reporter
    reporter = QualityReporter(output_directory="./quality_reports")
    
    # Simulate some analysis results (in practice, these would come from actual analysis)
    # For demonstration, we'll create mock data
    
    mock_completeness_results = [
        {
            'status': 'complete',
            'completeness_score': 95.0,
            'issues': [],
            'directory': '/path/to/sim1'
        },
        {
            'status': 'incomplete',
            'completeness_score': 45.0,
            'issues': ['missing_structure_file', 'simulation_not_completed'],
            'directory': '/path/to/sim2'
        },
        {
            'status': 'partial',
            'completeness_score': 75.0,
            'issues': ['trajectory_too_small'],
            'directory': '/path/to/sim3'
        }
    ]
    
    mock_validation_results = [
        {
            'status': 'valid',
            'validation_score': 90.0,
            'issues': [],
            'chain_analysis': {'protein_chains': 5, 'chain_count_normal': True}
        },
        {
            'status': 'problematic',
            'validation_score': 30.0,
            'issues': ['unusual_chain_count', 'no_protein_found'],
            'chain_analysis': {'protein_chains': 2, 'chain_count_normal': False}
        },
        {
            'status': 'acceptable',
            'validation_score': 70.0,
            'issues': ['unusual_chain_count'],
            'chain_analysis': {'protein_chains': 4, 'chain_count_normal': False}
        }
    ]
    
    mock_batch_data = {
        'total_batches': 3,
        'total_directories': 50,
        'unique_pdbs': 45,
        'cross_batch_duplicates': 5,
        'identification_rate': 90.0
    }
    
    print("Generating comprehensive quality report...")
    
    # Generate comprehensive report
    report = reporter.generate_comprehensive_report(
        completeness_results=mock_completeness_results,
        validation_results=mock_validation_results,
        batch_tracker_data=mock_batch_data,
        report_name="example_quality_report"
    )
    
    # Display key findings
    exec_summary = report['executive_summary']
    print(f"\nExecutive Summary:")
    print(f"  Overall Health Score: {exec_summary['overall_health_score']:.1f}/100")
    print(f"  Health Status: {exec_summary['health_status']}")
    print(f"  Success Rate: {exec_summary['successful_completion_rate']:.1f}%")
    print(f"  Validation Rate: {exec_summary['structure_validation_rate']:.1f}%")
    
    # Key recommendations
    recommendations = report['recommendations']
    if recommendations:
        print(f"\nTop Recommendations:")
        for i, rec in enumerate(recommendations[:3], 1):
            print(f"  {i}. [{rec['priority'].upper()}] {rec['title']}")
            print(f"     {rec['description']}")
    
    # Top issues
    issues = report['prioritized_issues']
    if issues:
        print(f"\nTop Issues:")
        for i, issue in enumerate(issues[:3], 1):
            print(f"  {i}. {issue['issue']} (Frequency: {issue['frequency']}, Severity: {issue['severity']})")
            print(f"     {issue['description']}")
    
    print(f"\nReport saved to: quality_reports/")
    
    return report


def main():
    """Run all quality analysis examples."""
    print("Immunex Quality Analysis Module Examples")
    print("=" * 60)
    
    try:
        # Run examples
        completeness_results = example_md_completeness_check()
        validation_results = example_structure_validation()
        batch_results, batch_stats = example_batch_tracking()
        comprehensive_report = example_comprehensive_reporting()
        
        print("\n" + "="*60)
        print("All examples completed successfully!")
        print("="*60)
        
        print("\nNext steps:")
        print("1. Adapt the file paths to your actual MD simulation directories")
        print("2. Adjust parameters (chain counts, file sizes, etc.) to match your system")
        print("3. Set up automated quality monitoring for your MD pipeline")
        print("4. Review generated reports and implement recommended improvements")
        
    except Exception as e:
        logger.error(f"Error running examples: {e}")
        print(f"\nNote: This example uses placeholder paths.")
        print(f"Replace '/path/to/...' with actual paths to run the analysis.")


if __name__ == "__main__":
    main()