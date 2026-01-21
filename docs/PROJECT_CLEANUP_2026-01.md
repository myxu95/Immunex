# Project Cleanup Report - 2026-01-20

## Summary
Cleaned up AfterMD project structure by removing redundant test scripts and consolidating batch processing scripts.

## Changes Made

### 1. Test Scripts Deletion (13 files)
All temporary test scripts removed:
- test_allostery_full.py
- test_allostery_quick.py
- test_batch_rmsd.py
- test_cdr_rmsd.py
- test_cdr_rmsf_curves.py
- test_clustered_heatmap.py
- test_contact_optimization.py
- test_discover_tasks.py
- test_mosaic_5tasks.py
- test_mosaic_fast.py
- test_mosaic_improved.py
- test_rmsd_hla_alpha.py
- test_single_mosaic.py

### 2. Deprecated Batch Scripts Deletion (4 files)
Consolidated multiple CDR RMSD versions, kept the latest (exact):
- batch_cdr_rmsd_phla_align.py (deleted)
- batch_cdr_rmsd_segid.py (deleted)
- batch_cdr_rmsd_smart.py (deleted)
- batch_cdr3_rmsd_simple.py (deleted)
- **batch_cdr_rmsd_exact.py** (retained - most recent)

### 3. Scripts Moved to development/ (28 files)

#### Case Studies (development/case_studies/)
- analyze_three_cases.py
- comparative_allostery_analysis.py

#### Utilities (development/utilities/)
- batch_rmsd_missing.py
- batch_process_raw_data_new.py
- calculate_rmsd.py
- compare_pdb_tasks.py
- export_index_chain_mapping.py
- extract_protein_chains.py
- find_representative_sample.py
- monitor_mosaic_batch.py
- move_unique_tasks.py
- regenerate_heatmaps.py
- rescue_tpr_files.py
- split_pdb_biopython.py

#### Plotting (development/plotting/)
- plot_cdr3_flexibility.py
- plot_cdr_interface_combined_boxplot.py
- plot_cdr_tcr_combined_boxplot.py
- plot_tcr_cdr3_comparison.py
- plot_tcr_overall_distribution.py

#### Clustering (development/clustering/)
- batch_clustering_analysis.py
- batch_mosaic_clustering.py
- cluster_contact_correlation.py
- extract_cluster_full_structures.py
- generate_bfactor_projection.py
- generate_pymol_cluster_visualization.py
- run_cdr3_clustering.py
- step1_cluster_static_angles.py

#### Visualization (development/visualization/)
- batch_pymol_visualization.py

#### Archived (development/archived_scripts/)
- run_comparative_analysis.py
- run_complete_rmsd_analysis.py
- run_full_cdr_rmsd_batch.py

## Final Structure

### scripts/ Directory (18 core scripts)
**Core Processing:**
- pbc_process.py - Single trajectory PBC processing
- batch_pbc_slurm.py - SLURM cluster batch processing
- batch_worker.py - Generic SLURM worker

**Quality Control:**
- md_quality_check.py - MD simulation quality checks
- md_task_organizer.py - Task organization
- md_workflow.py - Complete workflow orchestration

**Batch Analysis:**
- batch_allostery_analysis.py - Allostery correlation analysis
- batch_cdr_rmsd_exact.py - CDR RMSD calculation (consolidated)
- batch_cdr_rmsf.py - CDR RMSF analysis
- batch_contact_frequency.py - Contact frequency calculation
- batch_interface_rmsd.py - Interface RMSD analysis
- batch_phla_analysis.py - pHLA-specific analysis
- batch_rmsd_hla_alpha.py - HLA alpha chain RMSD
- batch_tcr_rmsd.py - TCR RMSD calculation
- batch_whole_protein_rmsf.py - Whole protein RMSF

**Analysis Pipelines:**
- analyze_chain_contacts.py - Chain contact analysis
- rmsf_analysis_pipeline.py - RMSF processing pipeline
- run_trajectory_analysis.py - Trajectory analysis runner

**Shell/TCL Scripts:**
- monitor_batch.sh - Batch job monitoring
- run_batch_cdr_with_env.sh - CDR batch with environment
- vmd_chain_contacts.tcl - VMD contact analysis

## Statistics
- **Before**: 66 Python scripts
- **After**: 18 Python scripts
- **Reduction**: 73% decrease
- **Deleted**: 17 files
- **Moved to development/**: 28 files

## Temporary Workspaces (Retained)
User chose to retain these for now:
- FEL_workspace/ (23GB)
- allosteric_workspace/ (1GB)

## Recommendations

### 1. Further Consolidation Opportunities
Consider merging similar batch analysis scripts:
- batch_cdr_rmsf.py + batch_whole_protein_rmsf.py → single RMSF tool
- batch_tcr_rmsd.py + batch_rmsd_hla_alpha.py → unified RMSD tool with --chain option

### 2. Documentation Updates Needed
- Update CLAUDE.md with new scripts/ structure
- Add scripts/ README explaining each core script
- Document development/ subdirectory organization

### 3. Future Cleanup Tasks
- Review development/archived_scripts/ for permanent deletion
- Clean up FEL_workspace/ and allosteric_workspace/ when analysis complete
- Convert useful development scripts to proper CLI tools in aftermd/

## Version Control
Remember to commit these changes:
```bash
git add -A
git commit -m "Major cleanup: Remove test scripts, consolidate batch tools, organize development scripts"
```
