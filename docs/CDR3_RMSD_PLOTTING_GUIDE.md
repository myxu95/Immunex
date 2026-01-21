# CDR3 RMSD Plotting Guide

This guide demonstrates how to use the enhanced PlotManager class for CDR3 RMSD visualization.

## Overview

The PlotManager now includes specialized methods for CDR3 RMSD analysis visualization:

1. **plot_cdr3_rmsd_summary()** - Comprehensive 4-panel summary
2. **plot_cdr3_stability_distribution()** - Stability category analysis
3. **plot_cdr3_rmsd_comparison()** - Multi-task trajectory comparison
4. **plot_cdr3_rmsd_heatmap()** - Statistical heatmap across tasks

## Quick Start

```python
from aftermd.utils.plotting import PlotManager
from pathlib import Path

# Initialize plotter
plotter = PlotManager(style='seaborn-v0_8-darkgrid')

# Define paths
summary_csv = "output/cdr3_analysis/cdr3_beta_rmsd_summary.csv"
trajectory_dir = "input/pbc_1000frames_2step"
output_dir = "output/cdr3_analysis/plots"
```

## Method Details

### 1. Comprehensive Summary Visualization

Generate a 4-panel figure with distribution, statistics, correlations, and top tasks:

```python
output_paths = plotter.plot_cdr3_rmsd_summary(
    summary_csv=summary_csv,
    output_dir=output_dir,
    show_plots=True  # Set False for batch processing
)
```

**Output**: 4-panel figure containing:
- **Panel 1**: RMSD distribution histogram with mean/median lines
- **Panel 2**: Box plots comparing Mean, Std, and Max RMSD
- **Panel 3**: Mean vs Std scatter plot (colored by Max RMSD)
- **Panel 4**: Top 20 highest RMSD tasks bar chart

### 2. Stability Distribution Analysis

Generate pie chart and bar chart showing stability categories:

```python
stability_paths = plotter.plot_cdr3_stability_distribution(
    summary_csv=summary_csv,
    output_dir=output_dir,
    show_plots=True
)
```

**Stability Categories**:
- Very Stable: < 0.10 nm
- Stable: 0.10-0.15 nm
- Moderate: 0.15-0.20 nm
- Flexible: 0.20-0.25 nm
- Highly Flexible: > 0.25 nm

### 3. Multi-Task RMSD Comparison

Compare RMSD trajectories for multiple tasks:

```python
# Compare stable tasks
stable_tasks = ['7nmg_1', '6am5_1', '2p5w_run1']
plotter.plot_cdr3_rmsd_comparison(
    summary_csv=summary_csv,
    task_names=stable_tasks,
    trajectory_dir=trajectory_dir,
    output_path="output/cdr3_rmsd_comparison_stable.png",
    show_stats=True  # Show mean RMSD lines
)

# Compare flexible tasks
flexible_tasks = ['4mji_1', '5c07_1', '5wkh_1']
plotter.plot_cdr3_rmsd_comparison(
    summary_csv=summary_csv,
    task_names=flexible_tasks,
    trajectory_dir=trajectory_dir,
    output_path="output/cdr3_rmsd_comparison_flexible.png",
    show_stats=True
)
```

### 4. Statistical Heatmap

Generate heatmap showing all RMSD statistics across tasks:

```python
# Top 30 most flexible tasks
plotter.plot_cdr3_rmsd_heatmap(
    summary_csv=summary_csv,
    output_path="output/cdr3_rmsd_heatmap_top30.png",
    top_n=30,
    sort_by='rmsd_mean_nm'  # Sort by mean RMSD
)

# All tasks
plotter.plot_cdr3_rmsd_heatmap(
    summary_csv=summary_csv,
    output_path="output/cdr3_rmsd_heatmap_all.png",
    top_n=None  # Show all tasks
)
```

**Heatmap includes**:
- Mean RMSD
- Std RMSD
- Min RMSD
- Max RMSD
- Median RMSD
- Range RMSD

### 5. Individual Task Analysis

For detailed analysis of single tasks:

```python
# Standard RMSD plot with moving average
plotter.plot_rmsd(
    rmsd_file="input/pbc_1000frames_2step/7nmg_1/cdr3_beta_rmsd.xvg",
    title="CDR3β RMSD - 7nmg_1",
    output_path="output/cdr3_rmsd_7nmg_1.png",
    show_stats=True,
    moving_average=100  # 100-frame moving average
)

# Convergence analysis
plotter.plot_rmsd_convergence(
    rmsd_file="input/pbc_1000frames_2step/7nmg_1/cdr3_beta_rmsd.xvg",
    window_sizes=[50, 100, 200, 500],
    title="CDR3β RMSD Convergence - 7nmg_1",
    output_path="output/cdr3_rmsd_convergence_7nmg_1.png"
)
```

## Complete Example Script

See `examples/plot_cdr3_rmsd_example.py` for a complete working example that generates all visualization types.

To run the example:
```bash
python examples/plot_cdr3_rmsd_example.py
```

## Output Files

All plots are saved as high-resolution PNG files (300 DPI) suitable for publication.

### Generated Files
- `cdr3_beta_rmsd_summary.png` - 4-panel comprehensive summary
- `cdr3_beta_stability_categories.png` - Stability distribution
- `cdr3_rmsd_comparison_stable.png` - Stable tasks comparison
- `cdr3_rmsd_comparison_flexible.png` - Flexible tasks comparison
- `cdr3_rmsd_heatmap_top30.png` - Statistical heatmap
- `cdr3_rmsd_{task}.png` - Individual task RMSD
- `cdr3_rmsd_convergence_{task}.png` - Convergence analysis

## Customization Options

### Style and Appearance

```python
# Initialize with different style
plotter = PlotManager(
    style='seaborn-v0_8-darkgrid',  # or 'ggplot', 'bmh', 'classic'
    figsize=(12, 8)  # Default figure size
)
```

### Color Palettes

PlotManager uses seaborn's "husl" palette by default. Colors are automatically assigned.

### Interactive Plots

For RMSD time series, you can generate interactive Plotly plots:

```python
plotter.plot_rmsd(
    rmsd_file="path/to/rmsd.xvg",
    interactive=True,  # Generate interactive HTML plot
    output_path="output/rmsd_interactive.html"
)
```

## Tips and Best Practices

1. **Batch Processing**: Set `show_plots=False` when generating multiple plots
2. **Memory Management**: Process tasks in batches for large datasets
3. **File Organization**: Use consistent output directories for easy access
4. **Publication Quality**: Default 300 DPI is suitable for most journals
5. **Comparison Plots**: Limit to 5-7 tasks per comparison for clarity

## Integration with Analysis Pipeline

```python
from aftermd.utils.plotting import PlotManager
from aftermd.analysis.cdr_manager import CDRManager

# Step 1: Run CDR3 recognition
cdr_manager = CDRManager()
batch_results = cdr_manager.batch_process_structures(pdb_dir, output_dir)

# Step 2: Calculate RMSD (see batch_cdr3_rmsd.py)
# ...

# Step 3: Visualize results
plotter = PlotManager()
plotter.plot_cdr3_rmsd_summary(summary_csv, output_dir)
plotter.plot_cdr3_stability_distribution(summary_csv, output_dir)
```

## Troubleshooting

### Common Issues

**Issue**: "ModuleNotFoundError: No module named 'aftermd'"
**Solution**: Install AfterMD in development mode:
```bash
pip install -e .
```

**Issue**: "File not found" errors
**Solution**: Ensure all paths are absolute or relative to working directory

**Issue**: Plots not displaying
**Solution**: Use `show_plots=True` or check if running in headless environment

**Issue**: Poor plot quality
**Solution**: Increase DPI in savefig or adjust figure size

## Further Information

- PlotManager API: See `aftermd/utils/plotting.py`
- RMSD Calculation: See `scripts/batch_cdr3_rmsd.py`
- CDR Recognition: See `aftermd/analysis/cdr_manager.py`
