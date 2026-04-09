#!/usr/bin/env python3
"""
Example script demonstrating CDR3 RMSD plotting functionality.
"""

from immunex.utils.plotting import PlotManager
from pathlib import Path

def main():
    """Demonstrate CDR3 RMSD plotting with PlotManager."""

    # Initialize PlotManager
    plotter = PlotManager(style='seaborn-v0_8-darkgrid', figsize=(12, 8))

    # Define paths
    base_dir = Path("/home/xumy/work/development/Immunex")
    summary_csv = base_dir / "output/cdr3_analysis/cdr3_beta_rmsd_summary.csv"
    trajectory_dir = base_dir / "input/pbc_1000frames_2step"
    output_dir = base_dir / "output/cdr3_analysis/plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("CDR3 RMSD Plotting Examples")
    print("=" * 80)

    # Example 1: Generate comprehensive summary plots
    print("\n1. Generating CDR3 RMSD summary visualization...")
    summary_paths = plotter.plot_cdr3_rmsd_summary(
        summary_csv=str(summary_csv),
        output_dir=str(output_dir),
        show_plots=False
    )
    print(f"   Summary plot saved to: {summary_paths['summary']}")

    # Example 2: Generate stability distribution plots
    print("\n2. Generating stability distribution visualization...")
    stability_paths = plotter.plot_cdr3_stability_distribution(
        summary_csv=str(summary_csv),
        output_dir=str(output_dir),
        show_plots=False
    )
    print(f"   Stability plot saved to: {stability_paths['stability']}")

    # Example 3: Compare RMSD trajectories for selected tasks
    print("\n3. Comparing RMSD trajectories for top stable and flexible tasks...")

    # Top 3 most stable tasks
    stable_tasks = ['7nmg_1', '6am5_1', '2p5w_run1']
    plotter.plot_cdr3_rmsd_comparison(
        summary_csv=str(summary_csv),
        task_names=stable_tasks,
        trajectory_dir=str(trajectory_dir),
        output_path=str(output_dir / "cdr3_rmsd_comparison_stable.png"),
        show_stats=True
    )

    # Top 3 most flexible tasks
    flexible_tasks = ['4mji_1', '5c07_1', '5wkh_1']
    plotter.plot_cdr3_rmsd_comparison(
        summary_csv=str(summary_csv),
        task_names=flexible_tasks,
        trajectory_dir=str(trajectory_dir),
        output_path=str(output_dir / "cdr3_rmsd_comparison_flexible.png"),
        show_stats=True
    )

    # Example 4: Generate heatmap for top 30 tasks
    print("\n4. Generating RMSD heatmap for top 30 flexible tasks...")
    plotter.plot_cdr3_rmsd_heatmap(
        summary_csv=str(summary_csv),
        output_path=str(output_dir / "cdr3_rmsd_heatmap_top30.png"),
        top_n=30,
        sort_by='rmsd_mean_nm'
    )

    # Example 5: Plot individual task RMSD with convergence analysis
    print("\n5. Generating detailed RMSD analysis for sample task...")
    sample_task = '7nmg_1'
    sample_rmsd_file = trajectory_dir / sample_task / "cdr3_beta_rmsd.xvg"

    if sample_rmsd_file.exists():
        # Standard RMSD plot
        plotter.plot_rmsd(
            rmsd_file=str(sample_rmsd_file),
            title=f"CDR3β RMSD - {sample_task}",
            output_path=str(output_dir / f"cdr3_rmsd_{sample_task}.png"),
            show_stats=True,
            moving_average=100
        )

        # Convergence analysis
        plotter.plot_rmsd_convergence(
            rmsd_file=str(sample_rmsd_file),
            window_sizes=[50, 100, 200, 500],
            title=f"CDR3β RMSD Convergence - {sample_task}",
            output_path=str(output_dir / f"cdr3_rmsd_convergence_{sample_task}.png")
        )

    print("\n" + "=" * 80)
    print("All plots generated successfully!")
    print(f"Output directory: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
