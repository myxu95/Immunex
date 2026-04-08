#!/bin/bash
#SBATCH --job-name=<job_name>
#SBATCH -p multi
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --gres=gpu:1

# Environment setup
export LD_LIBRARY_PATH=/public/software/lib/:$LD_LIBRARY_PATH
source /public/software/profile.d/apps_gromacs_2023.2.sh

# Activate conda environment
source /public/home/xmy/app/miniconda/bin/activate immunex

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

# Run Immunex batch processing
echo "Starting Immunex batch processing..."
{batch_command}

# Cleanup temporary files (optional)
if [ "$CLEANUP_TEMP_FILES" != "false" ]; then
    echo "Cleaning up temporary files..."
    find . -name "*.temp_*.xtc" -delete 2>/dev/null || true
    find . -name "*.temp_*.trr" -delete 2>/dev/null || true
    find . -name "#*.xtc.*#" -delete 2>/dev/null || true
    find . -name "#*.trr.*#" -delete 2>/dev/null || true
    find . -name "temp_*.ndx" -delete 2>/dev/null || true
    find . -name "*.backup.*" -delete 2>/dev/null || true
    echo "Cleanup completed"
fi

echo "End time: $(date)"
echo "Job completed successfully"