#!/bin/bash
#SBATCH --job-name=<job_name>
#SBATCH -p quick
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
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