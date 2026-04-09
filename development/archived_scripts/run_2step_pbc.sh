#!/bin/bash
#SBATCH -J pbc_2step_1000f
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o /public/home/xmy/dataset/workspace/pbc_1000frames_2step/logs/pbc_2step_%j.out
#SBATCH -e /public/home/xmy/dataset/workspace/pbc_1000frames_2step/logs/pbc_2step_%j.err

# 2-Step PBC Processing for 1000frames Trajectories
# SLURM批处理脚本

echo "=========================================="
echo "2-Step PBC Processing for 1000frames"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# 激活conda环境
source /public/home/xmy/app/miniconda/bin/activate immunex

# 切换到工作目录
cd /public/home/xmy/dataset/workspace/pbc_1000frames_2step

# 运行Python脚本
python3 process_2step_pbc.py

echo ""
echo "=========================================="
echo "Processing completed"
echo "End time: $(date)"
echo "=========================================="
