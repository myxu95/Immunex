#!/bin/bash
# Wrapper script to run batch CDR RMSD analysis with correct environment

# Set PATH to include aftermd conda environment
export PATH="/home/xumy/work/miniconda3/envs/aftermd/bin:$PATH"

# Set Python interpreter
PYTHON="/home/xumy/work/miniconda3/envs/aftermd/bin/python"

# Change to working directory
cd /home/xumy/work/development/AfterMD

# Run the analysis
$PYTHON scripts/batch_cdr_rmsd_smart.py "$@"
