#!/bin/bash
# Immunex Environment Setup Script
# This script creates a complete conda environment for Immunex

set -e  # Exit on error

echo "========================================================================="
echo "                  Immunex Environment Setup"
echo "========================================================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found!"
    echo "Please install Anaconda or Miniconda first."
    exit 1
fi

echo "Conda found: $(which conda)"
echo "Conda version: $(conda --version)"
echo ""

# Check if immunex environment exists
if conda env list | grep -q "^immunex "; then
    echo "Warning: immunex environment already exists!"
    read -p "Do you want to remove and recreate it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing immunex environment..."
        conda env remove -n immunex -y
    else
        echo "Updating existing environment..."
        conda env update -f environment.yml --prune
        echo ""
        echo "========================================================================="
        echo "Environment updated successfully!"
        echo "Activate with: conda activate immunex"
        echo "========================================================================="
        exit 0
    fi
fi

echo "Creating new immunex environment..."
echo "This may take 10-20 minutes..."
echo ""

# Create environment
conda env create -f environment.yml

echo ""
echo "========================================================================="
echo "Environment created successfully!"
echo ""
echo "To activate the environment, run:"
echo "  conda activate immunex"
echo ""
echo "To test the installation, run:"
echo "  conda activate immunex"
echo "  python -c 'import immunex; print(immunex.__version__)'"
echo "  python development/validation/test_environment.py"
echo "========================================================================="
