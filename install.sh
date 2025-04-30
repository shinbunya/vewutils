#!/bin/bash

# Check if environment name is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide a conda environment name"
    echo "Usage: ./install_conda.sh <environment_name>"
    exit 1
fi

ENV_NAME=$1

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed. Please install conda first."
    exit 1
fi

# Create conda environment with Python 3.11
echo "Creating conda environment '$ENV_NAME' with Python 3.11..."
conda create -n $ENV_NAME python=3.11 -y

# Activate the environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Install dependencies from environment.yml
echo "Installing dependencies from environment.yml..."
conda env update -f environment.yml

# Install unresolved dependencies for adcircpy
echo "Installing unresolved dependencies for adcircpy..."
pip install stormevents

# Install adcircpy from GitHub
echo "Installing adcircpy from a fork at GitHub..."
pip install --no-deps git+https://github.com/shinbunya/adcircpy.git

# Install erddapy
echo "Installing erddapy..."
pip install erddapy

# Install adcircutils in development mode without dependencies
echo "Installing adcircutils in development mode..."
pip install --no-deps -e .

echo "Installation complete! To activate the environment, run:"
echo "conda activate $ENV_NAME"
