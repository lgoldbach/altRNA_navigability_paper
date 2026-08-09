#!/bin/bash
#SBATCH --job-name=jupyter-notebook # Job name
#SBATCH --output=jupyter_%j.log     # Log file
#SBATCH --time=0-02:00:00           # Maximum run time
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=4           # Number of threads
#SBATCH --mem=20G                    # RAM size

micromamba activate evo_nav

# Set up a port for Jupyter
PORT="8888"

NDIR="./"

# Run Jupyter Notebook
jupyter-notebook --no-browser --ip=0.0.0.0 --port=$PORT --notebook-dir=$NDIR
