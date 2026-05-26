#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=100:00:00      # Adjust the time limit as needed (e.g., 24:00:00 for 24 hours)
#SBATCH --mem=42G            # Request more memory if needed
#SBATCH --job-name=sandder_mm
#SBATCH --output=sandder_mm-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hao.wu@vanderbilt.edu

module purge
module load gcc/11.2.0
module load r/4.4.0
module load openmpi/4.1.4-intel_20.2

export LD_LIBRARY_PATH=/nas/longleaf/rhel8/apps/gcc/11.2.0/lib64:$LD_LIBRARY_PATH

Rscript sand_derivation_RUN.R
