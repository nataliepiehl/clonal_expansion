#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name parameter_screen
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 120GB
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1042/GoyalLab/nat/clonal_expansion/logs/probability_curve_parameter_screen/%j.log
#SBATCH --verbose

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load matlab/r2022b

# Run job
matlab -batch /home/ncp2306/clonal_expansion/code/probability_curve_parameter_screen/data_generation.m