#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name parameter_screen
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 64GB
#SBATCH --time 4:00:00
#SBATCH --output /projects/p31666/nat/clonal_expansion/logs/probability_curve_parameter_screen/%j.log
#SBATCH --verbose

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load matlab/r2022b

# Run data generation
# matlab -nodisplay -nosplash -nodesktop -r "run('/projects/p31666/nat/clonal_expansion/code/hierarchal_neutral_competition/probability_curve_parameter_screen/data_generation.m'); exit;"

# Run visualization
matlab -nodisplay -nosplash -nodesktop -r "run('/projects/p31666/nat/clonal_expansion/code/hierarchal_neutral_competition/probability_curve_parameter_screen/visualization.m'); exit;"