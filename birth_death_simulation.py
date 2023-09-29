# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----          Clonal expansion in 2D adult cycling tissues              -----
# -----                                                                    -----
# -----                          Goyal Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 09-26-2023
# Written by: Natalie Piehl
# Summary: Generate a simple birth-death process simulation
#
# Basic process as described by https://doi.org/10.1119/PICUP.Exercise.BirthDeath
# 1. Calculate the mean rate at which events (any event) occur
# 2. Calculate the wait time until the next event (a random process)
# 3. Determine which event actually occurs (a second random process)
# 4. Adjust population and time parameters based on the outcomes of the event
# 5. Repeat
#
#-------------------------------------------------------------------------------
# Initialization

# Read in packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import os
from statistics import mean
from operator import add
# import dill

# Import visualization script
import birth_death_visualization as vis

# Set random seed
seed = 123
random.seed(seed)

# Define output folder
run_date = "09292023"
output_dir = os.path.join(os.getcwd(), "results", "birth_death_process", "coin_flip", run_date)
os.makedirs(output_dir, exist_ok=True)

#-------------------------------------------------------------------------------
# Define clone class
class Clone:
    # Add cell number attributes
    def __init__(self, num_stem, num_diff, time=0):
        self.num_stem = [num_stem]
        self.num_diff = [num_diff]
        self.time = [time]

    # Add time step update method
    def add_time_step(self, time_step):
        self.time.append(self.time[-1]+time_step)

    # Add self-renewal method
    def self_renew(self):
        self.num_stem.append(self.num_stem[-1] + 1)
        self.num_diff.append(self.num_diff[-1])

    # Add differentation method
    def differentiate(self):
        self.num_stem.append(self.num_stem[-1] - 1)
        self.num_diff.append(self.num_diff[-1] + 1)

#-------------------------------------------------------------------------------
# Define simulation step
def birth_death_simulation(clone, tstop, birth_rate, death_rate):

    # Calculate birth probability
    birth_prob = birth_rate / (birth_rate + death_rate)

    # Draw a time step
    while clone.time[-1] < tstop:
        # Add another time step
        time_step = np.random.exponential(scale=1/(clone.num_stem[-1]*(birth_rate+death_rate)), size=None)
        clone.add_time_step(time_step)

        # Draw random number from 0 to 1
        random_draw = random.random()

        # Self-renew or differentiate depending on results
        if random_draw <= birth_prob:
            clone.self_renew()
        else:
            clone.differentiate()

        if clone.num_stem[-1] == 0:
            return clone

    # Return the resulting clone
    return clone

#-------------------------------------------------------------------------------
# Define parameters and run simulation

# Define starting parameters
n0 = 10000
tstop = 100
birth_rate = 0.05 # right now just arbitrary, must be equal to death rate
death_rate = 0.05

# Define simulation runner
def run_simulation(n0, tstop, birth_rate, death_rate):
    # Initialize clones
    clones = [Clone(1,0) for _ in range(n0)]                                                                          

    # Run simulation on each clone
    expanded_clones = list()
    for clone in clones:
        expanded_clones.append(birth_death_simulation(clone, tstop, birth_rate, death_rate))

    return expanded_clones

# Run simulation
expanded_clones = run_simulation(n0, tstop, birth_rate, death_rate)

# # Testing
# clone = Clone(1,0)
# expanded_clone = birth_death_simulation(clone, tstop, birth_rate, death_rate)

# # Export clones object
# filename = output_dir + "/birth_death_expanded_clones.pkl"
# with open(filename, 'wb') as outp:
#     dill.dump(expanded_clones, outp)

#-------------------------------------------------------------------------------
# Run visualizations

# Generate frequency plots
vis.frequency_histogram(expanded_clones, "allcells", output_dir)
vis.frequency_histogram(expanded_clones, "stemcells", output_dir)

# Generate cumulative frequency plots
vis.cumulative_frequency_plot(expanded_clones, "allcells", output_dir)
vis.cumulative_frequency_plot(expanded_clones, "allcells", output_dir, active_only=True)

# Generate clone size group plot
vis.clone_size_over_time_plots(expanded_clones, output_dir)
