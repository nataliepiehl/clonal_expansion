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
import random
import os
from statistics import mean
import matplotlib.pyplot as plt

# Set random seed
seed = 123
random.seed(seed)

# Define output folder
output_dir = os.path.join(os.getcwd(), "results", "birth_death_process", "seed" + str(seed))
os.makedirs(output_dir, exist_ok=True)

#-------------------------------------------------------------------------------
# Define clone class
class Clone:
    # Add cell number attributes
    def __init__(self, num_stem, num_diff, time):
        self.num_stem = num_stem
        self.num_diff = num_diff
        self.time = time

    # Add self-renewal method
    def self_renew(self):
        self.num_stem += 1

    # Add differnentation method
    def differentiate(self):
        self.num_stem -= 1
        self.num_diff += 1

#-------------------------------------------------------------------------------
# Define simulation step
def birth_death_simulation(clone, tstop, birth_rate, death_rate):

    # Calculate birth probability
    birth_prob = birth_rate / (birth_rate + death_rate)

    # Define mean event rate
    mean_event_rate = clone.num_stem*birth_rate*death_rate

    # Add a time step
    time_step = np.random.exponential(scale=mean_event_rate, size=None)
    clone.time = clone.time + time_step

    # Draw a time step
    while clone.time < tstop and clone.num_stem > 0:
        # Draw random number from 0 to 1
        random_draw = random.random()

        # Self-renew or differentiate depending on results
        if random_draw <= birth_prob:
            clone.self_renew()
        else:
            clone.differentiate()

        # Add another time step
        time_step = np.random.exponential(scale=mean_event_rate, size=None)
        clone.time = clone.time + time_step

    # Return the resulting clone
    return clone

#-------------------------------------------------------------------------------
# Define parameters and run simulation

# Define starting parameters
n0 = 2000
t0 = 0
tstop = 100
birth_rate = 1 # right now just arbitrary, must be equal to death rate
death_rate = 1

# Define simulation runner
def run_simulation(n0, t0, tstop, birth_rate, death_rate):
    # Initialize clones
    clones = [Clone(1,0,t0) for i in range(n0)]

    # Run simulation on each clone
    expanded_clones = list()
    for clone in clones:
        expanded_clones.append(birth_death_simulation(clone, tstop, birth_rate, death_rate))

    return expanded_clones

# Run simulation
expanded_clones = run_simulation(n0, t0, tstop, birth_rate, death_rate)

#-------------------------------------------------------------------------------
# Produce cumulative proportion over scaled clone size plot

# Grab num_stem, num_diff, and total cell number from all clones
num_stem_list = [clone.num_stem for clone in expanded_clones]
num_diff_list = [clone.num_diff for clone in expanded_clones]
num_cells_list = num_stem_list + num_diff_list

# Mean-scale number of cells
mean_scaled_num_cells = [num_cells / mean(num_cells_list) for num_cells in num_cells_list]

# Sort values
mean_scaled_num_cells.sort()

# Define steps
steps = np.arange(len(mean_scaled_num_cells))[::-1] / len(mean_scaled_num_cells)

# Plot cumulative frequency
plt.step(mean_scaled_num_cells, steps)
plt.xlabel("Rescaled average, n/<n(t)>")
plt.ylabel("Cumulative frequency (%)")
plt.title("birth-death model")
plt.margins(x = 0.01, y = 0.01)
plt.rc('font', size=14)

# Save plot
plt.savefig(output_dir + "/birth_death_process_seed123.png", bbox_inches='tight', dpi=300)
plt.close()

# Save parameters?