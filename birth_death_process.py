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

# Set random seed
seed = 123
random.seed(seed)

# Define output folder
run_date = "09272023"
output_dir = os.path.join(os.getcwd(), "results", "birth_death_process", "coin_flip", run_date, "10000n0")
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

    # Add differentation method
    def differentiate(self):
        self.num_stem -= 1
        self.num_diff += 1

#-------------------------------------------------------------------------------
# Define simulation step
def birth_death_simulation(clone, tstop, birth_rate, death_rate):

    # Calculate birth probability
    birth_prob = birth_rate / (birth_rate + death_rate)

    # Add a time step
    time_step = np.random.exponential(scale=1/(clone.num_stem*(birth_rate+death_rate)), size=None)
    clone.time = clone.time + time_step

    # Draw a time step
    while clone.time < tstop:
        # Draw random number from 0 to 1
        random_draw = random.random()

        # Self-renew or differentiate depending on results
        if random_draw <= birth_prob:
            clone.self_renew()
        else:
            clone.differentiate()

        if clone.num_stem == 0:
            return clone

        # Add another time step
        time_step = np.random.exponential(scale=1/(clone.num_stem*(birth_rate+death_rate)), size=None)
        clone.time = clone.time + time_step

    # Return the resulting clone
    return clone

#-------------------------------------------------------------------------------
# Define parameters and run simulation

# Define starting parameters
n0 = 10000
t0 = 0
tstop = 100
birth_rate = 0.05 # right now just arbitrary, must be equal to death rate
death_rate = 0.05

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
# Produce cumulative proportion over scaled clone size step plot

# Define function to generate cumulative frequency vs rescaled average
def cumulative_frequency_plot(expanded_clones, tag, active_only = False):
    # Retain only clones with at least one stem cells (if active only)
    if active_only == True:
        expanded_clones = [clone for clone in expanded_clones if clone.num_stem > 0]
        active_tag = "_onlyActiveClones"
    else:
        active_tag = ""

    # Grab num_stem, num_diff, and total cell number from all clones
    num_stem_list = [clone.num_stem for clone in expanded_clones]
    num_diff_list = [clone.num_diff for clone in expanded_clones]
    num_cells_list = list(map(add, num_stem_list, num_diff_list))

    # Grab either all cells or only stem cells
    if tag == "allcells":
        cells_list = num_cells_list
    elif tag == "stemcells":
        cells_list = num_stem_list

    # Mean-scale number of cells
    mean_scaled_num_cells = [num_cells / mean(cells_list) for num_cells in cells_list]

    # Sort values
    mean_scaled_num_cells.sort()

    # Define steps
    steps = np.arange(len(mean_scaled_num_cells))[::-1] / len(mean_scaled_num_cells)

    # Plot cumulative frequency (simulation)
    plt.step(mean_scaled_num_cells, steps, "red", label = "Simulation")

    # Plot exponential decay (theoretical)
    exp_y = [np.exp(-x) for x in mean_scaled_num_cells]
    plt.plot(mean_scaled_num_cells, exp_y, "gray", label = "Theoretical: exp(-x)")

    # Add formatting
    fontsize = 16
    plt.xlabel("Rescaled average, n/<n>", fontsize = fontsize)
    plt.ylabel("Cumulative frequency (%)", fontsize = fontsize)
    plt.title("birth-death model: " + tag + active_tag, fontsize = fontsize)
    plt.legend(loc="upper right", fontsize = fontsize)
    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    plt.margins(x = 0.01, y = 0.01)

    # Save plot
    plt.savefig(output_dir + "/birth_death_process_cumulative_frequency_" + tag + active_tag + ".png", bbox_inches='tight', dpi=300)
    plt.close()

    return None

# Generate cumulative frequency plot
cumulative_frequency_plot(expanded_clones, "allcells")
cumulative_frequency_plot(expanded_clones, "allcells", active_only=True)

#-------------------------------------------------------------------------------
# Produce frequency over clone size histogram

# Define function to generate histogram of clone size frequency
def frequency_histogram(expanded_clones, tag):
    # Grab num_stem, num_diff, and total cell number from all clones
    num_stem_list = [clone.num_stem for clone in expanded_clones]
    num_diff_list = [clone.num_diff for clone in expanded_clones]
    num_cells_list = list(map(add, num_stem_list, num_diff_list))

    # Grab either all cells or only stem cells
    if tag == "allcells":
        cells_list = num_cells_list
    elif tag == "stemcells":
        cells_list = num_stem_list

    # Generate histogram
    plt.hist(cells_list, bins = max(cells_list),
             color = "red", alpha = 0.5, edgecolor = "black")

    # Add formatting
    fontsize = 16
    plt.xlabel("Clone size, n", fontsize = fontsize)
    plt.ylabel("Frequency", fontsize = fontsize)
    plt.title("birth-death model: " + tag, fontsize = fontsize)
    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    plt.margins(x = 0.01, y = 0.01)

    # Save plot
    plt.savefig(output_dir + "/birth_death_process_frequency_histogram_" + tag + ".png", bbox_inches='tight', dpi=300)
    plt.close()

    return None

# Generate frequency plot
frequency_histogram(expanded_clones, "allcells")
frequency_histogram(expanded_clones, "stemcells")