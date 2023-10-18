# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----          Clonal expansion in 2D adult cycling tissues              -----
# -----                                                                    -----
# -----                          Goyal Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 10-06-2023
# Written by: BingKan Xue (writer of UF reference), Natalie Piehl
# Summary: Generate a simple birth-death process simulation
#
# Copied with slight modifications from University of Florida
# reference: https://cmp.phys.ufl.edu/PHZ4710/files/unit3/birth-death.html
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

# Set random seed
seed = 123
random.seed(seed)

# Define output folder
run_date = "10062023"
output_dir = os.path.join(os.getcwd(), "results", "birth_death_process", "UF_version", run_date)
os.makedirs(output_dir, exist_ok=True)

#-------------------------------------------------------------------------------
# Define simulation class

class BirthDeath:
    """"
    Simulate the birth-death process using gillespie algorithm
    """
    def __init__(self, birth_rate, death_rate, N0 = 1):
        """
        initialize the population.
        inputs:
        birth_rate: float, birth rate per individual.
        death_rate: float, death rate per individual.
        N0: int, initial population size.
        """
        self.birth_rate = birth_rate
        self.death_rate = death_rate
        self.cells = N0
        # set initial time to always be 0
        self.time = 0.
        # alternative to always grabbing the last element in the cell num list
        # is to keep the cell num history as a separate variable
        self.cells_hist = [N0] 
        self.time_hist = [0.]

    def next_event(self):
        """
        generates both the waiting time and identity of next event

        This is an alternate approach to the method I did previously where we
        pulled a random number from 0 to 1 (uniform) and classified any results
        where this random number is less than probability of birth as birth,
        otherwise as death.
        outputs:
        tau: float, waiting time before next event
        event: int, 0 means birth and 1 means death
        """
        # Set total birth and death rates
        k_b = self.cells * self.birth_rate
        k_d = self.cells * self.death_rate

        # Draw a random number from an exponential distribution of each total rate
        # (this is the heart of the Gillespie algorithm, I believe)
        t_b = np.random.exponential(1/k_b)
        t_d = np.random.exponential(1/k_d)

        # Check which event occurs firts and set that as the identity
        if t_b < t_d:
            event = 0
            return t_b, event
        else:
            event = 1
            return t_d, event
        
    def run(self, T):
        """
        Run simulation until time T since the beginning.
        inputs:
        T: float, time since the beginning of the simulation.
        """
        while self.time < T:
            # Check if have at least one cell and exit simulation if not
            if self.cells == 0:
                break

            # Draw next event
            tau, event = self.next_event()

            # Upate time
            self.time += tau

            # Change cell number based on event
            if event == 0:
                self.cells += 1
            elif event == 1:
                self.cells -= 1
            
            # Record time and cell number history
            self.time_hist.append(self.time)
            self.cells_hist.append(self.cells)        


#-------------------------------------------------------------------------------
# Run the simulation

# Define parameters
beta = 2 # starting with equal birth and death rates
delta = 1
T = 10. 
trials = 100
bd_list = [] # list to save all simulations

# Run simulation for each trial
for i in range(trials):
    # Create the simulation
    bd1 = BirthDeath(beta, delta)

    # Run the simulation
    bd1.run(T)

    # Print status update
    print(f'current time = {bd1.time}, current population size = {bd1.cells}')

    # Save result to a list
    bd_list.append(bd1)

#-------------------------------------------------------------------------------
# Plot clone size over time

plt.figure()
for bd1 in bd_list:
    plt.plot(bd1.time_hist, bd1.cells_hist, drawstyle='steps-post')    # stochastic realizations
# plt.plot(t_array, N_analytic, 'k', linewidth=2, label='deterministic')    # deterministic solution from above
plt.yscale('log')
plt.xlabel('time')
plt.ylabel('clone size')
# plt.legend(loc='upper left')
# plt.show()

 # Save plot
plt.savefig(output_dir + "/birth_death_process_clone_size_over_time_birth2death1.png", bbox_inches='tight', dpi=300)
plt.close()

#-------------------------------------------------------------------------------
# Run the simulation

# Define parameters
beta = 1 # starting with equal birth and death rates
delta = 1
T = 10. 
trials = 100
bd_list = [] # list to save all simulations

# Run simulation for each trial
for i in range(trials):
    # Create the simulation
    bd1 = BirthDeath(beta, delta)

    # Run the simulation
    bd1.run(T)

    # Print status update
    print(f'current time = {bd1.time}, current population size = {bd1.cells}')

    # Save result to a list
    bd_list.append(bd1)

#-------------------------------------------------------------------------------
# Plot clone size over time

plt.figure()
for bd1 in bd_list:
    plt.plot(bd1.time_hist, bd1.cells_hist, drawstyle='steps-post')    # stochastic realizations
# plt.plot(t_array, N_analytic, 'k', linewidth=2, label='deterministic')    # deterministic solution from above
plt.yscale('log')
plt.xlabel('time')
plt.ylabel('clone size')
# plt.legend(loc='upper left')
# plt.show()

 # Save plot
plt.savefig(output_dir + "/birth_death_process_clone_size_over_time_birth1death1.png", bbox_inches='tight', dpi=300)
plt.close()

#-------------------------------------------------------------------------------
# Plot cumulative distribution function

# Grab final cell number at end of simulation
num_cells = [clone.cells for clone in bd_list]

# Sort values
num_cells.sort()

# Define steps
steps = np.arange(len(num_cells))[::-1] / len(num_cells)

# Plot cumulative frequency (simulation)
plt.step(num_cells, steps, "red", label = "Simulation")

# Plot exponential decay (theoretical)
exp_y = [np.exp(-x) for x in num_cells]
plt.plot(num_cells, exp_y, "gray", label = "Theoretical: exp(-x)")

# Add formatting
fontsize = 16
plt.xlabel("Clone size", fontsize = fontsize)
plt.ylabel("Cumulative frequency (%)", fontsize = fontsize)
plt.title("birth-death model", fontsize = fontsize)
plt.legend(loc="upper right", fontsize = fontsize)
plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
plt.margins(x = 0.01, y = 0.01)

# Save plot
plt.savefig(output_dir + "/birth_death_process_birth1death1_cumulative_distribution.png", bbox_inches='tight', dpi=300)
plt.close()