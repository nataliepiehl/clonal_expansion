# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----          Clonal expansion in 2D adult cycling tissues              -----
# -----                                                                    -----
# -----                          Goyal Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 09-28-2023
# Written by: Natalie Piehl
# Summary: Visualize results from simple birth-death process simulation
#
#-------------------------------------------------------------------------------
# Initialization

# Read in packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
from statistics import mean
from operator import add
import math
import seaborn as sns

# Set random seed
seed = 123
random.seed(seed)

#-------------------------------------------------------------------------------
# Produce frequency over clone size histogram

# Define function to generate histogram of clone size frequency
def frequency_histogram(expanded_clones, tag, output_dir):
    # Grab num_stem, num_diff, and total cell number from all clones
    num_stem_list = [clone.num_stem[-1] for clone in expanded_clones]
    num_diff_list = [clone.num_diff[-1] for clone in expanded_clones]
    num_cells_list = list(map(add, num_stem_list, num_diff_list))

    # Grab either all cells or only stem cells
    if tag == "allcells":
        cells_list = num_cells_list
    elif tag == "stemcells":
        cells_list = num_stem_list
    else:
        print("Unrecognizable tag, please use 'allcells' or 'stemcells'")
        return None

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

#-------------------------------------------------------------------------------
# Produce cumulative proportion over scaled clone size step plot

# Define function to generate cumulative frequency vs rescaled average
def cumulative_frequency_plot(expanded_clones, tag, output_dir, active_only = False):
    # Retain only clones with at least one stem cells (if active only)
    if active_only == True:
        expanded_clones = [clone for clone in expanded_clones if clone.num_stem[-1] > 0]
        active_tag = "_onlyActiveClones"
    else:
        active_tag = ""

    # Grab num_stem, num_diff, and total cell number from all clones
    num_stem_list = [clone.num_stem[-1] for clone in expanded_clones]
    num_diff_list = [clone.num_diff[-1] for clone in expanded_clones]
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

#-------------------------------------------------------------------------------
# Produce powers of 2 clone size distribution over time plot

# Define function to generate histogram of clone size frequency
def clone_size_over_time_plots(expanded_clones, output_dir):
    # Grab number of cells per clone for each timestep
    num_stem_list = [clone.num_stem for clone in expanded_clones]
    num_diff_list = [clone.num_diff for clone in expanded_clones]
    num_cells_list = [list(map(add, clone.num_stem, clone.num_diff)) for clone in expanded_clones]

    # Identify max timesteps taken by any clone
    timesteps = max([len(clone_timeline) for clone_timeline in num_cells_list])

    # Iterate over each clone timeline
    for i in range(len(num_cells_list)):
        # Append NAs to fill in to max timesteps
        for _ in range(timesteps - len(num_cells_list[i])):
            num_cells_list[i].append(np.nan)
            num_stem_list[i].append(np.nan)
            num_diff_list[i].append(np.nan)

    # Create dataframe
    cells_df = pd.DataFrame(num_cells_list)
    stem_df = pd.DataFrame(num_stem_list)
    diff_df = pd.DataFrame(num_diff_list)

    # Add column names
    column_names = {cells_df.columns.values[i] : "t"+str(i) for i in cells_df.columns.values}
    cells_df = cells_df.rename(columns = column_names)
    stem_df = stem_df.rename(columns = column_names)
    diff_df = diff_df.rename(columns = column_names)

    # Add cell type column
    cells_df['celltype'] = pd.Series(np.repeat("allcells", len(cells_df)))
    stem_df['celltype'] = pd.Series(np.repeat("stemcells", len(cells_df)))
    diff_df['celltype'] = pd.Series(np.repeat("diffcells", len(cells_df)))

    # Add clone ID column
    cells_df['cloneID'] = pd.Series(["clone"+str(i) for i in range(len(cells_df))])
    stem_df['cloneID'] = pd.Series(["clone"+str(i) for i in range(len(cells_df))])
    diff_df['cloneID'] = pd.Series(["clone"+str(i) for i in range(len(cells_df))])

    # Merge DFs together
    df = pd.concat([cells_df, stem_df, diff_df])

    # Move celltype and cloneID to front
    df = df[['celltype', 'cloneID'] + [col for col in df.columns if col not in ['celltype', 'cloneID']]]

    # Export table
    df.to_csv(output_dir + "/clone_trajectories_over_time.csv")

    # Isolate only all cells
    df_allonly = df[df['celltype'] == 'allcells']

    # Generate clone size bins per time point
    bins = [0, 2, 4, 8, 16, 32, 64, 128]
    timepoints = cells_df.columns.values[~cells_df.columns.isin(['celltype', 'cloneID'])]
    bin_columns = list()
    for timepoint in timepoints:
        bin_columns.append(pd.cut(df_allonly[timepoint], bins = bins))
    bin_frequencies = pd.concat(bin_columns, axis = 1)

    # Find proportions of each bin per timepoint
    bin_frequencies = bin_frequencies[timepoints].apply(pd.Series.value_counts)
    bin_frequencies = bin_frequencies.apply(lambda x: x / x.sum())

    # Plot proportions
    ax = bin_frequencies.T.plot()

    # Log-scale x axis
    plt.xscale("log")

    # Add formatting
    fontsize = 16
    plt.xlabel("Time step", fontsize = fontsize)
    plt.ylabel("Clones (%)", fontsize = fontsize)
    plt.title("birth-death model", fontsize = fontsize)
    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    plt.margins(x = 0.01, y = 0.01)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

    # Save plot
    plt.savefig(output_dir + "/birth_death_process_clone_bin_over_time.png", bbox_inches='tight', dpi=300)
    plt.close()
    
    #---------------------------------------------------------------------------
    # Average # cells over time

    # Produce average per celltype
    grouped = df.groupby('celltype')
    average_cells = grouped.mean()

    # Plot means
    ax = average_cells.T.plot()

    # Log-scale x axis
    plt.xscale("log")

    # Add formatting
    fontsize = 16
    plt.xlabel("Time step", fontsize = fontsize)
    plt.ylabel("Avg. cells per clone", fontsize = fontsize)
    plt.title("birth-death model", fontsize = fontsize)
    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    plt.margins(x = 0.01, y = 0.01)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

    # Save plot
    plt.savefig(output_dir + "/birth_death_process_cell_numbers_over_time.png", bbox_inches='tight', dpi=300)
    plt.close()

    #---------------------------------------------------------------------------
    # Scatter of clone size over timepoints

    # Isolate five equally spaced timepoints
    interval = math.floor(len(timepoints) / 5)
    timepoints_oi = ['t' + str(timepoint) for timepoint in [interval*i for i in range(1,5)]]
    timepoints_oi.append(timepoints[-1])

    # Isolate timepoints of interest
    timepoints_oi.append("cloneID")
    df_timepoints_oi = df_allonly.loc[:, timepoints_oi]

    # Pivot to long
    df_timepoints_oi_long = pd.wide_to_long(df_timepoints_oi, stubnames = 't', i = 'cloneID', j = 'timepoint').reset_index()
    df_timepoints_oi_long.columns = ['cloneID', 'timepoint', 'numcells']
    df_timepoints_oi_long['timepoint'] = ['t' + str(timepoint) for timepoint in df_timepoints_oi_long['timepoint']]

    # Remove rows with NAs (nonactive clones)
    df_timepoints_oi_long = df_timepoints_oi_long.dropna()

    # Generate swarmplot
    ax = sns.swarmplot(data=df_timepoints_oi_long, x="timepoint", y="numcells", s = 2)
        
    # Add formatting
    fontsize = 16
    plt.xlabel("Time step", fontsize = fontsize)
    plt.ylabel("Total cells per clone", fontsize = fontsize)
    plt.title("birth-death model", fontsize = fontsize)
    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    plt.margins(x = 0.01, y = 0.01)

    # Save plot
    plt.savefig(output_dir + "/birth_death_process_clone_size_over_time_scatter.png", bbox_inches='tight', dpi=300)
    plt.close()

    return None