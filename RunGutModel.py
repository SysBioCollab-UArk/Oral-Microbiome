from GutModel import model
from ReadGutlogo import read_Gutlogo
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import os
from scroll import ScrollableWindow
import csv
import glob
import warnings

# todo: figure out why our simulations always end with zero bacteria while most GutLogo simulations end
#  with an equilibrium level of bacteria

file_dir = 'GutLogo Simulations'
file_prefix = 'GutLogo plots_'  # 'GutLogo Populations'
file_suffix = 'csv'
files=glob.glob(os.path.join(file_dir, file_prefix + "*"))
n_files = len(files)
file_names = []

fig, axs = plt.subplots(nrows=n_files, ncols=2, sharex='all', sharey='all', figsize=[6.4, 4.8*n_files/2]) # [6.4, 4.8]
legend_created = False  # Track if the legend has been created

all_settings = {}  # dictionary to store the settings for all files
sim = ScipyOdeSimulator(model, verbose=False)
run_which = 'all' # [24, 34, 44]

n=0
for i in range(n_files):  # todo: figure out why the simulation stops at 6-->7
    file_name = '%s/%s%d.%s' % (file_dir, file_prefix, n, file_suffix)
    while not os.path.exists(file_name):
        n += 1
        file_name = '%s/%s%d.%s' % (file_dir, file_prefix, n, file_suffix)
    #####
    print(file_name)
    #####
    file_names.append(file_name)
    # read time courses and settings from gutlogo
    sim_steps, cell_counts, settings = read_Gutlogo(file_name)

    # ##### Plot GutLogo time courses #####
    for species in sorted(cell_counts.keys()):
        axs[i][0].plot(sim_steps, cell_counts[species], label=species)
    if i == n_files - 1:
        axs[i][0].set_xlabel('step')
    axs[i][0].set_ylabel('# of agents')
    axs[i][0].set_title("%s%d" % (file_prefix, n))
    axs[i][0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    # ##### Run PySB model with GutLogo settings #####

    # remove extraneous system settings
    param_names = [p.name for p in model.parameters]
    for key in [k for k in settings.keys()]:
        if key not in param_names:
            # print('discarding', key)
            warnings.warn("Discarding GutLogo parameter '%s'" % key)  # TODO: figure out why this isn't printing to the screen
            settings.pop(key)

    # save current settings in all_settings dictionary to output to csv file later
    for key in settings.keys():
        if i == 0:
            all_settings[key] = [settings[key]]
        else:
            all_settings[key] += [settings[key]]

    # run simulation
    if run_which == 'all' or n in run_which:
        t_span = np.array(sim_steps) * 60  # 1 sim_step = 60 seconds
        result = sim.run(tspan=t_span, param_values=settings)

        # plot pysb time courses
        for obs in sorted(['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']):
            axs[i][1].plot(t_span / 60, result.observables[obs], label=obs)
            if i == n_files - 1:
                axs[i][1].set_xlabel('time (min)')
            axs[i][1].set_ylabel('# of cells')
            axs[i][1].set_title("%s%d (PySB)" % (file_prefix, n))

    n += 1  # update file name index

# Create a subplot for the legend (invisible)
legend_ax = fig.add_subplot(111, frameon=False)
legend_ax.axis('off')  # Turn off the axis for the legend subplot

# Get handles and labels from the last plotted subplot (axs[-1][0])
handles, labels = axs[-1][0].get_legend_handles_labels()

# Create the legend in the "invisible" subplot
legend_ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0, 1.02), ncol=4)

# Adjust the layout to incorporate the legend
plt.tight_layout()

# Set legend_created to True so that it won't be created again
legend_created = True

# Create a scrollable window for the figure
a = ScrollableWindow(fig)
a.fig.tight_layout()
a.fig.savefig("GutModel.pdf", format="pdf")

# Output settings from all runs into a csv file
with open('GutLogoSettings.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(['parameter'] + [os.path.splitext(os.path.basename(file))[0] for file in file_names])
    for key in all_settings.keys():
        writer.writerow([key] + all_settings[key])
