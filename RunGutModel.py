from GutModel import model
from ReadGutlogo import read_Gutlogo
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import os
from scroll import ScrollableWindow

file_dir = 'GutLogo Simulations'
file_prefix = 'GutLogo Populations'
file_suffix = 'csv'
n_files = len(os.listdir(file_dir))

sim = ScipyOdeSimulator(model, verbose=True)

fig, axs = plt.subplots(nrows=n_files, ncols=2, sharex=False, sharey=True, figsize=[6.4, 4.8*n_files/2]) #6.4, 4.8 default
legend_created = False  # Track if the legend has been created

for i in range(n_files): #todo: figure out why the simulation stops at 6-->7
    file_name = '%s/%s%d.%s' % (file_dir, file_prefix, i, file_suffix)
    print(file_name)
    # read time courses and settings from gutlogo
    sim_steps, cell_counts, settings = read_Gutlogo(file_name)

    # plot gutlogo time courses
    for species in sorted(cell_counts.keys()):
        axs[i][0].plot(sim_steps, cell_counts[species], label=species)
    axs[i][0].set_xlabel('step')
    axs[i][0].set_ylabel('# of agents')

    # run pysb model with gutlogo settings

    # remove extraneous system settings
    param_names = [p.name for p in model.parameters]
    for key in [k for k in settings.keys()]:
        if key not in param_names:
            settings.pop(key)

    # run simulation
    t_span = np.array(sim_steps) * 60 # 1 sim_step = 60 seconds
    result = sim.run(tspan=t_span, param_values=settings)

    # plot pysb time courses
    for obs in sorted(['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']):
        axs[i][1].plot(t_span / 60, result.observables[obs], label=obs)
        axs[i][1].set_xlabel('time (min)')
        axs[i][1].set_ylabel('# of cells')

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

#todo: figure is good but pdf is bad--fix tightlayout and axis issue. Maybe the legend is messing it up.
# Comment out 47-50 and see what the tightlayout looks like. If it is still screwey look at documentationa
# and if it is still screwey
#https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tight_layout.html

# TODO (11/18/23):
#  1) Work on adding a single legend to the figure using the example code provided in 'single_legend.py'
#  2) Complete the example code in the `if __name__ == "__main__"' block of 'ReadGutlogo.py' to produce two plots,
#  one for GutLogo and one for PySB

# todo: create table with files on x and settings on the y using tabulate library https://www.geeksforgeeks.org/how-to-make-a-table-in-python/
# todo: make the biorender figure


