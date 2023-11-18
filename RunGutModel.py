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

for i in range(n_files):
    print(i)
    # read time courses and settings from gutlogo
    t_span, cell_counts, settings = read_Gutlogo('%s/%s%d.%s' % (file_dir, file_prefix, i, file_suffix))

    # plot gutlogo time courses
    for species in sorted(cell_counts.keys()):
        axs[i][0].plot(t_span, cell_counts[species], label=species)
    axs[i][0].set_xlabel('step')
    axs[i][0].set_ylabel('# of agents')

    # run pysb model with gutlogo settings

    # remove extraneous system settings
    param_names = [p.name for p in model.parameters]
    for key in [k for k in settings.keys()]:
        if key not in param_names:
            settings.pop(key)

    # run simulation
    result = sim.run(tspan=t_span, param_values=settings)

    # plot pysb time courses
    for obs in sorted(['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']):
        axs[i][1].plot(t_span, result.observables[obs], label=obs)
        axs[i][1].set_xlabel('time')
        axs[i][1].set_ylabel('# of cells')

    # Only create the legend for the first subplot
    if not legend_created:
        axs[i][0].legend(loc='center', bbox_to_anchor=(1.1, 1.2), ncol=4)
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
