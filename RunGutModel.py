from GutModel import model
from ReadGutlogo import read_Gutlogo
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import os

file_dir = '/Users/ahmedtolba/desktop/GutLogo Simulations'
file_prefix = 'GutLogo Populations'
file_suffix = 'csv'

n_files = len(os.listdir(file_dir))
fig, axs = plt.subplots(nrows=n_files, ncols=2, sharex=True, sharey=True, figsize=[6.4, 4.8*n_files/2]) #6.4, 4.8 default

for i in range(len(os.listdir(file_dir))):
    print(i)
    # read time courses and settings from gutlogo
    t_span, cell_counts, settings = read_Gutlogo('%s/%s%d.%s' % (file_dir, file_prefix, i, file_suffix))

    # plot gutlogo time courses
    #plt.figure()
    for species in sorted(cell_counts.keys()):
        axs[i][0].plot(t_span, cell_counts[species], label=species)
    axs[i][0].set_xlabel('Time (number of time steps)')
    axs[i][0].set_ylabel('Population (number of agents)')
#    axs[i][0].legend(loc=0)
    #plt.tight_layout()

    # run pysb model with gutlogo settings
    sim = ScipyOdeSimulator(model, t_span, verbose=True)

    # remove extraneous system settings
    param_names = [p.name for p in model.parameters]
    for key in [k for k in settings.keys()]:
        if key not in param_names:
            settings.pop(key)

    # run simulation
    result = sim.run(param_values=settings)

    #plt.figure()
    for obs in sorted(['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']):
        axs[i][1].plot(t_span, result.observables[obs], label=obs)
        axs[i][1].set_xlabel('time (number of time steps)')
        axs[i][1].set_ylabel('# of cells')
        # plt.yscale('log', base = 10)
#    axs[i][1].legend(loc=0)

    plt.tight_layout()

plt.show()

# todo: get the plots to display in a useful way, maybe break up into chunks of 5 at a time,
#  very automated and general way. Ex.) 14 --> 5, 5, 4
# todo: how to generate a legend for the whole figure --> one possible answer is
# todo: have a spreadsheet or something that lists what the differences are in each file to figure out how each of them
# are different
# https://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots
#  systematically go through and run all of them and compare our results to theirs


