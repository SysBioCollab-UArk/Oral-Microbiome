from GutModel import model
from ReadGutlogo import read_Gutlogo
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

# read time courses and settings from gutlogo
t_span, cell_counts, settings = read_Gutlogo('GutLogoPopulations_long.csv')

# plot gutlogo time courses
plt.figure()
for species in sorted(cell_counts.keys()):
    plt.plot(t_span, cell_counts[species], label=species)
plt.xlabel('Time (number of time steps)')
plt.ylabel('Population (number of agents)')
plt.legend(loc=0)
plt.tight_layout()

# run pysb model with gutlogo settings
sim = ScipyOdeSimulator(model, t_span, verbose=True)

# remove extraneous system settings
param_names = [p.name for p in model.parameters]
for key in [k for k in settings.keys()]:
    if key not in param_names:
        settings.pop(key)

# run simulation
result = sim.run(param_values=settings)

plt.figure()
for obs in sorted(['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']):
    plt.plot(t_span, result.observables[obs], label=obs)
    plt.xlabel('time (number of time steps)')
    plt.ylabel('# of cells')
    # plt.yscale('log', base = 10)
plt.legend(loc=0)
plt.tight_layout()

plt.show()

# todo: run a bunch of gutlogo simulations with different inputs; save the csv files somewhere. Eventually, we will
#  systematically go through and run all of them and compare our results to theirs
