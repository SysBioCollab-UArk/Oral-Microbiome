from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.util import alias_model_components
from pysb.core import as_complex_pattern
from GutModel import model
import numpy as np
import matplotlib.pyplot as plt

# print(model.parameters)
# quit()

alias_model_components()

Observable('Bacteroides_stuck', Bacteroides(stuck='s') + Bacteroides(stuck='p'))

tspan = np.arange(2001) * 60  # 1 sim_step = 60 seconds
sim = ScipyOdeSimulator(model, tspan)

# obs_to_plot = ['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']
obs_to_plot = ['Bacteroides_stuck']

perturbations = [
    {'k_Death_Age': 0, 'seedchance': 0},
    {'k_Death_Age': 0, 'seedchance': 0, 'unstuckchance': 0}
]

for param_values in perturbations:
    plt.figure(constrained_layout=True)
    result = sim.run(param_values = param_values)
    for obs in obs_to_plot:
        plt.plot(tspan / 60, result.observables[obs], lw=2, label=obs)
    plt.xlabel('Time (min)')
    plt.ylabel('# of cells')
    plt.legend(loc='best')

plt.show()
