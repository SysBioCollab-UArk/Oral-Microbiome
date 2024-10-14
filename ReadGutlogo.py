import csv


def read_Gutlogo(filename):
    t = []
    pops = [[] for _ in range(4)]
    initial_conditions = [0] * 4  # Initialize initial conditions with zeros

    with open(filename, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        sim_steps = []
        cells = []
        for i, line in enumerate(reader):
            # read settings keys from csv file
            if i == 5:
                settings = dict.fromkeys(line)
            # read settings values from csv file
            elif i == 6:
                for j, key in enumerate(settings.keys()):
                    if line[j].lower() == 'true':
                        settings[key] = True
                    else:
                        settings[key] = float(line[j])
                # replace dashes with underscores in settings names
                for key in [k for k in settings.keys()]:
                    new_key = key.replace('-', '_')
                    if new_key != key:
                        settings[new_key] = settings.pop(key)
            elif i == 18:
                cells = [name.strip('\"') for name in line if name != '']
                cell_counts = {}
                for name in cells:
                    cell_counts[name] = []
            elif i > 19:
                sim_steps.append(int(line[0]))
                offset = 0
                for name in cells:
                    cell_counts[name].append(int(line[1 + offset]))
                    offset += 4

    return sim_steps, cell_counts, settings


if __name__ == '__main__':
    from GutModel import model
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt
    import os

    # Read GutLogo data
    file=os.path.join("GutLogo Simulations","GutLogo plots_2.csv")
    steps_gutlogo, cell_counts_gutlogo, settings_gutlogo = read_Gutlogo(file)

    # Plot GutLogo data
    plt.figure(constrained_layout=True)
    for species in sorted(cell_counts_gutlogo.keys()):
        plt.plot(steps_gutlogo, cell_counts_gutlogo[species], label=species)
    plt.xlabel('Time Step')
    plt.ylabel('Population (number of agents)')
    plt.legend(loc=0)
    plt.title('GutLogo Simulation')


    # remove extraneous system settings
    param_names = [p.name for p in model.parameters]
    for key in [k for k in settings_gutlogo.keys()]:
        if key not in param_names:
            settings_gutlogo.pop(key)

    # Set up PySB simulator
    sim = ScipyOdeSimulator(model, verbose=True)

    # Run PySB simulation with GutLogo settings
    t_span = np.array(steps_gutlogo) * 60 # 1 sim_step = 60 seconds
    result = sim.run(tspan=t_span, param_values=settings_gutlogo)

    # Plot PySB data
    plt.figure(constrained_layout=True)
    for obs in sorted(['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']):
        plt.plot(t_span / 60, result.observables[obs], label=obs)
    plt.xlabel('Time (min)')
    plt.ylabel('# of cells')
    plt.legend(loc=0)
    plt.title('PySB Simulation')
    plt.show()

    # param_values = {'param1': 10, 'param2': 50,...} - key word for changing parameters of imported model
    # make  a parameter (it is currently a variable in gutmodel)
    # initials = {species_1: 100, species_2: 50,...}
    # output = sim.run(param_values = param_values, initials = initials)
    # use initials to change initials, parameters to change everything else
    # map lists of settings with model, do we need all of them?


#  testconst, reservefraction, randist, flowdist, plots-on?, absorption
