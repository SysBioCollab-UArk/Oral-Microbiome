import matplotlib.pyplot as plt
import csv


def read_Gutlogo(filename):
    t = []
    pops = [[] for _ in range(4)]
    initial_conditions = [0] * 4  # Initialize initial conditions with zeros

    with open(filename, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        time = []
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
                time.append(int(line[0]))
                offset = 0
                for name in cells:
                    cell_counts[name].append(int(line[1 + offset]))
                    offset += 4

    return time, cell_counts, settings


if __name__ == '__main__':
    from GutModel import model

    time, cell_counts, settings = read_Gutlogo('GutLogoPopulations_long.csv')

    for species in cell_counts.keys():
        plt.plot(time, cell_counts[species], label=species)
    plt.xlabel('Time (number of time steps)')
    plt.ylabel('Population (number of agents)')
    plt.legend(loc=0)

    # todo: now run a simulation with our code, using the settings returned with this code--put it here.
    #  we have imported our gutmodel (line 50) and we have read the csv (line 52). Now, we must modify our model
    #  to have the same system settings as the GutLogo model from netlogo

    for key in settings.keys():
        print(key, settings[key])

    # param_values = {'param1': 10, 'param2': 50,...} - key word for changing parameters of imported model
    # make  a parameter (it is currently a variable in gutmodel)
    # initials = {species_1: 100, species_2: 50,...}
    # output = sim.run(param_values = param_values, initials = initials)
    # use initials to change initials, parameters to change everything else
    # map lists of settings with model, do we need all of them?

    plt.show()

#  testconst, reservefraction, randist, flowdist, plots-on?, absorption
