import matplotlib.pyplot as plt
import csv


def read_Gutlogo(filename):
    t = []
    y_1 = []
    y_2 = []
    y_3 = []
    y_4 = []

    with open(filename, 'r', encoding='utf-8-sig') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        next(plots)  # skip header
        for row in plots:
            t.append(int(row[0]))
            y_1.append(int(row[1]))
            y_2.append(int(row[2]))
            y_3.append(int(row[3]))
            y_4.append(int(row[4]))

    labels = ['Closts', 'Bifidos', 'Desulfos', 'Bacteroides']
    pops = [y_1, y_2, y_3, y_4]

    return t, pops, labels


if __name__ == '__main__':
    x, pops, label = read_Gutlogo('GutLogo_Populations.csv')
    for y, l in zip(pops, label):
        plt.plot(x, y, label=l)
    plt.xlabel('Time (number of time steps)')
    plt.ylabel('Population (number of agents)')
    plt.title('Kinetics Model of Gut Microbiome')
    plt.legend(loc = 0)
    plt.show()
