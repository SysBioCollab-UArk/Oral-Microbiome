import matplotlib
import matplotlib.pyplot as plt
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from GutModel import model
from ReadGutlogo import read_Gutlogo
from pysb.simulator import ScipyOdeSimulator
import os

matplotlib.use('Qt5Agg')


class ScrollableWindow(QtWidgets.QMainWindow):
    def __init__(self, fig):
        self.qapp = QtWidgets.QApplication([])

        QtWidgets.QMainWindow.__init__(self)
        self.widget = QtWidgets.QWidget()
        self.setCentralWidget(self.widget)
        self.widget.setLayout(QtWidgets.QVBoxLayout())
        self.widget.layout().setContentsMargins(0, 0, 0, 0)
        self.widget.layout().setSpacing(0)

        self.fig = fig
        self.canvas = FigureCanvas(self.fig)
        self.canvas.draw()
        self.scroll = QtWidgets.QScrollArea(self.widget)
        self.scroll.setWidget(self.canvas)

        self.nav = NavigationToolbar(self.canvas, self.widget)
        self.widget.layout().addWidget(self.nav)
        self.widget.layout().addWidget(self.scroll)

        #self.show()
        #exit(self.qapp.exec_())


if __name__ == "__main__":
    # example of everything up above
    # Directory and file information
    file_dir = '/Users/ahmedtolba/desktop/GutLogo Simulations'
    file_prefix = 'GutLogo Populations'
    file_suffix = 'csv'

    n_files = len(os.listdir(file_dir))
    fig, axs = plt.subplots(nrows=n_files, ncols=2, sharex=False, sharey=True, figsize=[6.4, 4.8 * n_files / 2])  # 6.4, 4.8 default

    legend_created = False  # Track if the legend has been created

    for i in range(n_files):
        print(i)
        # Read time courses and settings from GutLogo
        t_span, cell_counts, settings = read_Gutlogo('%s/%s%d.%s' % (file_dir, file_prefix, i, file_suffix))

        # Plot gutlogo time courses
        for species in sorted(cell_counts.keys()):
            axs[i][0].plot(t_span, cell_counts[species], label=species)
        axs[i][0].set_xlabel('step')
        axs[i][0].set_ylabel('# of agents')

        # Run Pysb model with GutLogo settings
        sim = ScipyOdeSimulator(model, t_span, verbose=True)

        # Remove extraneous system settings
        param_names = [p.name for p in model.parameters]
        for key in [k for k in settings.keys()]:
            if key not in param_names:
                settings.pop(key)

        # Run simulation
        result = sim.run(param_values=settings)

        # Plot observables
        for obs in sorted(['Bact_tot', 'Clost_tot', 'Bifido_tot', 'Desulfo_tot']):
            axs[i][1].plot(t_span, result.observables[obs], label=obs)
            axs[i][1].set_xlabel('time')
            axs[i][1].set_ylabel('# of cells')

        # Only create the legend for the first subplot
        if not legend_created:
            axs[i][0].legend(loc='center', bbox_to_anchor=(1.1, 1.2), ncol=4)
            legend_created = True

    # ax[0].legend(bbox_to_anchor=(0., 1.02, 2.2, .102), loc=3,
    #                ncol=4, mode="expand", borderaxespad=0)

    # Create a scrollable window for the figure
    a = ScrollableWindow(fig)
