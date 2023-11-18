import matplotlib.pyplot as plt
import numpy as np

# Create some example data
x = np.linspace(0, 2 * np.pi, 100)
y1 = np.sin(x)
y2 = np.cos(x)

# Create a 2x2 grid of subplots
fig, axes = plt.subplots(nrows=2, ncols=2)

# Create subplots and store line objects for each plot
lines = []  # List to store line objects
labels = ['Data 1', 'Data 2']  # Labels for the legend

for i in range(2):
    for j in range(2):
        ax = axes[i, j]
        if i == 0 and j == 0:
            line, = ax.plot(x, y1, label=labels[0])
            lines.append(line)
        elif i == 1 and j == 1:
            line, = ax.plot(x, y2, label=labels[1])
            lines.append(line)
        else:
            ax.plot(x, y1)
            ax.plot(x, y2)

# Create a subplot for the legend (invisible)
legend_ax = fig.add_subplot(111, frameon=False)
legend_ax.axis('off')  # Turn off the axis for the legend subplot

# Create the legend in the "invisible" subplot
legend_ax.legend(lines, labels, loc='upper left', bbox_to_anchor=(0, 1.15), ncol=2)

plt.tight_layout()
plt.show()
