import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline

# Change here the title of the graph
graph_title = "Distryrylbenzene in vacuum"
# Change here the number of states to display
number_states = 2
# Change here the file root for the file
file_root = 'spectra_other'
# Change to true if absorbance, false if fluorescence
absorbance = True


input_file = file_root + ".output"

color_code = ['g', 'r', 'b', 'y', 'm']

data = np.loadtxt(input_file)

nm_index = 1
state_intensity_index_start = 2
state_intensity_index_end = 2 + number_states
x = data[:, nm_index]
ys = data[:, state_intensity_index_start:state_intensity_index_end]

for i in range(number_states):
    y = ys[:, i]
    color = color_code[i]
    label = 'S' + str(i+1)
    plt.plot(x, y, color=color, label=label)

plt.xlabel('Wavelength, nm')
if absorbance:
    ylabel = 'Normalized Absorbance'
else:
    ylabel = 'Normalized Fluorescence'

plt.ylabel(ylabel)
plt.title(graph_title)
plt.legend(loc=2)
plt.savefig(file_root + '.png')
plt.show()


