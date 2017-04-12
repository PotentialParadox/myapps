import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline

# Change here the title of the graph
graph_title = "Distryrylbenzene in CHCL3"
# Change here the number of states to display
number_states = 2
# Change here the file root for the file
file_root = 'spectra'
# Change to true if absorbance, false if fluorescence
absorbance = True


abs_input_file = file_root + "_abs.output"
flu_input_file = file_root + "_flu.output"

color_code = ['g', 'r', 'b', 'y', 'm']

data_abs = np.loadtxt(abs_input_file)
data_flu = np.loadtxt(flu_input_file)

nm_index = 1
state_intensity_index_start = 2
x_abs = data_abs[:, nm_index]
y_abs = data_abs[:, state_intensity_index_start]
x_flu = data_flu[:, nm_index]
y_flu = data_flu[:, state_intensity_index_start]

plt.plot(x_abs, y_abs, color='r', label='absorption')
plt.plot(x_flu, y_flu, color='g', label='fluorescence')

plt.xlabel('Wavelength, nm')

plt.title(graph_title)
plt.legend(loc=4)
plt.savefig(file_root + '.png')
plt.show()


