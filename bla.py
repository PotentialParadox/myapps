# Calculate the Bond Length Alteration of Using data provided by vmd's
# ANSII graph output
import numpy as np

root_name = 'dg'

data1 = np.loadtxt(root_name+'1.dat')
data2 = np.loadtxt(root_name+'2.dat')
data3 = np.loadtxt(root_name+'3.dat')
d1 = data1[:,1]
d2 = data2[:,1]
d3 = data3[:,1]

bl = ((d1 + d3) / 2) - d2
bla = np.average(bl)
print(bla)
