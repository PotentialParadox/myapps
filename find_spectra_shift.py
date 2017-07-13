'''
Print the spectra shift between two
spectra given the filename of the
file containing both spectra
and the units you prefer it to be in
'''
import sys
import find_shift
import numpy as np

FILENAME1 = sys.argv[1]
FILENAME2 = sys.argv[2]
UNIT = sys.argv[3]
DATA1 = np.loadtxt(FILENAME1)
DATA2 = np.loadtxt(FILENAME2)
DATA = DATA1
DATA[:, 3] = DATA2[:, 3]
print(find_shift.find_shift(DATA, UNIT))
