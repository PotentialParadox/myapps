'''
The callable python script
'''
import sys
import re
import numpy as np
import force
import my_utils

FILE1 = sys.argv[1]
FILE2 = sys.argv[2]
FORCE_STRING = my_utils.numpy_to_txt(str(force.find_forces(FILE1, FILE2)))
print(FORCE_STRING)
