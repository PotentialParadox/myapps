'''
Commonly used utilities
'''
import re
import numpy as np

def numpy_to_txt(numpy_array):
    numpy_to_txt = str(numpy_array)
    numpy_to_txt = re.sub("\[", ' ', numpy_to_txt)
    numpy_to_txt = re.sub("\]", ' ', numpy_to_txt)
    return numpy_to_txt
