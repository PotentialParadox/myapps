import numpy as np
import os
import sys

def print_average(averages):
    for avg in averages:
        line = ""
        for a in avg:
            line += "{0:18.8f}".format(a)
        print(line)

def my_empty(fin):
    filename = "{}/{}".format(os.getcwd(), fin)
    return os.path.getsize(filename) == 0

def main():
    nonempty_files = [f for f in sys.argv[1:] if not my_empty(f)]
    raw_data = [np.loadtxt(fin) for fin in nonempty_files]
    max_length = max([len(d) for d in raw_data])
    min_length = max_length
    filter_data = np.array([d[:min_length] for d in raw_data if len(d) >= min_length])
    averages = np.average(filter_data, axis=0)
    print_average(averages)

main()
