import numpy as np
import sys

def print_average(averages):
    for avg in averages:
        line = ""
        for a in avg:
            line += "{0:18.8f}".format(a)
        print(line)

def main():
    raw_data = [np.loadtxt(fin) for fin in sys.argv[1:]]
    max_length = max([len(d) for d in raw_data])
    min_length = max_length
    filter_data = np.array([d[:min_length] for d in raw_data if len(d) >= min_length])
    averages = np.average(filter_data, axis=0)
    print_average(averages)

main()
