import re
import argparse
import matplotlib.pyplot as plt
from python_scripts.libquantumcoeffs import *

def main(filename, title):
    filestring = open(filename, 'r').read()
    times = get_times(filestring)
    coeffs = get_coeffs(filestring)
    graph_coeffs(times, coeffs, title)

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", default="coefficient.out")
parser.add_argument("-t", "--title", default="State Coefficients")
args = parser.parse_args()

main(args.filename, args.title)
