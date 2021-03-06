'''
Functions used to compare amber output files
'''
import argparse
import pynasqm.amberout
import numpy as np
import matplotlib.pyplot as plt
from my_utils import numpy_to_txt

def compare_scf_energies(time_step, n_steps, file_1_stream, file_2_stream,
                           file_1_label="File1", file_2_label="File2", step=1):
    '''
    Plots a graph of the total enegies comparing two separate amber output files,
    time is in fs
    '''
    file_1_energies_string = pynasqm.amberout.find_scf_energies(file_1_stream)
    file_2_energies_string = pynasqm.amberout.find_scf_energies(file_2_stream)
    file_1_energies = np.fromstring(file_1_energies_string, sep="\n")[::step]
    file_2_energies = np.fromstring(file_2_energies_string, sep="\n")[::step]
    file_1 = open(file_1_label + ".energies", 'w')
    for i in file_1_energies:
        file_1.write(numpy_to_txt(i) + "\n")
    file_1.close()
    file_2 = open(file_2_label + ".energies", 'w')
    for i in file_2_energies:
        file_2.write(numpy_to_txt(i) + "\n")
    file_2.close()
    print(n_steps)
    print(step)
    x = [i*time_step for i in range(0, n_steps+1, step)]
    plt.plot(x, file_1_energies, label=file_1_label)
    plt.plot(x, file_2_energies, label=file_2_label)
    title = "Total energy of " + file_1_label + " and " + file_2_label
    plt.title(title)
    plt.ylabel("Energy (kcal/mol)")
    plt.xlabel("Time (ps)")
    plt.legend()
    plt.show()

parser = argparse.ArgumentParser()
parser.add_argument("--file1", help="the first file you want to parse")
parser.add_argument("--file2", help="the second file you want to parse")
parser.add_argument("--file1_label", help="the label you want to have"\
                    "your graph for the first file", default="file1")
parser.add_argument("--file2_label", help="the label you want to have your graph"\
                    " for the second file", default="file2")
parser.add_argument("--time_step", help="the time step measured in PS", type=float)
parser.add_argument("--n_steps", help="the number of steps of the trajectories", type=int)
parser.add_argument("--step", help="step size you want to use in your graph", type=int)
args = parser.parse_args()
FILE1STREAM = open(args.file1, 'r')
FILE2STREAM = open(args.file2, 'r')
compare_scf_energies(args.time_step, args.n_steps, FILE1STREAM, FILE2STREAM, args.file1_label,
                     args.file2_label, args.step)
