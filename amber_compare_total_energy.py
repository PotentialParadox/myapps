'''
Functions used to compare amber output files
'''
import argparse
import pynasqm.amberout
import matplotlib.pyplot as plt

def compare_total_energies(file_1_stream, file_2_stream,
                           file_1_label="File1", file_2_label="File2"):
    '''
    Plots a graph of the total enegies comparing two separate amber output files
    '''
    file_1_energies = pynasqm.amberout.find_total_energies(file_1_stream)
    file_2_energies = pynasqm.amberout.find_total_energies(file_2_stream)
    plt.plot(file_1_energies, label=file_1_label)
    plt.plot(file_2_energies, label=file_2_label)
    title = "Total energy of " + file_1_label + " and " + file_2_label
    plt.title(title)
    plt.show()

parser = argparse.ArgumentParser()
parser.add_argument("file1", help="The first file you want to parse for data")
parser.add_argument("file2", help="The second file you want to parse for data")
parser.add_argument("--file1_label", help="The label you want to give to data of the first file",
                    default="File1")
parser.add_argument("--file2_label", help="The label you want to give to data of the second file",
                    default="File2")
args=parser.parse_args()
FILE_1_STREAM = open(args.file1, 'r')
FILE_2_STREAM = open(args.file2, 'r')
compare_total_energies(FILE_1_STREAM, FILE_2_STREAM, args.file_1_label, args.file_2_label)
