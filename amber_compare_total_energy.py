'''
Functions used to compare amber output files
'''
import sys
import amber_out
import matplotlib.pyplot as plt

def compare_total_energies(file_1_stream, file_2_stream,
                           file_1_label="File1", file_2_label="File2"):
    '''
    Plots a graph of the total enegies comparing two separate amber output files
    '''
    file_1_energies = amber_out.find_total_energies(file_1_stream)
    file_2_energies = amber_out.find_total_energies(file_2_stream)
    plt.plot(file_1_energies, label=file_1_label)
    plt.plot(file_2_energies, label=file_2_label)
    title = "Total energy of " + file_1_label + " and " + file_2_label
    plt.title(title)
    plt.show()

FILE_1_FILENAME = str(sys.argv[1])
FILE_2_FILENAME = str(sys.argv[2])
try:
    LABEL_1 = str(sys.argv[3])
    LABEL_2 = str(sys.argv[4])
except IndexError:
    LABEL_1 = "File1"
    LABEL_2 = "File2"
FILE_1_STREAM = open(FILE_1_FILENAME, 'r')
FILE_2_STREAM = open(FILE_2_FILENAME, 'r')
compare_total_energies(FILE_1_STREAM, FILE_2_STREAM, LABEL_1, LABEL_2)
