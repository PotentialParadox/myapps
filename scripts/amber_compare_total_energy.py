'''
Functions used to compare amber output files
'''
import argparse
import pynasqm.amberout
import seaborn as sns
import matplotlib.pyplot as plt

def compare_total_energies(file_1_stream, file_2_stream,
                           file_1_label="File1", file_2_label="File2", timestep=1):
    '''
    Plots a graph of the total enegies comparing two separate amber output files
    '''
    file_1_energies = pynasqm.amberout.find_total_energies(file_1_stream)
    file_2_energies = pynasqm.amberout.find_total_energies(file_2_stream)

    sns.set()
    sns.set_style("white")
    sns.set_style("ticks")

    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    x = [x*timestep for x in range(len(file_1_energies))]
    ax.plot(x,file_1_energies, label=file_1_label, color=sns.xkcd_rgb["medium green"])
    p = ax.plot(x,file_2_energies, label=file_2_label, color=sns.xkcd_rgb["pale red"])
    ax.set_ylabel("Energy (kcal/mol)", color=p[0].get_color())
    ax.legend()

    ax2 = ax.twinx()
    p2 = ax2.plot(x,file_2_energies-file_1_energies, label=file_2_label, color=sns.xkcd_rgb["denim blue"])
    ax2.set_ylim(0,0.002)
    ax2.set_ylabel("Energy Difference (kcal/mol)", color=p2[0].get_color())

    plt.show()

parser = argparse.ArgumentParser()
parser.add_argument("file1", help="The first file you want to parse for data")
parser.add_argument("file2", help="The second file you want to parse for data")
parser.add_argument("--label1", help="The label you want to give to data of the first file",
                    default="File1")
parser.add_argument("--label2", help="The label you want to give to data of the second file",
                    default="File2")
parser.add_argument("--timestep", help="Timestep in ps", default=1, type=float)
args=parser.parse_args()
FILE_1_STREAM = open(args.file1, 'r')
FILE_2_STREAM = open(args.file2, 'r')
compare_total_energies(FILE_1_STREAM, FILE_2_STREAM, args.label1, args.label2, args.timestep)
