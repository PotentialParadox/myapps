import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from python_scripts.libdihedral import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--files", help="labels of the data", default=[1], nargs="+")
    parser.add_argument("--labels", help="file names", default=[1], nargs="+")
    parser.add_argument("--plot", help="print graph", action="store_true")
    parser.add_argument("--solvent", help="the solvent in your system", default="")
    parser.add_argument("--letter", help="the solvent in your system", default="")
    parser.add_argument("--ylims", help="limits", default=[], nargs="+", type=float)
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.plot:
        sns.set()
        sns.set_style("white")
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        for (filename, label) in zip(args.files, args.labels):
            dihs = np.load(filename)
            print("dihss", dihs.shape)
            dih = np.average(dihs[:args.n_trajs], axis=0)
            window_width = 1
            cumsum_vec = np.cumsum(np.insert(dih, 0, 0))
            ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
            plotter(ax, ma_vec, label, args.traj_time)
        if args.ylims != []:
            print(args.ylims)
            ax.set_ylim(args.ylims[0], args.ylims[1])
        ax.text(-0.15, 0.95, args.letter, transform=ax.transAxes,
                fontsize=15, fontweight='bold', va='top')
        plt.legend()
        fig.savefig("{}_dihedral".format(args.solvent))
        plt.show()
    else:
        dihs = getDihedrals(args.n_trajs, suffix, [[18, 17, 16, 15],
                                                   [16, 15, 14, 13]])
        print(dihs.shape)
        dihs_average = np.average(dihs, axis=1)
        print(dihs_average.shape)
        np.save("dihedral_{}.npy".format(suffix), dihs_average)



main()
