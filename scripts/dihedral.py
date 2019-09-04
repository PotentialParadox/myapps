import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from python_scripts.libdihedral import *
from scipy.optimize.minpack import curve_fit

def null_list():
    while True:
        yield ""

def exp_decay(t, A, K, C):
    return A*np.exp(-K*t)+C

def fit_exp_nonlinear(t, y):
    opt_parms, _ = curve_fit(exp_decay, t, y, maxfev=1000)
    A, K, C = opt_parms
    return A, K, C

def fit_exp_linear(t, y, C=0):
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    return A, K


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--window_width", help="larger number greater smoothness 1=raw", type=int,
                        default=1)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--abs_files", help="labels of the data", default=[1], nargs="+")
    parser.add_argument("--flu_files", help="labels of the data", default=[1], nargs="+")
    parser.add_argument("--labels", help="file names", default=[""], nargs="+")
    parser.add_argument("--plot", help="print graph", action="store_true")
    parser.add_argument("--solvent", help="the solvent in your system", default="")
    parser.add_argument("--letter", help="the solvent in your system", default=null_list, nargs="+")
    parser.add_argument("--ylims", help="limits", default=[], nargs="+", type=float)
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.plot:
        sns.set()
        sns.set_style("white")
        sns.set_style("ticks")
        fig, axs = plt.subplots(1,2)
        fig.set_size_inches(12, 5)
        window_width = args.window_width
        for (filename, label) in zip(args.abs_files, args.labels):
            dihs = np.load(filename)
            print("dihss", dihs.shape)
            dih = np.average(dihs[:args.n_trajs], axis=0)
            cumsum_vec = np.cumsum(np.insert(dih, 0, 0))
            ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
            plotter(axs[0], ma_vec, label, args.traj_time)
        for (filename, label) in zip(args.flu_files, args.labels):
            dihs = np.load(filename)
            dih = np.average(dihs[:args.n_trajs], axis=0)
            cumsum_vec = np.cumsum(np.insert(dih, 0, 0))
            ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
            plotter(axs[1], ma_vec, label, args.traj_time)

            # Non-linear fitting
            t = np.linspace(0, args.traj_time, len(ma_vec), endpoint=True)
            A, K, C = fit_exp_nonlinear(t, ma_vec)
            print("Time constant:", (1/K))
            fit_y = exp_decay(t, A, K, C)
            # axs[1].plot(t, fit_y,
            #             label='Fitted Function:\n $y = %0.2f e^{%0.2f t} + %0.2f$' % (A, K, C))

            # Linear fitting
            # half = int(len(ma_vec)/2)
            # C0 = np.average(ma_vec[half:])
            # print(C0)
            # A, K = fit_exp_linear(t, ma_vec, C0)
            # print("Linear Time constant:", (1/K))

            axs[1].legend(frameon=False)
        if args.ylims != []:
            print(args.ylims)
            axs[0].set_ylim(args.ylims[0], args.ylims[1])
            axs[1].set_ylim(args.ylims[0], args.ylims[1])
        axs[0].text(-0.15, 0.95, args.letter[0], transform=axs[0].transAxes,
                fontsize=15, fontweight='bold', va='top')
        axs[1].text(-0.15, 0.95, args.letter[1], transform=axs[1].transAxes,
                fontsize=15, fontweight='bold', va='top')
        fig.savefig("{}_dihedral".format(args.solvent), bbox_inches='tight')
        plt.show()
    else:
        dihs = getDihedrals(args.n_trajs, suffix, [[18, 17, 16, 15],
                                                   [16, 15, 14, 13]])
        print(dihs.shape)
        dihs_average = np.average(dihs, axis=1)
        print(dihs_average.shape)
        np.save("dihedral_{}.npy".format(suffix), dihs_average)



main()
