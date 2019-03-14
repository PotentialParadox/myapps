import numpy as np
import argparse
from python_scripts.libdihedral import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--plot", help="print graph", action="store_true")
    parser.add_argument("--solvent", help="the solvent in your system", default="")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.plot:
        dihs = np.load("dihedral_{}.npy".format(suffix))
        print("dihss", dihs.shape)
        dih = np.average(dihs[:args.n_trajs], axis=0)
        window_width = 1
        cumsum_vec = np.cumsum(np.insert(dih, 0, 0))
        ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
        plotter(ma_vec, suffix, args.traj_time, args.solvent)
    else:
        dihs = getDihedrals(args.n_trajs, suffix, [[22, 23, 24, 35]])
        print(dihs.shape)
        dihs_average = np.average(dihs, axis=1)
        print(dihs_average.shape)
        np.save("dihedral_{}.npy".format(suffix), dihs_average)



main()
