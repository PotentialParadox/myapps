import numpy as np
import argparse
from python_scripts.libbla import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--plot", help="print graph", action="store_true")
    parser.add_argument("--solvent", help="solvent used in calculation", default="")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.plot:
        dss = remove_failures(np.load("bla_{}.npy".format(suffix)))
        plotter(dss, suffix, args.traj_time, args.solvent)
    else:
        d1 = getDistances(args.n_trajs, suffix, 6, 7)
        d2 = getDistances(args.n_trajs, suffix, 7, 8)
        d3 = getDistances(args.n_trajs, suffix, 8, 9)
        bla = np.array([d1, d2, d3])
        np.save("bla_{}.npy".format(suffix), bla)

main()
