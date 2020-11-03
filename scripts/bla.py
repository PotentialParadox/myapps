import numpy as np
import argparse
from python_scripts.libbla import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory in PS", type=float)
    parser.add_argument("--dirname", help="directory name of the trajfiles eg. (qmexcited)")
    parser.add_argument("--plot", help="print graph", action="store_true")
    parser.add_argument("--solvent", help="solvent used in calculation", default="")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--near", action="store_true")
    parser.add_argument("--output", "-o", help="output figure file")
    args = parser.parse_args()
    suffix = args.dirname
    distance = "close_" if args.near else ""
    if args.plot:
        dss_s0 = remove_failures(np.load("bla_{}abs.npy".format(distance)))
        dss_s1 = remove_failures(np.load("bla_{}flu.npy".format(distance)))
        plotter(dss_s0[:, :args.n_trajs, ::1], dss_s1[:, :args.n_trajs, ::1],
                args.traj_time, args.solvent, args.near)
    else:
        # PPV3-NO2
        pairs = [(6,7), (7,8), (8,9)] # Near Pairs
        # pairs = [(17,16), (16,15), (15,14)] # Far Pairs

        # PPV3
        # pairs = [(24,23), (23,22), (22,21)]
        # pairs = [(4,7), (7,8), (8,9)]
        # pairs = [(17,16), (16,15), (15,12)]

        data = getDistances(args.n_trajs, suffix, pairs)
        d1r = data[:,0]
        d2r = data[:,1]
        d3r = data[:,2]
        d1 = np.hstack((np.ones((data.shape[0], 1)), d1r))
        d2 = np.hstack((np.ones((data.shape[0], 1))*2, d2r))
        d3 = np.hstack((np.ones((data.shape[0], 1))*3, d3r))
        bla = np.vstack((d1, d2, d3))
        if args.debug:
            print(bla.shape)
        # np.save("bla_{}.npy".format(suffix), bla) #Used for far
        np.save("bla_close_{}.npy".format(suffix), bla) #Used for far
        if args.output:
            np.savetxt(args.output, bla, delimiter=',')

main()
