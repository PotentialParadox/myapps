import numpy as np
import pytraj as pt
from functools import reduce
import matplotlib.pyplot as plt
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--plot", help="print graph", action="store_true")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.plot:
        dihs = np.load("dihedral_{}.npy".format(suffix))
        dih = np.average(dihs[:args.n_trajs], axis=0)
        window_width = 1
        cumsum_vec = np.cumsum(np.insert(dih, 0, 0))
        ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
        plotter(ma_vec, suffix, args.traj_time)
    else:
        dihs = getDihedrals(args.n_trajs, suffix, [[22, 23, 24, 35]])
        print(dihs.shape)
        dihs_average = np.average(dihs, axis=1)
        print(dihs_average.shape)
        np.save("dihedral_{}.npy".format(suffix), dihs_average)


def getDihedral(suffix, traj, atomss):
    traj = pt.load('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj), top='m1.prmtop')
    return [dihedralAbs(pt.dihedral(traj, '@{} @{} @{} @{}'.format(atoms[0],
                                                                   atoms[1],
                                                                   atoms[2],
                                                                   atoms[3])))
            for atoms in atomss]

def getDihedrals(nTrajs, suffix, atoms):
    dihs = [getDihedral(suffix, traj, atoms) for traj in range(1, nTrajs+1)]
    ll = max([len(x[0]) for x in dihs])
    return np.array([x for x in dihs if len(x[0]) == ll])

def dihedralAbs(dihs):
    return [min(abs(di), 180-abs(di)) for di in dihs]

def plotter(dihs, suffix, time):
    t = np.linspace(0, time, len(dihs), endpoint=True)
    plt.plot(t, dihs)
    suffix = 'S0' if suffix == 'abs' else 'S1'
    plt.title("Dihedral {}".format(suffix))
    plt.xlabel("time ps")
    plt.savefig("{}_dihedral".format(suffix))
    plt.show()
    print("Average Dihedral: {} Degrees".format(np.average(dihs)))


main()
