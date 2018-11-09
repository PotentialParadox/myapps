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
        plotter(dih, suffix, args.traj_time)
    else:
        dihs = getDihedrals(args.n_trajs, suffix, [[15, 16, 17, 22],
                                                    [11, 14, 15, 16],
                                                    [7, 8, 9, 12],
                                                    [2, 4, 7, 8]])
        dihs_average = reduceDihs(dihs)
        np.save("dihedral_{}.npy".format(suffix), dihs_average)

def reduceDihs(dihs):
    return np.true_divide(reduce(np.add, dihs), len(dihs))

def getDihedral(suffix, traj, atomss):
    traj = pt.load('nasqm_{}_{}.nc'.format(suffix, traj), top='m1.prmtop')
    return [dihedralAbs(pt.dihedral(traj, '@{} @{} @{} @{}'.format(atoms[0], atoms[1], atoms[2], atoms[3])))
            for atoms in atomss]

def getDihedrals(nTrajs, suffix, atoms):
    return [getDihedral(suffix, traj, atoms) for traj in range(1, nTrajs+1)]

def dihedralAbs(dihs):
    return [min(abs(di), abs(180-abs(di))) for di in dihs]

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
