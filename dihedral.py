import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--plot", help="print graph")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.plot:
        dihs = np.load("dihedral_{}.npy".format(suffix))
        plotter(dihs, suffix, args.time)
    else:
        dihs = np.average(getDihedrals(args.n_trajs, suffix, [16, 15, 14, 11]), axis=0)
        np.save("dihedral_{}.npy".format(suffix), dihs)

def getDihedral(traj, suffix, atoms):
    traj = pt.load('nasqm_{}_{}.nc'.format(suffix, traj), top='m1.prmtop')
    return dihedralAbs(pt.dihedral(traj, '@{} @{} @{} @{}'.format(atoms[0], atoms[1], atoms[2], atoms[3])))

def getDihedrals(nTrajs, suffix, atoms):
    return [getDihedral(traj, suffix, atoms) for traj in range(1, nTrajs+1)]

def dihedralAbs(dihs):
    return [min(abs(di), abs(180-abs(di))) for di in dihs]

def plotter(dihs, suffix, time):
    t = np.linspace(0, time, len(dihs), endpoint=True)
    plt.plot(t, dihs)
    plt.title("Dihedral {}".format(suffix))
    plt.xlabel("time ps")
    plt.show()
    print("Average Dihedral: {} Degrees".format(np.average(dihs)))


main()
