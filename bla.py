import pytraj as pt
import numpy as np
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
        blas = np.load("bla_{}.npy".format(suffix))
        blas = np.average(blas, axis=0)
        plotter(blas, suffix, args.traj_time)
    else:
        d1 = getDistances(args.n_trajs, suffix, 6, 7)
        d2 = getDistances(args.n_trajs, suffix, 7, 8)
        d3 = getDistances(args.n_trajs, suffix, 8, 9)
        bla = np.subtract(np.true_divide(np.add(d1, d3), 2), d2)
        np.save("bla_{}.npy".format(suffix), bla)

def getDistance(traj, suffix, atom1, atom2):
    traj = pt.load('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj), top='m1.prmtop')
    return pt.distance(traj, '@{} @{}'.format(atom1, atom2))

def getDistances(nTrajs, suffix, atom1, atom2):
    return [getDistance(traj, suffix, atom1, atom2) for traj in range(1, nTrajs+1)]

def plotter(bla, suffix, time):
    t = np.linspace(0, time, len(bla), endpoint=True)
    plt.plot(t, bla)
    if suffix == 'abs':
        suffix = 'S0'
    else:
        suffix = 'S1'
    plt.title("BLA {}".format(suffix))
    plt.xlabel("time ps")
    plt.savefig("{}_bla".format(suffix))
    plt.show()
    print("Average Bla: {}A".format(np.average(bla)))

main()
