'''
Calculates the dipoles of a md simultion set using a file containing charges
dustin tracy (dtracy.uf@gmail.com)
'''
import argparse
import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt
from my_constants import AE_TO_DEBYE

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--full", help="print entire dipole list", action="store_true")
    parser.add_argument("--mags", help="print list of dipole magnitudes", action="store_true")
    parser.add_argument("--average", help="print average dipole", action="store_true")
    parser.add_argument("--plot", help="plot dipole vs time", action="store_true")
    parser.add_argument("--parse", help="parse data from amber outs", action="store_true")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.parse:
        print_dipoless(get_dipoles(args.n_trajs, suffix), suffix)
    dipoless = load_dipoless(suffix)
    if args.plot:
        plot_dipoles(traj_average(dipoless), args.traj_time)

def print_dipoless(dipoless, suffix):
    np.save("dipoles_{}.npy".format(suffix), dipoless)

def get_dipoles(nTrajs, suffix):
    return convert_to_debye(np.array([get_solvent_dipoles(traj, suffix) for traj in range(1, nTrajs+1)]))

def convert_to_debye(xss):
    return np.array([[x*AE_TO_DEBYE for x in xs] for xs in xss])

def get_solvent_dipoles(traj, suffix):
    return pt.analysis.vector.dipole(solvents(get_traj(traj, suffix)))

def solvents(traj):
    return traj['!:1']

def get_traj(traj, suffix):
    return pt.load('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj), top='m1.prmtop')

def load_dipoless(suffix):
    return np.load("dipoles_{}.npy".format(suffix))

def plot_dipoles(dips, traj_time):
    t = np.linspace(0, traj_time, len(dips), endpoint=True)
    plt.plot(t, split_dips(dips, 0))
    plt.show()

def split_dips(dips, index):
    return [dip[index] for dip in dips]

def traj_average(dips):
    return np.average(dips, axis=0)

main()
