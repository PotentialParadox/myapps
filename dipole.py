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
    parser.add_argument("--molecule", help="the name of the molecule you are using", default="Molecule")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.parse:
        print_dipoless(get_dipoles(args.n_trajs, suffix), suffix)
    dipoless = load_dipoless(suffix)
    if args.plot:
        plot_dipoles(traj_average(dipoless), args.traj_time, args.molecule, suffix)

def print_dipoless(dipoless, suffix):
    np.save("dipoles_{}.npy".format(suffix), dipoless)

def get_dipoles(nTrajs, suffix):
    return convert_to_debye(completed(foreach_traj(get_solvent_dipoles, nTrajs, suffix)))

def completed(trajs):
    maxlen = max([len(x) for x in trajs])
    return [traj for traj in trajs if len(traj) == maxlen]

def foreach_traj(func, nTrajs, suffix):
    return np.array([func(traj, suffix) for traj in range(1, nTrajs+1)])

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

def plot_dipoles(dips, traj_time, molecule, suffix):
    t = np.linspace(0, traj_time, len(dips), endpoint=True)
    plot_xyzs(dips, t, molecule, suffix)
    plot_magnitudes(dips, t, molecule, suffix)
    plt.legend()
    plt.show()

def plot_xyzs(dips, t, molecule, suffix):
    plt.plot(t, split_dips(dips, 0), label='x')
    plt.plot(t, split_dips(dips, 1), label='y')
    plt.plot(t, split_dips(dips, 2), label='y')
    state = "S1" if suffix == 'flu' else "S0"
    plt.title(plot_title(molecule, state))
    plt.xlabel("Time Ps")

def plot_magnitudes(dips, t, molecule, suffix):
    mags = get_magnitudes(dips)
    plt.plot(t, mags, label='magnitude')

def get_magnitudes(dips):
    return [np.linalg.norm(dip) for dip in dips]

def plot_title(molecule, state):
    return "Dipole vs Time of {} in State {}".format(molecule, state)

def split_dips(dips, index):
    return [dip[index] for dip in dips]

def traj_average(dips):
    return np.average(dips, axis=0)

main()
