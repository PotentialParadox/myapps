'''
Calculates the dipoles of a md simultion set using a file containing charges
dustin tracy (dtrac.uf@gmail.com)
'''
import argparse
import pytraj as pt
import numpy as np
from my_constants import AE_TO_DEBYE

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--full", help="print entire dipole list", action="store_true")
    parser.add_argument("--mags", help="print list of dipole magnitudes", action="store_true")
    parser.add_argument("--average", help="print average dipole", action="store_true")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    dipoless = get_dipoles(args.n_trajs, suffix)
    if args.average:
        print(average_dipole(dipoless))

def average_dipole(dipoles):
    print(diMags(dipoles))
    return np.average(diMags(dipoles))

def diMags(dipoless):
    return np.array([[np.linalg.norm(x) for x in dipoles] for dipoles in dipoless])

def get_dipoles(nTrajs, suffix):
    return convert_to_debye(np.array([get_dipole(traj, suffix) for traj in range(1, nTrajs+1)]))

def get_dipole(traj, suffix):
    traj = pt.load('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj), top='m1.prmtop')
    return pt.analysis.vector.dipole(traj)

def convert_to_debye(xss):
    return np.array([[x*AE_TO_DEBYE for x in xs] for xs in xss])

main()
