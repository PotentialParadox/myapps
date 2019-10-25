'''
Calculates the dipoles of a md simultion set using a file containing charges
dustin tracy (dtracy.uf@gmail.com)
'''
import argparse
from python_scripts.libdipole import print_dipoless, load_dipoless, get_dipoles
from python_scripts.libdipole import plot_dipoles, traj_average, traj_std

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
        plot_dipoles(traj_average(dipoless), traj_std(dipoless), args.traj_time, args.molecule, suffix)

main()

