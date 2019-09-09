'''
Calculates the dot product between the dipoles of the closest solvents to the
line that spans the length of the solute
Dustin Tracy (dtracy.uf@gmail.com)
'''
import argparse
from python_scripts.libsolventdipolerotation import loadDipoleDots
from python_scripts.libsolventdipolerotation import dipoleDots
from python_scripts.libsolventdipolerotation import writeDipoleDots
from python_scripts.libsolventdipolerotation import plot_dipoles
from python_scripts.libdipole import traj_average
from python_scripts.libdipole import traj_std

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--endpoints", help="atomIDs of the edges of the molecule", type=int, nargs="+", default=[1])
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--plot", help="plot dipole vs time", action="store_true")
    parser.add_argument("--solvent", help="the name of the molecule you are using", default="Molecule")
    parser.add_argument("--parse", help="parse data from amber outs", action="store_true")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.parse:
        writeDipoleDots(dipoleDots(args.n_trajs, suffix, args.endpoints[0], args.endpoints[1]), suffix)
    dipdots = loadDipoleDots(suffix)[:args.n_trajs]
    print(dipdots.shape)
    if args.plot:
        plot_dipoles(traj_average(dipdots), traj_std(dipdots), args.traj_time, args.solvent, suffix)

main()
