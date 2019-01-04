'''
Calculates the dot product between the dipoles of the closest solvents to the
line that spans the length of the solute
Dustin Tracy (dtracy.uf@gmail.com)
'''
import argparse
import re
import matplotlib.pyplot as plt
from libdipole import completed, convert_to_debye, foreach_traj
from libdipole import traj_average, traj_std, plot_dipoles
from pynasqm.closestreader import ClosestReader
import pytraj as pt
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("traj_time", help="time of each trajectory", type=float)
    parser.add_argument("--endpoints", help="atomIDs of the edges of the molecule", type=int, nargs="+", default=[1])
    parser.add_argument("--flu", help="apply to fluorescence", action="store_true")
    parser.add_argument("--plot", help="plot dipole vs time", action="store_true")
    parser.add_argument("--molecule", help="the name of the molecule you are using", default="Molecule")
    parser.add_argument("--parse", help="parse data from amber outs", action="store_true")
    args = parser.parse_args()
    suffix = 'flu' if args.flu else 'abs'
    if args.parse:
        writeDipoleDots(dipoleDots(args.n_trajs, suffix, args.endpoints[0], args.endpoints[1]), suffix)
    dipdots = loadDipoleDots(suffix)
    print(dipdots.shape)
    if args.plot:
        plot_dipoles(traj_average(dipdots), traj_std(dipdots), args.traj_time, args.molecule, suffix)

def plot_dipoles(dips, dipsstd, traj_time, molecule, suffix):
    t = np.linspace(0, traj_time, len(dips), endpoint=True)
    low = np.subtract(dips, dipsstd)
    high = np.add(dips, dipsstd)
    plt.plot(t, dips)
    plt.fill_between(t, low, high, alpha=0.5)
    plt.xlabel("Time (ps)")
    plt.ylabel("Debye")
    state = "S$_1$" if suffix == 'flu' else 'S$_0$'
    plt.title(plot_title(molecule, state))
    plt.legend()
    plt.savefig(save_filename(molecule, state))
    plt.show()

def plot_title(molecule, state):
    return "Dipole vs Time of {} in State {}".format(molecule, state)

def save_filename(molecule, state):
    cleaned_molecule = "".join([m for m in molecule if m not in ('_', '$', ' ', '\'')])
    cleaned_state = "".join([s for s in state if s not in ('_', '$', ' ', '\'')])
    return cleaned_molecule + cleaned_state + "dots"

def writeDipoleDots(dipoless, suffix):
    np.save("dipoledots_{}.npy".format(suffix), dipoless)

def loadDipoleDots(suffix):
    return np.load("dipoledots_{}.npy".format(suffix))

def dipoleDots(nTrajs, suffix, index1, index2):
    '''
    Int -> String -> Int -> Int -> [Float][trajs][time]
    '''
    return convert_to_debye(completed(foreach_traj(averageResiduesDotvsTime, nTrajs, suffix, index1, index2)))

def soluteVector(traj, index1, index2):
    return np.array([frame[index2] - frame[index1] for frame in traj])

def residueMask(residueID):
    return ':{}'.format(residueID)

def residuesDotvsTime(trajID, suffix, index1, index2):
    return np.array([residueDotvsTime(trajectory(suffix, trajID), index1, index2, resID)
                     for resID in getClosests(suffix, trajID)])

def averageResiduesDotvsTime(trajID, suffix, index1, index2):
    '''
    Int -> String -> Int -> Int -> [Float][time]
    '''
    return np.average(residuesDotvsTime(trajID, suffix, index1, index2), axis=0)

def residueDipolevsTime(traj, resID):
    '''
    Int -> String -> Int -> Int -> [Float][solvent][time]
    '''
    return pt.analysis.vector.dipole(traj[residueMask(resID)])

def residueDotvsTime(traj, index1, index2, resID):
    return [np.dot(solute, solvent) for (solute, solvent) in
            zip(soluteVector(traj, index1, index2), residueDipolevsTime(traj, resID))]

def trajectory(suffix, index):
    return pt.iterload('{}/nasqm_{}_{}.nc'.format(index, suffix, index), top='{}/m1.prmtop'.format(index))

def closestOutput(index):
    return "{}/closest_{}.txt".format(index, index)

def getClosests(suffix, index):
    p_mask = re.compile("qmmask=':1,(.*)',")
    return [int(x) for x in
            re.findall(p_mask, open('{}/nasqm_{}_{}.in'.format(index, suffix, index), 'r').read())[0].split(',')]

main()
