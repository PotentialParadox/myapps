'''
Calculates the dot product between the dipoles of the closest solvents to the
line that spans the length of the solute
Dustin Tracy (dtracy.uf@gmail.com)
'''
import re
import matplotlib.pyplot as plt
import pytraj as pt
import numpy as np
from python_scripts.libdipole import completed, foreach_traj


def plot_dipoles(dips, dipsstd, traj_time, molecule, suffix):
    t = np.linspace(0, traj_time, len(dips), endpoint=True)
    low = np.subtract(dips, dipsstd)
    high = np.add(dips, dipsstd)
    plt.plot(t, dips)
    plt.fill_between(t, low, high, alpha=0.5)
    plt.xlabel("Time (ps)")
    state = "S$_1$" if suffix == 'flu' else 'S$_0$'
    plt.title(plot_title(molecule, state))
    plt.legend()
    plt.savefig(save_filename(molecule, state))
    plt.show()

def plot_title(molecule, state):
    return "Dipole vs Time in {} of State {}".format(molecule, state)

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
    return completed(foreach_traj(averageResiduesDotvsTime, nTrajs, suffix, index1, index2))

def soluteEndpointsvsTime(traj, index1, index2):
    return np.array([(frame[index1], frame[index2]) for frame in traj])

def residueMask(residueID):
    return ':{}'.format(residueID)

def residuesDotvsTime(trajID, suffix, index1, index2):
    return np.array([residueDotvsTime(trajectory(suffix, trajID), index1, index2, resID)
                     for resID in sharedClosests(trajID)])

def sharedClosests(trajID):
    closest_abs = getClosests('abs', trajID)
    closest_flu = getClosests('flu', trajID)
    print(shared(closest_abs, closest_flu))
    return shared(closest_abs, closest_flu)

def shared(a, b):
    return [x for x in a if x in b]

def averageResiduesDotvsTime(trajID, suffix, index1, index2):
    '''
    Int -> String -> Int -> Int -> [Float][time]
    '''
    return np.average(residuesDotvsTime(trajID, suffix, index1, index2), axis=0)

def normalizedResidueDipolevsTime(traj, resID):
    '''
    Int -> String -> Int -> Int -> [Float][solvent][time]
    '''
    dipole = pt.analysis.vector.dipole(traj[residueMask(resID)])
    return normalizeVectors(dipole)

def normalizeVectors(vecs):
    norms = np.array([np.linalg.norm(x) for x in vecs])
    return np.array([np.divide(dip, norm) for (dip, norm) in zip(vecs, norms)])

def residueDotvsTime(traj, index1, index2, resID):
    return [residueDot(soluteEndpoints, solventCenter, normalizedSolventDipole)
            for (soluteEndpoints, solventCenter, normalizedSolventDipole) in
            zip(soluteEndpointsvsTime(traj, index1, index2), centerofMassvsTime(traj, resID),
                normalizedResidueDipolevsTime(traj, resID))]

def centerofMassvsTime(traj, resID):
    '''
    [pytrajTraj] -> Int -> [Float]
    Given a pytraj trajectory and a residue ID return the center of mass of the residue vs time
    '''
    return pt.center_of_mass(traj[residueMask(resID)])

def residueDot(soluteEndpoints, solventCenter, normalizedSolventDipole):
    '''
    ([Float], [Float]) -> [Float] -> [Float] -> Float
    Given solute endpoints, the center of geometry of the solvent, and the normalized dipole
    of the solvent return the dot product between the dipole and the shortest line from the center of geometry
    of the solvent to the line created by the solute endpoints
    '''
    return np.dot(normalizedSolventDipole, normalizedShortestVec(solventCenter, soluteEndpoints))

def normalizedShortestVec(solventCenter, soluteEndpoints):
    '''
    [Float] -> ([Float], [Float]) -> [Float]
    Given the center of mass of the solvent and the endpoints of the solute return the direction of
    the shortest distance between the solvent and the line of the solute
    '''
    cross = np.cross(soluteVector(soluteEndpoints), solventCenter - soluteEndpoints[0])
    shortest_line = np.cross(cross, soluteVector(soluteEndpoints))
    return np.divide(shortest_line, np.linalg.norm(shortest_line))

def soluteVector(soluteEndpoints):
    return soluteEndpoints[1] - soluteEndpoints[0]

def trajectory(suffix, index):
    print('{}/nasqm_{}_{}.nc'.format(index, suffix, index))
    return pt.iterload('{}/nasqm_{}_{}.nc'.format(index, suffix, index), top='{}/m1.prmtop'.format(index))

def closestOutput(index):
    return "{}/closest_{}.txt".format(index, index)

def getClosests(suffix, index):
    p_mask = re.compile("qmmask=':1,(.*)',")
    return [int(x) for x in
            re.findall(p_mask, open('{}/nasqm_{}_{}.in'.format(index, suffix, index), 'r').read())[0].split(',')]
