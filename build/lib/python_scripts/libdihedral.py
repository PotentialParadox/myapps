import numpy as np
import pytraj as pt
from functools import reduce
import matplotlib.pyplot as plt

def getDihedral(suffix, traj, atomss):
    traj = pt.load('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj), top='m1.prmtop')
    return [dihedralAbs(pt.dihedral(traj, '@{} @{} @{} @{}'.format(atoms[0],
                                                                   atoms[1],
                                                                   atoms[2],
                                                                   atoms[3])))
            for atoms in atomss]

def getDihedrals(nTrajs, suffix, atoms):
    dihs = [getDihedral(suffix, traj, atoms) for traj in range(1, nTrajs+1)]
    ll = max([len(x[0]) for x in dihs])
    return np.array([x for x in dihs if len(x[0]) == ll])

def dihedralAbs(dihs):
    return [min(abs(di), 180-abs(di)) for di in dihs]

def plotter(dihs, suffix, time, solvent):
    t = np.linspace(0, time, len(dihs), endpoint=True)
    plt.plot(t, dihs)
    suffix = 'S0' if suffix == 'abs' else 'S1'
    plt.title("{} Dihedral {}".format(solvent, suffix))
    plt.xlabel("time ps")
    plt.savefig("{}_dihedral".format(suffix))
    plt.show()
    print("Average Dihedral: {} Degrees".format(np.average(dihs)))
