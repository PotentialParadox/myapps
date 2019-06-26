import os
import numpy as np
import pytraj as pt
import matplotlib.pyplot as plt

def getDihedral(suffix, traj, atomss):
    traj_file = '{0}/traj_{1}/nasqm_{0}_{1}.nc'.format(suffix, traj)
    print(traj_file)
    traj = pt.load(traj_file, top='m1.prmtop')
    return [dihedralAbs(pt.dihedral(traj, '@{} @{} @{} @{}'.format(atoms[0],
                                                                   atoms[1],
                                                                   atoms[2],
                                                                   atoms[3])))
            for atoms in atomss]

def finished(suffix, traj):
    filename = "{0}/traj_{1}/nasqm_{0}_{1}.nc".format(suffix, traj)
    return os.path.isfile(filename)

def getDihedrals(nTrajs, suffix, atoms):
    dihs = [getDihedral(suffix, traj, atoms) for traj in range(1, nTrajs+1)
            if finished(suffix, traj)]
    ll = max([len(x[0]) for x in dihs])
    return np.array([x for x in dihs if len(x[0]) == ll])

def dihedralAbs(dihs):
    return [min(abs(di), 180-abs(di)) for di in dihs]

def plotter(ax, dihs, label, time):
    plot_dihedral(ax, dihs, time, label)
    print("Average Dihedral: {} Degrees".format(np.average(dihs)))

def plot_dihedral(ax, dihs, time, label):
    t = np.linspace(0, time, len(dihs), endpoint=True)
    ax.plot(t, dihs, label=label)
    ax.set_xlabel("Time (ps)")
