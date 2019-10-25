'''
Calculates the dipoles of a md simultion set using a file containing charges
dustin tracy (dtracy.uf@gmail.com)
'''
import argparse
import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt
from python_scripts.libmyconstants import AE_TO_DEBYE

def print_dipoless(dipoless, suffix):
    np.save("dipoles_{}.npy".format(suffix), dipoless)

def get_dipoles(nTrajs, suffix):
    '''
    Int -> String -> [Float][traj][time][coordinates]
    '''
    return convert_to_debye(completed(foreach_traj(get_solvent_dipoles, nTrajs, suffix)))

def completed(trajs):
    maxlen = max([len(x) for x in trajs])
    return [traj for traj in trajs if len(traj) == maxlen]

def foreach_traj(func, nTrajs, *args):
    return np.array([func(traj, *args) for traj in range(1, nTrajs+1)])

def convert_to_debye(xss):
    return np.array([[x*AE_TO_DEBYE for x in xs] for xs in xss])

def get_solvent_dipoles(traj, suffix):
    '''
    Traj -> String -> [Float][time][coordinates]
    '''
    return pt.analysis.vector.dipole(solvents(get_traj(traj, suffix)))

def solvents(traj):
    return traj['!:1']

def get_traj(traj, suffix):
    return pt.load('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj), top='{}/m1.prmtop'.format(traj))

def load_dipoless(suffix):
    return np.load("dipoles_{}.npy".format(suffix))

def plot_dipoles(dips, dipsstd, traj_time, molecule, suffix):
    t = np.linspace(0, traj_time, len(dips), endpoint=True)
    plot_xyzs(dips, dipsstd, t)
    plot_magnitudes(dips, t)
    plt.xlabel("Time (ps)")
    plt.ylabel("Debye")
    state = "S$_1$" if suffix == 'flu' else 'S$_0$'
    plt.title(plot_title(molecule, state))
    plt.legend()
    plt.savefig(save_filename(molecule, state))
    plt.show()

def save_filename(molecule, state):
    cleaned_molecule = "".join([m for m in molecule if m not in ('_', '$', ' ', '\'')])
    cleaned_state = "".join([s for s in state if s not in ('_', '$', ' ', '\'')])
    return cleaned_molecule + cleaned_state

def plot_xyzs(dips, dipsstd, t):
    for i in range(3):
        plot_std(t, dips, dipsstd, i)

def plot_std(t, dips, dipsstd, component):
    labels = ['x-component', 'y-component', 'z-component']
    low = np.subtract(split_dips(dips, component), split_dips(dipsstd, component))
    high = np.add(split_dips(dips, component), split_dips(dipsstd, component))
    plt.fill_between(t, low, high, label=labels[component], color='C{}'.format(component), alpha=0.1)
    plt.plot(t, split_dips(dips, component), color='C{}'.format(component))

def plot_magnitudes(dips, t):
    mags = get_magnitudes(dips)
    plt.plot(t, mags, label='magnitude', color='r')

def get_magnitudes(dips):
    return [np.linalg.norm(dip) for dip in dips]

def plot_title(molecule, state):
    return "Dipole vs Time of {} in State {}".format(molecule, state)

def split_dips(dips, index):
    return [dip[index] for dip in dips]

def traj_average(dips):
    return np.average(dips, axis=0)

def traj_std(dips):
    return np.std(dips, axis=0) / np.sqrt(len(dips))
