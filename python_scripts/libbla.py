import pytraj as pt
import numpy as np
import os.path
import matplotlib.pyplot as plt
import seaborn as sns
from python_scripts.libmyconstants import HZ_TO_WAVENUMBER
from python_scripts.libmymath import myFourierTransform

def remove_failures(dss):
    nsteps = max([len(x) for x in dss[0]])
    result = []
    for distance in dss:
        result.append([traj for traj in distance if len(traj) == nsteps])
    return np.array(result)

def getDistance(traj, suffix, atom1, atom2):
    print('{}/traj_{}/nasqm_{}_{}.nc'.format(suffix, traj, suffix, traj))
    traj = pt.load('{}/traj_{}/nasqm_{}_{}.nc'.format(suffix, traj, suffix, traj), top='m1.prmtop')
    return pt.distance(traj, '@{} @{}'.format(atom1, atom2))

def finished(suffix, traj):
    filename = "{}/traj_{}/nasqm_{}_{}.nc".format(suffix, traj, suffix, traj)
    return os.path.isfile(filename)

def getDistances(nTrajs, suffix, atom1, atom2):
    return [getDistance(traj, suffix, atom1, atom2) for traj in range(1, nTrajs+1)
            if finished(suffix, traj)]

def plotter(dss, suffix, time, solvent):
    d1s = np.average(dss[0], axis=0)
    d2s = np.average(dss[1], axis=0)
    d3s = np.average(dss[2], axis=0)
    timestep = time / len(d1s) * 1e-12

    print_averages(d1s, d2s, d3s, timestep)

    sns.set()
    sns.set_style("white")
    sns.set_style("ticks")

    fig1, ax1 = plt.subplots()
    plot_fourier(ax1, timestep, d1s)

    fig2, ax2 = plt.subplots()
    plot_individuals(ax2, time, d1s, d2s, d3s)

    fig3, ax3 = plt.subplots()
    plot_bla(ax3, time, d1s, d2s, d3s, suffix, solvent)

    plt.show()

def plot_bla(ax, time, d1s, d2s, d3s, suffix, solvent):
    bla = np.subtract(np.true_divide(np.add(d1s, d3s), 2), d2s)
    t = np.linspace(0, time, len(bla), endpoint=True)
    ax.plot(t, bla)
    if suffix == 'abs':
        suffix = 'S0'
    else:
        suffix = 'S1'
    ax.set_title("{} BLA {}".format(solvent, suffix))
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Deformation ($\AA$)")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(direction='in')
    print("Average Bla: {}A".format(np.average(bla)))

def plot_fourier(ax, ts, yi):
    x, y = myFourierTransform(yi, ts)
    x = x * HZ_TO_WAVENUMBER
    ax.axes.get_yaxis().set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel("cm$^{-1}$")
    ax.tick_params(direction='in')
    x, y = truncate_to(5000, x, y)
    ax.plot(x[1:], y[1:])

def truncate_to(limit, xs, ys):
    filtered = [(x[0], x[1]) for x in zip(xs, ys) if x[0] <= limit]
    new_xs = [x[0] for x in filtered]
    new_ys = [x[1] for x in filtered]
    return new_xs, new_ys

def print_averages(d1s, d2s, d3s, timestep):
    d1 = np.average(d1s[-100:])
    d2 = np.average(d2s[-100:])
    d3 = np.average(d3s[-100:])
    print("d1: ", d1)
    print("d2: ", d2)
    print("d3: ", d3)
    print("d1+d2/2", (d1+d3)/2)
    print("bla", (d1+d3)/2 - d2)
    print("nsteps", len(d1s))
    print("timestep", timestep)
    print("steps", len(d1s))

def plot_individuals(ax, time, d1s, d2s, d3s):
    t = np.linspace(0, time, len(d1s), endpoint=True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(direction='in')
    ax.plot(t, d1s, label='d1')
    ax.plot(t, d2s, label='d2')
    ax.plot(t, d3s, label='d3')
    ax.legend(frameon=False, bbox_to_anchor=(1, 1.1))
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Bond Length ($\AA$)")
