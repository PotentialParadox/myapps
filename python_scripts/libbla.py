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

def calc_bla(d1, d2, d3):
    return np.subtract(np.true_divide(np.add(d1, d3), 2), d2)

def seperate_bondlengths(dss):
    return (np.average(dss[0], axis=0), np.average(dss[1], axis=0), np.average(dss[2], axis=0))

def plotter(dss_s0, dss_s1, time, solvent):
    (d1s_s0, d2s_s0, d3s_s0) = seperate_bondlengths(dss_s0)
    (d1s_s1, d2s_s1, d3s_s1) = seperate_bondlengths(dss_s1)
    timestep = time / len(d1s_s0) * 1e-12
    bla_s0 = calc_bla(d1s_s0, d2s_s0, d3s_s0)
    bla_s1 = calc_bla(d1s_s1, d2s_s1, d3s_s1)

    print("S0 info")
    print_averages(d1s_s0, d2s_s0, d3s_s0, timestep)
    print("S1 info")
    print_averages(d1s_s1, d2s_s1, d3s_s1, timestep)

    sns.set()
    sns.set_style("white")
    sns.set_style("ticks")

    fig1, ax1 = plt.subplots(2)
    plot_fourier(ax1[1], timestep, d1s_s0, "S0")
    plot_fourier(ax1[0], timestep, d1s_s1, "S1")
    ax1[1].set_xlabel("cm$^{-1}$")
    fig1.savefig("{}-bla-fourier.png".format(solvent))

    fig2, ax2 = plt.subplots(1, 2)
    fig2.set_size_inches(10, 5)
    allpoints = np.dstack((d1s_s0, d2s_s0, d3s_s0, d1s_s1, d2s_s1, d3s_s1))
    ylims = allpoints.min(), allpoints.max()
    print('ylims', ylims)
    print('allpoints', allpoints.shape)
    plot_individuals(ax2[0], time, d1s_s0, d2s_s0, d3s_s0, ylims)
    plot_individuals(ax2[1], time, d1s_s1, d2s_s1, d3s_s1, ylims)
    ax2[1].legend(frameon=False, bbox_to_anchor=(1, 1.1))
    fig2.savefig("{}-bla-individuals.png".format(solvent))

    fig3, ax3 = plt.subplots(1,2)
    fig3.set_size_inches(10, 5)
    allpoints = np.dstack((bla_s0, bla_s1))
    ylims = allpoints.min(), allpoints.max()
    plot_bla(ax3[0], time, bla_s0, "S0", ylims)
    plot_bla(ax3[1], time, bla_s1, "S1", ylims)
    fig3.savefig("{}-bla.png".format(solvent))

    plt.show()

def plot_bla(ax, time, bla, state, ylims):
    t = np.linspace(0, time, len(bla), endpoint=True)
    ax.plot(t, bla)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Deformation ($\AA$)")
    ax.set_ylim(ylims[0], ylims[1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(direction='in')
    print("Average Bla: {}A".format(np.average(bla)))

def plot_fourier(ax, ts, yi, state):
    x, y = myFourierTransform(yi, ts)
    x = x * HZ_TO_WAVENUMBER
    ax.axes.get_yaxis().set_ticks([])
    ax.set_ylabel(state)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
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

def plot_individuals(ax, time, d1s, d2s, d3s, ylims):
    t = np.linspace(0, time, len(d1s), endpoint=True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(direction='in')
    ax.plot(t, d1s, label='d1')
    ax.plot(t, d2s, label='d2')
    ax.plot(t, d3s, label='d3')
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Bond Length ($\AA$)")
