import pytraj as pt
import numpy as np
import os.path
import matplotlib.pyplot as plt
import seaborn as sns
from python_scripts.libmyconstants import HZ_TO_WAVENUMBER
# from python_scripts.libmymath import myFourierTransform

def remove_failures(dss):
    nsteps = max([len(d[2]) for d in dss])
    return np.array([d for d in dss if len(d[2]) == nsteps])

def getDistance(traj, suffix, pairs):
    filename = '{}/traj_{}/nasqm_{}_{}.nc'.format(suffix, traj, suffix, traj)
    print(filename)
    traj = pt.load(filename, top='m1.prmtop')
    return  [pt.distance(traj, '@{} @{}'.format(pair[0], pair[1]))
             for pair in pairs]

def finished(suffix, traj):
    filename = "{}/traj_{}/nasqm_{}_{}.nc".format(suffix, traj, suffix, traj)
    return os.path.isfile(filename)

def getDistances(nTrajs, suffix, pairs):
    dss = [getDistance(traj, suffix, pairs) for traj in range(1, nTrajs+1)
                     if finished(suffix, traj)]
    return remove_failures(dss)

def calc_bla(d1, d2, d3):
    return np.subtract(np.true_divide(np.add(d1, d3), 2), d2)

def seperate_bondlengths(dss):
    return (np.average(dss[0], axis=0), np.average(dss[1], axis=0), np.average(dss[2], axis=0))

def plotter(dss_s0, dss_s1, time, solvent, near):
    (d1s_s0, d2s_s0, d3s_s0) = seperate_bondlengths(dss_s0)
    (d1s_s1, d2s_s1, d3s_s1) = seperate_bondlengths(dss_s1)
    timestep = time / len(d1s_s0) * 1e-12
    bla_s0 = calc_bla(d1s_s0, d2s_s0, d3s_s0)
    bla_s1 = calc_bla(d1s_s1, d2s_s1, d3s_s1)

    print("S0 info")
    s0_averages = print_averages(d1s_s0, d2s_s0, d3s_s0, timestep)
    print("S1 info")
    print_averages(d1s_s1, d2s_s1, d3s_s1, timestep)

    sns.set()
    sns.set_style("white")
    sns.set_style("ticks")

    # fig1, ax1 = plt.subplots(2)
    # plot_fourier(ax1[1], timestep, d1s_s0, "S0")
    # plot_fourier(ax1[0], timestep, d1s_s1, "S1")
    # ax1[1].set_xlabel("cm$^{-1}$")
    # fig1.savefig("{}-bla-fourier.png".format(solvent), bbox_inches='tight')

    fig2, ax2 = plt.subplots(1,2)
    fig2.set_size_inches(10, 5)

    bla_allpoints = np.dstack((bla_s0, bla_s1))
    bla_ylims = bla_allpoints.min(), bla_allpoints.max()
    plot_bla(ax2[0], time, bla_s1, "S1", bla_ylims)
    ax2[0].axhline(s0_averages[3], ls="--")

    colors = sns.color_palette()
    ax2[1].axhline(s0_averages[0], ls="--", color=colors[0])
    ax2[1].axhline(s0_averages[1], ls="--", color=colors[1])
    ax2[1].axhline(s0_averages[2], ls="--", color=colors[2])
    indiv_allpoints = np.dstack((d1s_s0, d2s_s0, d3s_s0, d1s_s1, d2s_s1, d3s_s1))
    indiv_ylims = indiv_allpoints.min(), indiv_allpoints.max()
    plot_individuals(ax2[1], time, d1s_s1, d2s_s1, d3s_s1, indiv_ylims, near)
    ax2[1].legend(frameon=False, bbox_to_anchor=(1, 1.1))

    fig2.savefig("{}-bla.png".format(solvent), bbox_inches='tight')
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
    half = int(len(bla)/2)
    print("half",half)
    print("Average Bla: {}A".format(np.average(bla[half:])))

# def plot_fourier(ax, ts, yi, state):
#     x, y = myFourierTransform(yi, ts)
#     x = x * HZ_TO_WAVENUMBER
#     ax.axes.get_yaxis().set_ticks([])
#     ax.set_ylabel(state)
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     ax.tick_params(direction='in')
#     x, y = truncate_to(5000, x, y)
#     ax.plot(x[1:], y[1:])

def truncate_to(limit, xs, ys):
    filtered = [(x[0], x[1]) for x in zip(xs, ys) if x[0] <= limit]
    new_xs = [x[0] for x in filtered]
    new_ys = [x[1] for x in filtered]
    return new_xs, new_ys

def print_averages(d1s, d2s, d3s, timestep):
    half = int(len(d1s)/2)
    print(d1s[half:])
    d1 = np.average(d1s[half:])
    d2 = np.average(d2s[half:])
    d3 = np.average(d3s[half:])
    print("d1: ", d1)
    print("d2: ", d2)
    print("d3: ", d3)
    print("d1+d2/2", (d1+d3)/2)
    bla = (d1+d3)/2 - d2
    print("bla", bla)
    print("nsteps", len(d1s))
    print("timestep", timestep)
    print("steps", len(d1s))
    return d1, d2, d3, bla

def plot_individuals(ax, time, d1s, d2s, d3s, ylims, near):
    l1, l2, l3 = ('d1', 'd2', 'd3') if near else ('d4', 'd5', 'd6')
    t = np.linspace(0, time, len(d1s), endpoint=True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(direction='in')
    ax.plot(t, d1s, label=l1)
    ax.plot(t, d2s, label=l2)
    ax.plot(t, d3s, label=l3)
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Bond Length ($\AA$)")
