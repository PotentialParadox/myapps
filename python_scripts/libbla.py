import pytraj as pt
import numpy as np
import os.path
import matplotlib.pyplot as plt
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
    d1 = np.average(d1s[-100:])
    d2 = np.average(d2s[-100:])
    d3 = np.average(d3s[-100:])
    print("d1: ", d1)
    print("d2: ", d2)
    print("d3: ", d3)
    print("d1+d2/2", (d1+d3)/2)
    print("bla", (d1+d3)/2 - d2)
    print("nsteps", len(d1s))
    timestep = time / len(d1s) * 10e-6
    print("timestep", timestep)
    print("steps", len(d1s))
    # Testing
    # ts = 0.001
    # N = 1000
    # xs = np.linspace(0,N*ts,N)
    # ys = np.sin(50.0 * 2.0*np.pi*xs) + 0.5*np.sin(80.0 * 2.0*np.pi*xs)
    # End of Testing
    ys = d1s
    ts = timestep
    x, y = myFourierTransform(ys, ts)
    x = x * HZ_TO_WAVENUMBER
    plt.plot(x[1:], y[1:])
    plt.show()
    # Fourier Testing
    t = np.linspace(0, time, len(d1s), endpoint=True)
    plt.plot(t, d1s)
    plt.plot(t, d2s)
    plt.plot(t, d3s)
    plt.xlabel("Time (ps)")
    plt.ylabel("Bond Length (A)")
    plt.show()
    bla = np.subtract(np.true_divide(np.add(d1s, d3s), 2), d2s)
    t = np.linspace(0, time, len(bla), endpoint=True)
    plt.plot(t, bla)
    if suffix == 'abs':
        suffix = 'S0'
    else:
        suffix = 'S1'
    plt.title("{} BLA {}".format(solvent, suffix))
    plt.xlabel("Time (ps)")
    plt.ylabel("Deformation (A)")
    plt.savefig("{}_bla".format(suffix))
    plt.show()
    print("Average Bla: {}A".format(np.average(bla)))
