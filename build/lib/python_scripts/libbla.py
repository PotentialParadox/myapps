import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt

def remove_failures(dss):
    nsteps = max([len(x) for x in dss[0]])
    result = []
    for distance in dss:
        result.append([traj for traj in distance if len(traj) == nsteps])
    return np.array(result)

def getDistance(traj, suffix, atom1, atom2):
    print('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj))
    traj = pt.load('{}/nasqm_{}_{}.nc'.format(traj, suffix, traj), top='m1.prmtop')
    return pt.distance(traj, '@{} @{}'.format(atom1, atom2))

def getDistances(nTrajs, suffix, atom1, atom2):
    return [getDistance(traj, suffix, atom1, atom2) for traj in range(1, nTrajs+1)]

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
    t = np.linspace(0, time, len(d1s), endpoint=True)
    plt.plot(t, d1s)
    plt.plot(t, d2s)
    plt.plot(t, d3s)
    plt.show()
    bla = np.subtract(np.true_divide(np.add(d1s, d3s), 2), d2s)
    t = np.linspace(0, time, len(bla), endpoint=True)
    plt.plot(t, bla)
    if suffix == 'abs':
        suffix = 'S0'
    else:
        suffix = 'S1'
    plt.title("{} BLA {}".format(solvent, suffix))
    plt.xlabel("time ps")
    plt.savefig("{}_bla".format(suffix))
    plt.show()
    print("Average Bla: {}A".format(np.average(bla)))
