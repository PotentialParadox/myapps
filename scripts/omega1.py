import warnings
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.optimize.minpack import curve_fit
import seaborn as sns

def exp_decay(t, A, B, K1, K2, C):
    return A*np.exp(-K1*t) + B*np.exp(-K2*t) + C

def fit_exp_nonlinear(t, y, c_guess=2.4):
    opt_parms, _ = curve_fit(exp_decay, t, y, p0=[0.05, 0.24, 0.09, 3.49, c_guess])
    A, B, K1, K2, C = opt_parms
    return A, B, K1, K2, C

def fitting(ax, t, omega1s):
    A, B, K1, K2, C = fit_exp_nonlinear(t, omega1s)
    fit_y = exp_decay(t, A, B, K1, K2, C)
    # t_conv = find_convergence_point(t, fit_y, C)
    # print("Converged at {}ps".format(t_conv))
    # ax.axvline(x=t_conv, linestyle="dashed", linewidth=1)
    ax.plot(t, fit_y,
            label='Fitted Function:\n $y = %2.2f e^{-%2.2f t} + %2.2f e^{-%2.2f t} + %2.2f$' % (A, B, K1, K2, C))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("traj_time", help="time of each trajectory in PS", type=float)
    parser.add_argument("n_trajs", help="number of trajectories", type=int)
    parser.add_argument("file_out", help="where to save figure")
    parser.add_argument("--labels", help="labels of omegas", nargs="+")
    parser.add_argument("--ymin", help="minimum value for y in eV", type=float)
    parser.add_argument("--ymax", help="maximum value for y in eV", type=float)
    args = parser.parse_args()
    ylims = [args.ymin, args.ymax]
    omega1ss = np.array(list(map(smooth, get_data(args.n_trajs))))
    avg_omega1s = np.average(omega1ss, axis=0)
    labels = [""] * len(omega1ss) if args.labels is None else args.labels
    plot_omegass(avg_omega1s, ylims, labels, args.traj_time, args.file_out)

def read_file(filename):
    np.seterr(all='raise')
    warnings.filterwarnings('error')
    try:
        return np.loadtxt(filename)
    except UserWarning:
        return []

def get_data(n_trajs):
    data = np.array([read_file("traj_{}.out".format(traj)) for traj in range(1, n_trajs+1)])
    m = max([len(d) for d in data])
    return np.array([d for d in data if len(d) == m])

def plot_omegass(omega1ss, ylims, labels, time, fileout):
    fig1, ax1 = plt.subplots()
    ax1.set_xlabel("Time (ps)")
    ax1.set_ylabel("Energy (eV)")
    # ax1.set_ylim(ylims)
    plot_omegas(ax1, time, omega1ss, labels[0])
    # plotter = lambda data: plot_omegas(ax1, time, data[0], data[1])
    # list(map(plotter, zip(omega1ss, labels)))
    fig1.savefig(fileout, bbox_inches='tight')
    plt.show()


def smooth(ys):
    window_width = 1
    cumsum_vec = np.cumsum(np.insert(ys, 0, 0))
    ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    return ma_vec

def plot_omegas(ax, time, omega1s, label):
    t = np.linspace(0, time, len(omega1s), endpoint=True)
    ax.plot(t, omega1s, label=label)
    # fitting(ax, t, omega1s)
    # ax.legend()

def print_average(data, label):
    half = int(len(data)/2)
    average = np.average(data[half:])
    print("Last half average of {} is {}".format(label, average))
    stdev = np.std(data[half:])
    print("Last half stdev of {} is {}".format(label, stdev))
    return average

def find_convergence_point(ts, fit_y, C):
    for (t, y) in zip(ts, fit_y):
        if y - C <= 0.01:
            return t

main()
