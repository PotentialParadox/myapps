import numpy as np
import seaborn as sns
from PIL import Image
import matplotlib.pyplot as plt
from scipy.optimize.minpack import curve_fit
import io
import argparse

def exp_decay(t, A, K, C):
    return A*np.exp(-K*t) + C

def fit_exp_nonlinear(ax, t, y):
    dy = 10 - y[-1]
    opt_parms, _ = curve_fit(exp_decay, t, y+dy, p0=[0.5, 0.5, 10])
    A, K, Copt = opt_parms
    C = Copt-dy
    fit_y = exp_decay(t, A, K, C)
    ax.plot(t, fit_y,
            label='Fitted Function:\n $y = %2.2f e^{-%2.4f t} + %2.2f$' % (A, K, C))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="data separated by newlines")
    parser.add_argument("--title", default="")
    parser.add_argument("--ylims", nargs="+", type=float)
    parser.add_argument("--xlims", nargs="+", type=float)
    parser.add_argument("--ylabel")
    parser.add_argument("--xlabel")
    parser.add_argument("--xscale", help="Scale factor for x axis", type=float)
    parser.add_argument("--firstx", help="The first column is x", action="store_true")
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--fitexp", action="store_true")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    sns.set()
    sns.set_style("white")
    sns.set_style("ticks")

    data = np.loadtxt(args.file)
    [xmin, xmax] = [0, (len(data)-1)]  if args.xlims is None else args.xlims
    if args.xscale is not None:
        xmax=xmax*args.xscale
    xs = np.linspace(xmin, xmax, len(data))
    fig, ax = plt.subplots(1)
    ax.set_title(args.title)
    if args.ylims is not None:
        plt.ylim(args.ylims)
    startindex = 0
    if args.firstx:
        startindex = 1
        xs = data[:, 0]
        data = data[:, 1:]
    if args.labels is not None:
        for i in range(startindex, len(data[1, :])):
            plt.plot(xs, data[:, i], label=args.labels[i])
            plt.legend(frameon=False)
    else:
        plt.plot(xs, data)
    if args.fitexp:
        fit_exp_nonlinear(ax, xs, data)
        plt.legend(frameon=False)
    if args.xlabel is not None:
        plt.xlabel(args.xlabel)
    if args.ylabel is not None:
        plt.ylabel(args.ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    if args.output is not None:
        fig.savefig(args.output, bbox_inches='tight')
    plt.show()

main()
