import re
import argparse
import matplotlib.pyplot as plt

def main(filename, title):
    filestring = open(filename, 'r').read()
    times = get_times(filestring)
    coeffs = get_coeffs(filestring)
    graph_coeffs(times, coeffs, title)

def graph_coeffs(times, coeffs, title):
    x = times
    n_states = len(coeffs[0])
    for state in range(n_states):
        y = [coeff_block[state] for coeff_block in coeffs]
        plt.plot(x, y, label="state {}".format(state+1))
    plt.legend()
    plt.xlabel("Time (fs)")
    plt.title(title)
    plt.show()

def get_coeffs(filestring):
    '''
    returns [timestep][state coeff]
    '''
    filestring += '\n blank line'
    coeff_blocks  = filestring.split('$COEFF')[1:]
    coeffs = []
    for block in coeff_blocks:
        temp_coeffs = block.split('\n')
        temp_coeffs = temp_coeffs[1:-3]
        temp_coeffs = [float(x.split()[0]) for x in temp_coeffs]
        coeffs.append(temp_coeffs)
    return coeffs

def get_times(filestring):
    # find time =         (timevalue)
    p_time = re.compile("time.*(\s\d+\.\d+)")
    times = [float(x) for x in re.findall(p_time, filestring)]
    return times

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", default="coefficient.out")
parser.add_argument("-t", "--title", default="State Coefficients")
args = parser.parse_args()

main(args.filename, args.title)
