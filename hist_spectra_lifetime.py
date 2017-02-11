import numpy as np
from math import exp

# Change here to the number of excited states propagated
number_states = 2
# Change here if you want more points per gaussian ?
number_gauss = 100
# Change the parameter for full width at half max in eV
full_width_half_max = 0.01
# Change here of bins you wish to distribute over
number_bins = 1000
# Change here the minimum energy in (eV)
x_minimum = 2.0
# Change here the maximum energy in (eV)
x_maximum = 5
# Change here the input file
file_input = 'nasqm_excited_omegas.txt'
# Change here the output file
file_output = 'nasqm_excited_spectra.txt'


def get_fwhm_energy_value(energy, fwhm, gauss_index, n_gauss):
    return energy - (2.0 * fwhm) + (gauss_index * 4.0 * fwhm / n_gauss)


def calculate_exponential(energy, strength, sigma, x_value_of_half_width):
    w = strength
    u = energy
    s = sigma
    x = x_value_of_half_width
    return (w / u**2) * exp((-(x - u)**2) / (2 * s**2))


def print_line(fo, ev_energy, nm_energy, weighted_state_intensities, total_intensity):
    line = "{: 8.3f}{: 18.10f}".format(ev_energy, nm_energy)
    for intensity in weighted_state_intensities:
        line += "{: 18.10f}".format(intensity)
    line += "{: 18.10f}\n".format(total_intensity)
    fo.write(line)


def calculate_spectra(n_states, n_gauss, fwhm, n_bins, x_min, x_max, file_in, file_out):

    divisor = 0.02

    sigma = fwhm / 2.35482

    # Rows will correspond to time steps while columns 0, 1
    # are frequency and strengths
    data = np.loadtxt(file_in)

    strengths = np.zeros((n_states, n_bins))
    total_strength = 0

    for state in range(n_states):
        for time_step in range(int(len(data)/n_states)):
            energy = data[time_step + state][0]
            strength = data[time_step + state][1]
            for gauss_point in range(n_gauss):
                energy_at_lower_fwhm = get_fwhm_energy_value(energy, fwhm, gauss_point, n_gauss)
                normalized_strength = calculate_exponential(energy, strength, sigma, energy_at_lower_fwhm)
                x_index = int(energy_at_lower_fwhm/divisor) + 1
                strengths[state, x_index] += normalized_strength
                total_strength += normalized_strength

    fo = open(file_out, 'w')
    for x_index in range(int(x_min/divisor), int(x_max/divisor)):
        ev_energy = x_index * divisor - divisor / 2
        nm_energy = 1E7 / (ev_energy * 8.06554465E3)
        weighted_state_intensities = np.divide(strengths[:, x_index], total_strength)
        total_intensity = np.sum(strengths[:, x_index]) / total_strength
        print_line(fo, ev_energy, nm_energy, weighted_state_intensities, total_intensity)


calculate_spectra(number_states, number_gauss, full_width_half_max, number_bins, x_minimum, x_maximum,
                  file_input, file_output)

