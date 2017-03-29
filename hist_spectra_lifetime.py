from calculate_spectra import calculate_spectra

# Change here to the number of excited states propagated
number_states = 5
# Change here the file root name
root_name = 'spectra_abs'
# Change here if you want more points per gaussian ?
number_gauss = 100
# Change the parameter for full width at half max in eV
full_width_half_max = 0.05
# Change here of bins you wish to distribute over
number_bins = 1000
# Change here the minimum energy in (eV)
x_minimum = 2.5
# Change here the maximum energy in (eV)
x_maximum = 4.1
# Change here the input file
file_input = root_name +'.input'
# Change here the output file
file_output = root_name + '.output'


calculate_spectra(number_states, number_gauss, full_width_half_max, number_bins, x_minimum, x_maximum,
                  file_input, file_output)

