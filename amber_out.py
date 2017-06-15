'''
Functions to parse amber output
'''
import re
import io
import numpy as np

def find_dipoles(file_stream):
    '''
    Return a numpy array of dipoles
    '''
    file_string = file_stream.read()
    dipoles = re.compile(r"MM DIPOLE\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+.\d+\s+(\-?\d+\.\d+)")
    search_results = re.findall(dipoles, file_string)
    dipoles = np.zeros(len(search_results))
    for i, value in enumerate(search_results):
        dipoles[i] = float(value)
    return dipoles

def read_nasqm_excited_states(input_stream, n_states):
    '''
    Return a tupple of of energies and strenghts
    read from an amber output file
    '''
    p_omega = re.compile(r'Frequencies \(eV\) and Oscillator')
    p_float = re.compile(r'-?\d+\.\d+E?-?\d*')
    energies = []
    strengths = []
    for line in input_stream:
        if re.search(p_omega, line):
            input_stream.readline()
            for _ in range(n_states):
                line2 = input_stream.readline()
                search_results = re.findall(p_float, line2)
                energies.append(search_results[0])
                strengths.append(search_results[-1])
    return np.array(energies), np.array(strengths)

def create_spectra_string(output_stream, energies, strengths, n_states):
    '''
    Prints and returns a formated string using the appropriate inputs
    '''
    if not output_stream:
        output_stream = io.StringIO()
    n_steps = int(len(energies) / n_states)
    # NAESMD will run twice on the first iteration of a md simulation
    # we only want to count the second one. During singlepoint calculations
    # we don't need to worry about this.
    begin_step = 1 if n_steps > 1 else 0
    for step in range(begin_step, n_steps):
        for state in range(n_states):
            index = n_states * step + state
            output_stream.write("{: 24.14E}{: 24.14E}".format(float(energies[index]),
                                                              float(strengths[index])))
        output_stream.write('\n')
    return output_stream.getvalue()


def find_nasqm_excited_state(input_stream, output_stream=None, n_states=1):
    '''
    Write the firt n_states excited states energies and strengths to the output stream
    in eV
    '''
    energies, strengths = read_nasqm_excited_states(input_stream, n_states)
    return create_spectra_string(output_stream, energies, strengths, n_states)


def find_excited_energy(input_stream, output_stream=None, state=1):
    '''
    Write the total energies of the excited state in eV
    '''
    p_energy = re.compile('Total energies of excited states')
    p_float = re.compile(r'-?\d+\.\d+E?-?\d*')
    for line in input_stream:
        if re.search(p_energy, line):
            for state_value in range(state):
                line2 = input_stream.readline()
                search_results = re.findall(p_float, line2)
                if state_value == state - 1:
                    output_stream.write("{: 24.14E}".format(float(search_results[0])) + '\n')


def find_nasqm_transition_dipole(input_stream, output_stream=None):
    '''
    FIXME Write the transition dipoles to the output stream
    '''
    if not output_stream:
        output_stream = io.StringIO()
    p_energy_block = re.compile(r'Omega.*\n\s+\d\s+(-?\d+\.\d+E?-?\d*\s*){5}')
    p_float = re.compile(r'-?\d+\.\d+E?-?\d*')
    search_results = re.findall(p_energy_block, input_stream)
    dipole_array = []
    for i in search_results:
        other_results = re.findall(p_float, i)
        dipole_array.append(float(other_results[0]))
    for i, dipole in enumerate(dipole_array):
        if i % 2 != 0:
            output_stream.write(str(dipole) + '\n')
    return output_stream.getvalue()

def find_total_energies(input_stream):
    '''
    Find the total energies in the amber output file
    and return a numpy array
    '''
    p_energies = re.compile(r"Etot\s*=\s+(\-?\d+\.\d+)")
    file_string = input_stream.read()
    search_results = re.findall(p_energies, file_string)
    energy_list = []
    for result in search_results:
        energy_list.append(float(result))
    return np.array(energy_list[:-2])
