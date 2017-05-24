'''
Functions to parse amber output
'''
import re
import numpy as np

def find_dipoles(file):
    '''
    Retrun a numpy array of dipoles
    '''
    file_string = open(file, 'r').read()
    dipoles = re.compile(r"MM DIPOLE\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+.\d+\s+(\-?\d+\.\d+)")
    search_results = re.findall(dipoles, file_string)
    dipoles = np.zeros(len(search_results))
    for i, value in enumerate(search_results):
        dipoles[i] = float(value)
    return dipoles


def get_num_atoms_amber_out(nasqm_root):
    '''
    Retrun the number of atoms
    '''
    file_in = open(nasqm_root+'.out', 'r')
    p_natom = re.compile(r'NATOM\s*=*\s*(\d+)')
    for line in file_in:
        if re.search(p_natom, line):
            search_results = re.findall(p_natom, line)
            return int(search_results[0])


def find_nasqm_excited_state(input_stream, output_stream, n_states=1):
    '''
    Write the firt n_states excited states energies and strengths to the output stream
    FIXME What units?
    '''
    p_omega = re.compile(r'Frequencies \(eV\) and Oscillator')
    p_float = re.compile(r'-?\d+\.\d+E?-?\d*')
    energies = []
    strengths = []
    for line in input_stream:
        if re.search(p_omega, line):
            input_stream.readline()
            for state in range(n_states):
                line2 = input_stream.readline()
                search_results = re.findall(p_float, line2)
                energies.append(search_results[0])
                strengths.append(search_results[-1])
    n_steps = int(len(energies) / n_states)
    for step in range(n_steps):
        for state in range(n_states):
            index = n_states * step + state
            output_stream.write("{: 24.14E}{: 24.14E}".format(float(energies[index]),
                                                              float(strengths[index])))
        output_stream.write('\n')


def find_excited_energy(input_stream, output_stream, state):
    '''
    Write the total energies of the excited state
    FIXME what units?
    '''
    p_energy = re.compile(r'Total energies of excited states')
    p_float = re.compile(r'-?\d+\.\d+E?-?\d*')
    for line in input_stream:
        if re.search(p_energy, line):
            for state_value in range(state):
                line2 = input_stream.readline()
                search_results = re.findall(p_float, line2)
                if state_value == state - 1:
                    output_stream.write("{: 24.14E}".format(float(search_results[0])) + '\n')


def find_nasqm_transition_dipole(input_file, output_file, state=0):
    '''
    '''
    f = open(input_file, 'r').read()
    fo = open(output_file, 'w')
    p_energy_block = re.compile('Omega.*\n\s+\d\s+(-?\d+\.\d+E?-?\d*\s*){5}')
    p_float = '-?\d+\.\d+E?-?\d*'
    m = re.findall(p_energy_block, f)
    dipole_array = []
    for i in m:
        n = re.findall(p_float, i)
        dipole_array.append(float(n[0]))
    for i, d in enumerate(dipole_array):
        if i % 2 != 0:
            fo.write(str(d) + '\n')
