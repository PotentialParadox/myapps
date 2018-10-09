'''
Created by Dustin Tracy (2015)
To use, run an AMBER md simulation (qm/mm) with the printdipole flag set to 2.
Apply scrip to the output file.
The Height Width and Volume can be found in the restart file at the bottom
'''
from pynasqm.amberout import find_dipoles
import numpy as np
from my_math import quadratic_formula
from my_constants import E0, KB, DEBYE_TO_COULOMBMETER, ANGSTROM_TO_METER

def calculate_dielectric(file_name, V_meters, T_kelvon):
    '''
    Estimates the dielectric of a box of solvent using
    equation 10 of 'A systematic study of water models for molecular simulation:
    Derivation of water models optimized for use with a reaction field' by
    David van der Spoel with E(0) = E(r)
    '''
    file_stream = open(file_name, 'r')

    dipoles = find_dipoles(file_stream) # Debye
    dipoles = dipoles * DEBYE_TO_COULOMBMETER # C*m

    s_exp_dipoles = np.power(np.average(dipoles),2)

    dipoles_2 = np.power(dipoles, 2)

    expval_dipoles_2 = np.average(dipoles_2)
    stdev = expval_dipoles_2 - s_exp_dipoles
    possible_answer = 1 + (expval_dipoles_2 - s_exp_dipoles) / (3.0 * E0 * V_meters * KB * T_kelvon)
    print(possible_answer)


FILE_NAME = 'nasqm_ground.out'
TEMPERATURE = 295
HEIGHT = 35.4E-10 # Meter
WIDTH = 39.7E-10 # Meter
LENGTH = 39.3E-10 # Meter
VOLUME = HEIGHT * WIDTH * LENGTH
calculate_dielectric(FILE_NAME, VOLUME, TEMPERATURE)
