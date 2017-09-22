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

    dipoles = find_dipoles(file_stream)
    dipoles = dipoles * DEBYE_TO_COULOMBMETER # C*m

    std = np.std(dipoles)
    std2 = std**2

    avg = np.average(dipoles)
    avg2 = avg**2

    dipole2 = std2 + avg2 # debye^2

    # Right hand side denominator
    rhsd = 9 * E0 * V_meters * KB * T_kelvon # C^2 * N^-1 * m * J
    # Righ hand side
    rhs = dipole2 / rhsd # N * m * J^-1 = 1 Good!

    # The equation can be rewriten in the form ax^2 + bx + x with
    a = 2
    b = -(1 + 9 * rhs)
    c = -1
    possible_answer = quadratic_formula(a, b, c)

    print(possible_answer[0])


FILE_NAME = 'md_qmmm_amb.out'
TEMPERATURE = 300
HEIGHT = 36.62E-10 # Meter
WIDTH = 37.693E-10 # Meter
LENGTH = 45.684E-10 # Meter
VOLUME = HEIGHT * WIDTH * LENGTH
calculate_dielectric(FILE_NAME, VOLUME, TEMPERATURE)
