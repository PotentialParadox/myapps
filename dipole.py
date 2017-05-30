'''
Calculate the dipole of a simple text file
dustin tracy (dtrac.uf@gmail.com)
'''
import numpy as np
from my_constants import AMBERCHARGE_TO_FRANKLIN, ANGSTROM_TO_CM, FRANKLINCM_TO_DEBYE


def calculate_dipole(file_name):
    '''
    Calculate the dipole of a file with formated as amber_charge charge, X, Y, Z
    '''
    data = np.loadtxt(file_name)
    charges = data[:, 0] * AMBERCHARGE_TO_FRANKLIN
    coords = data[:, 1:] * ANGSTROM_TO_CM
    dipoles = np.dot(charges[:, None].T, coords) * FRANKLINCM_TO_DEBYE
    total = np.linalg.norm(dipoles)

    answer_string = "Dipole: X= {: 10.4f}, Y= {: 10.4f}, Z= {: 10.4f}"\
                                " Total= {: 10.4f}".format(dipoles[0][0], dipoles[0][1],
                                                           dipoles[0][2], total)

    print(answer_string)

FILE_NAME = 'm2_coord_charges.txt'
calculate_dipole(FILE_NAME)
