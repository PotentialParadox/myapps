'''
Calculates the dipoles of an md simultions using a file containing charges
dustin tracy (dtrac.uf@gmail.com)
'''
import argparse
import re
import numpy as np
from my_constants import AMBERCHARGE_TO_FRANKLIN, ANGSTROM_TO_CM, FRANKLINCM_TO_DEBYE

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parmtop_file", "-p", help="Parmtop FILE")
    parser.add_argument("--restart_file", "-r", help="Non-binary restart file")
    args = parser.parse_args()

    charges = find_charges(args.parmtop_file)
    coordinates = find_coordinates(args.restart_file)
    calculate_dipole(charges, coordinates)

def find_charges(parmtop_file):
    p_charges = re.compile("CHARGE")
    p_flags = re.compile("FLAG")
    charge_block = ""
    with open(parmtop_file, 'r') as filein:
        for line in filein:
            if re.search(p_charges, line):
                filein.readline()
                for morelines in filein:
                    if re.search(p_flags, morelines):
                        break
                    charge_block += morelines + "\n"
    charge_list = [float(x) for x in charge_block.split()]
    return np.array(charge_list)

def find_coordinates(restart_file):
    coordinate_block = ""
    number_atoms = None
    with open(restart_file, 'r') as filein:
        filein.readline()
        second_line = filein.readline()
        number_atoms = int(second_line.split()[0])
        coordinate_block = ""
        for line in filein:
            coordinate_block += line + "\n"
    n_coords = number_atoms * 3
    coord_list = [float(x) for x in coordinate_block.split()[:n_coords]]
    return reshape_coordinates(coord_list)

def reshape_coordinates(coordinates):
    n_coords = int(len(coordinates) / 3)
    coords = []
    for coord in range(n_coords):
        atom_coord = []
        for dim in range(3):
            index = 3*coord + dim
            atom_coord.append(coordinates[index])
        coords.append(atom_coord)
    return np.array(coords)

def calculate_dipole(amber_charges, coordinates):
    '''
    Calculate the dipole of a file with formated as amber_charge charge, X, Y, Z
    '''
    charges = amber_charges * AMBERCHARGE_TO_FRANKLIN
    coords = coordinates * ANGSTROM_TO_CM
    dipoles = np.dot(charges[:, None].T, coords) * FRANKLINCM_TO_DEBYE
    total = np.linalg.norm(dipoles)

    answer_string = "Dipole: X= {: 10.4f}, Y= {: 10.4f}, Z= {: 10.4f}"\
                                " Total= {: 10.4f}".format(dipoles[0][0], dipoles[0][1],
                                                           dipoles[0][2], total)

    print(answer_string)


main()
