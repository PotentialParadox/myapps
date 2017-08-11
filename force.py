'''
Module to find forces interaction between two files
'''
import numpy as np

def find_differences(coord1, coord2):
    '''
    Returns the differences between each row of
    of coord1 and each row of coord2
    '''
    differences = []
    for i in coord1:
        for j in coord2:
            differences.append(j - i)
    return np.array(differences)

def find_forces(file_1, file_2):
    '''
    Returns the force vectors between two molecules decribed in two separate files
    formated in rows x, y, z, charge
    '''
    data_1 = np.loadtxt(file_1)
    data_2 = np.loadtxt(file_2)
    if data_1.ndim == 1:
        data_1 = np.expand_dims(data_1, axis=0)
    if data_2.ndim == 1:
        data_2 = np.expand_dims(data_2, axis=0)
    coords_1 = data_1[:, 0:3]
    charge_1 = data_1[:, 3]
    coords_2 = data_2[:, 0:3]
    charge_2 = data_2[:, 3]
    differences = find_differences(coords_1, coords_2)
    differences_2 = differences**2
    diff_sum_square = np.sum(differences_2, axis=1)
    distance = np.sqrt(diff_sum_square)
    charges = (np.tensordot(charge_1, charge_2, axes=0))
    charges = np.reshape(charges, (charges.shape[0]*charges.shape[1]))
    force_mags = charges / diff_sum_square
    unit_forces = differences / distance[:, None]
    # print(unit_forces)
    n1 = len(coords_1)
    n2 = len(coords_2)
    forces = np.multiply(unit_forces, force_mags[:, None])
    forces_1 = np.zeros([n1, 3])
    for i in range(n1):
        forces_1[i] = np.sum(forces[i*n2:i*n2+n2], axis=0)
    return forces_1
