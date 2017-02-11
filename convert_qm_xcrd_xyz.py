import re

file_in = open('qm_xcrd.txt', 'r').readlines()
file_out = open('qm_xcrd.xyz', 'w')

array = []
for line in file_in:
    x = line.split()
    atom_type = ''
    if float(x[3]) > 0:
        atom_type = 'H'
    else:
        atom_type = 'O'
    array.append([atom_type, x[0], x[1], x[2]])

natoms = len(array)

file_out.write('{:>12}\n'.format(natoms))
file_out.write('MM atoms in QM cut\n')

coords = ''
for i in range(natoms):
    coords += "{:>2}{: 16.10f}{: 16.10f}{: 16.10f}\n".format(array[i][0], float(array[i][1]),
                                                             float(array[i][2]), float(array[i][3]))

file_out.write(coords)

