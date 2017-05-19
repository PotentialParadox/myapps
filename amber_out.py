import re
import math
import numpy as np

## Constants ##
e0 = 8.85418782E-12 #C^2*N-1*m-2
kb = 1.38064852E-23 #J*K-1
## Given parameters ##
file = 'md_qmm_amb.out'
T = 300 # K


def quadratic_formula(a, b, c):
    positive = (-b + math.sqrt(b**2 - 4*a*c))/2*a
    negative = (-b - math.sqrt(b**2 - 4*a*c))/2*a
    return (positive, negative)


def find_dipoles(file)
    file_string = open(file, 'r').read()
    dipoles = re.compile("MM DIPOLE\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+.\d+\s+(\-?\d+\.\d+)")
    m = re.findall(dipoles, file_string)
    dipoles = np.zeros(len(m))
    for i, value in enumerate(m):
        dipoles[i] = float(m[i])
    return dipoles

dipoles = find_dipoles(file)
dipoles = dipoles * 3.33564E-30 # C*m

std = np.std(dipoles)
std2 = std**2

avg = np.average(dipoles)
avg2 = avg**2

m2 = std2 + avg2 # debye^2

h = 29.175 # Angstrom
w = 30.179 # Angstrom
l = 36.468 # Angstrom
V = h*w*l * (10E-10)**3 # m^3

rhsd = 9 * e0 * V * kb * T # C^2 * N^-1 * m * J
rhs = m2 / rhsd # N * m * J^-1 = 1 Good!

a = 2
b = -(1 + 9 * rhs)
c = -1
possible_answer = quadratic_formula(a, b, c)


print(possible_answer[0])
