from functools import reduce
from subprocess import Popen, PIPE, STDOUT
from pynasqm.amberout import find_box
import numpy as np
import operator
from python_scripts.libmymath import quadratic_formula
from python_scripts.libmyconstants import E0, KB, AE_TO_COULOMBMETER, ANGSTROM_TO_METER, AE_TO_DEBYE
import pytraj as pt

def product(xs):
    return reduce(operator.mul, xs, 1)

def a_to_m(x):
    return x * 1E-10

def average_dipole(dipoles):
    return np.average(dipoles, axis=0)

def dipole_magnitudes(dipoles):
    return np.array([np.linalg.norm(x) for x in dipoles])

def get_dipoles(parmtop, trajectory_file):
    return pt.analysis.vector.dipole(trajectory_file, top=parmtop)

def convert_to_debye(xs):
    return [[y*AE_TO_COULOMBMETER for y in x] for x in xs]

def calculate_dielectric(dipoles, V_meters, T_kelvon):
    '''Estimates the dielectric of a box of solvent using "On the accurate
    calculation of the dielectric constant and the diffusion coefficient from
    molecular dynamics simulations: the case of SPC/E water Units are (debye,
    m^3, K)
    [dx,dy,dz] -> v -> t -> io ()
    '''
    s_exp_dipoles = np.dot(average_dipole(dipoles),average_dipole(dipoles))
    expval_dipoles_2 = np.power(np.average(dipole_magnitudes(dipoles)),2)
    possible_answer = 1 + (expval_dipoles_2 - s_exp_dipoles) / (3.0 * E0 * V_meters * KB * T_kelvon)
    print(possible_answer)
