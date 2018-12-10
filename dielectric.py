'''
Created by Dustin Tracy (2015)
To use, run an AMBER md simulation (qm/mm) with the printdipole flag set to 2.
Apply scrip to the output file.
The Height Width and Volume can be found in the restart file at the bottom
'''
from functools import reduce
from subprocess import Popen, PIPE, STDOUT
from pynasqm.amberout import find_box
import numpy as np
import argparse
import operator
from my_math import quadratic_formula
from my_constants import E0, KB, AE_TO_COULOMBMETER, ANGSTROM_TO_METER, AE_TO_DEBYE

def product(xs):
    return reduce(operator.mul, xs, 1)

def a_to_m(x):
    return x * 1E-10

def average_dipole(dipoles):
    return np.average(dipoles, axis=0)

def dipole_magnitudes(dipoles):
    return np.array([np.linalg.norm(x) for x in dipoles])

def get_dipoles(parmtop, trajectory_file):
    command = bytes("parm {}\n"\
                    "trajin {}\n"\
                    "vector ve1 (:*) dipole out dipole.out\n"\
                    "run\n"\
                    "quit\n".format(parmtop, trajectory_file), 'utf-8')
    proc = Popen(["cpptraj"], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    cpptraj_stdout = proc.communicate(input=command)[0]
    lines = open('dipole.out','r').readlines()
    dipole_file = "".join(lines[1::100])
    data = np.fromstring(dipole_file, sep=" ")
    data = np.reshape(data, (-1,7))
    return data[:,1:4]

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

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--amberout", help="Amber output file", default="nasqm_ground.out")
parser.add_argument("-x", "--trajfile",help="Amber's trajectory file", default="nasqm_ground.nc")
parser.add_argument("-p", "--prmtop", help="Amber's top file", default="m1.prmtop")
args = parser.parse_args()

DIPOLES = convert_to_debye(get_dipoles(args.prmtop, args.trajfile)) # Debye
volume = product([a_to_m(x) for x in find_box(open('nasqm_ground.out','r'))]) # m^3
temperature=300 # Kelvin
begin = 1
d = [calculate_dielectric(DIPOLES[:x], volume, temperature) for x in range(begin,len(DIPOLES))]
