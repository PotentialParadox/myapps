'''
Created by Dustin Tracy (2015)
To use, run an AMBER md simulation (qm/mm) with the printdipole flag set to 2.
Apply scrip to the output file.
The Height Width and Volume can be found in the restart file at the bottom
'''
import argparse
from python_scripts.libdielectric import calculate_dielectric
from python_scripts.libdielectric import product
from python_scripts.libdielectric import find_box
from python_scripts.libdielectric import get_dipoles
from python_scripts.libdielectric import convert_to_debye
from python_scripts.libdielectric import a_to_m


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
