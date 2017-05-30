"""
Functions to automate the use of NAESMD
"""
import shutil
import subprocess
import os
from inputceon import get_xyz_coordinates
from amber_out import find_nasqm_excited_state
from amber_out import find_nasqm_transition_dipole

def run_naesmd(root_name):
    '''
    The primary call to naesmd
    '''
    file_out = open(root_name+'.out', 'w')
    subprocess.run([str(os.getenv('NAESMDHOME') + '/sqmceonaesmd.exe'), 'input.ceon'],
                   stdout=file_out)
    shutil.copy('energy-ev.out', root_name+'_energy.out')
    shutil.copy('coords.xyz', root_name+'_coords.xyz')
    shutil.copy('temperature.out', root_name+'_temp.out')
    subprocess.run('rm coordinates*', shell=True)


def run_naesmd_snapshots(root_name, n_coordinates, input_ceon):
    '''
    Run multiple naesmd calls
    '''
    file_in = open(root_name+'_coords.xyz', 'r')
    out_file = root_name + '_snapshots.out'
    file_out = open(out_file, 'w')
    for _ in range(n_coordinates):
        coordinates = get_xyz_coordinates(file_in)
        input_ceon.set_input(coordinates=coordinates)
        subprocess.run([str(os.getenv('NAESMDHOME') + '/sqmceonaesmd.exe'), 'input.ceon'],
                       stdout=file_out)
        find_nasqm_excited_state(file_in, file_out)
        find_nasqm_transition_dipole(file_in, file_out)
