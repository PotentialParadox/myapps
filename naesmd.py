import sys
import time
import shutil
import re
import subprocess
import os
import numpy as np
from sed import sed_inplace, sed_global
from periodic_table import periodic_table


def get_xyz_coordinates(file_stream):
    number_atoms = int(file_stream.readline())
    file_stream.readline()
    coords = ""
    for i in range(number_atoms):
        l_coords = file_stream.readline().split()
        try:
            l_coords[0] = periodic_table(l_coords[0])
        except KeyError:
            l_coords[0] = l_coords[0]
        coords += "{:>2}{: 16.10f}{: 16.10f}{: 16.10f}\n".format(l_coords[0], float(l_coords[1]),
                                                                 float(l_coords[2]), float(l_coords[3]))
    return coords


class InputCeon:
    def __init__(self, n_steps=None, n_exc_states_propagate=None, n_steps_to_print=None, exc_state_init=None,
                 coordinates=None, verbosity=None, periodic=None):
        self.log = ''
        self.set_input(n_steps, n_exc_states_propagate, n_steps_to_print, exc_state_init,
                       coordinates, verbosity, periodic)

    def set_input(self, n_steps=None, n_exc_states_propagate=None, n_steps_to_print=None, exc_state_init=None,
                  coordinates=None, verbosity=None, periodic=None):
        if n_steps is not None:
            sed_inplace('input.ceon', 'n_class_steps=\d*', 'n_class_steps=' + str(n_steps))
            sed_inplace('md_qmmm_amb.in', 'nstlim\s*=\s*\d*', 'nstlim='+str(n_steps))
        if n_exc_states_propagate is not None:
            sed_inplace('input.ceon', 'n_exc_states_propagate=\d*', 'n_exc_states_propagate='
                        + str(n_exc_states_propagate))
        if n_steps_to_print is not None:
            sed_inplace('input.ceon', 'out_data_steps=\d*', 'out_data_steps=' + str(n_steps_to_print))
            sed_inplace('md_qmmm_amb.in', 'ntwx=\s*\d*', 'ntwx=' + str(n_steps_to_print))
        if coordinates is not None:
            sed_global('input.ceon', '&coord(.|\s|\n)*?&endcoord', '&coord\n' + coordinates + '&endcoord')
        if exc_state_init is not None:
            sed_inplace('input.ceon', 'exc_state_init=\d*', 'exc_state_init=' + str(exc_state_init))
        if verbosity is not None:
            sed_inplace('input.ceon', 'verbosity=\d*', 'verbosity=' + str(verbosity))
            sed_inplace('md_qmmm_amb.in', 'verbosity\s*=\s*\d*', 'verbosity=' + str(verbosity))
        if periodic is True:
            sed_inplace('md_qmmm_amb.in', 'ntb\s*=\s*\d+', 'ntb=2')
            sed_inplace('md_qmmm_amb.in', 'ntp\s*=\s*\d+', 'ntp=1')
        if periodic is False:
            sed_inplace('md_qmmm_amb.in', 'ntb\s*=\s*\d+', 'ntb=0')
            sed_inplace('md_qmmm_amb.in', 'ntp\s*=\s*\d+', 'ntp=0')
        if n_steps == 0:
            sed_inplace('md_qmmm_amb.in', 'irest\s*=\s*\d+', 'irest=0')
            sed_inplace('md_qmmm_amb.in', 'ntx\s*=\s*\d+', 'ntx=1')
        self.log += open('input.ceon', 'r').read()

    def write_log(self):
        open('input_ceon.log', 'w').write(self.log)


def run_naesmd(root_name):
    file_out = open(root_name+'.out', 'w')
    subprocess.run([str(os.getenv('NAESMDHOME') + '/sqmceonaesmd.exe'), 'input.ceon'], stdout=file_out)
    shutil.copy('energy-ev.out', root_name+'_energy.out')
    shutil.copy('coords.xyz', root_name+'_coords.xyz')
    shutil.copy('temperature.out', root_name+'_temp.out')
    subprocess.run('rm coordinates*', shell=True)


def run_nasqm(root_name, coordinate_file=None):
    cd_file = 'm1.inpcrd'
    if coordinate_file is not None:
        cd_file = coordinate_file
    subprocess.run(['sander', '-O', '-i', 'md_qmmm_amb.in', '-o', root_name+'.out', '-c', cd_file, '-p',
                    'm1.prmtop', '-r', root_name+'.rst', '-x', root_name+'.nc'])


def run_naesmd_snapshots(root_name, n_coordinates, input_ceon):
    file_in = open(root_name+'_coords.xyz', 'r')
    out_file = root_name + '_snapshots.out'
    file_out = open(out_file, 'w')
    for i in range(n_coordinates):
        coordinates = get_xyz_coordinates(file_in)
        input_ceon.set_input(coordinates=coordinates)
        subprocess.run([str(os.getenv('NAESMDHOME') + '/sqmceonaesmd.exe'), 'input.ceon'], stdout=file_out)
        find_nasqm_excited_state(out_file, root_name+'_omegas')
        find_nasqm_transition_dipole(out_file, root_name+'_dipoles')


def find_nasqm_excited_state(input_file, output_file, n_states=1):
    f = open(input_file, 'r')
    out_file = output_file
    fo = open(out_file, 'w')
    p_omega = re.compile('Frequencies \(eV\) and Oscillator')
    p_float = '-?\d+\.\d+E?-?\d*'
    energies = []
    strengths = []
    for line in f:
        if re.search(p_omega, line):
            f.readline()
            for state in range(n_states):
                line2 = f.readline()
                m = re.findall(p_float, line2)
                energies.append(m[0])
                strengths.append(m[-1])

    n_steps = int(len(energies) / n_states)
    for step in range(n_steps):
        for state in range(n_states):
            index = n_states * step + state
            fo.write("{: 24.14E}{: 24.14E}".format(float(energies[index]), float(strengths[index])))
        fo.write('\n')


def find_nasqm_transition_dipole(input_file, output_file, state=0):
    f = open(input_file, 'r').read()
    fo = open(output_file, 'w')
    p_energy_block = re.compile('Omega.*\n\s+\d\s+(-?\d+\.\d+E?-?\d*\s*){5}')
    p_float = '-?\d+\.\d+E?-?\d*'
    m = re.findall(p_energy_block, f)
    dipole_array = []
    for i in m:
        n = re.findall(p_float, i)
        dipole_array.append(float(n[0]))
    for i, d in enumerate(dipole_array):
        if i % 2 != 0:
            fo.write(str(d) + '\n')


def create_convert_to_crd():
    ctc = 'parm m1.prmtop\n' \
          'trajin nasqm_ground.nc 1 last\n' \
          'center :1\n' \
          'trajout output.crd restart\n' \
          'run\n' \
          'quit'
    open('convert_to_crd.in', 'w').write(ctc)


def get_num_atoms_amber_out(nasqm_root):
    file_in = open(nasqm_root+'.out', 'r')
    p_natom = re.compile('NATOM\s*=*\s*(\d+)')
    for line in file_in:
        if re.search(p_natom, line):
            m = re.findall(p_natom, line)
            return int(m[0])


def get_atom_types_prmtop(nasqm_root):
    number_atoms = get_num_atoms_amber_out(nasqm_root)
    n_lines = int(number_atoms / 20)
    remainder = number_atoms % 20
    if remainder > 0:
        n_lines += 1
    file_in = open('m1.prmtop', 'r')
    p_atom_flag = re.compile('%FLAG ATOM_NAME')
    p_atom_type = re.compile('([A-Z]+)\d+')
    a_types = []
    for line in file_in:
        if re.search(p_atom_flag, line):
            file_in.readline()
            for _ in range(n_lines):
                line2 = file_in.readline()
                m = re.findall(p_atom_type, line2)
                a_types.extend(m)
    return a_types


def create_nasqm_restart(nasqm_root):
    create_convert_to_crd()
    sed_inplace('convert_to_crd.in', 'trajin .*\.nc', 'trajin '+nasqm_root+'.nc')
    subprocess.run(['cpptraj', '-i', 'convert_to_crd.in'])


def run_nasqm_snapshots(nasqm_root, n_coordinates):
    create_nasqm_restart(nasqm_root)
    out_file = 'nasqm_ground_snapshots.out'
    temp_file = 'nasqm_ground_snap_temp'
    file_out = open(out_file, 'w')
    for i in range(n_coordinates):
        coordinate_file = 'output.crd.'+str(i+1)
        run_nasqm(temp_file, coordinate_file=coordinate_file)
        temp_in = open('nasqm_ground_snap_temp.out', 'r').read()
        file_out.write(temp_in)
    subprocess.run('rm nasqm_ground_snap_temp*', shell=True)
    find_nasqm_excited_state(out_file, nasqm_root+'_omegas.txt')
    find_nasqm_transition_dipole(out_file, nasqm_root+'_dipoles.txt')
    subprocess.run('rm output.crd.*', shell=True)


def set_inpcrd(coordinates):
    temp_list = coordinates.split()
    c_list = []
    for i, c in enumerate(temp_list):
        if i % 4 != 0:
            c_list.append(float(c))
    n_atoms = int(len(c_list) / 3)
    inpcrd = "MOL\n"
    inpcrd += "    " + str(n_atoms) + "\n"
    n_full_lines = int(n_atoms / 2)
    for i in range(n_full_lines):
        inpcrd += '{: 12.7f}{: 12.7f}{: 12.7f}{: 12.7f}{: 12.7f}{: 12.7f}'.format(c_list[i*6], c_list[i*6+1],
                                                                                  c_list[i*6+2], c_list[i*6+3],
                                                                                  c_list[i*6+4], c_list[i*6+5])
        inpcrd += '\n'
    open('m1.inpcrd', 'w').write(inpcrd)


def main(n_steps):

    is_qmmm = True
    is_box = True

    # Copy inputs
    input_ceon_bac = open('input.ceon', 'r').read()
    md_qmmm_amb = open('md_qmmm_amb.in', 'r').read()
    m1_inpcrd = open('m1.inpcrd', 'r').read()

    n_steps_to_print = 50
    n_coordinates = int(n_steps / n_steps_to_print)

    start_time = time.time()

    # Run the Ground State
    n_exc_states_propagate = 0
    exc_state_init = 0
    verbosity = 0
    input_ceon = InputCeon(n_steps, n_exc_states_propagate, n_steps_to_print, exc_state_init, verbosity=verbosity)

    # NASQM
    coordinate_file = None
    if is_qmmm:
        coordinate_file = 'm1_md2.rst'
    # run_nasqm('nasqm_ground', coordinate_file=coordinate_file)
    # NAESMD
    # run_naesmd('naesmd_ground')

    # Now we want to take the geometry snapshots and run single-point calculations
    # on these to get the S1 excited states
    n_exc_states_propagate = 20
    n_steps = 0
    n_steps_to_print = 1
    exc_state_init = 0
    verbosity = 3
    input_ceon.set_input(n_steps, n_exc_states_propagate, n_steps_to_print, exc_state_init, verbosity=verbosity)
    # run_naesmd_snapshots('naesmd_ground', n_coordinates, input_ceon)
    run_nasqm_snapshots('nasqm_ground', n_coordinates)

    # Restore Original Inputs
    open('input.ceon', 'w').write(input_ceon_bac)
    open('md_qmmm_amb.in', 'w').write(md_qmmm_amb)
    open('m1.inpcrd', 'w').write(m1_inpcrd)
    input_ceon.write_log()

    end_time = time.time()
    print("Job finished in %s seconds" % (end_time - start_time))

# n_steps = int(sys.argv[1])
# main(n_steps)
find_nasqm_excited_state('nasqm_ground_snapshots.out', 'nasqm_excited_omegas.txt', n_states=2)

