import sys
import time
import shutil
import re
import subprocess
import os
import numpy as np
from sed import sed_inplace, sed_global
from periodic_table import periodic_table
from amber import run_amber_parallel, run_hpc_trajectories, create_snapshot_slurm_script


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
                  coordinates=None, verbosity=None, periodic=None, time_step=None):
        if n_steps is not None:
            sed_inplace('input.ceon', 'n_class_steps=\d+', 'n_class_steps=' + str(n_steps))
            sed_inplace('md_qmmm_amb.in', 'nstlim\s*=\s*\d+\.?\d*', 'nstlim='+str(n_steps))
        if n_exc_states_propagate is not None:
            sed_inplace('input.ceon', 'n_exc_states_propagate=\d*', 'n_exc_states_propagate='
                        + str(n_exc_states_propagate))
        if n_steps_to_print is not None:
            sed_inplace('input.ceon', 'out_data_steps=\d*', 'out_data_steps=' + str(n_steps_to_print))
            sed_inplace('md_qmmm_amb.in', 'ntwx=\s*\d*', 'ntwx=' + str(n_steps_to_print))
            sed_inplace('md_qmmm_amb.in', 'ntpr=\s*\d*', 'ntpr=' + str(n_steps_to_print))
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
        if time_step is not None:
            sed_inplace('md_qmmm_amb.in', 'dt=\s*\d+\.?\d*', 'dt=' +str(time_step/1000))
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


def run_nasqm(root_name, coordinate_file=None, pmemd_available=False):
    cd_file = 'm1.inpcrd'
    if coordinate_file is not None:
        cd_file = coordinate_file
    amber = 'sander'
    if pmemd_available:
        amber = 'pmemd.cuda'
    subprocess.run([amber, '-O', '-i', 'md_qmmm_amb.in', '-o', root_name+'.out', '-c', cd_file, '-p',
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


def find_nasqm_excited_state(input_stream, output_stream, n_states=1):
    f = input_stream
    fo = output_stream
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


def find_excited_energy(input_stream, output_stream, state):
    p_energy = re.compile('Total energies of excited states')
    p_float = '-?\d+\.\d+E?-?\d*'
    for line in input_stream:
        if re.search(p_energy, line):
            for s in range(state):
                line2 = input_stream.readline()
                m = re.findall(p_float, line2)
                if s == state - 1:
                    output_stream.write("{: 24.14E}".format(float(m[0])) + '\n')


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


def create_restarts(input, output, step=None):
    ctc = 'parm m1.prmtop\n'
    if step is not None:
        ctc += 'trajin ' + input + '.nc' + ' 1 last ' + str(step) + "\n"
    else:
        ctc += 'trajin ' + input + '.nc' + ' 1 last\n'
    ctc += 'center :1\n' \
           'trajout ' + output + ' restart\n' \
           'image\n' \
           'run\n' \
           'quit'
    open('convert_to_crd.in', 'w').write(ctc)
    subprocess.run(['cpptraj', '-i', 'convert_to_crd.in'])


def run_ground_state_snapshots(nasqm_root, output_root, n_coordinates, n_snapshots, is_hpc, pmemd_available, ppn):
    restart_step = int(n_coordinates / n_snapshots)
    create_restarts(input=nasqm_root, output='ground_snap', step=restart_step)
    snap_restarts = []
    snap_trajectories = []
    if n_snapshots == 1:
        snap_restarts.append("ground_snap")
        snap_trajectories.append(output_root + '1')
    else:
        for i in range(n_snapshots):
            snap_restarts.append("ground_snap."+str(i+1))
            snap_trajectories.append(output_root + str(i + 1))
    if is_hpc:
        run_hpc_trajectories(n_trajectories=n_snapshots, n_processor_per_node=ppn,
                             root_name=output_root)
    else:
        run_amber_parallel(pmemd_available, snap_trajectories, snap_restarts, number_processors=ppn)


def run_abs_snapshots(output_root, n_trajectories, n_frames, is_hpc):
    pmemd_available = False # We require ESMD
    for i in range(n_trajectories):
        amber_restart = 'nasqm_abs_' + str(i+1)
        create_restarts(input=amber_restart, output=amber_restart)
    # We will now have files that look like nasqm_abs_[trajectory]_[frame]
    # Lets run it
    snap_singles = []
    snap_restarts = []
    for traj in range(n_trajectories):
        for frame in range(n_frames):
            snap_singles.append("nasqm_abs_" + str(traj+1) + "_" + str(frame+1))
            snap_restarts.append("nasqm_abs_" + str(traj+1) + "." + str(frame+1))
    if is_hpc:
        create_snapshot_slurm_script(script_file_name='run_abs_snapshot.sbatch', n_trajectories=n_trajectories,
                                     n_frames=n_frames, root_name=output_root, crd_file=output_root)
        print("Please now submit run_abs_snapshots.sbatch, then run abs collection")
        sys.exit('Slurm submission exception')
    else:
        run_amber_parallel(pmemd_available, snap_singles, snap_restarts, number_processors=8)


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


def accumulate_flu_spectra(n_trajectories, n_states=1):
    output_stream = open('spectra_flu.input', 'w')
    for i in range(n_trajectories):
        amber_outfile = 'nasqm_flu_' + str(i+1) + ".out"
        input_stream = open(amber_outfile, 'r')
        find_nasqm_excited_state(input_stream, output_stream, n_states=n_states)
    output_stream.close()
    # We may also want the average omega_1s over time
    average_omegas_time = open('omega_1_time.txt', 'w')
    data = np.loadtxt('spectra_flu.input')
    n_rows_per_trajectory = int(data.shape[0] / n_trajectories)
    for i in range(n_rows_per_trajectory):
        omega = np.average(data[i::n_rows_per_trajectory, 0])
        average_omegas_time.write(str(omega) + '\n')
    average_omegas_time.close()
    # FIXME Need to split this up
    output_stream = open('nasqm_flu_energies.txt', 'w')
    for i in range(n_trajectories):
        amber_outfile = 'nasqm_flu_' + str(i+1) + ".out"
        input_stream = open(amber_outfile, 'r')
        find_excited_energy(input_stream, output_stream, n_states)
    output_stream.close()
    average_energies_time = open('nasqm_flu_energy_time.txt', 'w')
    data = np.loadtxt('nasqm_flu_energies.txt')
    subprocess.run(['rm', 'nasqm_flu_energies.txt'])
    n_rows_per_trajectory = int(data.shape[0] / n_trajectories)
    for i in range(n_rows_per_trajectory):
        e = np.average(data[i::n_rows_per_trajectory])
        average_energies_time.write(str(e) + '\n')


def accumulate_abs_spectra(n_trajectories, n_frames, n_states=20):
    output_stream = open('spectra_abs.input', 'w')
    for traj in range(n_trajectories):
        for frame in range(n_frames):
            amber_out = 'nasqm_abs_' + str(traj+1) + '_' + str(frame+1) + '.out'
            input_stream = open(amber_out, 'r')
            find_nasqm_excited_state(input_stream, output_stream, n_states)
    output_stream.close()


def clean_up_abs(n_trajectories, n_frame):
    base_name = 'nasqm_abs_'
    for i in range(n_trajectories):
        for j in range(n_frame):
            traj = i + 1
            frame= j + 1
            subprocess.run(['rm', base_name + str(traj) + '_' + str(frame) + '.out'])
            subprocess.run(['rm', base_name + str(traj) + '_' + str(frame) + '.nc'])
            subprocess.run(['rm', base_name + str(traj) + '_' + str(frame) + '.rst'])
            subprocess.run(['rm', base_name + str(traj) + '.' + str(frame)])



def main():
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Begin Inputs
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    is_qmmm = True
    is_hpc = False
    processor_per_node = 8
    run_ground_dynamics = False
    run_absorption_trajectories = False
    run_absorption_snapshots = True
    run_absorption_collection = True
    run_excited_state = False
    run_fluorescence_collection = False

    # Change here the number of snapshots you wish to take
    # from the initial ground state trajectory to run the
    # further ground state dynamics
    n_snapshots_gs = 8

    # Change here the number of snapshots you wish to take
    # from the initial ground state trajectory to run the
    # new excited state dynamics
    n_snapshots_ex = 8

    # Change here the time step that will be shared by
    # each trajectory
    time_step = 0.2  # fs

    # Change here the runtime of the initial ground state MD
    ground_state_run_time = 300 # ps

    # Change here the runtime for the the trajectories
    # used to create calculated the absorption
    abs_run_time = 1 # ps

    # Change here the runtime for the the trajectories
    # used to create calculated the fluorescence
    exc_run_time = 1 # ps

    # Change here the number of excited states you 
    # with to have in the CIS calculation
    n_exc_states_propagate_ex_param = 1

    # Change here the initial state
    exc_state_init_ex_param = 1

    # PMEMD NOT YET WORKING SET TO FALSE
    pmemd_available = False

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # End Inputs
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Copy inputs
    input_ceon_bac = open('input.ceon', 'r').read()
    md_qmmm_amb = open('md_qmmm_amb.in', 'r').read()
    m1_inpcrd = open('m1.inpcrd', 'r').read()

    # Create the input_ceon object
    input_ceon = InputCeon()

    n_steps_gs = int(ground_state_run_time / time_step * 1000)
    n_steps_to_print_gs = 50
    n_frames_gs = int(n_steps_gs / n_steps_to_print_gs)

    n_steps_abs = int(abs_run_time / time_step * 1000)
    # We will do absorption calculation on all 
    # steps printed out, so 1 would do absorption
    # for each step during the run_abs_snapshot step
    n_steps_to_print_abs = 50 
    n_frames_abs = int(n_steps_abs / n_steps_to_print_abs)

    n_steps_exc = int(exc_run_time / time_step * 1000)
    n_steps_to_print_exc = 1
    n_frames_exc = int(n_steps_exc / n_steps_to_print_exc)

    start_time = time.time()

    # Run the Ground State
    if run_ground_dynamics:
        print("!!!!!!!!!!!!!!!!!!!! Running Ground Dynamics !!!!!!!!!!!!!!!!!!!!")
        n_exc_states_propagate = 0
        exc_state_init = 0
        verbosity = 0
        input_ceon.set_input(n_steps_gs, n_exc_states_propagate, n_steps_to_print_gs, exc_state_init, verbosity=verbosity,
                             time_step=time_step)
        coordinate_file = None
        if is_qmmm:
            coordinate_file = 'm1_md2.rst'
        run_nasqm('nasqm_ground', coordinate_file=coordinate_file, pmemd_available=pmemd_available)
    if run_absorption_trajectories:
        print("!!!!!!!!!!!!!!!!!!!! Running Absorbance Trajectories !!!!!!!!!!!!!!!!!!!!")
        # Now we want to take the original trajectory snapshots and run more trajectories
        # from those at the ground state
        n_exc_states_propagate = 0
        exc_state_init = 0
        verbosity = 0
        input_ceon.set_input(n_steps_abs, n_exc_states_propagate, n_steps_to_print_abs, exc_state_init, verbosity=verbosity,
                             time_step=time_step)
        run_ground_state_snapshots('nasqm_ground', 'nasqm_abs_', n_frames_gs, n_snapshots_gs, is_hpc,
                                   pmemd_available=pmemd_available, ppn=processor_per_node)
    if run_absorption_snapshots:
        print("!!!!!!!!!!!!!!!!!!!! Running Absorbance Snapshots !!!!!!!!!!!!!!!!!!!!")
        # Once the ground state trajectory files are made, we need
        # to calculate snapshots the Si - S0 energies
        n_exc_states_propagate = 20
        exc_state_init = 0
        verbosity = 3
        n_steps = 0
        input_ceon.set_input(n_steps, n_exc_states_propagate, n_steps_to_print_abs, exc_state_init, verbosity=verbosity,
                             time_step=time_step)
        run_abs_snapshots(output_root='nasqm_abs_', n_trajectories=n_snapshots_gs, n_frames=n_frames_abs, is_hpc=is_hpc)
    if run_absorption_collection:
        print("!!!!!!!!!!!!!!!!!!!! Parsing Absorbance !!!!!!!!!!!!!!!!!!!!")
        accumulate_abs_spectra(n_trajectories=n_snapshots_gs, n_frames=n_frames_abs)
        clean_up_abs(n_snapshots_gs, n_frames_abs)
    if run_excited_state:
        print("!!!!!!!!!!!!!!!!!!!! Running Excited States !!!!!!!!!!!!!!!!!!!!")
        # We take the original trajectory snapshots and run further trajectories
        # from those at the excited state
        verbosity = 3
        exc_state_init = exc_state_init_ex_param
        n_exc_states_propagate = n_exc_states_propagate_ex_param
        input_ceon.set_input(n_steps_exc, n_exc_states_propagate, n_steps_to_print_exc, exc_state_init,
                             verbosity=verbosity, time_step=time_step)
        run_ground_state_snapshots('nasqm_ground', 'nasqm_flu_', n_frames_gs, n_snapshots_ex, is_hpc,
                                   pmemd_available=False, ppn=processor_per_node)
    if run_fluorescence_collection:
        print("!!!!!!!!!!!!!!!!!!!! Parsing Fluorescences !!!!!!!!!!!!!!!!!!!!")
        exc_state_init = exc_state_init_ex_param
        accumulate_flu_spectra(n_trajectories=n_snapshots_ex, n_states=exc_state_init)

    # Restore Original Inputs
    if not is_hpc:
        open('input.ceon', 'w').write(input_ceon_bac)
        open('md_qmmm_amb.in', 'w').write(md_qmmm_amb)
        open('m1.inpcrd', 'w').write(m1_inpcrd)
    input_ceon.write_log()

    end_time = time.time()
    print("Job finished in %s seconds" % (end_time - start_time))

main()
# input_stream = open('naesmd.out', 'r')
# output_stream = open('naesmd_energies.txt', 'w')
# input_stream = open('nasqm_flu_16.out', 'r')
# output_stream = open('nasqm_energies16.txt', 'w')
# find_excited_energy(input_stream, output_stream, state=1)
