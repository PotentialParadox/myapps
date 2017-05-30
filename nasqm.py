'''
Run NASQM
created by Dustin Tracy (dtracy.uf@gmail.com)
This program is used to automate NASQM job creations.
You'll find the parameters to change in the file nasqm_user_input.py
'''
import sys
import time
import subprocess
from amber import run_amber_parallel
from inputceon import InputCeon
from nasqm_write import (accumulate_abs_spectra, write_omega_vs_time,
                         write_nasqm_flu_energie, write_spectra_flu_input)
from nasqm_user_input import UserInput


def run_nasqm(root_name, coordinate_file=None, pmemd_available=False):
    '''
    Command line command to call nasqm
    '''
    cd_file = 'm1.inpcrd'
    if coordinate_file is not None:
        cd_file = coordinate_file
    amber = 'sander'
    if pmemd_available:
        amber = 'pmemd.cuda'
    subprocess.run([amber, '-O', '-i', 'md_qmmm_amb.in', '-o', root_name+'.out', '-c',
                    cd_file, '-p', 'm1.prmtop', '-r', root_name+'.rst', '-x', root_name+'.nc'])


def create_restarts(amber_input, output, step=None):
    '''
    Create amber restart files using cpptraj
    '''
    ctc = 'parm m1.prmtop\n'
    if step is not None:
        ctc += 'trajin ' + amber_input + '.nc' + ' 1 last ' + str(step) + "\n"
    else:
        ctc += 'trajin ' + amber_input + '.nc' + ' 1 last\n'
    ctc += 'center :1\n' \
           'trajout ' + output + ' restart\n' \
           'image\n' \
           'run\n' \
           'quit'
    open('convert_to_crd.in', 'w').write(ctc)
    subprocess.run(['cpptraj', '-i', 'convert_to_crd.in'])


def run_simulation_from_trajectory(nasqm_root, output_root, n_coordinates, n_snapshots,
                                   is_hpc, ppn):
    '''
    Run n_snapshots simulations using nasqm_root as the basis for the generation of the
    inital geometries. This will output data to output_root+str(i). Restart_step is the
    number of steps between the snapshots of the trajectory you are using as your geometries
    generator. Generally you will want to make the restart step to be the number of prints from
    this trajectory divided by the number of tajectories you wish to perform
    '''
    restart_step = int(n_coordinates / n_snapshots)
    create_restarts(amber_input=nasqm_root, output='ground_snap', step=restart_step)
    snap_restarts = []
    snap_trajectories = []
    if n_snapshots == 1:
        snap_restarts.append("ground_snap")
        snap_trajectories.append(output_root + '1')
    else:
        for i in range(n_snapshots):
            snap_restarts.append("ground_snap."+str(i+1))
            snap_trajectories.append(output_root + str(i + 1))
    # if is_hpc:
    #     run_hpc_trajectories(n_trajectories=n_snapshots, n_processor_per_node=ppn,
    #                          root_name=output_root)
    # else:
    pmemd_available = False
    run_amber_parallel(pmemd_available, snap_trajectories, snap_restarts, number_processors=ppn)


def run_abs_snapshots(output_root, n_trajectories, n_frames, is_hpc):
    '''
    Run snapshots on the nasqm_abs_* ground state trajectory runs
    '''
    pmemd_available = False # We require ESMD
    for i in range(n_trajectories):
        amber_restart = 'nasqm_abs_' + str(i+1)
        create_restarts(amber_input=amber_restart, output=amber_restart)
    # We will now have files that look like nasqm_abs_[trajectory]_[frame]
    # Lets run it
    snap_singles = []
    snap_restarts = []
    for traj in range(n_trajectories):
        for frame in range(n_frames):
            snap_singles.append("nasqm_abs_" + str(traj+1) + "_" + str(frame+1))
            snap_restarts.append("nasqm_abs_" + str(traj+1) + "." + str(frame+1))
    # if is_hpc:
    #     create_snapshot_slurm_script(script_file_name='run_abs_snapshot.sbatch',
    #                                  n_trajectories=n_trajectories, n_frames=n_frames,
    #                                  root_name=output_root, crd_file=output_root)
    #     print("Please now submit run_abs_snapshots.sbatch, then run abs collection")
    #     sys.exit('Slurm submission exception')
    # else:
    run_amber_parallel(pmemd_available, snap_singles, snap_restarts, number_processors=8)




def clean_up_abs(is_tully, n_trajectories, n_frame):
    '''
    Removes the files created by the absorption routine
    '''
    base_name = 'nasqm_abs_'
    if is_tully:
        for i in range(n_trajectories):
            for j in range(n_frame):
                traj = i + 1
                frame = j + 1
                subprocess.run(['rm', base_name + str(traj) + '_' + str(frame) + '.out'])
                subprocess.run(['rm', base_name + str(traj) + '_' + str(frame) + '.nc'])
                subprocess.run(['rm', base_name + str(traj) + '_' + str(frame) + '.rst'])
                subprocess.run(['rm', base_name + str(traj) + '.' + str(frame)])
    else:
        subprocess.run('rm nasqm_abs_*', shell=True)
        subprocess.run('rm ground_snap*', shell=True)


def run_ground_state_dynamics(input_ceon, user_input):
    '''
    Run the ground state trajectory that will be used to generate initial geometries
    for future calculations
    '''
    input_ceon.set_n_steps(user_input.n_steps_gs)
    input_ceon.set_exc_state_propagate(0)
    input_ceon.set_n_steps_to_print(user_input.n_steps_to_print_gs)
    input_ceon.set_exc_state_init(0)
    input_ceon.set_verbosity(0)
    input_ceon.set_time_step(user_input.time_step)
    input_ceon.set_random_velocities(False)
    coordinate_file = None
    if user_input.is_qmmm:
        coordinate_file = 'm1_md2.rst'
    run_nasqm('nasqm_ground', coordinate_file=coordinate_file, pmemd_available=False)


def run_absorption_trajectories(input_ceon, user_input):
    '''
    Now we want to take the original trajectory snapshots and run more trajectories
    using random velocities to make them different from each other
    '''
    print("!!!!!!!!!!!!!!!!!!!! Running Absorbance Trajectories !!!!!!!!!!!!!!!!!!!!")
    input_ceon.set_n_steps(user_input.n_steps_abs)
    input_ceon.set_exc_state_propagate(0)
    input_ceon.set_n_steps_to_print(user_input.n_steps_to_print_abs)
    input_ceon.set_exc_state_init(0)
    input_ceon.set_verbosity(0)
    input_ceon.set_time_step(user_input.time_step)
    input_ceon.set_random_velocities(True)
    run_simulation_from_trajectory('nasqm_ground', 'nasqm_abs_', user_input.n_frames_gs,
                                   user_input.n_snapshots_gs, user_input.is_hpc,
                                   ppn=user_input.processors_per_node)


def run_absorption_snapshots(input_ceon, user_input):
    '''
    Once the ground state trajectory file are made, we need can calculate the absorption
    snapshots.
    '''
    print("!!!!!!!!!!!!!!!!!!!! Running Absorbance Snapshots !!!!!!!!!!!!!!!!!!!!")
    input_ceon.set_n_steps(0)
    input_ceon.set_exc_state_propagate(user_input.n_abs_exc)
    input_ceon.set_n_steps_to_print(user_input.n_steps_to_print_abs)
    input_ceon.set_exc_state_init(0)
    input_ceon.set_verbosity(3)
    input_ceon.set_time_step(user_input.time_step)
    if user_input.is_tully:
        run_abs_snapshots(output_root='nasqm_abs_', n_trajectories=user_input.n_snapshots_gs,
                          n_frames=user_input.n_frames_abs, is_hpc=user_input.is_hpc)
    else:
        run_simulation_from_trajectory('nasqm_ground', 'nasqm_abs_', user_input.n_frames_gs,
                                       user_input.n_snapshots_gs, user_input.is_hpc,
                                       ppn=user_input.processors_per_node)


def run_absorption_collection(user_input):
    '''
    Parse the output data from amber for absorption energies and create a spectra_abs.input
    file
    '''
    print("!!!!!!!!!!!!!!!!!!!! Parsing Absorbance !!!!!!!!!!!!!!!!!!!!")
    accumulate_abs_spectra(user_input.is_tully, n_snapshots_gs=user_input.n_snapshots_gs,
                           n_frames=user_input.n_frames_abs, n_states=user_input.n_abs_exc)
    # clean_up_abs(user_input.is_tully, user_input.n_snapshots_gs, user_input.n_frames_abs)


def run_excited_state_trajectories(input_ceon, user_input):
    '''
    We take the original trajectory snapshots and run further trajectories
    from those at the excited state
    '''
    print("!!!!!!!!!!!!!!!!!!!! Running Excited States !!!!!!!!!!!!!!!!!!!!")
    input_ceon.set_n_steps(user_input.n_steps_exc)
    input_ceon.set_exc_state_propagate(user_input.n_exc_states_propagate_ex_param)
    input_ceon.set_n_steps_to_print(user_input.n_steps_to_print_exc)
    input_ceon.set_exc_state_init(user_input.exc_state_init_ex_param)
    input_ceon.set_verbosity(3)
    input_ceon.set_time_step(user_input.time_step)
    input_ceon.set_random_velocities(False)
    run_simulation_from_trajectory('nasqm_ground', 'nasqm_flu_', user_input.n_frames_gs,
                                   user_input.n_snapshots_ex, user_input.is_hpc,
                                   ppn=user_input.processors_per_node)


def run_fluorescence_collection(user_input):
    '''
    Parse the output data from amber for fluorescene energies and create a spectra_flu.input
    file
    '''
    print("!!!!!!!!!!!!!!!!!!!! Parsing Fluorescences !!!!!!!!!!!!!!!!!!!!")
    exc_state_init = user_input.exc_state_init_ex_param
    write_spectra_flu_input(user_input)
    write_omega_vs_time(n_trajectories=user_input.n_snapshots_ex, n_states=exc_state_init)
    write_nasqm_flu_energie(n_trajectories=user_input.n_snapshots_ex, n_states=exc_state_init)


def main():
    '''
    The primary nasqm automation function call. All changable parameters can be
    found in nasqm_user_input.py
    '''

    user_input = UserInput()

    # Copy inputs
    input_ceon_bac = open('input.ceon', 'r').read()
    md_qmmm_amb = open('md_qmmm_amb.in', 'r').read()
    m1_inpcrd = open('m1.inpcrd', 'r').read()

    # Create the input_ceon object
    input_ceon = InputCeon(amber_input='md_qmmm_amb.in')

    start_time = time.time()

    if user_input.run_ground_state_dynamics:
        run_ground_state_dynamics(input_ceon, user_input)
    if user_input.run_absorption_trajectories:
        run_absorption_trajectories(input_ceon, user_input)
    if user_input.run_absorption_snapshots:
        run_absorption_snapshots(input_ceon, user_input)
    if user_input.run_absorption_collection:
        run_absorption_collection(user_input)
    if user_input.run_excited_state_trajectories:
        run_excited_state_trajectories(input_ceon, user_input)
    if user_input.run_fluorescence_collection:
        run_fluorescence_collection(user_input)

    # Restore Original Inputs
    if not user_input.is_hpc:
        open('input.ceon', 'w').write(input_ceon_bac)
        open('md_qmmm_amb.in', 'w').write(md_qmmm_amb)
        open('m1.inpcrd', 'w').write(m1_inpcrd)
    input_ceon.write_log()

    end_time = time.time()
    print("Job finished in %s seconds" % (end_time - start_time))

main()

# stripper(n_ground_snaps=1000, n_atoms=48)
# input_stream = open('naesmd.out', 'r')
# output_stream = open('naesmd_energies.txt', 'w')
# input_stream = open('nasqm_flu_16.out', 'r')
# output_stream = open('nasqm_energies16.txt', 'w')
# find_excited_energy(input_stream, output_stream, state=1)
