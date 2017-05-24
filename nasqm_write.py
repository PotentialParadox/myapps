"""
Functions that write outputs of the NASQM script
"""
import subprocess
import numpy as np
from amber_out import find_nasqm_excited_state, find_excited_energy

def accumulate_flu_spectra(n_trajectories, n_states=1, time_delay=0):
    """
    Create the spectra_flu.input file using the nasqm_flu_* files
    """
    output_stream = open('spectra_flu.input', 'w')
    for i in range(n_trajectories):
        amber_outfile = 'nasqm_flu_' + str(i+1) + ".out"
        input_stream = open(amber_outfile, 'r')
        find_nasqm_excited_state(input_stream, output_stream, n_states=n_states)
    output_stream.close()


def write_omega_vs_time(n_trajectories, n_states=1):
    '''
    Reads data from spectra_flu.input and writes the omegas vs time
    to omega_1_time.txt
    '''
    average_omegas_time = open('omega_1_time.txt', 'w')
    data = np.loadtxt('spectra_flu.input')
    n_rows_per_trajectory = int(data.shape[0] / n_trajectories)
    for i in range(n_rows_per_trajectory):
        omega = np.average(data[i::n_rows_per_trajectory, 0])
        average_omegas_time.write(str(omega) + '\n')
    average_omegas_time.close()


def write_nasqm_flu_energie(n_trajectories, n_states=1):
    '''
    Reads the data from the amber output files and writes the data
    to nasqm_flue_energies.txt
    '''
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
        energy = np.average(data[i::n_rows_per_trajectory])
        average_energies_time.write(str(energy) + '\n')


def accumulate_abs_spectra(is_tully, n_snapshots_gs, n_frames, n_states=20):
    '''
    Reads the data from the amber output files and writes the data to spectra_abs.input
    '''
    output_stream = open('spectra_abs.input', 'w')
    if is_tully:
        for traj in range(n_snapshots_gs):
            for frame in range(n_frames):
                amber_out = 'nasqm_abs_' + str(traj+1) + '_' + str(frame+1) + '.out'
                input_stream = open(amber_out, 'r')
                find_nasqm_excited_state(input_stream, output_stream, n_states)
                input_stream.close()
    else:
        for snap in range(n_snapshots_gs):
            amber_out = 'nasqm_abs_' + str(snap+1) + '.out'
            input_stream = open(amber_out, 'r')
            find_nasqm_excited_state(input_stream, output_stream, n_states)
            input_stream.close()
    output_stream.close()
