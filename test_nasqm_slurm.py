'''
A unit tester for nasqm_slurm.py
Things were tested on Hypergator 2 at the University of Florida
'''
import pytest
from nasqm_user_input import UserInput
import nasqm_slurm

@pytest.fixture
def userinput():
    '''
    Create a test user input
    '''
    user_input = UserInput()
    user_input.email = "dtracy.uf@gmail.com"
    user_input.email_options = 2
    user_input.number_nodes = 1
    user_input.processors_per_node = 16
    user_input.memory_per_node = "2000mb"
    user_input.walltime = "00:01:00"
    user_input.max_jobs = 4
    return user_input


def test_create_slurm_header(userinput):
    '''
    Test the simple create slurm header function, could
    break if someone messes with the user input class
    '''
    slurm_header = nasqm_slurm.create_slurm_header(userinput)
    comparison = {'n_nodes': 1, 'ppn': 16, 'email_options': 2,
                  'max_jobs': 4, 'email': 'dtracy.uf@gmail.com',
                  'walltime': '00:01:00', 'memory': '2000mb'}
    assert slurm_header == comparison


def test_build_command():
    '''
    Tests the build command
    '''
    restart_root = "ground_snap"
    output_root = "nasqm_abs_"
    n_trajectories = 16
    result = nasqm_slurm.build_command(restart_root, output_root, n_trajectories)
    test = open("tests/nasqm_slurm_build_command.txt", 'r').read()
    assert result == test


def test_slurm_trajectory_file_1(userinput):
    '''
    Tests to see if slurm trajectory_file is capable of running
    one trajectory
    '''
    restart_root = "ground_snap"
    output_root = "nasqm_abs_"
    title = "MyJob"
    n_trajectories = 1
    result = nasqm_slurm.slurm_trajectory_files(userinput, restart_root,
                                               output_root, title, n_trajectories)
    test = open("tests/nasqm_slurm_1.txt", 'r').read()
    assert result == (None, test)


def test_slurm_trajectory_file_16(userinput):
    '''
    Tests to see if slurm trajectory_file is capable of running
    one whole trajectory
    '''
    restart_root = "ground_snap"
    output_root = "nasqm_abs_"
    title = "MyJob"
    n_trajectories = 16
    result = nasqm_slurm.slurm_trajectory_files(userinput, restart_root,
                                               output_root, title, n_trajectories)
    test = open("tests/nasqm_slurm_16.txt", 'r').read()
    open("result.txt", 'w').write(result[0])
    assert result == (test, None)


def test_slurm_trajectory_file_33(userinput):
    '''
    Tests to see if slurm trajectory_file is capable of running
    multiple whole trajectories with remainder
    '''
    restart_root = "ground_snap"
    output_root = "nasqm_abs_"
    title = "MyJob"
    n_trajectories = 33
    result = nasqm_slurm.slurm_trajectory_files(userinput, restart_root,
                                               output_root, title, n_trajectories)
    test_0 = open("tests/nasqm_slurm_32.txt", 'r').read()
    test_1 = open("tests/nasqm_slurm_1.txt", 'r').read()
    open("result.txt", 'w').write(result[1])
    assert result == (test_0, test_1)
