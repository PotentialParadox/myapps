'''
Units tests for the cpptraj wrappers for nasqm
'''
import pytest
import numpy as np
import nasqm_cpptraj
import nasqm_user_input

@pytest.fixture
def userinput():
    '''
    Create a test user input
    '''
    user_input = nasqm_user_input.UserInput()
    return user_input

def test_closest_script_3(userinput):
    '''
    Test to see if program is capable of building a script
    to include the nearest three solvent molecules
    '''
    userinput.number_nearest_solvents = 3
    snap_id = 1
    result = nasqm_cpptraj.closest_script(userinput, snap_id)
    test = open('tests/nasqm_cpptraj_nearest_three.txt', 'r').read()
    assert result == test

def test_closest_script_4():
    '''
    Test to see if program is capable of building a script
    to include the nearest four solvent molecules
    '''
    userinput.number_nearest_solvents = 4
    snap_id = 1
    result = nasqm_cpptraj.closest_script(userinput, snap_id)
    test = open('tests/nasqm_cpptraj_nearest_four.txt', 'r').read()
    assert result == test

def test_closest_script_5_snap_2():
    '''
    Test to see if it can build a script to include
    the nearest 5 but from the 2nd ground snapshot
    '''
    userinput.number_nearest_solvents = 5
    snap_id = 2
    result = nasqm_cpptraj.closest_script(userinput, snap_id)
    test = open('tests/nasqm_cpptraj_nearest_five.txt', 'r').read()
    assert result == test

def test_read_closest_3():
    '''
    Test to determine if the read closest can read
    the output from the closest three script
    '''
    closest_steam = open('tests/nasqm_cpptraj_closest_3.txt', 'r')
    result = nasqm_cpptraj.read_closest(closest_steam)
    np.testing.assert_array_equal(result, np.array([400, 456, 152]))

def test_read_closest_4():
    '''
    Test to determine if the read closest can read
    the output from the closest 4 script
    '''
    closest_steam = open('tests/nasqm_cpptraj_closest_4.txt', 'r')
    result = nasqm_cpptraj.read_closest(closest_steam)
    np.testing.assert_array_equal(result, np.array([420, 560, 252, 397]))
