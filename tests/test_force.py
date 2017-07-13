'''
Unit tests for force
'''
import os
import pytest
import force
import numpy as np

def setup_module(module):
    '''
    Switch to test directory
    '''
    os.chdir("tests")

def teardown_module(module):
    '''
    Return to main directory
    '''
    os.chdir("..")

def test_find_force_1_1():
    '''
    Test the force between two atoms
    '''
    file_1 = "force_1_1.txt"
    file_2 = "force_2_1.txt"
    result = force.find_forces(file_1, file_2)
    answer = np.array([[0.00044788, -0.00079723, -0.00026873]])
    np.testing.assert_almost_equal(result, answer, 8)

def test_find_force_1_2():
    '''
    Test the force interaction between one atom and two others
    '''
    file_1 = "force_1_1.txt"
    file_2 = "force_2_2.txt"
    result = force.find_forces(file_1, file_2)
    answer = np.array([[0.00079267, -0.00205524, -0.00040851]])
    np.testing.assert_almost_equal(result, answer, 8)

def test_find_force_2_1():
    '''
    Test the force interaction between two atoms of interest and
    one solvent atom
    '''
    file_1 = "force_1_2.txt"
    file_2 = "force_2_1.txt"
    result = force.find_forces(file_1, file_2)
    answer = np.array([[0.00044788, -0.00079723, -0.00026873],
                       [0.00052826, -0.0008254, -0.00033016]])
    np.testing.assert_almost_equal(result, answer, 8)

def test_find_force_2_2():
    '''
    Test the force interaction between two atoms of interest and
    two solvent atoms
    '''
    file_1 = "force_1_2.txt"
    file_2 = "force_2_2.txt"
    result = force.find_forces(file_1, file_2)
    answer = np.array([[0.00079267, -0.00205524, -0.00040851],
                       [0.0009055, -0.00212957, -0.00049183]])
    np.testing.assert_almost_equal(result, answer, 8)
