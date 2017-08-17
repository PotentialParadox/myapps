'''
Unit tests for amber_out
'''
import io
import os
import pytest
import amber_out
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


def test_find_excited_energies():
    '''
    Only tests the first state
    FIXME need test for multiple states
    '''
    test_string = "   1     2.72992400735139        -3.72988788959956       -0.503715618713691"\
                "0.101289783571188E-01     14.1658956897201\n"\
                "Total energy of the ground state (eV,AU)\n" \
                "0  -4012.58104634772       -147.459579686974\n"\
                "Total energies of excited states (eV,AU)\n"\
                "1  -4009.85112234037       -147.359256866811\n"\
                "QMMM:"

    input_stream = io.StringIO(test_string)
    output_stream = io.StringIO()
    amber_out.find_excited_energies(input_stream, output_stream, [1])
    output_string = output_stream.getvalue()
    output_stream.close()
    assert output_string == '   -4.00985112234037E+03\n'


def test_find_nasqm_excited_state_1():
    '''
    Tests to see if we can find the omega and total oscillator strength of multiple states
    '''
    input_stream = open("nasqm_flu_1.out")
    result = amber_out.find_nasqm_excited_state(input_stream)
    assert result == "    2.90923255131416E+00    9.08120295476811E-01\n" \
        "    2.90923255131440E+00    9.08120295476795E-01\n" \
        "    2.90576054156170E+00    9.12245042622513E-01\n" \
        "    2.89718863282344E+00    9.17508693665317E-01\n" \
        "    2.88542766911918E+00    9.23660450221756E-01\n"


def test_find_dipole():
    '''
    Make sure that we can get the mm dipoles
    '''
    test_string = "Ewald error estimate:   0.1964E-03" \
                  "-------------------------------------------------"\
                  "-----------------------------\n"\
                  "\n"\
                  "                  X        Y        Z     TOTAL  \n"\
                  "  MM DIPOLE    -6.487   81.075    4.482   81.458\n"\
                  "                  X        Y        Z     TOTAL  \n"\
                  "  QM DIPOLE       NaN      NaN      NaN      NaN\n"\
                  " QMMM: No. QMMM Pairs per QM atom:          289\n"\
                  "Ewald error estimate:   0.1964E-03" \
                  "-------------------------------------------------"\
                  "-----------------------------\n"\
                  "\n"\
                  "                  X        Y        Z     TOTAL  \n"\
                  "  MM DIPOLE    -6.487   81.075    4.482   87.879\n"\
                  "                  X        Y        Z     TOTAL  \n"\
                  "  QM DIPOLE       NaN      NaN      NaN      NaN\n"\
                  " QMMM: No. QMMM Pairs per QM atom:          289\n"
    input_stream = io.StringIO(test_string)
    result = amber_out.find_dipoles(input_stream)
    np.testing.assert_array_equal(result, np.array([81.458, 87.879]))


def test_find_total_energies():
    '''
    Tests the search for the total energies
    '''
    input_stream = open("nasqm_flu_1.out")
    result = amber_out.find_total_energies(input_stream)
    np.testing.assert_array_equal(result, np.array([167.2595, 167.3586, 167.9726,
                                                    168.5343, 168.3066]))

def test_find_mulliken_ground():
    '''
    Tests the search for the ground mulliken charges
    '''
    input_stream = open("mulliken.out")
    state = 0
    result = amber_out.find_mulliken(input_stream, state)
    answer = np.array([-0.138, -0.119, -0.123, -0.043, -0.122, -0.120,
                       -0.114, -0.120, -0.031, 0.132, 0.140, 0.131, 0.134,
                       0.132, 0.121, 0.122, -0.121, -0.116, -0.118, -0.122,
                       -0.031, -0.119, -0.114, -0.037, 0.134, 0.128, 0.142,
                       0.119, 0.113, 0.116, -0.115, -0.132, -0.124, -0.132,
                       -0.122, 0.138, 0.130, 0.139, 0.126, 0.135])
    np.testing.assert_array_equal(result, answer)

def test_find_mulliken_first():
    '''
    Tests the search for the ground mulliken charges
    '''
    input_stream = open("mulliken.out")
    state = 1
    result = amber_out.find_mulliken(input_stream, state)
    answer = np.array([-0.138, -0.121, -0.122, -0.041, -0.122,
                       -0.120, -0.114, -0.113, -0.038, 0.133, 0.140,
                       0.132, 0.134, 0.132, 0.122, 0.122, -0.121,
                       -0.120, -0.122, -0.119, -0.040, -0.110, -0.114,
                       -0.035, 0.133, 0.126, 0.141, 0.117, 0.112, 0.118,
                       -0.115, -0.131, -0.123, -0.131, -0.122, 0.139,
                       0.131, 0.140, 0.127, 0.135])
    np.testing.assert_array_equal(result, answer)


def test_find_number_excited_states():
    input_stream = open("nexmd.out", 'r')
    result = amber_out.find_number_excited_states(input_stream)
    answer = 3
    assert result == answer
