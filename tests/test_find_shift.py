import find_shift
import numpy as np

def test_find_ev_max1():
    '''
    Find Max in the Middle
    '''
    test_array = np.array([[1.1, 810.0, 0.0], [1.2, 799.0, 2.0], [1.3, 789.0, 1.0]])
    result = find_shift.find_ev_max(test_array)
    assert result == 1.2

def test_find_ev_max2():
    '''
    Find Max on the edge
    '''
    test_array = np.array([[1.1, 810.0, 3.0], [1.2, 799.0, 2.0], [1.3, 789.0, 1.0]])
    result = find_shift.find_ev_max(test_array)
    assert result == 1.1

def test_find_nm_max():
    '''
    Find Max in the Middle
    '''
    test_array = np.array([[1.1, 810.0, 0.0], [1.2, 799.0, 2.0], [1.3, 789.0, 1.0]])
    result = find_shift.find_nm_max(test_array)
    assert result == 799.0

def test_find_ev_shift1():
    '''
    Find the eV shift between two spectra
    '''
    test_array = np.array([[1.53, 810.0, 0, 1],
                           [1.55, 799.8, 1, 2],
                           [1.57, 789.7, 2, 1],
                           [1.59, 779.7, 3, 0]])
    unit = "ev"
    result = find_shift.find_shift(test_array, unit)
    assert result == -0.04


def test_find_ev_shift2():
    '''
    Find the eV shift between two spectra
    '''
    test_array = np.array([[1.53, 810.0, 3, 1],
                           [1.55, 799.8, 4, 2],
                           [1.57, 789.7, 2, 1],
                           [1.59, 779.7, 3, 0]])
    unit = "ev"
    result = find_shift.find_shift(test_array, unit)
    assert result == 0


def test_find_ev_shift3():
    '''
    Find the eV shift between two spectra
    '''
    test_array = np.array([[1.53, 810.0, 3, 1],
                           [1.55, 799.8, 4, 2],
                           [1.57, 789.7, 2, 8],
                           [1.59, 779.7, 3, 0]])
    unit = "ev"
    result = find_shift.find_shift(test_array, unit)
    assert result == 0.02


def test_find_nm_shift1():
    '''
    Find the eV shift between two spectra
    '''
    test_array = np.array([[1.53, 810.0, 0, 1],
                           [1.55, 799.8, 1, 2],
                           [1.57, 789.7, 2, 1],
                           [1.59, 779.7, 3, 0]])
    unit = "nm"
    result = find_shift.find_shift(test_array, unit)
    assert result == 20.1
