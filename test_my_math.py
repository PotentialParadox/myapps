'''
Unit tests for my_math
'''
import my_math

def test_quadratic_formula():
    '''
    Test reversed engineered
    '''
    result = my_math.quadratic_formula(2, -8, 6)
    assert result == (3.0, 1.0)
