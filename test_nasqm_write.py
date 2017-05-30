'''
Unit Tests for nasqm_write
'''
import nasqm_write
import numpy as np

def test_accumulate_flu_spectra():
    '''
    Tests 2 small nasqm_flu trajectories
    '''
    results = nasqm_write.accumulate_flu_spectra(n_trajectories=2)
    assert results == '    2.90923255131416E+00    9.08120295476811E-01\n'\
                      '    2.90923255131440E+00    9.08120295476795E-01\n'\
                      '    2.90576054156170E+00    9.12245042622513E-01\n'\
                      '    2.89718863282344E+00    9.17508693665317E-01\n'\
                      '    2.88542766911918E+00    9.23660450221756E-01\n'\
                      '    3.02432904582440E+00    7.59972440955212E-01\n'\
                      '    3.02432904608962E+00    7.59972440957533E-01\n'\
                      '    3.02035138173248E+00    7.69148450276112E-01\n'\
                      '    3.01028876491316E+00    7.78144025689657E-01\n'\
                      '    2.99529803383703E+00    7.85693947577494E-01\n'


def test_strip_timedelay_flu_spectra():
    '''
    Test to see if the timedelay strip routine removes the
    approriate data
    '''
    test_string = nasqm_write.accumulate_flu_spectra(n_trajectories=2)
    results = nasqm_write.strip_timedelay_flu_spectra(test_string, 2, 0.5, 1)
    assert results == '    2.90576054156170E+00    9.12245042622513E-01\n'\
                      '    2.89718863282344E+00    9.17508693665317E-01\n'\
                      '    2.88542766911918E+00    9.23660450221756E-01\n'\
                      '    3.02035138173248E+00    7.69148450276112E-01\n'\
                      '    3.01028876491316E+00    7.78144025689657E-01\n'\
                      '    2.99529803383703E+00    7.85693947577494E-01\n'


def test_numpy_to_spectra_string():
    '''
    Tests whether numpy to string converts properly using data from the
    previous test
    '''
    data = np.array([[2.90923255131416, 9.08120295476811E-01],
                     [2.90923255131440, 9.08120295476795E-01]])
    results = nasqm_write.numpy_to_specta_string(data)
    assert results == '    2.90923255131416E+00    9.08120295476811E-01\n'\
                      '    2.90923255131440E+00    9.08120295476795E-01\n'
