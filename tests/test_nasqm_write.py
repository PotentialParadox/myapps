'''
Unit Tests for nasqm_write
'''
import os
import nasqm_write
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


def test_accumulate_abs_spectra_1():
    '''
    Test to see if capable of reading one trajectory, one frame, 3 states
    '''
    n_trajectories = 1
    n_frames = 1
    n_states = 3
    result = nasqm_write.accumulate_abs_spectra(n_trajectories, n_frames, n_states)
    assert result == "    2.89512919738290E+00    9.02843432872959E-01    " \
        "3.24727908968971E+00    2.43418572091824E-02    3.35012834746067E+00" \
        "    5.87942644022577E-04\n"

def test_accumulate_abs_spectra_2():
    '''
    Test to see if capable of reading one trajectory, two frames, two states
    '''
    n_trajectories = 1
    n_frames = 2
    n_states = 2
    result = nasqm_write.accumulate_abs_spectra(n_trajectories, n_frames, n_states)
    assert result == "    2.89512919738290E+00    9.02843432872959E-01    3.24727908968971E+00" \
        "    2.43418572091824E-02\n" \
        "    2.72735213287108E+00    9.69081788117023E-01    3.21395494210826E+00" \
        "    2.17477133910216E-04\n"

def test_accumulate_abs_spectra_3():
    '''
    Test to see if capable of reading two trajectory, two frames, one states
    '''
    n_trajectories = 2
    n_frames = 2
    n_states = 1
    result = nasqm_write.accumulate_abs_spectra(n_trajectories, n_frames, n_states)
    assert result == "    2.89512919738290E+00    9.02843432872959E-01\n" \
        "    2.72735213287108E+00    9.69081788117023E-01\n" \
        "    2.91123799741646E+00    6.66345055839272E-01\n" \
        "    2.92816920322448E+00    6.11902594708590E-01\n"