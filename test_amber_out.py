'''
Unit tests for amber_out
'''
import io
import amber_out
import numpy as np

def test_find_excited_energy():
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
    amber_out.find_excited_energy(input_stream, output_stream, 1)
    output_string = output_stream.getvalue()
    output_stream.close()
    assert output_string == '   -4.00985112234037E+03\n'

def test_find_nasqm_excited_state():
    '''
    Tests to see if we can find the omega and total oscillator strength of multiple states
    '''
    test_string = "5 +++    3.772738023318104     0.23E-05 0.10E-04\n"\
                  "6 +++    3.904739638623995     0.39E-06 0.10E-04\n"\
                  "-------------------------------------------------\n"\
                  "\n"\
                  "Frequencies (eV) and Oscillator strengths (unitless)\n"\
                  "       Omega            fx              fy              fz          ftotal\n"\
                  " 1     3.04211901717783        0.814746171892645        0.209220955980193E-01  "\
                  "  0.103966885259747E-03    0.835772234375924\n"\
                  " 2     3.32932626530115        0.347146719311554E-02    0.129044151463899E-03  "\
                  "  0.381772659210805E-04    0.363868861050052E-02\n"\
                  "   3     3.69744792800296        0.155341793108174E-02    0.808029139455998E-03"\
                  "  0.709446925075805E-03    0.307089399561354E-02\n"
    input_stream = io.StringIO(test_string)
    output_stream = io.StringIO()
    amber_out.find_nasqm_excited_state(input_stream, output_stream, 2)
    output_string = output_stream.getvalue()
    output_stream.close()
    assert output_string == '    3.04211901717783E+00    8.35772234375924E-01'\
        '    3.32932626530115E+00    3.63868861050052E-03\n'



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
