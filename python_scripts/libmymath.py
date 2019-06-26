'''
My personal math functions
'''
import math
import scipy.fftpack
import numpy as np
import matplotlib.pyplot as plt

def quadratic_formula(a, b, c):
    '''
    returns a tuple with + and - versions of the quadratic formula
    '''
    positive = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
    negative = (-b - math.sqrt(b**2 - 4*a*c))/(2*a)
    return (positive, negative)

def myFourierTransform(data, timestep):
    nsamples = len(data)
    normalized_data = data - np.average(data)
    y = scipy.fftpack.fft(normalized_data)
    x = np.linspace(0.0, 1.0/(2.0*timestep), nsamples/2)
    y = 2.0/nsamples * np.abs(y[:nsamples//2])
    return x, y

def test():
    N = 600
    # sample spacing
    T = 1.0 / 800.0
    x = np.linspace(0.0, N*T, N)
    y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
    print("y:", len(y))
    print("x:", len(x))
    xf, yf = myFourierTransform(y, T)
    _, ax = plt.subplots()
    ax.plot(xf, yf)
    plt.show()
