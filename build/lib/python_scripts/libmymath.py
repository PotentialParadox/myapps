'''
My personal math functions
'''
import math

def quadratic_formula(a, b, c):
    '''
    returns a tuple with + and - versions of the quadratic formula
    '''
    positive = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
    negative = (-b - math.sqrt(b**2 - 4*a*c))/(2*a)
    return (positive, negative)
