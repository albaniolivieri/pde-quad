from sympy import *
from sympy import Derivative as D
from functools import reduce

def get_order(set_derivs):
    max_order = 0
    for deriv in set_derivs:
        max_order = reduce(max, [der.args[1][1] for der in deriv.atoms(D)], max_order)
    return max_order
    