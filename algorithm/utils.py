from sympy import *
from sympy import Derivative as D
from functools import reduce

def get_order(set_derivs):
    max_order = 0
    for deriv in set_derivs:
        max_order = reduce(max, [der.args[1][1] for der in deriv.atoms(D)], max_order)
    return max_order

def reduction(pol1, pol2):
    coef = pol1[1].coeff_monomial(pol2[2])
    if coef != 0:
        new_pol = pol1[1] - coef * pol2[1]
        return (pol1[0] - coef * pol2[0], new_pol, new_pol.LM())
    else: 
        return pol1
    
def reduction_sparse(pol1, pol2):
    if pol2[1] != 0:
        coef = pol1[1].coeff(pol2[2])
        if coef != 0:
            new_pol = pol1[1] - coef * pol2[1]
            return (pol1[0] - coef * pol2[0], new_pol, new_pol.leading_monom())
    return pol1