from functools import reduce
from sympy import *
from sympy import Derivative as D
from .PolySys import *
from .utils import  diff_dict


def test_quadratization(func_eq, new_vars: list, n_diff: int): 
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    
    # check if there's a fraction in V2

    poly_syst = PolySys(func_eq, n_diff, x_var, new_vars)
    
    # convert back before returning
    return poly_syst.get_quad()

