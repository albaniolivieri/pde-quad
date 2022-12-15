from sympy import *
from sympy import Derivative as D
from .check_quad import get_dics, build_ring
from .utils import get_order
from .branch_and_bound import PolynomialSystem, bnb 
from .var_selection import *

def quadratize(func_eq, n_diff, sort_fun=by_order_degree, nvars_bound=10):
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    
    max_order = get_order([expr for _, expr in func_eq])
    
    poly_syms, eqs_pol = build_ring(func_eq, n_diff, x_var, max_order)
    dics = get_dics(func_eq, poly_syms, eqs_pol, n_diff, max_order)
    
    poly_syst = PolynomialSystem(dics, poly_syms, n_diff, x_var, eqs_pol)
    quad = bnb([], nvars_bound, poly_syst, sort_fun)
    
    return quad
    

    