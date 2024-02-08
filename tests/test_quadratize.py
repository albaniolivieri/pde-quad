import sys
from sympy import symbols
from sympy import Derivative as D
from test_check_quad import test_quad

sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.utils import ring_to_expr, revert_frac_decomp 

def test_try_quadratize(func_eq, n_diff, sort_fun):
    quad_prop, frac_vars = quadratize(func_eq, n_diff, sort_fun)
    
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    quad_prop_expr = list(map(lambda x: ring_to_expr(x)[0], quad_prop))
    quad_prop = revert_frac_decomp(quad_prop_expr, frac_vars, quad=False)
    der_subs = []
    for fun, _ in func_eq:
        der_subs += [(symbols(f'{fun.name}_{x_var}{i}'), D(fun, x_var, i)) 
                for i in range(n_diff+2)] + [(symbols(fun.name), fun)] 
    quad_prop_expr = [expr.subs(der_subs) for expr in quad_prop_expr]
    
    if test_quad(func_eq, quad_prop_expr, n_diff): 
        return True
    return False

