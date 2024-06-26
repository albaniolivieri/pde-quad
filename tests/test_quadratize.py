import sys
from sympy import symbols
from sympy import Derivative as D
from test_check_quad import test_quad

sys.path.append("..")

from algorithm.quadratize import quadratize
from algorithm.quadratize_nearest_neighbor import quadratize as quad_nearest_neighbor
from algorithm.utils import ring_to_expr

def test_try_quadratize(func_eq, n_diff, sort_fun, nvars_bound=5, max_der_order=None, search_alg='b&b'):
    if search_alg == 'near_neighbor':
        quad_prop, frac_vars, _ = quad_nearest_neighbor(func_eq, n_diff, sort_fun)
    elif search_alg == 'b&b': 
        quad_prop, frac_vars, _ = quadratize(
        func_eq, n_diff, sort_fun=sort_fun, nvars_bound=nvars_bound, max_der_order=max_der_order)
    
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    if not quad_prop and not frac_vars:
        print('Quadratization not found')
        return False
    quad_prop_expr = list(map(lambda x: ring_to_expr(x), quad_prop))
    frac_vars = list(map(lambda x: ring_to_expr(x), frac_vars))
    
    if test_quad(func_eq, quad_prop_expr, n_diff, frac_vars, vars_from_alg=True): 
        return True
    return False

