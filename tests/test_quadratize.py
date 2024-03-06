import sys
from sympy import symbols
from sympy import Derivative as D
from test_check_quad import test_quad

sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.utils import ring_to_expr

def test_try_quadratize(func_eq, n_diff, sort_fun, nvars_bound=5):
    quad_prop, frac_vars, _ = quadratize(func_eq, n_diff, sort_fun, nvars_bound)
    
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    quad_prop_expr = list(map(lambda x: ring_to_expr(x), quad_prop))
    frac_vars = list(map(lambda x: ring_to_expr(x), frac_vars))
    der_subs = []
    for fun, _ in func_eq:
        der_subs += [(symbols(f'{fun.name}_{x_var}{i}'), D(fun, x_var, i)) 
                for i in range(n_diff+2)] + [(symbols(fun.name), fun)] 
    quad_prop_expr = [expr.subs(der_subs) for expr in quad_prop_expr]
    
    frac_vars = [(q, 1/expr.subs(der_subs)) for q, expr in frac_vars]
    # frac_vars_subs = [(symbols(f'q_{i}'), 1/frac_vars[i]) for i in range(len(frac_vars))]
    # frac_vars = [1/var.subs(frac_vars_subs) for var in [expr.subs(der_subs) for expr in frac_vars]]
    
    if test_quad(func_eq, quad_prop_expr, n_diff, frac_vars): 
        return True
    return False

