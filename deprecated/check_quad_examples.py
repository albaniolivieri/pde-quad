from sympy import *
from sympy import Derivative as D
from ..algorithm.quadratization import is_a_quadratization
from ..algorithm.utils import get_order
from functools import reduce

def get_quadratization(func_eq, new_vars: list, n_diff: int):
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    
    vars_t = [(f'w_{i}', new_vars[i]) for i in range(len(new_vars))] 
    deriv_t = differentiate_t(func_eq, vars_t) + [(symbols(eqs[0].name + '_t'), eqs[1]) for eqs in func_eq] 
    
    quad_vars = differentiate_x(x_var, new_vars, n_diff)
    max_order = max(get_order([der for _, der in deriv_t]), get_order([der for _, der in quad_vars]))

    # ordering matters !
    refac = []
    for fun in undef_fun:
        refac += [(D(fun, x_var, i), symbols(f'{fun.name}_{x_var}{i}')) 
                  for i in range(max_order, 0, -1)] + [(fun, symbols(fun.name))]
    
    poly_vars = reduce(lambda name1, name2: name1 + name2[1].name + ',', refac, '')
    ring_polys = ring(poly_vars, QQ) 
    
    V = [(name, ring_polys[0].ring_new(exprs.subs(refac))) for name, exprs in quad_vars]
     
    V.extend([(var, ring_polys[0].ring_new(var)) for _, var in refac + [(1, 1)]])  
     
    deriv_t = list(map(lambda exprs: (exprs[0], ring_polys[0].ring_new(exprs[1].subs(refac))), deriv_t))

    return is_a_quadratization(V, deriv_t)

def differentiate_t(funcs_eqs, new_vars):
    deriv_t = []
    refac = [(D(deriv[0], symbols('t')), deriv[1]) for deriv in funcs_eqs]
    for i in range(len(new_vars)):
        wt = D(new_vars[i][1], symbols('t')).doit().subs(refac)
        deriv_t.append((symbols(f'{new_vars[i][0]}t'), wt.doit()))
    return deriv_t

def differentiate_x(var_indep, new_vars, n):
    quad_vars = []
    for i in range(len(new_vars)):
        quad_vars.extend([(symbols(f'w_{i}{var_indep}{j}'), D(new_vars[i], var_indep, j).doit()) 
                          for j in range(1, n + 1)] + [(symbols(f'w_{i}'), new_vars[i])])  
    return quad_vars