from sympy import *
from sympy import Derivative as D
from .quadratization import is_quadratization
from .utils import get_order, diff_dict
from functools import reduce

def get_quadratization(func_eq, new_vars: list, n_diff: int):
    
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    
    new_vars_pol, poly_syms, eqs_pol = build_ring(func_eq, n_diff, new_vars, x_var)
    
    dic_t, dic_x = get_dics(func_eq, poly_syms, eqs_pol, n_diff)

    return get_quad(dic_t, dic_x, new_vars_pol, eqs_pol, n_diff, x_var, poly_syms)

def build_ring(func_eq, order, new_vars, var_indep):
    # build ring for polys    
    refac = []
    for fun, expr in func_eq:
        refac += [(D(fun, var_indep, i), symbols(f'{fun.name}_{var_indep}{i}')) 
                  for i in range(get_order([expr]), 0, -1)] + [(fun, symbols(fun.name))] 
    
    poly_vars = []
    for fun, eq in func_eq:
        poly_vars.append(f'{fun.name}')
        poly_vars.extend([f'{fun.name}_{var_indep}{i}' for i in range(1, get_order([eq]) + order + 1)])
        
    R, pol_sym = xring(poly_vars, QQ)
    print('ring:', R, pol_sym)
    
    vars_pol = [(symbols(f'w_{i}'), R.ring_new(new_vars[i].subs(refac))) for i in range(len(new_vars))]
    expr_pol = [(symbols(f'{fun.name}_t'), R.ring_new(eq.subs(refac))) for fun, eq in func_eq]
    
    return vars_pol, pol_sym, expr_pol

def get_dics(func_eq, symb, eqs_pol, order): 
    # produce dict for derivatives 
    dic_x = {}
    dic_t = {}
    
    offset = 0
    for _, eq in func_eq:
        max_order = get_order([eq]) + order
        print('max_order', max_order)
        for i in range(max_order): 
            dic_x[symb[i + offset]] = symb[i + offset + 1]
        offset += max_order + 1
    
    offset = 0
    k = 0
    for _, eq in func_eq:
        max_order = get_order([eq]) + order
        for i in range(max_order): 
            if i != 0:
                dic_t[symb[i + offset]] = diff_dict(dic_t[symb[i - 1 + offset]], dic_x)
            else: 
                dic_t[symb[offset]] = eqs_pol[k][1]
        offset += max_order + 1
        k += 1
        
    return dic_t, dic_x

def get_quad(dic_t, dic_x, new_vars, eqs_pol, order, var_indep, poly_syms):
    # get derivs
    # build V 
    print('dic_t', dic_t)
    print('dic_x', dic_x)
    new_vars_t, new_vars_x = differentiate_dict(dic_t, dic_x, new_vars, order, var_indep)
    
    deriv_t = new_vars_t + eqs_pol
        
    V = [(1, new_vars[0][1].ring(1))] + [(symbols(f'{sym}'), sym) for sym in poly_syms] + new_vars + new_vars_x
  
    return is_quadratization(V, deriv_t)

def differentiate_dict(dic_t, dic_x, new_vars, order, var_indep):
    deriv_t = []
    deriv_x = []
    
    vars_name = [(f'w_{i}', new_vars[i]) for i in range(len(new_vars))] 
    
    for name, expr in vars_name:
        deriv_t.append((symbols(f'{name}t'), diff_dict(expr[1], dic_t)))
    
    for name, expr in vars_name:
        for i in range(1, order + 1):
            deriv_x.append((symbols(f'{name}{var_indep}{i}'), diff_dict(expr[1], dic_x, i)))
            
    return deriv_t, deriv_x

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