from functools import reduce
from sympy import *
from sympy import Derivative as D
from .quadratization import is_quadratization
from .utils import get_order, diff_dict


def get_quadratization(func_eq, new_vars: list, n_diff: int): 
    undef_fun = [symbol for symbol, _ in func_eq] 
    x_var = [symbol for symbol in undef_fun[0].free_symbols if symbol != symbols('t')].pop()
    
    max_order = get_order([expr for _, expr in func_eq])
    
    new_vars_pol, poly_syms, eqs_pol = build_ring(func_eq, n_diff, x_var, max_order, new_vars)
    dic_t, dic_x = get_dics(func_eq, poly_syms, eqs_pol, n_diff, max_order)

    return get_quad(dic_t, dic_x, new_vars_pol, eqs_pol, n_diff, x_var, poly_syms)

def build_ring(func_eq, order, var_indep, max_order, new_vars=None):
    refac = []
    for fun, _ in func_eq:
        refac += [(D(fun, var_indep, i), symbols(f'{fun.name}_{var_indep}{i}')) 
                  for i in range(max_order, 0, -1)] + [(fun, symbols(fun.name))] 
    
    poly_vars = []
    for fun, _ in func_eq:
        poly_vars.append(f'{fun.name}')
        poly_vars.extend([f'{fun.name}_{var_indep}{i}' for i in range(1, max_order + order + 1)])
        
    R, pol_sym = xring(poly_vars, QQ)
    
    expr_pol = [(symbols(f'{fun.name}_t'), R.ring_new(eq.subs(refac))) for fun, eq in func_eq]
    
    if new_vars != None: 
        vars_pol = [R.ring_new(new_vars[i].subs(refac)) for i in range(len(new_vars))]
        return vars_pol, pol_sym, expr_pol 
    
    return pol_sym, expr_pol

def get_dics(func_eq, symb, eqs_pol, order, max_order): 
    dic_x = {}
    dic_t = {}
    
    der_order = max_order + order
    for j in range(len(func_eq)):
        for i in range(der_order): 
            dic_x[symb[i + (der_order+1)*j]] = symb[i + (der_order+1)*j + 1]
    
    for k in range(len(func_eq)):
        for i in range(der_order - 1): 
            if i != 0:
                dic_t[symb[i + (der_order+1)*k]] = diff_dict(dic_t[symb[i - 1 + (der_order+1)*k]], dic_x)
            else: 
                dic_t[symb[(der_order+1)*k]] = eqs_pol[k][1]
        
    return dic_t, dic_x

def get_quad(dic_t, dic_x, new_vars, eqs_pol, order, var_indep, poly_syms):   
    new_vars_named = [(symbols(f'w_{i}'), pol) for i, pol in enumerate(new_vars)] 
    new_vars_t, new_vars_x = differentiate_dict(dic_t, dic_x, new_vars_named, order, var_indep)  
    deriv_t = new_vars_t + eqs_pol   
    V = [(1, list(dic_t.keys())[0].ring(1))] + [(symbols(f'{sym}'), sym) for sym in poly_syms] + new_vars_named + new_vars_x
    return is_quadratization(V, deriv_t)

def differentiate_dict(dic_t, dic_x, new_vars, order, var_indep):
    deriv_t = []
    deriv_x = []
    
    for name, expr in new_vars:
        deriv_t.append((symbols(f'{name}t'), diff_dict(expr, dic_t)))
    
    for name, expr in new_vars:
        for i in range(1, order + 1):
            deriv_x.append((symbols(f'{name}{var_indep}{i}'), diff_dict(expr, dic_x, i)))
            
    return deriv_t, deriv_x