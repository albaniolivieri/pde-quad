from quadratization import is_a_quadratization
from sympy import *
from sympy import Derivative as D
import numpy as np

def get_quadratization(vars_func: tuple, var_indep, new_var, n_diff):
    eq_t = vars_func[1]
    wt = D(new_var, t).doit().subs(D(Function(vars_func[0][0])(var_indep[0], var_indep[1]), var_indep[0]), eq_t) 
    
    refac = list(
        (D(Function(vars_func[0][0])(var_indep[0], var_indep[1]), var_indep[1], len(vars_func[0][i])-1), symbols(vars_func[0][i]))
        for i in range(1, len(vars_func[0])))
    refac.append((Function(vars_func[0][0])(var_indep[0], var_indep[1]), symbols(vars_func[0][0])))
    
    deriv_x = list((symbols(f'wx{i}'), D(new_var, var_indep[1], i).doit().subs(refac)) for i in range(1, n_diff+1))
        
    new_vars = deriv_x + [(symbols('w'), new_var.subs(refac))]
    vars = list(symbols(var) for var in vars_func[0])
    
    V = list((name, poly(exprs, vars)) for name, exprs in new_vars)
    V.extend(list((var, poly(var, vars)) for var in vars) + [(1, poly(1, vars))])
    
    deriv_t = [(symbols(f'{vars_func[0][0]}t'), poly(eq_t.subs(refac), vars)), (symbols('wt'), poly(wt.subs(refac), vars))]
    
    return is_a_quadratization(V, deriv_t)

#test 
t, x = symbols('t x')
u = Function('u')(t,x)
ut = u**2*D(u, x) + u
w0 = u**2
get_quadratization((['u', 'ux'], ut), (t,x), w0, 1)