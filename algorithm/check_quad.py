from quadratization import *
from sympy import *
import numpy as np

def get_quadratization(vars_func: tuple[list, sp.Eq], new_var, n_diff=1):
    ux = vars_func[1]
    w0 = symbols('w0')
    if poly(ux.rhs.subs(new_var, w0)).total_degree() <= 2: 
        diff_t = diff(new_var, t)
        diff_x = dict()
        for i in range(1, n_diff+1):
            diff_x[f"deriv_{i}"] = diff(new_var, i)
        V = list(diff_x.values()).append(new_var) + vars_func[0].append(1)
        return is_a_quadratization(V, diff_t)
    else: 
        print("It is not possible to get a quadratization with this variable")
        return False
    
#test 
u, t, x = symbols('u t x')


get_quadratization(([u, ux], ))