from quadratization import is_a_quadratization
from sympy import *
from sympy import Derivative as D
import numpy as np

def get_quadratization(vars_func: tuple, new_vars: list, n_diff: int):
    ut = vars_func[1]
    undef_fun = vars_func[0] 
    sec_indep = list(vars_func[0].free_symbols)
    sec_indep.remove(symbols('t'))
    
    refac = [(D(undef_fun, sec_indep[0], i), symbols(f'{str(undef_fun)[0]}_{sec_indep[0]}{i}')) 
             for i in range(n_diff+2, 0, -1)] + [(undef_fun, symbols(str(undef_fun)[0]))]
    
    poly_vars = [name for _, name in refac]
    deriv_t = []
    quad_vars = []
    
    for i in range(len(new_vars)):
        wt = D(new_vars[i], symbols('t')).doit().subs(D(undef_fun, symbols('t')), ut)
        quad_vars.extend([(symbols(f'w_{i}{sec_indep[0]}{j}'), D(new_vars[i], sec_indep[0], j).doit().subs(refac)) 
                          for j in range(1, n_diff+1)] + [(symbols(f'w_{i}'), new_vars[i].subs(refac))])  
        deriv_t.append((symbols(f'w_{i}t'), poly(wt.doit().subs(refac), poly_vars)))
    
    V = [(name, poly(exprs, poly_vars)) for name, exprs in quad_vars] 
    V.extend([(var, poly(var, poly_vars)) for var in poly_vars] + [(1, poly(1, poly_vars))])    
    deriv_t.extend([(symbols(str(undef_fun)[0] + '_t'), poly(ut.subs(refac), poly_vars))])

    return is_a_quadratization(V, deriv_t)

#tests
t, x = symbols('t x')
u = Function('u')(t,x)

# u_t = u**2 * ux + u
# w = u**2
# w_t = 2u * (u**2 * ux + u)
ut1 = u**2*D(u, x) + u
w01 = u**2
#get_quadratization((u, ut1), [w01], 1)

# u_t = u**2 * uxx
# w = u**2
# w_t = 2u * (u**2 * uxx)
ut2 = u**2*D(u, x, 2)
w02 = u**2
#get_quadratization((u, ut2), [w02], 2)

# u_t = u * (ux**2 + u * uxx)
# w = u**2
# w_t = 2u**2 * (ux**2 + u * uxx)
ut3 = u * (D(u, x)**2 + u * D(u, x, 2))
w03 = u**2
#get_quadratization((u, ut3), [w03], 2)

# u_t = u * (2 ux * uxx + u * uxxx + 1)
# w = u**2
# w_t = 2 * u**2 * (2 ux * uxx + u * uxxx + 1)
ut4 = u * (3 * D(u, x) * D(u, x, 2) + u * D(u, x, 3) + 1)
w04 = u**2
#get_quadratization((u, ut4), [w04], 3)

ut5 = u**3 * D(u, x, 3)
get_quadratization((u, ut5), [u**3, u * D(u, x)**2], 3)
