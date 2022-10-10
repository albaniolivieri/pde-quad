from sympy import *
from sympy import Derivative as D
from quadratization import is_a_quadratization
from utils import get_order

def get_quadratization(func_eq: list[tuple], new_vars: list, n_diff: int):
    undef_fun = [symbol for symbol, _ in func_eq] 
    
    for fun_name in undef_fun: 
        indep = fun_name.free_symbols
        indep.remove(symbols('t'))
        x_var = indep.pop()
    
    vars_t = [(f'w_{i}', new_vars[i]) for i in range(len(new_vars))] 
    deriv_t = differentiate_t(func_eq, vars_t) + [(symbols(str(eqs[0]).split("(")[0] + '_t'), eqs[1]) for eqs in func_eq] 
    
    quad_vars = differentiate_x(x_var, new_vars, n_diff)
    
    max_order = max(get_order([der for _, der in deriv_t]), get_order([der for _, der in quad_vars]))

    # ordering matters !
    refac = []
    for fun in undef_fun:
        refac += [(D(fun, x_var, i), symbols(f'{str(fun).split("(")[0]}_{x_var}{i}')) 
                  for i in range(max_order, 0, -1)] + [(fun, symbols(str(fun).split("(")[0]))]

    poly_vars = [name for _, name in refac] 
    V = [(name, poly(exprs.subs(refac), poly_vars)) for name, exprs in quad_vars]
    V.extend([(var, poly(var, poly_vars)) for var in poly_vars + [1]])   
     
    deriv_t = list(map(lambda exprs: (exprs[0], poly(exprs[1].subs(refac), poly_vars)), deriv_t))

    return is_a_quadratization(V, deriv_t)

def differentiate_t(func_eq, new_vars):
    deriv_t = []
    refac = [(D(deriv[0], symbols('t')), deriv[1]) for deriv in func_eq]
    for i in range(len(new_vars)):
        wt = D(new_vars[i][1], symbols('t')).doit().subs(refac)
        deriv_t.append((symbols(f'{new_vars[i][0]}t'), wt.doit()))
    return deriv_t

def differentiate_x(func_var, new_vars, n):
    quad_vars = []
    for i in range(len(new_vars)):
        quad_vars.extend([(symbols(f'w_{i}{func_var}{j}'), D(new_vars[i], func_var, j).doit()) 
                          for j in range(1, n+1)] + [(symbols(f'w_{i}'), new_vars[i])])  
    return quad_vars

#tests
t, x = symbols('t x')
u = Function('u')(t,x)
u1 = Function('u1')(t,x)

# u_t = u**2 * ux + u
# w = u**2
# w_t = 2u * (u**2 * ux + u)
ut1 = u**2*D(u, x) + u
w01 = u**2
#get_quadratization([(u, ut1)], [w01], 1)

# u_t = u**2 * uxx
# w = u**2
# w_t = 2u * (u**2 * uxx)
ut2 = u**2*D(u, x, 2)
w02 = u**2
#get_quadratization([(u, ut2)], [w02], 2)

# u_t = u * (ux**2 + u * uxx)
# w = u**2
# w_t = 2u**2 * (ux**2 + u * uxx)
ut3 = u * (D(u, x)**2 + u * D(u, x, 2))
w03 = u**2
#get_quadratization([(u, ut3)], [w03], 2)

# u_t = u * (2 ux * uxx + u * uxxx + 1)
# w = u**2
# w_t = 2 * u**2 * (2 ux * uxx + u * uxxx + 1)
ut4 = u * (3 * D(u, x) * D(u, x, 2) + u * D(u, x, 3) + 1)
w04 = u**2
#get_quadratization([(u, ut4)], [w04], 3)

ut5 = u**3 * D(u, x, 3)
#get_quadratization([(u, ut5)], [u**3, u * D(u, x)**2], 3)

#get_quadratization([(u, ut5)], [u**3], 5)

u1t = u1**3 * D(u1, x, 1)

get_quadratization([(u, ut5), (u1, u1t)], [u**3, u * D(u, x)**2, u1**3], 3)