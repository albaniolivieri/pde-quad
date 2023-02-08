from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm import check_quad as quad
from algorithm.utils import get_order

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

def test_quad(func_eq, new_vars: list, n_diff: int):
    x_var = [symbol for symbol in func_eq[0][0].free_symbols if symbol != symbols('t')].pop()
    var_dic = [(symbols(f'w_{i}'), new_vars[i]) for i in range(len(new_vars))] 
    quad_vars = differentiate_x(x_var, new_vars, n_diff)
    undef_fun = [symbol for symbol, _ in func_eq]
    deriv_t = differentiate_t(func_eq, var_dic) + [(symbols(eqs[0].name + '_t'), eqs[1]) for eqs in func_eq] 
    max_order = max(get_order([der for _, der in deriv_t]), get_order([der for _, der in quad_vars]))
    
    refac = []
    for fun in undef_fun:
        refac += [(symbols(f'{fun.name}_{x_var}{i}'), D(fun, x_var, i)) 
                  for i in range(max_order, 0, -1)] + [(symbols(fun.name), fun)]
    refac += quad_vars
    exprs_orig = [expr for _, expr in deriv_t]
    results = quad.get_quadratization(func_eq, new_vars, n_diff)
    if not results[0]: return False 
    
    for i in range(len(exprs_orig)):
        print('passed eq:', results[1][i])
        if simplify(exprs_orig[i]) - simplify(results[1][i].rhs.subs(refac)) != 0:
            print('Test failed: expressions are not equal')
            print('equation: ', results[1][i])
            print('Original expression: ', simplify(exprs_orig[i]))
            print('Quad expression: ', simplify(results[1][i].rhs.subs(refac)))
            return False
    return True

t, x = symbols('t x')
u = Function('u')(t,x)

tests = []

# u_t = u**2 * ux + u
# w = u**2
# w_t = 2u * (u**2 * ux + u)
ut1 = u**2*D(u, x) + u
w01 = u**2
tests.append(test_quad([(u, ut1)], [w01], 1))

ut2 = u**2*D(u, x, 2)
w02 = u**2
tests.append(test_quad([(u, ut2)], [w02], 2))

ut3 = u * (D(u, x)**2 + u * D(u, x, 2))
w03 = u**2
tests.append(test_quad([(u, ut3)], [w03], 2))

ut4 = u * (3 * D(u, x) * D(u, x, 2) + u * D(u, x, 3) + 1)
w04 = u**2
tests.append(test_quad([(u, ut4)], [w04], 3))

#Dym 
ut5 = u**3 * D(u, x, 3)
tests.append(test_quad([(u, ut5)], [u**3, u * D(u, x)**2], 3))

u1 = Function('u1')(t,x)
u1t = u1**3 * D(u1, x, 1)
ut5 = u**3 * D(u, x, 3)
tests.append(test_quad([(u, ut5), (u1, u1t)], [u**3, u * D(u, x)**2, u1**3], 3))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))


    
    