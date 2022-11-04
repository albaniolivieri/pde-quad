from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm import check_quad_sparse as quad
from algorithm.utils import get_order


def test_quad(func_eq, new_vars: list, n_diff: int):
    x_var = [symbol for symbol in func_eq[0][0].free_symbols if symbol != symbols('t')].pop()
    var_dic = [(symbols(f'w_{i}'), new_vars[i]) for i in range(len(new_vars))] 
    quad_vars = quad.differentiate_x(x_var, new_vars, n_diff)
    undef_fun = [symbol for symbol, _ in func_eq]
    deriv_t = quad.differentiate_t(func_eq, var_dic) + [(symbols(eqs[0].name + '_t'), eqs[1]) for eqs in func_eq] 
    max_order = max(get_order([der for _, der in deriv_t]), get_order([der for _, der in quad_vars]))
    
    refac = []
    for fun in undef_fun:
        refac += [(symbols(f'{fun.name}_{x_var}{i}'), D(fun, x_var, i)) 
                  for i in range(max_order, 0, -1)] + [(symbols(fun.name), fun)]
    refac += quad_vars
    exprs_orig = [expr for _, expr in deriv_t]
    results = quad.get_quadratization(func_eq, new_vars, n_diff)
    
    for i in range(len(exprs_orig)):
        if simplify(exprs_orig[i]) - simplify(results[i].rhs.subs(refac)) != 0:
            print('Test failed: expressions are not equal')
            return False
    return True

t, x = symbols('t x')
u = Function('u')(t,x)

# u_t = u**2 * ux + u
# w = u**2
# w_t = 2u * (u**2 * ux + u)
ut1 = u**2*D(u, x) + u
w01 = u**2

print(test_quad([(u, ut1)], [w01], 1))

ut2 = u**2*D(u, x, 2)
w02 = u**2
print(test_quad([(u, ut2)], [w02], 2))

ut3 = u * (D(u, x)**2 + u * D(u, x, 2))
w03 = u**2
print(test_quad([(u, ut3)], [w03], 2))

ut4 = u * (3 * D(u, x) * D(u, x, 2) + u * D(u, x, 3) + 1)
w04 = u**2
print(test_quad([(u, ut4)], [w04], 3))

ut5 = u**3 * D(u, x, 3)
print(test_quad([(u, ut5)], [u**3, u * D(u, x)**2], 3))



    
    