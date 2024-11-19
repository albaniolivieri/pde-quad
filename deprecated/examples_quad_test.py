# In this file we test verifying a quadratization for some toy and practical PDE examples 
from sympy import symbols, Function
import time
from sympy import Derivative as D
from deprecated.quadratize_test import test_try_quadratize
from qupde.var_selection import by_order_degree, by_fun, by_degree_order

# examples

t, x = symbols('t x')
u = Function('u')(t, x)
v = Function('v')(t, x)

tests = []

search_alg_to_test = 'bnb' # 'near_neighbor' or 'bnb'

# tests
ut1_1 = u**2 * D(u, x, 2) + 2
vt1_2 = D(v, x, 2)/u**3 + u
tests.append(test_try_quadratize([(u, ut1_1), (v, vt1_2)], 3, by_order_degree, search_alg=search_alg_to_test))

vt = D(v, x, 1) * u - 2 * D(v, x, 1)
ut = - D(v, x, 1) * u**3 - 2 * D(v, x, 1) * u**2
tests.append(test_try_quadratize([(v, vt), (u, ut)], 2, by_fun, max_der_order=2, search_alg=search_alg_to_test))

ut2 = u**3 * D(u, x, 3)
tests.append(test_try_quadratize([(u, ut2)], 3, by_order_degree, search_alg=search_alg_to_test))

ut3 = D(u, x)**3 + u**3
tests.append(test_try_quadratize([(u, ut3)], 3, by_order_degree, search_alg=search_alg_to_test))

ut4 = D(u, x)**4
tests.append(test_try_quadratize([(u, ut4)], 3, by_degree_order, max_der_order=2, search_alg=search_alg_to_test))
  
ut5 = D(u, x)**3 * u
tests.append(test_try_quadratize([(u, ut5)], 3, by_order_degree, max_der_order=3, search_alg=search_alg_to_test))

ut6 = D(u, x)**3 
tests.append(test_try_quadratize([(u, ut6)], 3, by_order_degree, max_der_order=3, search_alg=search_alg_to_test))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))     
