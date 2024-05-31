# In this file we test verifying a quadratization for some toy and practical PDE examples 
from sympy import symbols, Function
import time
from sympy import Derivative as D
from test_quadratize import test_try_quadratize
from algorithm.var_selection import by_order_degree, by_fun, by_degree_order

# examples

t, x = symbols('t x')
u = Function('u')(t, x)
v = Function('v')(t, x)

tests = []

search_alg_to_test = 'near_neighbor' # 'near_neighbor' or 'bnb'

t_i = time.time()
# tests
ut1_1 = u**2 * D(u, x, 2) + 2
vt1_2 = D(v, x, 2)/u**3 + 1.7
# print('\nu**2*D(u, x, 2) + 1/3')
tests.append(test_try_quadratize([(u, ut1_1), (v, vt1_2)], 3, by_order_degree, search_alg=search_alg_to_test))
# print('time elapsed: ', time.time() - t_i)

v = Function('v')(t,x)
u = Function('u')(t,x)
vt = D(v, x, 1) * u - 2 * D(v, x, 1)
ut = - D(v, x, 1) * u**3 - 2 * D(v, x, 1) * u**2
print('\nvt = D(v, x, 1) * u - 2 * D(v, x, 1)', '\nut = - D(v, x, 1) * u**3 - 2 * D(v, x, 1) * u**2')
tests.append(test_try_quadratize([(v, ut), (u, ut)], 2, by_fun, search_alg=search_alg_to_test))
# tests.append(test_try_quadratize([(v, ut), (u, ut)], 2, by_fun, nvars_bound=6, max_der_order=10))

ut2 = u**3 * D(u, x, 3)
print('\nu**3 * D(u, x, 3)')
tests.append(test_try_quadratize([(u, ut2)], 3, by_order_degree, search_alg=search_alg_to_test))

ut3 = D(u, x)**3 + u**3
print('\nD(u,x)**3 + u**3')
tests.append(test_try_quadratize([(u, ut3)], 3, by_order_degree, search_alg=search_alg_to_test))

ut4 = D(u, x)**4
print('\nD(u,x)**3 + u**3')
tests.append(test_try_quadratize([(u, ut4)], 3, by_degree_order, search_alg=search_alg_to_test))
  
ut5 = D(u, x)**3 * u
print('\nD(u,x)**3 * u')
tests.append(test_try_quadratize([(u, ut5)], 3, by_order_degree, max_der_order=3, search_alg=search_alg_to_test))

ut6 = D(u, x)**3 
print('\nD(u, x)**3')
tests.append(test_try_quadratize([(u, ut6)], 3, by_order_degree, max_der_order=3, search_alg=search_alg_to_test))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))     
