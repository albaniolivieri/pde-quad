# In this file we test verifying a quadratization for some toy and practical PDE examples 
from sympy import symbols, Function
from sympy import Derivative as D
from test_quadratize import test_try_quadratize
from algorithm.var_selection import by_order_degree, by_fun, by_degree_order

# examples

t, x = symbols('t x')
u = Function('u')(t, x)

tests = []

# tests
# ut1 = u**2 * D(u, x, 2) + 1/3
# print('\nu**2*D(u, x, 2) + 1/3')
# tests.append(test_try_quadratize([(u, ut1)], 3, by_order_degree, search_alg='b&b'))

# v = Function('v')(t,x)
# u = Function('u')(t,x)
# vt = D(v, x, 1) * u - 2 * D(v, x, 1)
# ut = - D(v, x, 1) * u**3 - 2 * D(v, x, 1) * u**2
# print('\nvt = D(v, x, 1) * u - 2 * D(v, x, 1)', '\nut = - D(v, x, 1) * u**3 - 2 * D(v, x, 1) * u**2')
# tests.append(test_try_quadratize([(v, ut), (u, ut)], 2, by_fun, nvars_bound=6, max_der_order=10))

# ut2 = u**3 * D(u, x, 3)
# print('\nu**3 * D(u, x, 3)')
# tests.append(test_try_quadratize([(u, ut2)], 3, by_order_degree, search_alg='b&b'))

# ut3 = D(u, x)**3 + u**3
# print('\nD(u,x)**3 + u**3')
# tests.append(test_try_quadratize([(u, ut3)], 3, by_order_degree, search_alg='b&b'))

ut4 = D(u, x)**4
print('\nD(u,x)**3 + u**3')
tests.append(test_try_quadratize([(u, ut4)], 3, by_degree_order, search_alg='near_neighbor'))
  
# ut5 = D(u, x)**3 * u
# print('\nD(u,x)**3 * u')
# tests.append(test_try_quadratize([(u, ut5)], 3, by_order_degree, max_der_order=3, search_alg='b&b'))

# ut6 = D(u, x)**3 
# print('\nD(u, x)**3')
# tests.append(test_try_quadratize([(u, ut6)], 3, by_order_degree, max_der_order=3, search_alg='b&b'))

# # Summary
# print('\nTests passed: ', tests.count(True))
# print('Tests failed: ', tests.count(False))     
