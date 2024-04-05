# This file contains examples of PDEs with fractions
from sympy import symbols, Function
from sympy import Derivative as D
from test_quadratize import test_try_quadratize 
from algorithm.var_selection import by_order_degree, by_fun, by_degree_order

t, x = symbols('t x')
u = Function('u')(t,x)
v = Function('v')(t,x)
y = Function('y')(t,x)
u1 = symbols('u1')
tests = []

# Solar wind
ut = 7*D(u, x)/u - 5*D(u, x)
tests.append(test_try_quadratize([(u, ut)], 3, by_order_degree))

ut1 = 1 / (5*(u + 1)) 
tests.append(test_try_quadratize([(u, ut1)], 1, by_order_degree))

ut2 = 1 / (u + 1)**2
tests.append(test_try_quadratize([(u, ut2)], 3, by_order_degree))

ut3 = 1 / (u**2 + 1)
tests.append(test_try_quadratize([(u, ut3)], 3, by_order_degree))

ut4 = D(u, x) / (u + 1)
tests.append(test_try_quadratize([(u, ut4)], 3, by_order_degree))

ut5 = 1 / (u + 1)**2 + 1 / (u - 1)
tests.append(test_try_quadratize([(u, ut5)], 4, by_fun, nvars_bound=3, max_der_order=10))

# toy example from fitz-hugh-nagamo
# handle how to check if the results are equal
# v_t = -D(v, x, 2)/v - v**2 - v + 5
# y_t = v/y - y + 5
# tests.append(test_try_quadratize([(v, v_t), (y, y_t)], 3, by_fun, nvars_bound=3, max_der_order=10))

#Hard example
# ut6 = 1/(u + 1) + 1/u + D(u,x)*v
# vt1 = 1/(u * v) + 1/u 
# tests.append(test_try_quadratize([(u, ut6), (v, vt1)], 3, by_fun, nvars_bound=10))

#obs: we need to figure out how to check that results are equal
# ut5 = (D(u, x) + u**2) / ((u + 1)**2) + 1 / u
# tests.append(test_try_quadratize([(u, ut5)], 2, by_order_degree, nvars_bound=4, max_der_order=10))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))    
