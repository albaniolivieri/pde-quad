# This file contains examples of PDEs with fractions
from sympy import symbols, Function
from sympy import Derivative as D
from test_quadratize import test_try_quadratize 
from algorithm.var_selection import by_order_degree, by_fun, by_degree_order

t, x = symbols('t x')
u = Function('u')(t,x)
u1 = symbols('u1')
tests = []

# Solar wind
ut = 7*D(u, x)/u - 5*D(u, x)
tests.append(test_try_quadratize([(u, ut)], 3, by_order_degree))

ut1 = 1 / (5*(u + 1)) 
tests.append(test_try_quadratize([(u, ut1)], 3, by_order_degree))

ut2 = 1 / (u + 1)**2
tests.append(test_try_quadratize([(u, ut2)], 3, by_order_degree))

ut3 = 1 / (u**2 + 1)
tests.append(test_try_quadratize([(u, ut3)], 3, by_order_degree))

ut4 = D(u, x) / (u + 1)
tests.append(test_try_quadratize([(u, ut4)], 3, by_order_degree))

# Hard example
# ut5 = (D(u,x) + u**2) / ((u+1)**2) + 1/u
# tests.append(test_try_quadratize([(u, ut5)], 3, by_order_degree))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))    