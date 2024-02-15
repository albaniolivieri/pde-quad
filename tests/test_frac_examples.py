# This file contains examples of PDEs with fractions
from sympy import symbols, Function
from sympy import Derivative as D
from test_check_quad import test_quad 

t, x = symbols('t x')
u = Function('u')(t,x)

tests = []

# Solar wind
ut = 7*D(u, x)/u - 5*D(u, x)
tests.append(test_quad([(u, ut)], [D(u, x)/u, D(u,x)/u**2], 3, [1/u]))
tests.append(test_quad([(u, ut)], [1/u**2], 4, [1/u]))
tests.append(test_quad([(u, ut)], [], 4, [1/u]))

# Euler equations
rho = Function('rho')(t,x)
u = Function('u')(t,x)
p = Function('p')(t,x)

rho_t = -u * D(rho, x) - rho * D(u, x)
u_t = -u * D(u, x) - D(p, x) / rho
p_t = -p * D(u, x) - u * D(p, x)

tests.append(test_quad([(u, u_t), (rho, rho_t), (p, p_t)], [], 3, [1/rho]))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))    