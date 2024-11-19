# This file contains examples of PDEs with fractions
from sympy import symbols, Function
from sympy import Derivative as D
from deprecated.check_quad_test import test_quad 

t, x = symbols('t x')
u = Function('u')(t,x)
q_0 = symbols('q_0')
omega = symbols('omega', constant=True)

tests = []

# Solar wind modified
# ut = 7*D(u, x)/u - 5*D(u, x)

# tests.append(test_quad([(u, ut)], [D(u, x)/u, D(u,x)/u**2], 1, [(q_0, 1/u)]))
# tests.append(test_quad([(u, ut)], [1/u**2], 4, [(q_0, 1/u)]))
# tests.append(test_quad([(u, ut)], [], 4, [(q_0, 1/u)]))

# Solar wind 
ut = omega * D(u, x)/u

tests.append(test_quad([(u, ut)], [D(u,x)/u**2], 1, [(q_0, 1/u)]))
tests.append(test_quad([(u, ut)], [1/u**2], 4, [(q_0, 1/u)]))
tests.append(test_quad([(u, ut)], [], 4, [(q_0, 1/u)]))

# Euler equations
# rho = Function('rho')(t,x)
# u = Function('u')(t,x)
# p = Function('p')(t,x)

# rho_t = -u * D(rho, x) - rho * D(u, x)
# u_t = -u * D(u, x) - D(p, x) / rho
# p_t = -p * D(u, x) - u * D(p, x)

# tests.append(test_quad([(u, u_t), (rho, rho_t), (p, p_t)], new_vars = [], n_diff = 3, frac_vars = [(q_0, 1/rho)]))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))    