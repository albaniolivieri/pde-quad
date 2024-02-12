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
tests.append(test_quad([(u, ut)], [1/u, 1/u**2], 4, [1/u]))

# Summary
print('\nTests passed: ', tests.count(True))
print('Tests failed: ', tests.count(False))    