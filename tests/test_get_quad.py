from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.main import quadratize

# tests
t, x = symbols('t x')
u = Function('u')(t,x)
ut = u**2*D(u, x, 2)

# print(quadratize([(u, ut)], 2))

ut2 = u**3 * D(u, x, 3)
print(quadratize([(u, ut2)], 3))
