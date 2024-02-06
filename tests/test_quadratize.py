import sys
import time

from sympy import *
from sympy import Derivative as D

sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import by_order_degree

# Gleb: would be great to have:
#  - more test cases
#  - automatic check that the produced thing is a quadratization (as in test_check_quad.py)
#  - summary at the end

ti = time.time()
# tests
t, x = symbols('t x')
u = Function('u')(t,x)
ut = u**2*D(u, x, 2)
print(quadratize([(u, ut)], 3, by_order_degree))

r, p  = symbols('r p')
v = Function('v')(r,p)
#u = Function('u')(r,p)
#vr = D(v, p, 1) * u - 2 * D(v, p, 1)
#ur = - D(v, p, 1) * u**3 - 2 * D(v, p, 1) * u**2
#print(quadratize([(u, ur), (v, vr)], 4))

ut2 = u**3 * D(u, x, 3)
ut3 = D(u, x)**3
# print(quadratize([(u, ut3)], 2))
print(time.time() - ti)

