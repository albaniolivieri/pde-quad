import sys
import time

from sympy import *
from sympy import Derivative as D

sys.path.append("..")
from algorithm.main import quadratize

# tests
t, x = symbols('t x')
u = Function('u')(t,x)
ut = u**2*D(u, x, 2)
#print(quadratize([(u, ut)], 3))

r, p  = symbols('r p')
v = Function('v')(r,p)
#u = Function('u')(r,p)
#vr = D(v, p, 1) * u - 2 * D(v, p, 1)
#ur = - D(v, p, 1) * u**3 - 2 * D(v, p, 1) * u**2
#print(quadratize([(u, ur), (v, vr)], 4))

ut2 = u**3 * D(u, x, 3)
ti = time.time()
print(quadratize([(u, ut2)], 3))
print(time.time() - ti)

