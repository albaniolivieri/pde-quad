from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
v = Function('v')(t,x)
u = Function('u')(t,x)

u_t = D(u, x, 2) + D(v, x, 2) + 1 - u + u**2 * v
v_t = D(v, x, 2) + D(u, x, 2) + 1 - u**2 * v

ti = time.time()
print(quadratize([(v, v_t), (u, u_t)], 5, by_fun, 3))
print('time', time.time() - ti)