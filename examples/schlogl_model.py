from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)

u_t = D(u, x, 2) - (u - 1)*(u - 2)*(u - 3)

ti = time.time()
print(quadratize([(u, u_t)], 5, by_fun, 3))
print('time', time.time() - ti)