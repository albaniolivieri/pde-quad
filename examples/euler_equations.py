from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
rho = Function('rho')(t,x)
u = Function('u')(t,x)
p = Function('p')(t,x)

rho_t = -u * D(rho, x) - rho * D(u, x)
u_t = -u * D(u, x) - D(p, x) / rho
p_t = -p * D(u, x) - u * D(p, x)

ti = time.time()
print(quadratize([(rho, rho_t), (u, u_t), (p, p_t)], 5, by_fun, 3))
print('time', time.time() - ti)