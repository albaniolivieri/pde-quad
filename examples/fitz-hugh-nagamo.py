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
y = Function('y')(t,x)

v_t = -D(v, x, 2) - v**3 + v**2 -v + 0.05
y_t = v - y + 0.05

ti = time.time()
print(quadratize([(v, v_t), (y, y_t)], 5, by_fun, 3))
print('time', time.time() - ti)