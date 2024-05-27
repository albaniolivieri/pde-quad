from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
import math
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

#test
t, x = symbols('t x')
v = Function('v')(t,x)
u = Function('u')(t,x)

v_t = D(v, x) / v - round((math.pi), 5) * D(v, x)
u_t = 3.4 * D(u, x) * D(v, x)
print(quadratize([(v, v_t), (u, u_t)], 4, sort_fun=by_fun))