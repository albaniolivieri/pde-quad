from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

#test
t, x = symbols('t x')
v = Function('v')(t,x)

v_t = D(v, x) / v - D(v, x)
print(quadratize([(v, v_t)], 4, sort_fun=by_fun))