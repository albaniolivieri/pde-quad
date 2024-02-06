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

v_t = 7 * D(v, x) / v - 5 * D(v, x)
quadratize([(v, v_t)], 4, sort_fun=by_fun)