from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.main import quadratize
from algorithm.var_selection import *

#tests
t, x = symbols('t x')
v = Function('v')(t,x)
vinv = Function('vinv')(t, x)

v_t = 7 * D(v, x) * vinv - 5 * D(v, x)
quadratize([(v, v_t), (vinv, -v_t * vinv**2)], 4, sort_fun=by_fun)
