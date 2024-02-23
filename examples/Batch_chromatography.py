from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.quadratize import *
t, x = symbols('t x')
c = Function('c')(t,x)
q = Function('q')(t,x)

ct = - D(c, x) + c**2 * D(c, x, 2)
qt = q**2 - q
quadratize([(c, ct), (q, qt)], 4, sort_fun=by_fun, nvars_bound=3, max_order=3)