from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.check_quad import *
t, x = symbols('t x')
c = Function('c')(t,x)
q = Function('q')(t,x)

ct = -D(q, t) - D(c, x) + c**2 * D(c, x, 2)
qt = q**2 - q
get_quadratization([(c, ct), (q, qt)], [c**2], 4)