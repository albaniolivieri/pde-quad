from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.check_quad import *
t, x = symbols('t x')
u = Function('u')(t,x)

ut5 = u**3 * D(u, x, 3)
get_quadratization([(u, ut5)], [u**3, u**2*D(u, x)**3], 6)

