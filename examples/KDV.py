from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.check_quad import *
t, x = symbols('t x')
u = Function('u')(t,x)

ut = - D(u, x, 3) - 6 * u**2 * D(u,x)
get_quadratization([(u, ut)], [u**2], 1)