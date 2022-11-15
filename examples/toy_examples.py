from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm import check_quad as quad

#tests
t, x = symbols('t x')
u = Function('u')(t,x)

# u_t = u**2 * ux + u
# w = u**2
# w_t = 2u * (u**2 * ux + u)
ut1 = u**2*D(u, x) + u
w01 = u**2
quad.get_quadratization([(u, ut1)], [w01], 1)

# u_t = u**2 * uxx
# w = u**2
# w_t = 2u * (u**2 * uxx)
ut2 = u**2*D(u, x, 2)
w02 = u**2
quad.get_quadratization([(u, ut2)], [w02], 2)

# u_t = u * (ux**2 + u * uxx)
# w = u**2
# w_t = 2u**2 * (ux**2 + u * uxx)
ut3 = u * (D(u, x)**2 + u * D(u, x, 2))
w03 = u**2
quad.get_quadratization([(u, ut3)], [w03], 2)

# u_t = u * (2 ux * uxx + u * uxxx + 1)
# w = u**2
# w_t = 2 * u**2 * (2 ux * uxx + u * uxxx + 1)
ut4 = u * (3 * D(u, x) * D(u, x, 2) + u * D(u, x, 3) + 1)
w04 = u**2
quad.get_quadratization([(u, ut4)], [w04], 3)