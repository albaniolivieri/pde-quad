from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import by_fun, by_order_degree, by_degree_order, by_fun2

t, s = symbols('t s')
psi = Function('psi')(t,s)
theta = Function('theta')(t,s)
y_0 = Function('y_0')(t,s)
y_1 = Function('y_1')(t,s)

psi_t = D(psi, s, 2) - D(psi, s) - psi * y_0
theta_t = D(theta, s, 2) - D(theta, s) - theta - 1 + psi * y_0
y_0_t = y_0 * y_1**2 * theta_t
y_1_t = -y_1**2 * theta_t

print(quadratize([(psi, psi_t), (theta, theta_t), (y_0, y_0_t), (y_1, y_1_t)], 3, sort_fun=by_fun, nvars_bound=3, max_order=3))
