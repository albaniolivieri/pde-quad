from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.check_quad_examples import *
t, s = symbols('t s')
psi = Function('psi')(t,s)
theta = Function('theta')(t,s)
y_0 = Function('y_0')(t,s)
y_1 = Function('y_1')(t,s)

psi_t = D(psi, s, 2) - D(psi, s) - psi * y_0
theta_t = D(theta, s, 2) - D(theta, s) - theta - 1 + psi * y_0
y_0_t = y_0 * y_1**2 * theta_t
y_1_t = -y_1**2 * theta_t

get_quadratization([(psi, psi_t), (theta, theta_t), (y_0, y_0_t), (y_1, y_1_t)], [y_1**2, y_0 * y_1, y_0 * y_1**2, psi * y_0 * y_1, psi * y_0 * y_1**2, D(theta, s) * y_1], 5) 