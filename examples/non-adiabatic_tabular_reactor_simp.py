from sympy import *
from sympy import Derivative as D
import sys
sys.path.append("..")
from algorithm.main import quadratize
t, s = symbols('t s')
psi = Function('psi')(t,s)
theta = Function('theta')(t,s)
y_0 = Function('y_0')(t,s)
y_1 = Function('y_1')(t,s)

psi_t = D(psi, s, 2) - D(psi, s) - psi * (theta**3 + 6 * theta**2 + theta + 1)
theta_t = D(theta, s, 2) - D(theta, s) - theta - 1 + psi * (theta**3 + 6 * theta**2 + theta + 1)


quadratize([(psi, psi_t), (theta, theta_t)], 5) 