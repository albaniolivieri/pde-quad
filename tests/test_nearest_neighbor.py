from sympy import *
from sympy import Derivative as D
import sys
import time
sys.path.append("..")
from algorithm.quadratize_nearest_neighbor import quadratize
from algorithm.var_selection import *


# r, phi = symbols('r phi')
# omega = symbols('omega', constant=True)
# v = Function('v')(r,phi)
# ti = time.time()
# v_r = (omega*D(v, phi)) / v
# print(quadratize([(v, v_r)], n_diff=4, sort_fun=by_fun, first_indep=r), 'time:', time.time() - ti)

# t, x = symbols('t x')
# u = Function('u')(t,x)
# ti = time.time()
# ut = - D(u, x, 3) + 6 * u**2 * D(u,x)
# print(quadratize([(u, ut)], n_diff=2, sort_fun=by_fun), 'time:', time.time() - ti)

# t, x = symbols('t x')
# beta = symbols('beta', constant=True)
# delta = symbols('delta', constant=True)
# gamma = symbols('gamma', constant=True)
# u = Function('u')(t,x)
# v = Function('v')(t,x)
# ti = time.time()
# ut = D(u, x) + (u*v - beta)*u + delta
# vt = D(v, x) - u**2*v + u + gamma*u
# print(quadratize([(u, ut), (v, vt)], n_diff=2, sort_fun=by_fun), 'time:', time.time() - ti)

# t, x = symbols('t x')
# a, b, gamma = symbols('a b gamma', constant=True)
# v = Function('v')(t,x)
# u = Function('u')(t,x)
# ti = time.time()
# u_t = D(u, x, 2) + D(v, x, 2) + gamma*a - u + u**2 * v
# v_t = D(v, x, 2) + D(u, x, 2) + gamma*b - u**2 * v
# print(quadratize([(v, v_t), (u, u_t)], 4, by_fun), 'time:', time.time() - ti)

t, s = symbols('t s')
psi = Function('psi')(t,s)
theta = Function('theta')(t,s)
y_0 = Function('y_0')(t,s)
y_1 = Function('y_1')(t,s)
ti = time.time()
psi_t = D(psi, s, 2) - D(psi, s) - psi * y_0
theta_t = D(theta, s, 2) - D(theta, s) - theta - 1 + psi * y_0
y_0_t = y_0 * y_1**2 * theta_t
y_1_t = -y_1**2 * theta_t

print(quadratize([(psi, psi_t), (theta, theta_t), (y_0, y_0_t), (y_1, y_1_t)], 3, sort_fun=by_fun), 
      'time:', time.time() - ti)
