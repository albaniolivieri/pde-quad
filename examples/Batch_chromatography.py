from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
c = Function('c')(t,x)
q = Function('q')(t,x)
epsilon = symbols('epsilon', constant=True)
Pe = symbols('Pe', constant=True)
H1 = symbols('H1', constant=True)
K = symbols('K', constant=True)
L = symbols('L', constant=True)
Q = symbols('Q', constant=True)
Ac = symbols('Ac', constant=True)
ki = symbols('ki', constant=True)

qt = (L/epsilon*Ac) * ki * (2 * (H1*c)/(1 + K*c**2) - q)
ct = - ((1 - epsilon)/epsilon) * qt - D(c, x) + (1/Pe) * D(c, x, 2)

quadratize([(c, ct), (q, qt)], 4, sort_fun=by_fun, search_alg='bnb')