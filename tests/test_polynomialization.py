import sys
from sympy import symbols, Function, sin, exp
from sympy import Derivative as D

sys.path.append("..")

from algorithm.polynomialization import polynomialize

x,t = symbols('x t')
u = Function('u')


#ut = sin(u(t,x))
#ut = u(t,x)**3.4
ut = sin(exp(u(t,x)))

print(polynomialize([(u(t,x), ut)]))