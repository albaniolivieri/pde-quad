from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

#test
r, phi = symbols('r phi')
omega = symbols('omega', constant=True)
v = Function('v')(r,phi)

# v_r = (omega*D(v, phi)) / v

v_r = (1/4)*(D(v, phi)) / (omega*v)

# print(quadratize([(v, v_r)], n_diff=4, sort_fun=by_fun, first_indep=r))

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(v, v_r)], n_diff=4, sort_fun=heur, first_indep=r))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))
print('averages', avg)
print('standard deviations', std)