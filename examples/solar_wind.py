from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

r, phi = symbols('r phi')
omega = symbols('omega', constant=True)
v = Function('v')(r,phi)

v_r = (omega*D(v, phi)) / v

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(v, v_r)], n_diff=2, sort_fun=heur, first_indep=r, search_alg='bnb'))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))
    
print('averages', avg)
print('standard deviations', std)