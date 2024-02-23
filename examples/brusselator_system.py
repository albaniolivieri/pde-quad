from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
beta = symbols('beta', constant=True)
delta = symbols('delta', constant=True)
gamma = symbols('gamma', constant=True)
u = Function('u')(t,x)
v = Function('v')(t,x)

ut = D(u, x) + (u*v - beta)*u + delta
vt = D(v, x) - u**2*v + u + gamma*u

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, ut), (v, vt)], 2, heur))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)