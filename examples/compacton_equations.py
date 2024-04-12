from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize_nearest_neighbor import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)

u_t = -5 * u**4 * D(u, x) - 6 * D(u, x)**3 - 12 * u * D(u, x) * D(u, x, 2) \
    - 6 * u * D(u, x) * D(u, x, 2) - 3 * u**2 * D(u, x, 3)

ti = time.time()
print(quadratize([(u, u_t)], 5, by_degree_order))
print('time', time.time() - ti)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(5):
        print(heur)
        ti = time.time()
        print(quadratize([(u, u_t)], 3, heur))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)