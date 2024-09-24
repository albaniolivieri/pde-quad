from sympy import *
from sympy import Derivative as D
import statistics
import sys
import time
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)

u_t = u**3 * D(u, x, 3)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times= []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, u_t)], 3, heur, search_alg = 'nn'))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)