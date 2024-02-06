from sympy import *
from sympy import Derivative as D
import statistics
import sys
import time
sys.path.append("..")
from algorithm.quadratize import *
from algorithm.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)

ut = u**3 * D(u, x, 3)

#quadratize([(u, ut5)], 3)

funcs = [by_order_degree] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, ut)], 3, heur))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)