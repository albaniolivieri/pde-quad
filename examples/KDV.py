from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)

ut = - D(u, x, 3) - 6 * u**2 * D(u,x)

funcs = [by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, ut)], 2, heur))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)