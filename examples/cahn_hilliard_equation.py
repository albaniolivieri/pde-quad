from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from qupde.quadratize import quadratize
from qupde.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)
epsilon = symbols('epsilon', constant=True)

u_t = 3 * u**2 * D(u, x) - D(u, x) - epsilon * D(u, x, 2)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times= []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, u_t)], 3, heur, search_alg='nn'))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)
