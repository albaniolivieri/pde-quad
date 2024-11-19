from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from qupde.quadratize import quadratize
from qupde.var_selection import *

t, x = symbols('t x')
d_1, d_2, a, b = symbols('d_1 d_2 a b', constant=True)
l = symbols('lambda', constant=True)
u = Function('u')(t,x)
v = Function('v')(t,x)

u_t = d_1 * D(u, x) + l * (1 - (b + 1) * u + b * u**2 * v)
v_t = d_2 * D(v, x) + l * a**2 * (u - u**2 * v)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times= []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, u_t), (v, v_t)], 3, heur, search_alg='bnb'))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)