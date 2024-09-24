from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
a, b, gamma, d_uv, d_vu, d_u, d_v = symbols('a b gamma d_uv d_vu d_u d_v', constant=True)
v = Function('v')(t,x)
u = Function('u')(t,x)

u_t = d_u*D(u, x, 2) + d_uv * D(v, x, 2) + gamma*(a - u + u**2 * v)
v_t = d_v*D(v, x, 2) + d_vu * D(u, x, 2) + gamma*(b - u**2 * v)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times= []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, u_t), (v, v_t)], 3, heur, search_alg='nn'))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)
