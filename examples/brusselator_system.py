from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from algorithm.quadratize import quadratize as quad_bb
from algorithm.quadratize_nearest_neighbor import quadratize as quad_nn
from algorithm.var_selection import *

t, x = symbols('t x')
d_1, d_2, a, b = symbols('d_1 d_2 a b', constant=True)
l = symbols('lambda', constant=True)
u = Function('u')(t,x)
v = Function('v')(t,x)

ut = d_1 * D(u, x) + l * (1 - (b + 1) * u + b * u**2 * v)
vt = d_2 * D(v, x) + l * a**2 * (u - u**2 * v)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg_bb = []
std_bb = []
avg_nn = []
std_nn = []

for heur in funcs: 
    times_bb = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quad_bb([(u, ut), (v, vt)], 3, heur))
        times_bb.append(time.time() - ti) 
    avg_bb.append(statistics.mean(times_bb))
    std_bb.append(statistics.stdev(times_bb))
    times_nn = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quad_nn([(u, ut), (v, vt)], 3, heur))
        times_nn.append(time.time() - ti) 
    avg_nn.append(statistics.mean(times_nn))
    std_nn.append(statistics.stdev(times_nn))

print('averages branch-and-bound', avg_bb)
print('standard deviations branch-and-bound', std_bb)

print('averages nearest neighbor', avg_nn)
print('standard deviations nearest neighbor', std_nn)