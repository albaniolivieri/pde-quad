from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize as quad_bb
from algorithm.quadratize_nearest_neighbor import quadratize as quad_nn
from algorithm.var_selection import *

t, x = symbols('t x')
rho = Function('rho')(t,x)
u = Function('u')(t,x)
p = Function('p')(t,x)

rho_t = -u * D(rho, x) - rho * D(u, x)
u_t = -u * D(u, x) - D(p, x) / rho
p_t = -p * D(u, x) - u * D(p, x)

# ti = time.time()
# print(quadratize([(rho, rho_t), (u, u_t), (p, p_t)], 5, by_fun, 3))
# print('time', time.time() - ti)

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
        print(quad_bb([(rho, rho_t), (u, u_t), (p, p_t)], 2, heur))
        times_bb.append(time.time() - ti) 
    avg_bb.append(statistics.mean(times_bb))
    std_bb.append(statistics.stdev(times_bb))
    times_nn = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quad_nn([(rho, rho_t), (u, u_t), (p, p_t)], 2, heur))
        times_nn.append(time.time() - ti) 
    avg_nn.append(statistics.mean(times_nn))
    std_nn.append(statistics.stdev(times_nn))

print('averages branch-and-bound', avg_bb)
print('standard deviations branch-and-bound', std_bb)

print('averages nearest neighbor', avg_nn)
print('standard deviations nearest neighbor', std_nn)