from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
a, b, gamma = symbols('a b gamma', constant=True)
v = Function('v')(t,x)
u = Function('u')(t,x)

u_t = D(u, x, 2) + D(v, x, 2) + gamma*a - u + u**2 * v
v_t = D(v, x, 2) + D(u, x, 2) + gamma*b - (1/2)*u**2 * v

# u_t = D(u, x, 2) + D(v, x, 2) + 1 - u + u**2 * v
# v_t = D(v, x, 2) + D(u, x, 2) + 1 - (1/2)*u**2 * v

# ti = time.time()
# print(quadratize([(v, v_t), (u, u_t)], 5, by_fun))
# print('time', time.time() - ti)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(v, v_t), (u, u_t)], 3, heur))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)