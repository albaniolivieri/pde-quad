from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)
v_1, v_2, v_3 = symbols('v_1 v_2 v_3', constant=True)

u_t = D(u, x, 2) - (u - v_1)*(u - v_2)*(u - v_3)

# ti = time.time()
# print(quadratize([(u, u_t)], 5, by_fun, 3))
# print('time', time.time() - ti)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, u_t)], 2, heur))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)