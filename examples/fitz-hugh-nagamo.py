from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, x = symbols('t x')
v = Function('v')(t,x)
y = Function('y')(t,x)
epsilon, h, gamma, r = symbols('epsilon h gamma r', constant=True)

v_t = epsilon * D(v, x, 2) - v**3 + 1.1*v**2 - 0.1*v - y + r
y_t = h * v - gamma * y + r

# ti = time.time()
# print(quadratize([(v, v_t), (y, y_t)], 5, by_fun, 3))
# print('time', time.time() - ti)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(v, v_t), (y, y_t)], 5, heur))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)

