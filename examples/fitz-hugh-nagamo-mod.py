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

v_t = -D(v, x, 2)/v - v**2 -v + 5
y_t = v/y - y + 5

# ti = time.time()
# print(quadratize([(v, v_t), (y, y_t)], 5, by_fun, 3))
# print('time', time.time() - ti)

funcs = [by_order_degree] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(2):
        print(heur)
        ti = time.time()
        print(quadratize([(v, v_t), (y, y_t)], 5, heur, nvars_bound=11))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)

