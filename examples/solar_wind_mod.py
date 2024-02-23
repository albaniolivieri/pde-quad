from sympy import *
from sympy import Derivative as D
import sys
import time
import statistics
sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

#tests
t, x = symbols('t x')
v = Function('v')(t,x)
vinv = Function('vinv')(t, x)

v_t = 7 * D(v, x) * vinv - 5 * D(v, x)
# quadratize([(v, v_t), (vinv, -v_t * vinv**2)], 4, sort_fun=by_fun)


funcs = [by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times = []
    for i in range(2):
        print(heur)
        ti = time.time()
        quadratize([(v, v_t), (vinv, -v_t * vinv**2)], 4, sort_fun=heur)
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)