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
v_1, v_2, v_3 = symbols('v_1 v_2 v_3', constant=True)
k = symbols('k', constant=True)

u_t = D(u, x, 2) - k*(u - v_1)*(u - v_2)*(u - v_3)

funcs = [by_order_degree, by_degree_order, by_fun] 
avg = []
std = []

for heur in funcs: 
    times= []
    for i in range(10):
        print(heur)
        ti = time.time()
        print(quadratize([(u, u_t)], 2, heur, search_alg = 'bnb'))
        times.append(time.time() - ti) 
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print('averages', avg)
print('standard deviations', std)
