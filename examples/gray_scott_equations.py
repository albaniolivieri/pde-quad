from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys
sys.path.append("..")
from qupde.quadratize import quadratize
from qupde.var_selection import *

t, x = symbols('t x')
e_1, e_2, F, k = symbols('epsilon_1 epsilon_2 F k', constant=True)
v = Function('v')(t,x)
u = Function('u')(t,x)

u_t = e_1 * D(u, x) - u * v**2 + F*(1 - u)
v_t = e_2 * D(v, x) + u * v**2 - (F + k) * v

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