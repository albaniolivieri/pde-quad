from sympy import *
from sympy import Derivative as D
import statistics
import sys
import time
sys.path.append("..")
from algorithm.quadratize import quadratize as quad_bb
from algorithm.quadratize_nearest_neighbor import quadratize as quad_nn
from algorithm.var_selection import *

t, x = symbols('t x')
u = Function('u')(t,x)

ut = u**3 * D(u, x, 3)

print(quad_nn([(u, ut)], 3, by_order_degree))

funcs = [by_order_degree, by_degree_order, by_fun] 
avg_bb = []
std_bb = []
avg_nn = []
std_nn = []

# for heur in funcs: 
#     times_bb = []
#     for i in range(2):
#         print(heur)
#         ti = time.time()
#         print(quad_bb([(u, ut)], 3, heur))
#         times_bb.append(time.time() - ti) 
#     avg_bb.append(statistics.mean(times_bb))
#     std_bb.append(statistics.stdev(times_bb))
#     times_nn = []
#     for i in range(2):
#         print(heur)
#         ti = time.time()
#         print(quad_nn([(u, ut)], 3, heur))
#         times_nn.append(time.time() - ti) 
#     avg_nn.append(statistics.mean(times_nn))
#     std_nn.append(statistics.stdev(times_nn))

print('averages branch-and-bound', avg_bb)
print('standard deviations branch-and-bound', std_bb)

print('averages nearest neighbor', avg_nn)
print('standard deviations nearest neighbor', std_nn)