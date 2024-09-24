from sympy import *
from sympy import Derivative as D
import time
import statistics
import sys

sys.path.append("..")
from algorithm.quadratize import quadratize
from algorithm.var_selection import *

t, s = symbols("t s")
psi = Function("psi")(t, s)
theta = Function("theta")(t, s)
Pe = symbols("Pe", constant=True)
B = symbols("B", constant=True)
D_cte = symbols("D", constant=True)
beta = symbols("beta", constant=True)
theta_ref = symbols("theta_ref", constant=True)
c_0, c_1, c_2, c_3 = symbols("c_0 c_1 c_2 c_3", constant=True)

psi_t = (
    (1 / Pe) * D(psi, s, 2)
    - D(psi, s)
    - D_cte * (psi * (c_3 * theta**3 + c_2 * theta**2 + c_1 * theta + c_0))
)
theta_t = (
    (1 / Pe) * D(theta, s, 2)
    - D(theta, s)
    - beta * (theta - theta_ref)
    + B * D_cte * psi * (c_3 * theta**3 + c_2 * theta**2 + c_1 * theta + c_0)
)

funcs = [by_order_degree, by_degree_order, by_fun]
avg = []
std = []

for heur in funcs:
    times = []
    for i in range(2):
        ti = time.time()
        quadratize(
            [(psi, psi_t), (theta, theta_t)],
            n_diff=2,
            nvars_bound=5,
            sort_fun=by_degree_order,
            max_der_order=3,
            search_alg="bnb",
        )
        times.append(time.time() - ti)
    avg.append(statistics.mean(times))
    std.append(statistics.stdev(times))

print("averages", avg)
print("standard deviations", std)
