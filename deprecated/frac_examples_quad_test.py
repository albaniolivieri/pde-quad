# This file contains examples of PDEs with fractions
from sympy import symbols, Function
from sympy import Derivative as D
import time
import math
from deprecated.quadratize_test import test_try_quadratize
from qupde.mon_heuristics import by_order_degree, by_fun, by_degree_order

t, x = symbols("t x")
omega = symbols("omega", constant=True)
u = Function("u")(t, x)
v = Function("v")(t, x)
y = Function("y")(t, x)
u1 = symbols("u1")
tests = []

# ut = 7.15666*D(u, x)/u - 5.677*D(u, x)
tests.append(test_try_quadratize([(u, ut)], 3, by_order_degree))

v_t = D(v, x) / v - round((math.pi), 5) * D(v, x)
u_t = 3.4 * D(u, x) * D(v, x)
tests.append(test_try_quadratize([(v, v_t), (u, u_t)], 4, by_fun))

ut1 = 1 / (5 * (u + 1))
tests.append(test_try_quadratize([(u, ut1)], 1, by_order_degree))

ut2 = 1 / (0.6 * u + 1) ** 2
tests.append(test_try_quadratize([(u, ut2)], 3, by_order_degree))

ut3 = 1 / (u**2 + 1)
tests.append(test_try_quadratize([(u, ut3)], 3, by_order_degree))

ut4 = D(u, x) / (u + 1)
tests.append(test_try_quadratize([(u, ut4)], 3, by_order_degree))

ut5 = 1 / (u + 1) ** 2 + 1 / (u - 1)
tests.append(
    test_try_quadratize([(u, ut5)], 4, by_fun, nvars_bound=3, max_der_order=10)
)

ut6 = 1 / (omega * u**2) + 0.5 / (u - 1)
tests.append(
    test_try_quadratize([(u, ut6)], 4, by_fun, nvars_bound=3, max_der_order=10)
)

ut8 = 1 / ((u**2) - 0.5 * u + 1)
tests.append(
    test_try_quadratize([(u, ut8)], 4, by_fun, nvars_bound=3, max_der_order=10)
)

# toy example from fitz-hugh-nagumo (hard example: finds a quadratization of order 4)
v_t = -D(v, x, 2) / v - v**2 - v + 5
y_t = v / y - y + 5
tests.append(
    test_try_quadratize(
        [(v, v_t), (y, y_t)],
        2,
        by_degree_order,
        nvars_bound=10,
        max_der_order=10,
        search_alg="near_neighbor",
    )
)

# Hard example (finds a quadratization of order 4)
ut6 = 1 / (u + 1) + 1 / u + D(u, x) * v
vt1 = 1 / (u * v) + 1 / u
tests.append(
    test_try_quadratize(
        [(u, ut6), (v, vt1)],
        n_diff=3,
        sort_fun=by_fun,
        nvars_bound=10,
        max_der_order=10,
        search_alg="near_neighbor",
    )
)

ut7 = (D(u, x) + u**2) / (3.5 * (u + 1) ** 2) + 1 / u
tests.append(
    test_try_quadratize(
        [(u, ut7)],
        2,
        by_order_degree,
        nvars_bound=5,
        max_der_order=10,
        search_alg="near_neighbor",
    )
)

# Summary
print("\nTests passed: ", tests.count(True))
print("Tests failed: ", tests.count(False))
