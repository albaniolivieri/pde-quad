import sympy as sp
from qbee import *


if __name__ == '__main__':
    u, ux, uxx, uxxx = functions("u, ux, uxx, uxxx")
    system = [
        (u, ux**3),
        (ux, 3 * ux**2 * uxx),
        (uxx, 6 * ux * uxx**2 + 3 * ux**2 * uxxx)
    ]

    quadr_system = polynomialize_and_quadratize(system, input_der_orders={uxxx: 0})
    if quadr_system:
        print(quadr_system)