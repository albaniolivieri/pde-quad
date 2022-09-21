import sympy as sp
from qbee import *


if __name__ == '__main__':
    u, ux, uxx, uxxx = functions("u, ux, uxx, uxxx")
    system = [
        (u, ux**2 * u),
        (ux, 2 * ux * uxx * u + ux**3),
        (uxx, 5 * ux**2 * uxx + 2 * uxx**2 * u + 2 * ux * uxxx * u)
    ]

    quadr_system = polynomialize_and_quadratize(system, input_der_orders={uxxx: 0})
    if quadr_system:
        print(quadr_system)