import sympy as sp
from qbee import *


if __name__ == '__main__':
    u, ux, uxx, uxxx = functions("u, ux, uxx, uxxx")
    system = [
        (u, ux * u**2 + u**3),
        (ux, 2 * u * ux**2 + u**2 * uxx + 3 * u**2 * ux),
        (uxx, 2 * ux**3 + 6 * ux * uxx * u + u**2 * uxxx + 3 * u**2 * uxx + 6 * u * ux**2)
    ]

    quadr_system = polynomialize_and_quadratize(system, input_der_orders={uxxx: 0})
    if quadr_system:
        print(quadr_system)